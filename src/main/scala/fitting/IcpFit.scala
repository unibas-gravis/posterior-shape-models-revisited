package fitting

import breeze.linalg.{DenseMatrix, DenseVector}
import norms.L2norm
import scalismo.common.PointId
import scalismo.geometry.{EuclideanVector, EuclideanVector3D, Point, Point3D, _3D}
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.{MultivariateNormalDistribution, PointDistributionModel}
import scalismo.transformations.{RigidTransformation, TranslationAfterRotationSpace3D}
import scalismo.{ModelRecenter, ModelUtils}
import utility.{MathHelp, Tchange}

object IcpFit {
  /**
   *
   * @param sampling
   * @param optAlignment
   * @param maxIterations
   * @param sigma2
   * @param sigma2InitFactor
   * @param stopTol
   * @param perpendicularFactor
   * @param init careful -> used after landmark alignment
   * @param initPose careful -> this is applied after the ldms alignment
   * @tparam A
   */
  case class IcpSettings[A <: SamplingStrategy](sampling: A, optAlignment: Boolean = true, maxIterations: Int = 25, sigma2: Double = 1.0, sigma2InitFactor: Double = 1000.0, stopTol: Double = 1e-2, perpendicularFactor: Double = 100.0, init: Option[DenseVector[Double]] = None, initPose: Option[RigidTransformation[_3D]] = None)

  /**
   * applies multiple non rigid icp optimizations.
   */
  def apply[A <: SamplingStrategy : Manifest](modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[A], ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty): (TriangleMesh[_3D], Int) = {
    val it = getIterator(modelOrig, target, settings, ldms)
    val results = it.take(settings.maxIterations).toIndexedSeq
    (results.last._1, results.count(_._2==false)) //TODO think about returning -1 if not converged. also see onlineIt
  }

  def getIterator[A <: SamplingStrategy : Manifest](modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[A], ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty): Iterator[(TriangleMesh[_3D], Boolean)] = {
    settings.sampling match {
      case _: SamplingStrategyNormals => iteratorSurfaceAware(modelOrig, target, settings.asInstanceOf[IcpSettings[SamplingStrategyNormals]], ldms)
      case _: SamplingStrategy => iterator(modelOrig, target, settings.asInstanceOf[IcpSettings[SamplingStrategy]], ldms)
    }
  }

  def iteratorWithIndex(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[SamplingStrategy], ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty): Iterator[(TriangleMesh[_3D], Boolean, Int)] = {
    val l2 = L2norm[_3D]()

    // TODO remove variable
    var model = if (ldms.nonEmpty) { //with ldms given, rigidly align once with ldms
      val transform = ModelUtils.alignLdmsToLdms(ldms.map(t => (modelOrig.mean.pointSet.point(t._1), t._2)))._2
      modelOrig.transform(transform) //TODO also one regression?
    } else modelOrig
    val initshape = settings.init.getOrElse(DenseVector.zeros[Double](model.rank))
    val init = model.instance(initshape).transform(settings.initPose.getOrElse(TranslationAfterRotationSpace3D.apply(Point3D.origin).identityTransformation))
    Iterator.iterate((init, false, 0))(t => {
      val (currentMesh, converged, i) = t
      if (converged) (t._1,t._2,t._3+1) else {
        val obs = settings.sampling.establishCorrespondence(currentMesh, target)
        val current = Tchange.getDef(currentMesh.pointSet, model.mean)
        val correctedObs = if (settings.optAlignment) {
          //TODO only an approximation as the weights are not regarded for non uniform strategies
          val curCor = obs.map(t => (t._1, t._2)) ++ (if (ldms.nonEmpty) ldms else IndexedSeq.empty)
          val tf = ModelUtils.alignShapeToLdms(currentMesh, curCor)._2
          model = model.transform(tf)
          settings.sampling.establishCorrespondence(currentMesh.transform(tf), target)
        } else obs

        val cursigma2 = getCurrentSigma2(i, settings.maxIterations, settings.sigma2, settings.sigma2 * settings.sigma2InitFactor)
        val dv = DenseVector.zeros[Double](3)
        val dcov = DenseMatrix.eye[Double](3) * cursigma2
        val td = correctedObs.map(t => (t._1, t._2, MultivariateNormalDistribution(dv, dcov * t._3)))
        val tdldm = td ++ (if (ldms.nonEmpty) {
          val mvn = Tchange.getZeroMeanNoise(settings.sigma2)
          ldms.map(t => (t._1,t._2,mvn))
        } else IndexedSeq.empty)
        val pmesh = model.posterior(tdldm).mean
        //val pmesh = FitScripts.getSurfaceAwarePosterior(model, t)
        val next = Tchange.getDef(pmesh.pointSet, model.mean)
        val change = l2.norm2Vector(next.zip(current).map(t => t._1 - t._2))
        if (change < settings.stopTol) (pmesh, true, i + 1) else (pmesh, false, i + 1)
      }
    })
  }

  private def iterator(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[SamplingStrategy], ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty): Iterator[(TriangleMesh[_3D], Boolean)] = {
    iteratorWithIndex(modelOrig, target, settings, ldms).map(t => (t._1, t._2))
  }

  def fitWithOnlineAlignment(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[SamplingStrategy], withRotation: Boolean = false, mainAxis: EuclideanVector[_3D], alignmentSettings: Option[SamplingStrategyUniform] = Option(TargetSamplingUnique()), ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty): (TriangleMesh[_3D], Int) = {
    val results = fitWithOnlineAlignmentIterator(modelOrig, target, settings, withRotation, mainAxis, alignmentSettings, ldms).take(settings.maxIterations).toIndexedSeq
    (results.last._1, results.count(_._2==false))
  }

  def fitWithOnlineAlignmentIteratorWithIndex(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[SamplingStrategy], withRotation: Boolean = false, mainAxis: EuclideanVector[_3D], alignmentSettings: Option[SamplingStrategyUniform] = Option(TargetSamplingUnique()), ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty): Iterator[(TriangleMesh[_3D], Boolean, Int)] = {
    val l2 = L2norm[_3D]()
    val axis = MathHelp.listOrthogonalVectors(mainAxis)

    //TODO remove variables
    var (model, posedModel) = if (ldms.nonEmpty) { //with ldms given, rigidly align once with ldms
      val transform = ModelUtils.alignLdmsToLdms(ldms.map(t => (modelOrig.mean.pointSet.point(t._1), t._2)))._2
      val m = modelOrig.transform(transform) //TODO also one regression?
      (m,m)
    } else (modelOrig, modelOrig)
    val initshape = settings.init.getOrElse(DenseVector.zeros[Double](model.rank))
    val init = model.instance(initshape).transform(settings.initPose.getOrElse(TranslationAfterRotationSpace3D.apply(Point3D.origin).identityTransformation))
    Iterator.iterate((init, false, 0))(t => {
      val (currentMesh, converged, i) = t
      if (converged) (t._1, t._2, t._3 + 1) else {
        val aobs = alignmentSettings match {
          case Some(a) => a.establishCorrespondenceUniform(currentMesh, target)
          case None => settings.sampling.establishCorrespondence(currentMesh, target).map(t => (t._1, t._2)) //TODO now ideal if non uniform strategy provided -> should just remove this freedom
        }
        model = if (withRotation) { //basis align model
          ModelRecenter.completeAlignPdm(posedModel, aobs.map(_._1), axis)
        } else ModelRecenter.recenterPdm(posedModel, aobs.map(_._1))
        val obs = settings.sampling.establishCorrespondence(currentMesh, target)
        val current = Tchange.getDef(currentMesh.pointSet, model.mean)
        val correctedObs = if (settings.optAlignment) { //pose align model
          val curCor = obs.map(t => (t._1, t._2)) ++ (if (ldms.nonEmpty) ldms else IndexedSeq.empty)
          val tf = ModelUtils.alignShapeToLdms(currentMesh, curCor)._2
          model = model.transform(tf)
          posedModel = posedModel.transform(tf)
          settings.sampling.establishCorrespondence(currentMesh.transform(tf), target)
        } else obs
        val cursigma2 = getCurrentSigma2(i, settings.maxIterations, settings.sigma2, settings.sigma2*settings.sigma2InitFactor)
        val dv = DenseVector.zeros[Double](3)
        val dcov = DenseMatrix.eye[Double](3) * cursigma2
        val td = correctedObs.map(t => (t._1, t._2, MultivariateNormalDistribution(dv, dcov * t._3)))
        val tdldm = td ++ (if (ldms.nonEmpty) {
          val mvn = Tchange.getZeroMeanNoise(settings.sigma2)
          ldms.map(t => (t._1,t._2,mvn))
        } else IndexedSeq.empty)
        val pmesh = model.posterior(tdldm).mean
        val next = Tchange.getDef(pmesh.pointSet, model.mean)
        val change = l2.norm2Vector(next.zip(current).map(t => t._1 - t._2))
        if (change < settings.stopTol) (pmesh, true, i + 1) else (pmesh, false, i + 1)
      }
    })
  }

  def fitWithOnlineAlignmentIterator(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[SamplingStrategy], withRotation: Boolean = false, mainAxis: EuclideanVector[_3D], alignmentSettings: Option[SamplingStrategyUniform] = Option(TargetSamplingUnique()), ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty): Iterator[(TriangleMesh[_3D], Boolean)] = {
    fitWithOnlineAlignmentIteratorWithIndex(modelOrig, target, settings, withRotation, mainAxis, alignmentSettings, ldms).map(t=>(t._1, t._2))
  }

  def iteratorSurfaceAware(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[SamplingStrategyNormals], ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty): Iterator[(TriangleMesh[_3D], Boolean)] = {
    fitSurfaceAwareIteratorWithIndex(modelOrig, target, settings, ldms=ldms).map(t => (t._1, t._2))
  }

  //TODO think about better way of having lookaheadRegression that provides a better path to show performance
  def fitSurfaceAwareIteratorWithIndex(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[SamplingStrategyNormals], mainAxis: EuclideanVector[_3D]=EuclideanVector3D.unitZ, useRealingment: Boolean=true, withRotation: Boolean = true, alignmentSettings: Option[SamplingStrategyUniform] = Option(TargetSamplingUnique()), ldms: IndexedSeq[(PointId, Point[_3D])]=IndexedSeq.empty, lookAheadRegression: Boolean=false): Iterator[(TriangleMesh[_3D], Boolean, Int)] = {
    val l2 = L2norm[_3D]()
    val aset = alignmentSettings.getOrElse(settings.sampling)
    val axis = MathHelp.listOrthogonalVectors(mainAxis)

    def getObs(m: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D], EuclideanVector[_3D], Double)] = {
      settings.sampling.establishCorrespondenceNormal(m, target)
    }

    //TODO remove variables
    var (model, posedModel) = if (ldms.nonEmpty) { //with ldms given, rigidly align once with ldms
      val transform = ModelUtils.alignLdmsToLdms(ldms.map(t => (modelOrig.mean.pointSet.point(t._1), t._2)))._2
      val m = modelOrig.transform(transform) //TODO also one regression?
      (m,m)
    } else (modelOrig, modelOrig)
    val initshape = settings.init.getOrElse(DenseVector.zeros[Double](model.rank))
    val init = model.instance(initshape).transform(settings.initPose.getOrElse(TranslationAfterRotationSpace3D.apply(Point3D.origin).identityTransformation))
    val it = Iterator.iterate((init, false, 0, modelOrig.mean))(t => {
      val (currentMesh, converged, i, currentMeshLookAhead) = t
      if (converged) (currentMesh, converged, i+1, currentMeshLookAhead) else {
        model = if (useRealingment) {
          val aobs = aset.establishCorrespondence(currentMesh, target)
          if (withRotation) { //basis align model
            ModelRecenter.completeAlignPdm(posedModel, aobs.map(_._1), axis)
          } else ModelRecenter.recenterPdm(posedModel, aobs.map(_._1))
        } else model
        val obs = getObs(currentMesh)
        val current = Tchange.getDef(currentMesh.pointSet, model.mean)
        if (settings.optAlignment) { //pose align model
          val curCor = obs.map(t => (t._1, t._2)) ++ (if (ldms.nonEmpty) ldms else IndexedSeq.empty)
          val tf = ModelUtils.alignShapeToLdms(currentMesh, curCor)._2
          model = model.transform(tf)
          posedModel = if (useRealingment) posedModel.transform(tf) else posedModel
        }
        val cursigma2 = getCurrentSigma2(i, settings.maxIterations, settings.sigma2, settings.sigma2*settings.sigma2InitFactor)
        val pmesh = {
          val td = FitScripts.getSurfaceAwareTd(model, obs, cursigma2, perpendicularFactor = settings.perpendicularFactor)
          val tdldm = td ++ (if (ldms.nonEmpty) {
            val mvn = Tchange.getZeroMeanNoise(settings.sigma2)
            ldms.map(t => (t._1,t._2,mvn))
          } else IndexedSeq.empty)
          model.posterior(tdldm).mean
        }
        val next = Tchange.getDef(pmesh.pointSet, model.mean)
        val change = l2.norm2Vector(next.zip(current).map(t => t._1 - t._2))
        if (lookAheadRegression && math.abs(cursigma2-settings.sigma2) > 1e-2) { //perform a regression with these observations but with a quickly converging sigma2
          val pmeshLookAhead = {
            val cs = settings.sigma2 //getCurrentSigma2(i, 10, settings.sigma2, settings.sigma2*settings.sigma2InitFactor,mode=1)
            val td = FitScripts.getSurfaceAwareTd(model, obs, cs)
            model.posterior(td).mean
          }
          if (change < settings.stopTol) (pmesh, true, i + 1, pmeshLookAhead) else (pmesh, false, i + 1, pmeshLookAhead)
        } else if (change < settings.stopTol) (pmesh, true, i + 1, pmesh) else (pmesh, false, i + 1, pmesh)
      }
    })
    it.map(t => (t._4, t._2, t._3))  //t._4 is t._1 if no lookahead is given
  }


  def getCurrentSigma2(it: Int, maxit: Int, finalS2: Double, robustS2: Double, mode: Int=0): Double = {
    val p = math.min(it.toDouble / maxit, 1.0) //maxit for until ends2 reached.
    mode match {
      case 0 =>
        val finalExp = math.log(finalS2)
        val robustExp = math.log(robustS2)
        val exp = p * finalExp + (1.0 - p) * robustExp
        val sigma2 = math.exp(exp)
        sigma2
      case 1 =>
        val sigma2 = p * finalS2 + (1.0 - p) * robustS2
        sigma2
    }
  }

}
