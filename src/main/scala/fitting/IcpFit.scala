package fitting

import norms.L2norm
import scalismo.common.PointId
import scalismo.geometry.{EuclideanVector, EuclideanVector3D, Point, _3D}
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.PointDistributionModel
import scalismo.{ModelRecenter, ModelUtils}
import utility.{MathHelp, Tchange}

object IcpFit {

  case class IcpSettings(sampling: SamplingStrategy, optAlignment: Boolean = true, maxIterations: Int = 25, sigma2: Double = 1.0, sigma2InitFactor: Double = 1000.0, stopTol: Double = 1e-2)

  /**
   * applies multiple non rigid icp optimizations.
   */
  def apply(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings): (TriangleMesh[_3D], Int) = {
    val l2 = L2norm[_3D]()
    def getObs(m: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = settings.sampling.establishCorrespondence(m, target)
    var model = modelOrig
    var current = model.mean.pointSet.pointSequence.map(_=>EuclideanVector3D.zero)

    (0 until settings.maxIterations).foreach(i => {
      val currentMesh = if (i == 0) model.mean else Tchange.meshPlusEv(model.mean, current)
      val obs = getObs(currentMesh)
      val correctedObs = if (settings.optAlignment) {
        val tf = ModelUtils.alignShapeToLdms(currentMesh, obs)._2
        model = model.transform(tf)
        getObs(currentMesh.transform(tf))
      } else obs

      val cursigma2 = {
        val p = i.toDouble / settings.maxIterations
        p*settings.sigma2 + (1.0-p)*settings.sigma2*settings.sigma2InitFactor
      }
      val pmesh = model.posterior(correctedObs, cursigma2).mean
      val next = Tchange.getDef(pmesh.pointSet, model.mean)
      val change = l2.norm2Vector(next.zip(current).map(t=>t._1-t._2))
      if (change < settings.stopTol) return (pmesh, i)
      current = next
    })

    val lastMesh = Tchange.meshPlusEv(model.mean, current)
    val lastObs = getObs(lastMesh)
    val result = model.posterior(lastObs, settings.sigma2).mean
    (result, settings.maxIterations)
  }

  def fitWithOnlineAlignment(modelOrig: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings, withRotation: Boolean = false, mainAxis: EuclideanVector[_3D], alignmentSettings: Option[SamplingStrategy] = Option(TargetSamplingUnique())): (TriangleMesh[_3D], Int) = {
    val l2 = L2norm[_3D]()
    val aset = alignmentSettings.getOrElse(settings.sampling)
    val axis = MathHelp.listOrthogonalVectors(mainAxis)
    def getObs(m: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = settings.sampling.establishCorrespondence(m, target)
    var model = modelOrig
    var posedModel = modelOrig
    var current = model.mean.pointSet.pointSequence.map(_=>EuclideanVector3D.zero)
    (0 until settings.maxIterations).foreach(i => {
      val currentMesh = if (i == 0) model.mean else Tchange.meshPlusEv(model.mean, current)
      val aobs = aset.establishCorrespondence(currentMesh, target)
      model = if (withRotation){ //basis align model
        ModelRecenter.completeAlignPdm(posedModel, aobs.map(_._1), axis)
      } else ModelRecenter.recenter(posedModel, aobs.map(_._1))
      val obs = if (alignmentSettings.isDefined) getObs(currentMesh) else aobs
      val correctedObs = if (settings.optAlignment) { //pose align model
        val tf = ModelUtils.alignShapeToLdms(currentMesh, obs)._2
        model = model.transform(tf)
        posedModel = posedModel.transform(tf)
        getObs(currentMesh.transform(tf))
      } else obs
      val cursigma2 = {
        val p = i.toDouble / settings.maxIterations
        p * settings.sigma2 + (1.0 - p) * settings.sigma2 * settings.sigma2InitFactor
      }
      val pmesh = model.posterior(correctedObs, cursigma2).mean
      val next = Tchange.getDef(pmesh.pointSet, model.mean)
      val change = l2.normVector(next.zip(current).map(t=>t._1-t._2))
      if (change < settings.stopTol) return (pmesh, i)
      current = next
    })

    val lastMesh = Tchange.meshPlusEv(model.mean, current)
    val lastObs = getObs(lastMesh)
    val result = model.posterior(lastObs, settings.sigma2).mean
    (result, settings.maxIterations)
  }
}
