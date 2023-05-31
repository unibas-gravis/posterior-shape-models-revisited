package scalismo

import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common._
import scalismo.common.interpolation.{FieldInterpolator, NearestNeighborInterpolator3D}
import scalismo.geometry.{EuclideanVector, EuclideanVector3D, Point, _3D}
import scalismo.kernels.MatrixValuedPDKernel
import scalismo.mesh.TriangleMesh
import scalismo.numerics.PivotedCholesky
import scalismo.numerics.PivotedCholesky.RelativeTolerance
import scalismo.registration.LandmarkRegistration
import scalismo.statisticalmodel.DiscreteLowRankGaussianProcess.{Eigenpair, KLBasis}
import scalismo.statisticalmodel._
import scalismo.statisticalmodel.dataset.DataCollection
import scalismo.transformations.{RigidTransformation, Translation}
import scalismo.utils.Random
import utility.Tchange.getMean
import utility.{MathHelp, Tchange}

object ModelUtils {

  /**
   * does full alignment. therefore needs at least 3 points
   */
  def pcaModel(ref: TriangleMesh[_3D], dataToUse: IndexedSeq[TriangleMesh[_3D]]): StatisticalMeshModel = {
    val model = StatisticalMeshModel.createUsingPCA(
      DataCollection.fromTriangleMesh3DSequence(ref, dataToUse),
      PivotedCholesky.RelativeTolerance(0)
    ).get
    model.truncate(model.gp.variance.toArray.count(_ > 1E-10))
  }

  def pcaAlignIdModel(ref: TriangleMesh[_3D], dataToUse: IndexedSeq[TriangleMesh[_3D]], alignIds: IndexedSeq[PointId], rotation: Boolean = true): StatisticalMeshModel = {
    val rp = getMean(ref,Option(alignIds))
    val aligned = dataToUse.map(mesh => alignShape(mesh,ref,Option(alignIds),rotation,Option(rp))._1)
    val dgp = getDgpSafe(ref, aligned)

    StatisticalMeshModel.apply(ref, dgp)
  }

  /**
   * building model by aligning to subset of pids
   */
  def pcaModelGpa(shapes: IndexedSeq[TriangleMesh[_3D]], alignIds: IndexedSeq[PointId]): StatisticalMeshModel = {
    val aligned = alignShapesGpa(shapes, Option(alignIds))
    pcaModel(aligned._2, aligned._1.map(_._1))
  }

  /**
   * building model by aligning to subset of pids with noisy observations.
   */
  def pcaModelGpaTildeLdm(shapes: IndexedSeq[TriangleMesh[_3D]], alignIds: IndexedSeq[PointId], tildeSigma: Double)(implicit rnd: Random): StatisticalMeshModel = {
    val data = shapes.map(shape => (shape, alignIds.map(pid => {
      val rv = EuclideanVector3D(rnd.scalaRandom.nextGaussian(),rnd.scalaRandom.nextGaussian(),rnd.scalaRandom.nextGaussian())
      shape.operations.closestPointOnSurface(shape.pointSet.point(pid) + rv * tildeSigma).point
    })))
    val aligned = alignShapesGpaLdms(data)
    pcaModel(aligned._2, aligned._1.map(_._1))
  }

  /**
   * same as pcaAlignIdModel but uses free points for the alignment
   */
  def pcaAlignPointModel(ref: (TriangleMesh[_3D], IndexedSeq[Point[_3D]]), dataToUse: IndexedSeq[(TriangleMesh[_3D], IndexedSeq[Point[_3D]])], rotation: Boolean = true): StatisticalMeshModel = {
    require(dataToUse.forall(_._2.length == ref._2.length), "all meshes need the same number of landmark observations")
    val aT = dataToUse.map(t => alignLdmsToLdms(t._2.zip(ref._2), rotation)._2)
    val dgp = getDgpSafe(ref._1, dataToUse.zip(aT).map(t => t._1._1.transform(t._2)))

    StatisticalMeshModel.apply(ref._1, dgp)
  }

  def alignShape(move: TriangleMesh[_3D], to: TriangleMesh[_3D], pidsOpt: Option[IndexedSeq[PointId]] = None, rotation: Boolean = true, rcenter: Option[Point[_3D]] = None): (TriangleMesh[_3D], RigidTransformation[_3D]) = {
    val pids = pidsOpt.getOrElse(move.pointSet.pointIds.toIndexedSeq)
    require(!rotation || pids.length>=3, "for rotational alignment, at least 3 ids are required")
    require(pids.nonEmpty, "for translational alignment, at least one id is required")

    val ldms = pids.map(pid => (move.pointSet.point(pid), to.pointSet.point(pid)))
    val t = if (rotation){
      val rp = rcenter.getOrElse(getMean(move,Option(pids)))
      LandmarkRegistration.rigid3DLandmarkRegistration(ldms,rp)
    } else {
      val movemean = getMean(move,Option(pids))
      val tomean = getMean(to,Option(pids))
      Translation[_3D](tomean-movemean)
    }

    (move.transform(t),t)
  }

  /**
   * aligns shapes to mean shape with only translation
   */
  def alignShapesT(shapes: IndexedSeq[TriangleMesh[_3D]], pidsOpt: Option[IndexedSeq[PointId]]): (IndexedSeq[(TriangleMesh[_3D], RigidTransformation[_3D])], TriangleMesh[_3D]) = {
    val mean = Tchange.getMean(shapes)
    (shapes.map(s => alignShape(s, mean, pidsOpt, rotation = false)), mean)
  }

  /**
   * aligns shapes to mean shape with rigid alignment for a number of iterations. also returns transformation from initial
   * shape to new shape
   */
  def alignShapesGpa(shapes: IndexedSeq[TriangleMesh[_3D]], pidsOpt: Option[IndexedSeq[PointId]] = None, iterations: Int = 4): (IndexedSeq[(TriangleMesh[_3D], RigidTransformation[_3D])], TriangleMesh[_3D]) = {
    val resShapes = (1 to iterations).foldLeft(shapes)((shapes,_) => {
      val mean = Tchange.getMean(shapes)
      val meanrp = Tchange.getMean(mean, pidsOpt)
      shapes.map(s => { //rigidTs are not composed as TbeforeR is difficult to compose efficiently due to built in derivatives etc
        alignShape(s, mean, pidsOpt, rotation = true, Option(meanrp))._1 //but it is possible if also chaining rotationcenters in the operations -> addto Tchange
      })
    })
    val res = shapes.zip(resShapes).map(t => alignShape(t._1, t._2, rotation = true)) //to recover final transformation as we're not composing
    (res, Tchange.getMean(res.map(_._1)))
  }

  /**
   * same as normal Gpa but uses free ldms
   */
  def alignShapesGpaLdms(shapes: IndexedSeq[(TriangleMesh[_3D], IndexedSeq[Point[_3D]])], iterations: Int = 4): (IndexedSeq[(TriangleMesh[_3D], IndexedSeq[Point[_3D]], RigidTransformation[_3D])], TriangleMesh[_3D]) = {
    val resShapes = (1 to iterations).foldLeft(shapes.map(_._2))((shapes,_) => {
      val mean = Tchange.getMeanPs(shapes)
      shapes.map(ps => {
        alignLdmsToLdms(ps.zip(mean), rotation = true)._1
      })
    })
    val resPs = shapes.map(_._2).zip(resShapes).map(t => alignLdmsToLdms(t._1.zip(t._2), rotation = true))
    val res = resPs.zip(shapes).map(t => (t._2._1.transform(t._1._2), t._1._1, t._1._2))
    (res, Tchange.getMean(res.map(_._1)))
  }

  def alignShapeToLdms(move: TriangleMesh[_3D], obs: IndexedSeq[(PointId,Point[_3D])], rotation: Boolean = true, rcenter: Option[Point[_3D]] = None): (TriangleMesh[_3D], RigidTransformation[_3D]) = {
    require(!rotation || obs.length>=3, "for rotational alignment, at least 3 ids are required")
    require(obs.nonEmpty, "for translational alignment, at least one id is required")
    val pids = Option(obs.map(_._1))

    val ldms = obs.map(t => (move.pointSet.point(t._1), t._2))
    val t = if (rotation){
      val rp = rcenter.getOrElse(getMean(move,pids))
      LandmarkRegistration.rigid3DLandmarkRegistration(ldms, rp)
    } else {
      val movemean = getMean(move,pids)
      val tomean = Tchange.getMean(obs.map(_._2))
      Translation[_3D](tomean-movemean)
    }

    (move.transform(t),t)
  }

  def alignLdmsToLdms(obs: IndexedSeq[(Point[_3D], Point[_3D])], rotation: Boolean = true, rcenter: Option[Point[_3D]] = None): (IndexedSeq[Point[_3D]], RigidTransformation[_3D]) = {
    require(!rotation || obs.length>=3, "for rotational alignment, at least 3 obs are required")
    require(obs.nonEmpty, "for translational alignment, at least one obs is required")


    val t = if (rotation){
      val rp = rcenter.getOrElse(getMean(obs.map(_._1)))
      LandmarkRegistration.rigid3DLandmarkRegistration(obs, rp)
    } else {
      val movemean = getMean(obs.map(_._1))
      val tomean = Tchange.getMean(obs.map(_._2))
      Translation[_3D](tomean-movemean)
    }

    (obs.map(pp => t.f(pp._1)), t)
  }


  case class ModelData(ref: TriangleMesh[_3D], mean: TriangleMesh[_3D], meanf: DenseVector[Double], basis: DenseMatrix[Double], eigv: DenseVector[Double])
  /**
   * returns mean and variance for !pre-aligned! meshes that are in correspondence.
   * the covariance is denoted in basismatrix and variance vector form.
   * this function should only be used if the number of meshes is sufficiently small.
   */
  def getModelDataLowLevel(ref: TriangleMesh[_3D], meshes: IndexedSeq[TriangleMesh[_3D]]): ModelData = {
    val mean = Tchange.getMean(meshes)
    val meanf = Tchange.seqToDenseVec(Tchange.getDef(mean.pointSet,ref))
    val defs = meshes.map(mesh => Tchange.getDef(mesh.pointSet, mean))
    val (vecs,vecnorms) = {
      val vecsUnnorm = defs.map(vs => Tchange.seqToDenseVec(vs))
      val vecnorms = vecsUnnorm.map(v => breeze.linalg.norm(v))
      (vecsUnnorm.zip(vecnorms).map(t => (1.0/t._2)*t._1), vecnorms)
    }
    val dscale = breeze.linalg.diag(DenseVector(vecnorms.map(math.sqrt).toArray))
    val basis = DenseVector.horzcat(vecs:_*)
    val c = basis.t * basis
    val scaledc = dscale*c*dscale
    val dec = breeze.linalg.eigSym(scaledc)
    val relrank = dec.eigenvalues.data.count(d => d > 1E-10)
    val ceigv = DenseVector(dec.eigenvalues.data.reverse.take(relrank))
    val cvecs = dec.eigenvectors(::,dec.eigenvectors.cols-1 to dec.eigenvectors.cols-relrank by -1)

    val newBasisUn = basis*cvecs
    val basisNorms = breeze.linalg.norm(newBasisUn,breeze.linalg.Axis._0).t
    val newBasis = newBasisUn*MathHelp.spdiag(basisNorms.map(d => 1.0/d))
    val newEigv = ceigv.*:*(basisNorms.*:*(basisNorms))
    throw new NotImplementedError("implementation error")
    ModelData(ref, mean, meanf, newBasis, newEigv)
  }

  def getDgpSafe(ref: TriangleMesh[_3D], meshes: IndexedSeq[TriangleMesh[_3D]], relTol: Double = 1e-6): DiscreteLowRankGaussianProcess[_3D, TriangleMesh, EuclideanVector[_3D]] = {
    val meanMesh = Tchange.getMean(meshes)
    val meanfunction = (p:Point[_3D]) => {
      val clp = ref.pointSet.findClosestPoint(p)
      meanMesh.pointSet.point(clp.id) - clp.point
    }
    val meanf: Field[_3D, EuclideanVector[_3D]] = Field.apply(EuclideanSpace3D, meanfunction)
    val kernel: MatrixValuedPDKernel[_3D] = new MatrixValuedPDKernel[_3D]() {
//      val km:DenseMatrix[Double] = {
//        val kernelMat = meshes.foldLeft(DenseMatrix.zeros[Double](ref.pointSet.numberOfPoints*3,ref.pointSet.numberOfPoints*3))((mat, mesh) => {
//          val v = Tchange.seqToDenseVec(Tchange.getDef(mesh.pointSet,meanMesh))
//          mat += v*v.t
//        })
//        kernelMat /= meshes.length.toDouble
//      }
      val kdata = meshes.map(_.pointSet.pointSequence.zip(meanMesh.pointSet.pointSequence).map(t => t._1-t._2).map(_.toBreezeVector)).transpose
      override protected def k(x: Point[_3D], y: Point[_3D]): DenseMatrix[Double] = {
        val xid = ref.pointSet.findClosestPoint(x).id.id
        val yid = ref.pointSet.findClosestPoint(y).id.id
//        km(xid*3 until xid*3+3, yid*3 until yid*3+3)
        val xps = kdata(xid)
        val yps = kdata(yid)
        val km = xps.zip(yps).map(t => t._1*t._2.t).reduce(_+_)
        km / xps.length.toDouble
      }
      override def outputDim: Int = 3
      override def domain: Domain[_3D] = meanf.domain
    }
    val gp = GaussianProcess.apply(meanf, kernel)
    implicit val vectorizer = EuclideanVector.Vector_3DVectorizer
    val interpolator: FieldInterpolator[_3D, TriangleMesh, EuclideanVector[_3D]] = NearestNeighborInterpolator3D()
    val lgp = LowRankGaussianProcess.approximateGPCholesky[_3D, TriangleMesh, EuclideanVector[_3D]](ref,gp,relTol, interpolator)
    val dgp = lgp.discretize(ref)

    dgp
  }

  def getModelData(ref: TriangleMesh[_3D], meshes: IndexedSeq[TriangleMesh[_3D]]): (DiscreteField[_3D, TriangleMesh, EuclideanVector[_3D]], KLBasis[_3D, TriangleMesh, EuclideanVector[_3D]]) = {
    val data = getModelDataLowLevel(ref,meshes)
    val df = DiscreteField.apply(ref, Tchange.dvToSeqEv3D(data.meanf))
    (df, (0 until data.basis.cols).map(i => Eigenpair(data.eigv(i), DiscreteField(ref, Tchange.dvToSeqEv3D(data.basis(::,i))))))
  }

  /**
   * creates from a model a new model by sampling a few shapes. uses correspondence for the new model
   */
  def buildSampleModel(model: PointDistributionModel[_3D, TriangleMesh], numSamples: Int, relTol: Double = 0.0)(implicit rnd: Random): PointDistributionModel[_3D, TriangleMesh] = {
    val samples = (1 to numSamples).map(_ => model.sample()).map(shape => DiscreteField3D.apply(model.reference,Tchange.getDef(shape.pointSet,model.reference)).interpolate(NearestNeighborInterpolator3D()))
    PointDistributionModel.createUsingPCA(model.reference,samples,RelativeTolerance(relTol))
  }

  //both transforms are cheap
  def pdmToSsm(model: PointDistributionModel[_3D, TriangleMesh]): StatisticalMeshModel = StatisticalMeshModel.apply(model.reference,model.gp)
  def ssmToPdm(model: StatisticalMeshModel): PointDistributionModel[_3D, TriangleMesh] = PointDistributionModel.apply(model.gp)

}
