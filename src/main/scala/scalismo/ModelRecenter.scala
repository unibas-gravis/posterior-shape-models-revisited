package scalismo

import breeze.linalg.{Axis, DenseMatrix, norm, sum}
import norms.L2norm
import scalismo.common._
import scalismo.common.interpolation.NearestNeighborInterpolator
import scalismo.geometry.EuclideanVector.Vector_3DVectorizer
import scalismo.geometry.{EuclideanVector, NDSpace, Point, _3D}
import scalismo.mesh.TriangleMesh
import scalismo.numerics.PivotedCholesky
import scalismo.numerics.PivotedCholesky.RelativeTolerance
import scalismo.statisticalmodel.DiscreteLowRankGaussianProcess.Eigenpair
import scalismo.statisticalmodel.{DiscreteLowRankGaussianProcess, PointDistributionModel, StatisticalMeshModel}
import utility.{MathHelp, Tchange}

object ModelRecenter {

    /**
     * recenters the model on the given points by using translation.
     * has undefined behaviour for constant eigenfunctions
     */
    def recenter(ssm: StatisticalMeshModel, centering: IndexedSeq[PointId]): StatisticalMeshModel = {
      //calc mean of centering set
      val centint = centering.map(_.id*3)
      val vectorsx = sum(DenseMatrix.tabulate(centint.length,ssm.rank)((x,y) => ssm.gp.basisMatrix(centint(x),y)),Axis._0)*(1.0/centint.length)
      val vectorsy = sum(DenseMatrix.tabulate(centint.length,ssm.rank)((x,y) => ssm.gp.basisMatrix(centint(x)+1,y)),Axis._0)*(1.0/centint.length)
      val vectorsz = sum(DenseMatrix.tabulate(centint.length,ssm.rank)((x,y) => ssm.gp.basisMatrix(centint(x)+2,y)),Axis._0)*(1.0/centint.length)
      //adjust existing eigenfunctions
      val eigfunc = {
        val base = ssm.gp.basisMatrix.copy
        val vx = vectorsx.t.toDenseVector.toDenseMatrix
        val vy = vectorsy.t.toDenseVector.toDenseMatrix
        val vz = vectorsz.t.toDenseVector.toDenseMatrix
        val translation = DenseMatrix.vertcat(vx,vy,vz)
        (0 until base.rows by 3).foreach(i => base(i until i+3,::):-=translation)
        base
      }
      //calculate eigenfunction norms
      val funcnorms = norm(eigfunc,Axis._0).t //TODO if desired filter for funcs too small -> replace with (0.0, normalized translation as eigfunction or original one or truncate)
      //norm eigenfunctions as inplace operation
      (0 until ssm.rank).map(i => eigfunc(::,i):/=funcnorms(i)) //TODO renorming is not strictly necessary -> non normed eigfuncs to avoid /0
      //change eigenvalue to reflect the normed eigenfunction
      val eigvCorrected = ssm.gp.variance*:*funcnorms*:*funcnorms

      //TODO rediagonalization for more than sampling. gram better
      val dgp = DiscreteLowRankGaussianProcess(ssm.gp._domain, ssm.gp.meanVector, eigvCorrected, eigfunc)
      val apprxeig = PivotedCholesky.computeApproximateEig(dgp.interpolate(NearestNeighborInterpolator()).cov, ssm.gp._domain.pointSet.pointSequence, RelativeTolerance(1e-6))
      StatisticalMeshModel(ssm.referenceMesh,DiscreteLowRankGaussianProcess(ssm.gp._domain, ssm.gp.meanVector, apprxeig._2, apprxeig._1))
    }

  /**
   * recenters the model on the given points by using translation.
   * has undefined behaviour for constant eigenfunctions
   */
  def recenter[D: NDSpace, PointRepr[D] <: DiscreteDomain[D]](pdm: PointDistributionModel[D, PointRepr], centering: IndexedSeq[PointId])(implicit warp: DomainWarp[D, DiscreteDomain], vectorizer: Vectorizer[EuclideanVector[D]]): PointDistributionModel[D, DiscreteDomain] = {
    //  def recenter(pdm: StatisticalMeshModel, centering: IndexedSeq[PointId]): StatisticalMeshModel = {
    val dim = pdm.gp.outputDim
    //calc mean of centering set
    val centint = centering.map(_.id * dim)
    val vectors = (0 until dim).map(i => sum(DenseMatrix.tabulate(centint.length, pdm.rank)((x, y) => pdm.gp.basisMatrix(centint(x) + i, y)), Axis._0) * (1.0 / centint.length))
    //adjust existing eigenfunctions
    val eigfunc = {
      val base = pdm.gp.basisMatrix.copy
      val v = vectors.map(_.t.toDenseMatrix)
      val translation = DenseMatrix.vertcat(v: _*)
      (0 until base.rows by dim).foreach(i => base(i until i + dim, ::) :-= translation)
      base
    }
    //calculate eigenfunction norms
    val funcnorms = norm(eigfunc, Axis._0).t //TODO if desired filter for funcs too small -> replace with (0.0, normalized translation as eigfunction or original one or truncate)
    //norm eigenfunctions as inplace operation
    (0 until pdm.rank).map(i => eigfunc(::, i) :/= funcnorms(i))
    //change eigenvalue to reflect the normed eigenfunction
    val eigvCorrected = pdm.gp.variance *:* funcnorms *:* funcnorms

    //TODO fast for pure sampling. unify with other approaches
    PointDistributionModel(pdm.gp._domain, pdm.gp.meanVector, eigvCorrected, eigfunc)
  }

  /**
   * recenters the model on the given points by using translation.
   * has undefined behaviour for constant eigenfunctions
   */
  def recenter(pdm: PointDistributionModel[_3D, TriangleMesh], centering: IndexedSeq[PointId]): PointDistributionModel[_3D, TriangleMesh] = {
    //  def recenter(pdm: StatisticalMeshModel, centering: IndexedSeq[PointId]): StatisticalMeshModel = {
    val dim = pdm.gp.outputDim
    //calc mean of centering set
    val centint = centering.map(_.id * dim)
    val vectors = (0 until dim).map(i => sum(DenseMatrix.tabulate(centint.length, pdm.rank)((x, y) => pdm.gp.basisMatrix(centint(x) + i, y)), Axis._0) * (1.0 / centint.length))
    //adjust existing eigenfunctions
    val eigfunc = {
      val base = pdm.gp.basisMatrix.copy
      val v = vectors.map(_.t.toDenseMatrix)
      val translation = DenseMatrix.vertcat(v: _*)
      (0 until base.rows by dim).foreach(i => base(i until i + dim, ::) :-= translation)
      base
    }
    //calculate eigenfunction norms
    val funcnorms = norm(eigfunc, Axis._0).t //TODO if desired filter for funcs too small -> replace with (0.0, normalized translation as eigfunction or original one or truncate)
    //norm eigenfunctions as inplace operation
    (0 until pdm.rank).map(i => eigfunc(::, i) :/= funcnorms(i))
    //change eigenvalue to reflect the normed eigenfunction
    val eigvCorrected = pdm.gp.variance *:* funcnorms *:* funcnorms

    //TODO rediagonalization for more than sampling. gram better
    val res = PointDistributionModel(pdm.gp._domain, pdm.gp.meanVector, eigvCorrected, eigfunc)
    val apprxeig = PivotedCholesky.computeApproximateEig(res.gp.interpolate(NearestNeighborInterpolator()).cov, res.gp.domain.pointSet.pointSequence, RelativeTolerance(0.0001))
    PointDistributionModel(pdm.gp._domain, pdm.gp.meanVector, apprxeig._2, apprxeig._1)
  }

  /**
   * recenters the model on the given points by using rotation only. rotates around mean of centering domain points.
   * requires at least three points in the domain. assumes translation alignment to points
   */
  def rerotate(ssm: StatisticalMeshModel, pids: IndexedSeq[PointId], axis: IndexedSeq[EuclideanVector[_3D]], rpoint: Option[Point[_3D]] = None): StatisticalMeshModel = {
    require(pids.length >= 3, "needs >=3 points")
    val l2 = L2norm[_3D]()
    //calc rotation eigenfunctions
    val ids = pids.map(_.id)
    val rp = rpoint.getOrElse(Tchange.getMean(ssm.mean, Option(pids)))
    val rotfuncs = axis.map(ax => {
      val rf = getRf(ssm, ax, rp)
      val norm = l2.normVector(ids.map(rf.apply))
      rf.map(_/norm) //norm the rotations to have a norm of 1 on the observed pids to not changed the later dot product
    })

    val M = DenseMatrix.apply(rotfuncs.map(rf => {
      Tchange.seqToDenseVec[_3D](rf)
    }):_*).t
    val MX = DenseMatrix.apply(rotfuncs.map(rf => {
      Tchange.seqToDenseVec[_3D](ids.map(rf))
    }):_*).t
    val Minv = breeze.linalg.pinv(MX)
    val correctedBasis = ssm.gp.klBasis.map(kl => {
      val (klv, klvX) = {
        (Tchange.seqToDenseVec[_3D](kl.eigenfunction.data), Tchange.seqToDenseVec[_3D](ids.map(kl.eigenfunction.data)))
      }
      val c = Minv*klvX
      val adjdv = klv - M*c
      val adjf = Tchange.dvToSeqEv3D(adjdv)
      val (normF, norm) = MathHelp.normalize(adjf)
      Eigenpair(kl.eigenvalue*(norm*norm), DiscreteField(ssm.gp.domain, normF))
    })


    val dgp = DiscreteLowRankGaussianProcess(ssm.gp.mean,correctedBasis)
    //TODO rediagonalization for more than sampling. gram better
    val apprxeig = PivotedCholesky.computeApproximateEig(dgp.interpolate(NearestNeighborInterpolator()).cov, ssm.gp._domain.pointSet.pointSequence, RelativeTolerance(1e-6))
    StatisticalMeshModel(ssm.referenceMesh,DiscreteLowRankGaussianProcess(ssm.gp._domain, ssm.gp.meanVector, apprxeig._2, apprxeig._1))
  }

  def completeAlign(ssm: StatisticalMeshModel, ids:IndexedSeq[PointId], axis:IndexedSeq[EuclideanVector[_3D]], rp: Option[Point[_3D]]=None): StatisticalMeshModel = {
    val recentered = recenter(ssm,ids)
    rerotate(recentered,ids,axis,rp)
  }

  def completeAlignPdm(pdm: PointDistributionModel[_3D, TriangleMesh], ids:IndexedSeq[PointId], axis:IndexedSeq[EuclideanVector[_3D]], rp: Option[Point[_3D]]=None): PointDistributionModel[_3D, TriangleMesh] = {
    val ssm = ModelUtils.pdmToSsm(pdm)
    val recentered = recenter(ssm,ids)
    val ssmres = rerotate(recentered,ids,axis,rp)
    ModelUtils.ssmToPdm(ssmres)
  }

  def getRf(model: StatisticalMeshModel, axis: EuclideanVector[_3D], rotPoint: Point[_3D]): IndexedSeq[EuclideanVector[_3D]] = {
    val l2norm = L2norm[_3D]()
    val eigFuncUnnormalized: IndexedSeq[EuclideanVector[_3D]] = model.mean.pointSet.points.map(p => {
      val vectFromRotAxis = MathHelp.projectToPlane(p - rotPoint, axis)
      axis.crossproduct(vectFromRotAxis)
    }).toIndexedSeq
    val norm = l2norm.normVector(eigFuncUnnormalized)
    val eigFunc = eigFuncUnnormalized.map(_ * (1.0 / norm))
    eigFunc
  }
}
