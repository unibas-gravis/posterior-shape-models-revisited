package scalismo

import breeze.linalg.{*, Axis, DenseMatrix, DenseVector, norm, sum}
import norms.L2norm
import scalismo.common._
import scalismo.geometry.EuclideanVector.Vector_3DVectorizer
import scalismo.geometry.{EuclideanVector, NDSpace, Point, _3D}
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.DiscreteLowRankGaussianProcess.Eigenpair
import scalismo.statisticalmodel.{DiscreteLowRankGaussianProcess, PointDistributionModel, StatisticalMeshModel}
import utility.{MathHelp, Tchange}

object ModelRecenter {

  /**
   * recenters the model on the given points by using translation.
   * has undefined behaviour for constant eigenfunctions
   */
  def recenterSsm(ssm: StatisticalMeshModel, centering: IndexedSeq[PointId], rediag: Boolean=true): StatisticalMeshModel = {
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

    val dgp = DiscreteLowRankGaussianProcess(ssm.gp._domain, ssm.gp.meanVector, eigvCorrected, eigfunc)
    if (rediag){
      val apprxeig = rediagonalizeGram(eigfunc, eigvCorrected) //about 50 times faster than complete rebuilding for femurs
      StatisticalMeshModel(ssm.referenceMesh,DiscreteLowRankGaussianProcess(ssm.gp._domain, ssm.gp.meanVector, apprxeig._2, apprxeig._1))
    } else  StatisticalMeshModel(ssm.referenceMesh,dgp)
  }

  /**
   * recenters the model on the given points by using translation.
   * has undefined behaviour for constant eigenfunctions
   */
  def recenterPdmGen[D: NDSpace, PointRepr[D] <: DiscreteDomain[D]](pdm: PointDistributionModel[D, PointRepr], centering: IndexedSeq[PointId], rediag: Boolean=true)(implicit warp: DomainWarp[D, DiscreteDomain], vectorizer: Vectorizer[EuclideanVector[D]]): PointDistributionModel[D, DiscreteDomain] = {
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

    val res: PointDistributionModel[D, DiscreteDomain] = PointDistributionModel(pdm.gp._domain, pdm.gp.meanVector, eigvCorrected, eigfunc)
    if (rediag) {
      val apprxeig = rediagonalizeGram(eigfunc, eigvCorrected)
      PointDistributionModel(pdm.gp._domain, pdm.gp.meanVector, apprxeig._2, apprxeig._1)
    } else res
  }

  /**
   * recenters the model on the given points by using translation.
   * has undefined behaviour for constant eigenfunctions
   */
  def recenterPdm(pdm: PointDistributionModel[_3D, TriangleMesh], centering: IndexedSeq[PointId], rediag: Boolean=true): PointDistributionModel[_3D, TriangleMesh] = {
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

    val res = PointDistributionModel(pdm.gp._domain, pdm.gp.meanVector, eigvCorrected, eigfunc)
    if (rediag) {
      val apprxeig = rediagonalizeGram(eigfunc, eigvCorrected)
      PointDistributionModel(pdm.gp._domain, pdm.gp.meanVector, apprxeig._2, apprxeig._1)
    } else res
  }

  /**
   * recenters the model on the given points by using rotation only. rotates around mean of centering domain points.
   * requires at least three points in the domain. assumes translation alignment to points
   */
  def rerotate(ssm: StatisticalMeshModel, pids: IndexedSeq[PointId], axis: IndexedSeq[EuclideanVector[_3D]], rpoint: Option[Point[_3D]] = None, rediag: Boolean=true): StatisticalMeshModel = {
    require(pids.length >= 3, "provided not enough points")
    val l2 = L2norm[_3D]()
    //calc rotation eigenfunctions
    val ids = pids.map(_.id)
    val rp = rpoint.getOrElse(Tchange.getMean(ssm.mean, Option(pids)))
    val rotfuncs = axis.map(ax => {
      val rf = RegressionHelp.getRf(ssm, ax, rp)
      val norm = l2.normVector(ids.map(rf.apply))
      rf.map(_/norm) //norm the rotations to have a norm of 1 on the observed pids to not changed the later dot product
    })

    val M = DenseMatrix.apply(rotfuncs.map(rf => {
      Tchange.seqToDenseVec[_3D](rf)
    }):_*).t
    val MX = DenseMatrix.apply(rotfuncs.map(rf => {
      Tchange.seqToDenseVec[_3D](ids.map(rf))
    }):_*).t
    val Minv = MathHelp.pseudoInverse(MX)._1  //same as breeze.linalg.pinv(MX) but numerically more efficient
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
    if (rediag) {
      val apprxeig = rediagonalizeGram(dgp.basisMatrix, dgp.variance)
      StatisticalMeshModel(ssm.referenceMesh,DiscreteLowRankGaussianProcess(ssm.gp._domain, ssm.gp.meanVector, apprxeig._2, apprxeig._1))
    } else StatisticalMeshModel(ssm.referenceMesh,dgp)
  }

  def completeAlign(ssm: StatisticalMeshModel, ids:IndexedSeq[PointId], axis:IndexedSeq[EuclideanVector[_3D]], rp: Option[Point[_3D]]=None, rediag: Boolean=true): StatisticalMeshModel = {
    val recentered = recenterSsm(ssm,ids,rediag)
    rerotate(recentered,ids,axis,rp,rediag)
  }

  def completeAlignPdm(pdm: PointDistributionModel[_3D, TriangleMesh], ids:IndexedSeq[PointId], axis:IndexedSeq[EuclideanVector[_3D]], rp: Option[Point[_3D]]=None, rediag: Boolean=true): PointDistributionModel[_3D, TriangleMesh] = {
    val ssm = ModelUtils.pdmToSsm(pdm)
    val recentered = recenterSsm(ssm,ids,rediag)
    val ssmres = rerotate(recentered,ids,axis,rp,rediag)
    ModelUtils.ssmToPdm(ssmres)
  }

  /**
   * assumes the basis nxr is non orthogonal, s is the variance (squared [eigen]scalars) of that basis.
   * returns a orthonormal basis with the adjusted variance. assumes no zero variance values
   */
  def rediagonalizeGram(basis: DenseMatrix[Double], s: DenseVector[Double]): (DenseMatrix[Double], DenseVector[Double]) = {
    val l = basis(*, ::) * breeze.numerics.sqrt(s)  //basis * breeze.linalg.diag(breeze.numerics.sqrt(s))
    val gram = l.t * l
    val svd = breeze.linalg.svd(gram)
    val newbasis = l * (svd.U(*, ::) * (1.0/breeze.numerics.sqrt(svd.S)))  //l * svd.U * breeze.linalg.diag(1.0/breeze.numerics.sqrt(svd.S))
    (newbasis, svd.S)
  }

}
