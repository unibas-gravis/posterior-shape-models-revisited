package utility

import breeze.linalg.svd.DenseSVD
import breeze.linalg.{CSCMatrix, DenseMatrix, DenseVector, det, diag, eig, svd}
import norms.L2norm
import scalismo.geometry._
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.PointDistributionModel

object MathHelp {

  def testOrtogonality(f1: IndexedSeq[EuclideanVector[_3D]], f2: IndexedSeq[EuclideanVector[_3D]]): Double = {
    f1.zip(f2).map(t => t._1.dot(t._2)).sum
  }

  def denseToSparse(dm: DenseMatrix[Double], tolerance: Double = 1E-10): CSCMatrix[Double] = {
    CSCMatrix.tabulate(dm.rows, dm.cols)((x,y) => if (dm(x,y) > tolerance) dm(x,y) else 0.0)
  }

  /**
   * returns the marginal variance of a given vector for the model. projection and scaling with the kl basis
   */
  def getMarginalVariance(model: PointDistributionModel[_3D, TriangleMesh], vec: DenseVector[Double]): Double = {
    throw new NotImplementedError()
    //TODO check if it makes sense
    println("test " + model.gp.klBasis.map(kl => math.pow(Tchange.seqToDenseVec(kl.eigenfunction.data).dot(vec),2)).sum)
    model.gp.klBasis.map(kl => kl.eigenvalue * math.pow(Tchange.seqToDenseVec(kl.eigenfunction.data).dot(vec),2)).sum
  }

  /**
   * pseudoinverse with svd and eigenvalue cutoff. rather expensive
   */
  def pseudoInverse(mat: DenseMatrix[Double], svdecomp: Option[DenseSVD]=None, cutoff: Double = 1E-10): (DenseMatrix[Double], DenseSVD) = {
    val sv = svdecomp.getOrElse(svd(mat))
    val isv = invertSVD(sv)
    (reconstruct(isv), isv)
  }

  def invertSVD(sv: DenseSVD, cutoff: Double = 1E-10): DenseSVD = {
    sv.copy(
      leftVectors = sv.rightVectors.t,
      singularValues = sv.singularValues.map(d => if(d>cutoff)1.0/d else 0.0),
      rightVectors = sv.leftVectors.t
    )
  }

  def invertWithSVD(mat: DenseMatrix[Double], sv: DenseSVD, cutoff: Double = 1E-10): (DenseMatrix[Double],DenseSVD) = {
    pseudoInverse(mat,Option(sv))
  }

  def reconstruct(sv: DenseSVD): DenseMatrix[Double] = {
    sv.leftVectors * diag(sv.singularValues) * sv.rightVectors
  }

  /**
   * replaces the cut eigenvalues with a 1 in the determinante calculation
   */
  def determinante(mat: DenseMatrix[Double], svdecomp: Option[DenseSVD]=None, cutoff: Double = 1E-10): Double = {
    if (svdecomp.isDefined){
      svdecomp.get.singularValues.data.foldLeft(1.0)((p,d) => if(d>cutoff) p*d else p)
    } else {
      val e = eig(mat)
      e.eigenvalues.data.foldLeft(1.0)((p,d) => if(d>cutoff) p*d else p)
    }
  }

  def logdeterminante(mat: DenseMatrix[Double], svdecomp: Option[DenseSVD]=None, cutoff: Double = 1E-10): Double = {
    if (svdecomp.isDefined){
      svdecomp.get.singularValues.data.foldLeft(0.0)((p,d) => if(d>cutoff) p+math.log(d) else p)
    } else {
      val e = eig(mat)
      e.eigenvalues.data.foldLeft(0.0)((p,d) => if(d>cutoff) p+math.log(d) else p)
    }
  }

  /**
   * calculates the kl distance. rather expensive due to the needed inversion. it is not a real distance and not symmetric.
   */
  def klMetric(mat1: DenseMatrix[Double], mat2: DenseMatrix[Double], svde1: Option[DenseSVD]=None, svde2: Option[DenseSVD]=None): Double = {
    val sv1 = svde1.getOrElse(svd(mat1))
    val sv2 = svde2.getOrElse(svd(mat2))
    val (mat1inv,isvd1) = pseudoInverse(mat1, Option(sv1))
    val invprod = mat1inv*mat2
    val tr = breeze.linalg.trace(invprod)
    val d = logdeterminante(mat1,Option(sv1))-logdeterminante(mat2,Option(sv2))
    0.5*(tr-1+d)
  }

//  def logklMetric(mat1: DenseMatrix[Double], mat2: DenseMatrix[Double], svde1: Option[DenseSVD]=None, svde2: Option[DenseSVD]=None): Double = {
//    val sv1 = svde1.getOrElse(svd(mat1))
//    val sv2 = svde1.getOrElse(svd(mat2))
//    val mat1inv = pseudoInverse(mat1, Option(sv1))
//    val logtr = diag(mat1inv*mat2).toArray.map(math.log).sum
//    //    val d = math.log(determinante(mat2,Option(sv2))/determinante(mat1,Option(sv1)))
//    val d = math.log(logdeterminante(mat2,Option(sv2))-logdeterminante(mat1,Option(sv1)))
//    res
//  }

  /**
   * simplified by assuming the same mean.
   */
  def bhattacharyyaDistance(mat1: DenseMatrix[Double], mat2: DenseMatrix[Double], svde1: Option[DenseSVD]=None, svde2: Option[DenseSVD]=None): Double = {
    val mat = 0.5*(mat1+mat2)
//    val d1 = determinante(mat1,svde1)
//    val d2 = determinante(mat2,svde2)
//    val d = determinante(mat)
//    0.5*math.log(d / math.sqrt(d1*d2))
    val d1 = logdeterminante(mat1,svde1)
    val d2 = logdeterminante(mat2,svde2)
    val d = logdeterminante(mat)
    0.5 * (d - 0.5*(d1+d2))
  }

  /**
   * expensive as the mahalanobis distance is calculated using a matrix inverse
   */
  def bhattacharyyaDistance(vec1: DenseVector[Double], mat1: DenseMatrix[Double], vec2: DenseVector[Double], mat2: DenseMatrix[Double]): Double = {
    val mat = 0.5*(mat1+mat2)
    val d1 = 0.125 * ((vec1-vec2).t*pseudoInverse(mat)._1*(vec1-vec2))
    val d2 = 0.5*math.log(det(mat) / math.sqrt(det(mat1)*det(mat2)))
    d1+d2
  }

  def correlationMatrixDistance(mat1: DenseMatrix[Double], mat2: DenseMatrix[Double]): Double = {
    val tr = breeze.linalg.trace(mat1*mat2)
    val d = frobeniusNorm(mat1)*frobeniusNorm(mat2)
    1.0-tr/d
  }

  def frobeniusNorm(mat: DenseMatrix[Double]): Double = {
//    math.sqrt(mat.data.foldLeft(0.0)((s,d) => s+d*d))
    math.sqrt(breeze.linalg.sum(mat*:*mat))
  }

  def spdiag(v: DenseVector[Double]): CSCMatrix[Double] = {
    val cm = CSCMatrix.zeros[Double](v.length,v.length)
    cm.reserve(v.length)
    (0 until v.length).foreach(i => cm.update(i,i,v(i)))
    cm
  }

  /**
   * projects to the plane defined by origin and normal n. n is assumed to be normalized
   */
  def projectToPlane(p: Point[_3D], po:Point[_3D], n: EuclideanVector[_3D]): Point[_3D] = po+projectToPlane(p-po, n)
  def projectToPlane(v: EuclideanVector[_3D], n: EuclideanVector[_3D]): EuclideanVector[_3D] = v - n*n.dot(v)//.normalize)

  case class AngleInfo(angle:Double, axis: EuclideanVector[_3D])
  /**
   * returns in angle and rotation axis if 3 dimensional. axis is normalized
   */
  def getAngleInfo(v1: EuclideanVector[_3D], v2: EuclideanVector[_3D]): AngleInfo = {
    val a = v1.normalize
    val b = v2.normalize
    AngleInfo(math.acos(a.dot(b)), a.crossproduct(b))
  }
  def getInAngle[D](v1: EuclideanVector[D], v2: EuclideanVector[D]): Double = {
    val dot = v1.normalize.dot(v2.normalize)
    if (dot < -1.0) math.Pi else if (dot > 1.0) 0.0 else math.acos(dot)
  }

  /**
   * returns two random orthogonal vectors to v1. all three are pairwise orthogonal
   */
  def getOrthogonalVectors(v: EuclideanVector[_3D]): (EuclideanVector[_3D],EuclideanVector[_3D]) = {
    val vb = v.toBreezeVector
    val m = svd(vb*vb.t) //TODO implement a more targeted approach
    val v1 = m.leftVectors(::,1).toDenseMatrix.toDenseVector
    val v2 = m.leftVectors(::,2).toDenseMatrix.toDenseVector
    (EuclideanVector.fromBreezeVector[_3D](v1), EuclideanVector.fromBreezeVector[_3D](v2))
  }

  def listOrthogonalVectors(v: EuclideanVector[_3D]): IndexedSeq[EuclideanVector[_3D]] = {
    val vs = getOrthogonalVectors(v)
    IndexedSeq(v,vs._1,vs._2)
  }

  /**
   * performs the edge check for all three edges.
   * also returns the int for which edge the point originates (0=12, 1=13, 2=23)
   */
  def intersectThreeEdges(t1:Point[_2D], t2:Point[_2D], t3:Point[_2D]): (Double,Int) = {
    Seq(checkEdge(t1,t2),checkEdge(t1,t3),checkEdge(t2,t3)).zipWithIndex.maxBy(_._1)
  }
  /**
   * performs the edge check for two edges 12,13 with the reason that 23 should not be considered to save time.
   */
  def intersectTwoEdges(t1:Point[_2D], t2:Point[_2D], t3:Point[_2D]): (Double,Int) = {
    Seq(checkEdge(t1,t2),checkEdge(t1,t3)).zipWithIndex.maxBy(_._1)
  }
  /**
   * returns the max x value of x-axis intersection of the given edge. returns -Inf if no intersections
   * exist.
   */
  def checkEdge(p1:Point[_2D], p2:Point[_2D]): Double = (p1.x,p1.y,p2.x,p2.y) match { //for speed changed to case class
    case (x1 ,0.0,x2 ,0.0) => math.max(x1,x2)
    case (_  ,_  ,x2 ,0.0) => x2
    case (x1 ,0.0,_  ,_  ) => x1
    case (x1 ,y1 ,x2 ,y2 ) => if (y1*y2<0.0){
      val k = -y1/(y2-y1)
      x1+k*(x2-x1)
    } else Double.NegativeInfinity
  }

  /**
   * t1 and t2 are on the xy plane and correspond to transformed triangle p1,p2,p3.
   * returns t3~p3 in xy plane. p3 is chosen to be on x+ side while the t points should be chosen such that the normal
   * under the same transformation is (0,0,1).
   */
  def projectLastPoint2D(t1:Point[_2D], t2:Point[_2D], p1:Point[_3D], p2:Point[_3D], p3:Point[_3D]): Point[_2D] = {
    val d12 = (p2-p1).norm
    val d13 = (p3-p1).norm
    val d23 = (p3-p2).norm
    val x = (d12*d12+d13*d13-d23*d23)/(2.0*d12)
    val y = math.sqrt(d13*d13-x*x)
    val v = (t2-t1).normalize
    //TODO handle v.y==0. probably handled in edge detection
    if (v.y==0) println("problematic v.y value in projectLastPoint2D")
    val n = if (v.y>0.0) EuclideanVector2D(v.y,-v.x) else EuclideanVector2D(-v.y,v.x)
    t1 + v*x + n*y
  }

  /**
   * shift t_i by -p then rotate around origin such that Rn = (0,0,1) and Rv = (1,0,0).
   * returns the resulting points expressed in two dims
   */
  def planeShift(p:Point[_3D], v:EuclideanVector[_3D], n:EuclideanVector[_3D], t1:Point[_3D], t2:Point[_3D], t3:Point[_3D]): (Point[_2D],Point[_2D],Point[_2D]) = {
    val fn = EuclideanVector3D(0.0,0.0,1.0)
    val fv = EuclideanVector3D(1.0,0.0,0.0)
    val r1 = rotMatrixFromAxisAngle(fn.crossproduct(n),-getInAngle(fn,n))
    val brv = r1*v.toBreezeVector
    val av = EuclideanVector3D(brv(0),brv(1),brv(2))
    val r2 = rotMatrixFromAxisAngle(fv.crossproduct(av),-getInAngle(fv,av))
    val r = r2*r1
    val p1 = r*(t1-p).toBreezeVector
    val p2 = r*(t2-p).toBreezeVector
    val p3 = r*(t3-p).toBreezeVector
    (Point2D(p1(0),p1(1)),Point2D(p2(0),p2(1)),Point2D(p3(0),p3(1)))
  }

  /**
   * rotmatrix to align a to b
   */
  def rotMatrixToAlignVectors(a:EuclideanVector[_3D], b:EuclideanVector[_3D]): DenseMatrix[Double] = {
    val info = getAngleInfo(a,b)
    rotMatrixFromAxisAngle(info.axis,-info.angle)
  }

  /**
   * returns the resulting rotation matrix from turning around the axis 'axis' by the amount 'angle'
   */
  def rotMatrixFromAxisAngle(axis: EuclideanVector[_3D], angle: Double): DenseMatrix[Double] = {
    val a = axis.normalize
    val sin = math.sin(angle)
    val cos = math.cos(angle)
    DenseMatrix(
      (cos+a.x*a.x*(1-cos),     a.x*a.y*(1-cos)-a.z*sin, a.x*a.z*(1-cos)+a.y*sin),
      (a.y*a.x*(1-cos)+a.z*sin, cos+a.y*a.y*(1-cos),     a.y*a.z*(1-cos)-a.x*sin),
      (a.z*a.x*(1-cos)-a.y*sin, a.z*a.y*(1-cos)+a.x*sin, cos+a.z*a.z*(1-cos))
    )
  }

  /**
   * returns the dot product of vectors which are not yet in vectorized form. requires same length vectors
   */
  def dot[D: NDSpace](v1: IndexedSeq[EuclideanVector[D]], v2: IndexedSeq[EuclideanVector[D]]): Double = {
    require(v1.length==v2.length)
    var i = 0
    var s = 0.0
    while(i < v1.length){
      s += v1(i).dot(v2(i))
      i+=1
    }
    s
  }

  /**
   * returns the normalized norm and the previous norm (that was used to divide the vector).
   */
  def normalize[D: NDSpace](v: IndexedSeq[EuclideanVector[D]]): (IndexedSeq[EuclideanVector[D]], Double) = {
    val l2 = L2norm[D]()
    val norm = l2.normVector(v)
    require(norm>1e-11, "provided vector too small. there might be numerical instability")
    (v.map(_/norm), norm)
  }

  def rotX(angle: Double): DenseMatrix[Double] = {
    val sin = math.sin(angle)
    val cos = math.cos(angle)
    DenseMatrix(
      (1.0,0.0,0.0),
      (0.0,cos,-sin),
      (0.0,sin,cos)
    )
  }
  def rotY(angle: Double): DenseMatrix[Double] = {
    val sin = math.sin(angle)
    val cos = math.cos(angle)
    DenseMatrix(
      (cos,0.0,sin),
      (0.0,1.0,0.0),
      (-sin,0.0,cos)
    )
  }
  def rotZ(angle: Double): DenseMatrix[Double] = {
    val sin = math.sin(angle)
    val cos = math.cos(angle)
    DenseMatrix(
      (cos,-sin,0.0),
      (sin,cos,0.0),
      (0.0,0.0,1.0)
    )
  }
}
