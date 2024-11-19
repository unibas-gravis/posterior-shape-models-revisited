package utility

import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common.{PointId, UnstructuredPoints}
import scalismo.geometry._
import scalismo.mesh._
import scalismo.mesh.boundingSpheres.{ClosestPointInTriangle, ClosestPointIsVertex, ClosestPointOnLine, ClosestPointWithType}
import scalismo.statisticalmodel.{MultivariateNormalDistribution, StatisticalMeshModel}

import scala.collection.mutable.ArrayBuffer

object Tchange {
  def mergeTwoSorted[A](d1: IndexedSeq[(Int,A)], d2: IndexedSeq[(Int,A)]): IndexedSeq[A] = {
    val buf = ArrayBuffer[A]()
    var ind1 = 0
    var ind2 = 0
    (0 until d1.length+d2.length).foreach(_ => (ind1, ind2) match {
      case (i1, i2) if i2 == d2.length =>
        buf += d1(i1)._2
        ind1 = i1 + 1
      case (i1, i2) if i1 == d1.length =>
        buf += d2(i2)._2
        ind2 = i2 + 1
      case (i1, i2) =>
        val id1 = d1(i1)
        val id2 = d2(i2)
        if (id1._1<id2._1) {
          buf += id1._2
          ind1 = i1 + 1
        } else {
          buf += id2._2
          ind2 = i2 + 1
        }
    })
    buf.toIndexedSeq
  }
  def seqToVec3D(data: IndexedSeq[Double]): IndexedSeq[EuclideanVector[_3D]] = {
    require(data.length % 3 == 0)
    data.grouped(3).map(t => EuclideanVector3D(t(0),t(1),t(2))).toIndexedSeq
  }
  def seqToDenseVec(data: IndexedSeq[Double],vLength: Int): IndexedSeq[DenseVector[Double]] = {
    require(data.length % vLength == 0)
    data.grouped(vLength).map(t => DenseVector(t.toArray)).toIndexedSeq
  }
  def dvToSeqEv3D(data: DenseVector[Double]): IndexedSeq[EuclideanVector[_3D]] = seqToVec3D(data.toArray.toIndexedSeq)
  def seqToDenseVec[D](data: IndexedSeq[EuclideanVector[D]]): DenseVector[Double] = DenseVector(data.flatMap(_.toArray):_*)
  def seqDvToEv[D: NDSpace](data: IndexedSeq[DenseVector[Double]]): IndexedSeq[EuclideanVector[D]] = data.map(dv => EuclideanVector.fromBreezeVector[D](dv))
  def vecToSeq[D](data: IndexedSeq[EuclideanVector[D]]): IndexedSeq[Double] = data.flatMap(_.toArray)
  def matToVec3D(data: IndexedSeq[IndexedSeq[Double]]): IndexedSeq[IndexedSeq[EuclideanVector[_3D]]] = data.map(seqToVec3D)
  def matToString(mat: DenseMatrix[Double]): String = if (mat.isTranspose) matToString(mat.t) else {
    val sb = mat.toArray.zipWithIndex.foldLeft(new StringBuilder(mat.rows*mat.cols*16,"["))((sb, t) => sb.append(if((t._2+1)%mat.cols==0) t._1+";" else t._1+","))
    sb.dropRight(1).append("]").toString()
  }
  def meshPlusV(mesh: TriangleMesh[_3D], data: IndexedSeq[Double]): TriangleMesh[_3D] = meshPlusEv(mesh, seqToVec3D(data))
  def meshPlusEv(mesh: TriangleMesh[_3D], data: IndexedSeq[EuclideanVector[_3D]]): TriangleMesh[_3D] = {
    mesh.copy(pointSet = UnstructuredPoints(mesh.pointSet.points.toIndexedSeq.zip(data).map(t => t._1+t._2)))
  }
  def meshPlusDv(mesh: TriangleMesh[_3D], data: IndexedSeq[DenseVector[Double]]): TriangleMesh[_3D] = meshPlusEv(mesh, seqDvToEv(data))
  def vecProj[D](a:EuclideanVector[D], b:EuclideanVector[D]): EuclideanVector[D] = b*(a.dot(b)/b.norm2)
  def vecSide[D](a:EuclideanVector[D], b:EuclideanVector[D]): Double = a.dot(b.normalize)
  def projectDiff(target: TriangleMesh[_3D], instance: DenseVector[Double], model: StatisticalMeshModel): TriangleMesh[_3D] = projectDiff(target, model.instance(instance),model)
  def projectDiff(target: TriangleMesh[_3D], instance: TriangleMesh[_3D], model: StatisticalMeshModel): TriangleMesh[_3D] = model.project(meshPlusEv(instance, getDef(target,instance)))
  def handleUcorSurface[D](surface: SurfacePointProperty[D], points: IndexedSeq[ClosestPointWithType]): IndexedSeq[D] = {
    points.map{
      case ClosestPointIsVertex(_,_,pid) => surface.apply(pid)
      case ClosestPointOnLine(_,_,pids,bc) => surface.interpolator.blend(surface.apply(pids._1),surface.apply(pids._2),bc)
      case ClosestPointInTriangle(_,_,tid,bc) => surface.onSurface(tid,bc)
      case _ => throw new NotImplementedError("no default case defined")
    }
  }
  def changeToVertexClps(mesh: TriangleMesh[_3D], points: IndexedSeq[ClosestPointWithType]): IndexedSeq[ClosestPointWithType] = {
    points.map{ //transforms all the points to ClosestPointInTriangle form. needed for symmetric proposals.
      case ClosestPointIsVertex(p,d,pid) =>
        val tr = mesh.triangulation.adjacentTrianglesForPoint(pid).head
        val i = mesh.triangulation.triangle(tr).pointIds.zipWithIndex.find(t => t._1.id==pid.id).get._2
        ClosestPointInTriangle(p,d,tr,BarycentricCoordinates(if(i==0)1.0 else 0.0,if(i==1)1.0 else 0.0,if(i==2)1.0 else 0.0))
      case ClosestPointOnLine(p,d,pids,bc) =>
        val tr = mesh.triangulation.adjacentTrianglesForPoint(pids._1).find(t => mesh.triangulation.triangle(t).pointIds.contains(pids._2)).get
        val bcVals = mesh.triangulation.triangle(tr).pointIds.map(pid => if(pid.id==pids._1.id) 1.0-bc else if (pid.id==pids._2.id) bc else 0.0)
        ClosestPointInTriangle(p,d,tr,BarycentricCoordinates(bcVals(0),bcVals(1),bcVals(2)))
      case p => p
    }
  }
  def applyClpsOnMesh(mesh: TriangleMesh[_3D], points: IndexedSeq[ClosestPointWithType]): IndexedSeq[Point[_3D]] = {
    points.map { //assumes the clps follow the triangulation of the given mesh
      case ClosestPointIsVertex(_, _, pid) => mesh.pointSet.point(pid)
      case ClosestPointOnLine(_, _, pids, bc) =>
        val p1 = mesh.pointSet.point(pids._1)
        val p2 = mesh.pointSet.point(pids._2)
        p1 + (p2-p1) * bc
      case ClosestPointInTriangle(_, _, tid, bcs) =>
        val points = mesh.triangulation.triangle(tid).pointIds.map(mesh.pointSet.point)
        val (p1, p2, p3) = (points(0), points(1), points(2))
        ((p1.toVector * bcs.a) + (p2.toVector * bcs.b) + (p3.toVector * bcs.c)).toPoint
    }
  }
  def handleBarycentricWithTid(mesh: TriangleMesh[_3D], tid: TriangleId, b: BarycentricCoordinates): Point[_3D] = {
    val tr = mesh.triangulation.triangle(tid)
    b.interpolateProperty(mesh.pointSet.point(tr.ptId1),mesh.pointSet.point(tr.ptId2),mesh.pointSet.point(tr.ptId3))
  }
  def getZeroMeanNoise(sigma2:Double = 1.0): MultivariateNormalDistribution = {
    val zm = DenseVector.zeros[Double](3)
    val eye = sigma2 * DenseMatrix.eye[Double](3)
    MultivariateNormalDistribution(zm,eye)
  }

  //uses closest point on surface to infer correspondence on the 'to' object
  def getDef(to: TriangleMesh[_3D], from: TriangleMesh[_3D]): IndexedSeq[EuclideanVector[_3D]] = getDef(to,from.pointSet.points)
  def getDef(to: TriangleMesh[_3D], from: UnstructuredPoints[_3D]): IndexedSeq[EuclideanVector[_3D]] = getDef(to,from.points)
  def getDef(to: TriangleMesh[_3D], from: IndexedSeq[Point[_3D]]): IndexedSeq[EuclideanVector[_3D]] = getDef(to,from.toIterator)
  def getDef(to: TriangleMesh[_3D], from: Iterator[Point[_3D]]): IndexedSeq[EuclideanVector[_3D]] = from.map(fp => to.operations.closestPointOnSurface(fp).point-fp).toIndexedSeq
  //assumes correspondence
  def getDef(to: UnstructuredPoints[_3D], from: TriangleMesh[_3D]): IndexedSeq[EuclideanVector[_3D]] = getDef(to,from.pointSet.points)
  def getDef(to: UnstructuredPoints[_3D], from: UnstructuredPoints[_3D]): IndexedSeq[EuclideanVector[_3D]] = getDef(to,from.points)
  def getDef(to: UnstructuredPoints[_3D], from: IndexedSeq[Point[_3D]]): IndexedSeq[EuclideanVector[_3D]] = getDef(to,from.toIterator)
  def getDef(to: UnstructuredPoints[_3D], from: Iterator[Point[_3D]]): IndexedSeq[EuclideanVector[_3D]] = to.points.zip(from).map(t => t._1-t._2).toIndexedSeq
  //list of unique ids on to that are mapped to. unsorted
  def getObsPids(to: TriangleMesh[_3D], from: TriangleMesh[_3D]): IndexedSeq[PointId] = to.pointSet.pointSequence.map(p => from.pointSet.findClosestPoint(p).id.id).distinct.map(PointId)
  def getMean(mesh: TriangleMesh[_3D], ids: Option[IndexedSeq[PointId]] = None): Point[_3D] = {
    val pids = ids.getOrElse(mesh.pointSet.pointIds.toIndexedSeq)
    pids.map(mesh.pointSet.point).reduce(_+_.toVector).map(_*(1.0/pids.length))
  }
  //assumes correspondence between all meshes and takes the trianglelist from the first mesh
  def getMean(meshes: IndexedSeq[TriangleMesh[_3D]]): TriangleMesh[_3D] = {
    val points = meshes.map(_.pointSet.pointSequence).transpose.map(_.reduce(_+_.toVector).map(_/meshes.length))
    TriangleMesh3D.apply(points, meshes.head.triangulation)
  }
  def getMeanPs(pointSeqs: IndexedSeq[IndexedSeq[Point[_3D]]]): IndexedSeq[Point[_3D]] = {
    pointSeqs.transpose.map(_.reduce(_+_.toVector).map(_/pointSeqs.length))
  }
  def getMean(points: IndexedSeq[Point[_3D]]): Point[_3D] = {
    points.map(_.toVector).reduce(_+_).map(_/points.length).toPoint
  }
  def getMeanVar(data: IndexedSeq[Double]): MeanVar[Double] = {
    require(data.length >= 2)
    val mean = data.sum/data.length
    MeanVar[Double](mean, data.map(d => (d-mean) * (d-mean)).sum / data.length)
  }
  def getMeanCovDV(data: IndexedSeq[DenseVector[Double]]): MeanCovDV = {
    require(data.length >= 2)
    val mean: DenseVector[Double] = data.reduce(_+_) / data.length.toDouble
    val variance: DenseMatrix[Double] = data.foldLeft(DenseMatrix.zeros[Double](mean.length, mean.length))((s, v) =>
      s + ((v-mean) * (v-mean).t) / data.length.toDouble
    )
    MeanCovDV(mean, variance)
  }

  def getMeanVarDV(data: IndexedSeq[DenseVector[Double]]): MeanVarDV = {
    require(data.length >= 2)
    val mean: DenseVector[Double] = data.reduce(_ + _) / data.length.toDouble
    val variance: DenseVector[Double] = data.foldLeft(DenseVector.zeros[Double](mean.length))((s, v) =>
      s + (v-mean)*:*(v-mean) / data.length.toDouble //*:* is elementwise multiplication
    )
    MeanVarDV(mean, variance)
  }
  case class MeanVar[T](mean: T, variance: T) {def m: T = mean; def v: T = variance}
  case class MeanCovDV(mean: DenseVector[Double], variance: DenseMatrix[Double]) {def m: DenseVector[Double] = mean; def v: DenseMatrix[Double] = variance}
  case class MeanVarDV(mean: DenseVector[Double], variance: DenseVector[Double]) {def m: DenseVector[Double] = mean; def v: DenseVector[Double] = variance}

  def reverseMapping[A](data: IndexedSeq[A], f: (A) => Option[A]): (A) => Option[A] = {
    val map = data.flatMap(a => f(a).map(b => (b, a))).toMap
    (a: A) => map.get(a)
  }
}
