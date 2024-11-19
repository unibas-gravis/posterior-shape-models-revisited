package utility

import scalismo.common.PointId
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.mesh.TriangleMesh

object MeshUtils {

  /**
   * sorts points based on distance to growP
   */
  def getPartialMeshPoint(mesh: TriangleMesh[_3D], growP: Point[_3D], ratio: Double, ratioForInterval: Boolean = true): (TriangleMesh[_3D], IndexedSeq[PointId], PointId => Option[PointId]) = {
    cutTo(mesh, (p,_) => (p-growP).norm2, ratio, ratioForInterval)
  }

  /**
   * sorts points based on dot product with growV
   */
  def getPartialMeshVector(mesh: TriangleMesh[_3D], growV: EuclideanVector[_3D], ratio: Double, ratioForInterval: Boolean = true): (TriangleMesh[_3D], IndexedSeq[PointId], PointId => Option[PointId]) = {
    cutTo(mesh, (p,_) => p.toVector.dot(growV), ratio, ratioForInterval)
  }

  /**
   * cuts given mesh down to provided point ids
   */
  def makePartial(mesh: TriangleMesh[_3D], obsset: Set[Int]): TriangleMesh[_3D] =
    mesh.operations.maskPoints(pid => obsset.contains(pid.id)).transformedMesh

  /**
   * returns (new mesh, old ids that are kept, a lazy mapping from old to new ids)
   * ratioForInterval if true ratio is relative to length; if false number of points
   */
  def cutTo(mesh: TriangleMesh[_3D], order: (Point[_3D],PointId)=>Double, ratio: Double, ratioForInterval: Boolean): (TriangleMesh[_3D], IndexedSeq[PointId], PointId => Option[PointId]) = {
    val ids = if (ratioForInterval) {
      val ordered = mesh.pointSet.pointsWithId.toIndexedSeq.map(t => (t, order(t._1,t._2))).sortBy(_._2)
      val h = (ordered.last._2-ordered.head._2)*ratio+ordered.head._2
      ordered.takeWhile(t => t._2 < h).map(_._1._2)
    } else {
      mesh.pointSet.pointsWithId.toIndexedSeq.sortBy(t => order(t._1,t._2))
        .take((mesh.pointSet.numberOfPoints*ratio).toInt).map(_._2)
    }
    val idset = ids.map(_.id).toSet
    val tm = mesh.operations.maskPoints(pid => idset.contains(pid.id)).transformedMesh
    (tm, ids, pid => {
      if (idset.contains(pid.id)) {
        Option(tm.pointSet.findClosestPoint(mesh.pointSet.point(pid)).id)
      } else None
    })
  }

  /**
   * iteratively cut the points of the mesh that have undefined/nan vertex normals.
   */
  def cutPointsWithUndefinedNormals(mesh: TriangleMesh[_3D], maxIter: Int=10): TriangleMesh[_3D] = {
    (0 until maxIter).foldLeft((mesh, false))((t,_) => {
      val (mesh, converged) = t
      if (converged) t else {
        val pids = mesh.pointSet.pointIds.filter(pid => mesh.vertexNormals.apply(pid).toArray.exists(d => d.isNaN)).map(_.id).toSet
        if (pids.isEmpty) (mesh, true) else {
          val umesh = mesh.operations.maskPoints(pid => !pids.contains(pid.id)).transformedMesh
          (umesh,false)
        }
      }
    })._1
  }

}
