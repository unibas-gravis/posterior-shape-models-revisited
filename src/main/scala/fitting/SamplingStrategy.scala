package fitting

import scalismo.common.PointId
import scalismo.geometry.{Point, _3D}
import scalismo.mesh.TriangleMesh

trait SamplingStrategy {
  def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])]
}
case class ModelSampling() extends SamplingStrategy {
  override def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = {
    from.pointSet.pointsWithId.map(t => (t._2, to.operations.closestPointOnSurface(t._1).point)).toIndexedSeq
  }
}
case class TargetSampling() extends SamplingStrategy {
  override def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = {
    to.pointSet.pointSequence.map(tp => {
      val clp = from.pointSet.findClosestPoint(tp)
      (clp.id,tp)
    })
  }
}

case class TargetSamplingUnique() extends SamplingStrategy {
  override def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = {
    to.pointSet.pointSequence.map(tp => {
      val clp = from.pointSet.findClosestPoint(tp)
      (clp,tp) //similar to TargetSampling but choosing minimum based on distance
    }).groupBy(_._1.id.id).map(t => t._2.minBy(cp => (cp._2-cp._1.point).norm)).map(t => (t._1.id,t._2)).toIndexedSeq
  }
}

/**
 * useful for posterior calculation for partial targets
 */
case class BidirectionalSamplingFromTarget() extends SamplingStrategy {
  override def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = {
    to.pointSet.pointSequence.map(tp => {
      val clp = from.pointSet.findClosestPoint(tp)
      clp.id
    }).distinct.map(pid => (pid, to.operations.closestPointOnSurface(from.pointSet.point(pid)).point))
  }
}

/**
 * useful for posterior calculation for partial targets where the targets are not clean or noisy.
 */
case class BidirectionalSamplingFromOrigin() extends SamplingStrategy {
  override def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = {
    from.pointSet.pointsWithId.flatMap(ppid => {
      val clp = from.operations.closestPointOnSurface(ppid._1).point
      if (from.pointSet.findClosestPoint(clp).id.id == ppid._2.id) Option((ppid._2,clp)) else None
    }).toIndexedSeq
  }
}
