package fitting

import scalismo.common.PointId
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.mesh.TriangleMesh
import scalismo.mesh.boundingSpheres.{ClosestPointInTriangle, ClosestPointIsVertex, ClosestPointOnLine}
import utility.{MathHelp, Tchange}


trait SamplingStrategy {
  /**
   * returns the PointIds of the 'from' pointcloud and associated point on the 'to' mesh and a double to scale the certainty of that observation
   */
  def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D], Double)] = establishCorrespondenceBackground(from, to)._1
  def establishCorrespondenceBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], Double)], Int)
  def cutof:Double = 1e-2
  /**
   * transforms the observation strength to a uncertainty factor, also removes almost useless observations.
   * best used with a flatmap. given doubles should be positive
   */
  def observationStrengthToUncertainty(d:Double): Option[Double] = {
    if (d < cutof) None else Some(1.0/d)
  }
  def filterObservationStrengthToUncertainty(d:Double): Boolean = d >= cutof
}

/**
 * Similar to TargetSampling. But here we find the closest point on the surface, and then split the observations to the
 * vertices and weight them based on the barycentric coordinates. should not be used to find domain
 * 0-n to 1
 */
case class TargetSamplingSplit() extends SamplingStrategy {
  override def establishCorrespondenceBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], Double)], Int) = {
    val cor = to.pointSet.pointSequence.flatMap(tp => {
      val clp = from.operations.closestPointOnSurface(tp)
      val obs = clp match {
        case ClosestPointIsVertex(_, _, pid) => IndexedSeq((pid, tp, 1.0))
        case ClosestPointOnLine(_, _, pids, bc) => IndexedSeq((pids._1, tp, bc),(pids._2, tp, 1.0-bc))
        case ClosestPointInTriangle(_, _, tid, bc) =>
          val triangle = from.triangulation.triangle(tid)
          IndexedSeq((triangle.ptId1, tp, bc.a),(triangle.ptId2, tp, bc.b),(triangle.ptId3, tp, bc.c))
        case _ => IndexedSeq.empty
      }
      obs.filter(t => filterObservationStrengthToUncertainty(t._3))
    })
    (cor, to.pointSet.numberOfPoints - cor.length) //bkg always 0
  }
}

/**
 * adds surface normals to the observations.
 */
trait SamplingStrategyNormals extends SamplingStrategy {
  override def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D], Double)] = establishCorrespondenceNormal(from, to).map(t => (t._1, t._2, t._4))
  def establishCorrespondenceNormal(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D], EuclideanVector[_3D], Double)] = establishCorrespondenceNormalBackground(from, to)._1
  def establishCorrespondenceNormalBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], EuclideanVector[_3D], Double)], Int)
  override def establishCorrespondenceBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], Double)], Int) = {
    val res = establishCorrespondenceNormalBackground(from, to)
    (res._1.map(t => (t._1,t._2,t._4)), res._2)
  }
}

/**
 * extends the provided strategy by only looking up the vertex normal at the pointid location. more an approximation of
 * the strategies that use the normal at the 'to' location.
 */
case class NormalSamplingSimpleExtension(strategy: SamplingStrategy) extends SamplingStrategyNormals {
  override def establishCorrespondenceNormalBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], EuclideanVector[_3D], Double)], Int) = {
    val (corr, bkg) = strategy.establishCorrespondenceBackground(from, to)
    (corr.map(t => (t._1,t._2,from.vertexNormals.apply(t._1), t._3)), bkg)
  }
}

/**
 * similar to NormalSamplingSimpleExtension but takes the normal of the 'to' shape.
 * This is a lot more expensive due to more clp operations; it is rarely useful
 */
case class NormalSamplingTargetExtension(strategy: SamplingStrategy) extends SamplingStrategyNormals {
  override def establishCorrespondenceNormalBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], EuclideanVector[_3D], Double)], Int) = {
    val (corr, bkg) = strategy.establishCorrespondenceBackground(from, to)
    val targetNormals = Tchange.handleUcorSurface(to.vertexNormals, corr.map(t => to.operations.closestPointOnSurface(t._2)))
    (corr.zip(targetNormals).map(t => (t._1._1,t._1._2, t._2, t._1._3)), bkg)
  }
}

/**
 * establishes correspondence by sampling over to triangles and checking normal direction for intersection.
 * uses no compatibility criteria
 * 0-n to 1vertex
 */
case class TargetNormalSampling() extends SamplingStrategyNormals {
  val samp = TargetCompatNormalSampling(Some((_:Double) => true))
  override def establishCorrespondenceNormalBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], EuclideanVector[_3D], Double)], Int) = {
    samp.establishCorrespondenceNormalBackground(from, to)
  }
}

/**
 * establishes correspondence by sampling over to triangles and checking normal direction for intersection.
 * also filters found intersections with a compatibility criteria on the norm
 * 0-n to 1vertex
 */
case class TargetCompatNormalSampling(compatiblef: Option[Double=>Boolean]=Some(_>0.5)) extends SamplingStrategyNormals {
  val usecf = compatiblef.isDefined
  val cf: Double=>Boolean = compatiblef.orNull
  override def establishCorrespondenceNormalBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], EuclideanVector[_3D], Double)], Int) = {
    val obs = to.triangulation.triangleIds.flatMap(tid => {
      val tr = to.triangulation.triangle(tid)
      val center = Tchange.getMean(tr.pointIds.map(to.pointSet.point))
      val normal = to.cellNormals.apply(tid)
      val clpsAll = from.operations.getIntersectionPointsOnSurface(center, normal)
      val clps = if (usecf) clpsAll.filter(t => cf(normal.dot(from.cellNormals(t._1)))) else clpsAll
      val distsTups = clps.map(t => {
        val pids = from.triangulation.triangle(t._1).pointIds
        val ps = pids.map(from.pointSet.point)
        val p = ps.zip(IndexedSeq(t._2.a,t._2.b,t._2.c)).map(t => t._1.toVector.*(t._2)).reduce(_+_).toPoint
        (pids, ps, p, (p-center).norm)
      })
      val dattups = clps.zip(distsTups)
      if (dattups.nonEmpty){
        IndexedSeq(dattups.minBy(_._2._4)).flatMap(clpInfo => {
          val pids = clpInfo._2._1
          val ps = clpInfo._2._2
          val pcenter = clpInfo._2._3
          val difs = ps.map(_-pcenter)
          val obsps = difs.map(d => center + d)
          val weights = IndexedSeq(clpInfo._1._2.a, clpInfo._1._2.b, clpInfo._1._2.c).map(observationStrengthToUncertainty)
          (0 until 3).flatMap(i => weights(i) match {
            case Some(w) => Some(pids(i), obsps(i), normal, w)
            case _ => None
          })
        })
      } else None
    })
    (obs, obs.map(_._4).sum.toInt)
  }
}

trait SamplingStrategyUniform extends SamplingStrategy {
  /**
   * returns the PointIds of the 'from' pointcloud and associated point on the 'to' mesh
   */
  def establishCorrespondenceUniform(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = establishCorrespondenceUniformBackground(from, to)._1

  override def establishCorrespondence(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D], Double)] = establishCorrespondenceUniform(from, to).map(t => (t._1, t._2, 1.0))
  def establishCorrespondenceUniformBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D])], Int)

  override def establishCorrespondenceBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D], Double)], Int) = {
    val res = establishCorrespondenceUniformBackground(from, to)
    (res._1.map(t => (t._1,t._2,1.0)), res._2)
  }
}

/**
 * every point on the 'from' mesh is associated with the closes point on the surface of the 'to' mesh.
 * 1 to 0-n
 */
case class ModelSampling() extends SamplingStrategyUniform {
  override def establishCorrespondenceUniformBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D])], Int) = {
    val cor = from.pointSet.pointsWithId.map(t => (t._2, to.operations.closestPointOnSurface(t._1).point)).toIndexedSeq
    (cor, from.pointSet.numberOfPoints - cor.length) //bkg always 0
  }
}

/**
 * every point on the 'from' point cloud that is a closest point of a point on the 'to' mesh is associated with the starting
 * 'to' point.
 * 0-n to 1
 */
case class TargetSampling() extends SamplingStrategyUniform {
  override def establishCorrespondenceUniformBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D])], Int) = {
    val cor = to.pointSet.pointSequence.map(tp => {
      val clp = from.pointSet.findClosestPoint(tp)
      (clp.id,tp)
    })
    (cor, to.pointSet.numberOfPoints - cor.length) //bkg always 0
  }
}

/**
 * same as TargetSampling but filters afterwards the mappings and chooses the closest one.
 * 0-1 to 1
 */
case class TargetSamplingUnique() extends SamplingStrategyUniform {
  override def establishCorrespondenceUniformBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D])], Int) = {
    val cor = to.pointSet.pointSequence.map(tp => {
      val clp = from.pointSet.findClosestPoint(tp)
      (clp,tp) //similar to TargetSampling but choosing minimum based on distance
    }).groupBy(_._1.id.id).map(t => t._2.minBy(cp => (cp._2-cp._1.point).norm)).map(t => (t._1.id,t._2)).toIndexedSeq
    (cor, to.pointSet.numberOfPoints - cor.length)
  }
}

/**
 * useful for posterior calculation for partial targets. Uses TargetSampling to identify candidates for modelSampling.
 * 0-1 to 0-1
 */
case class BidirectionalSamplingFromTarget() extends SamplingStrategyUniform {
  override def establishCorrespondenceUniformBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D])], Int) = {
    val cor = to.pointSet.pointSequence.map(tp => {
      val clp = from.pointSet.findClosestPoint(tp)
      clp.id
    }).distinct.map(pid => (pid, to.operations.closestPointOnSurface(from.pointSet.point(pid)).point))
    (cor, to.pointSet.numberOfPoints - cor.length)
  }
}

/**
 * useful for posterior calculation for partial targets where the targets are not clean or noisy. Starts with
 * ModelSampling and then checks if the backward direction is consistent.
 * 0-1 to 0-1
 */
case class BidirectionalSamplingFromOrigin() extends SamplingStrategyUniform {
  override def establishCorrespondenceUniformBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D])], Int) = {
    val cor = from.pointSet.pointsWithId.flatMap(ppid => {
      val clp = to.operations.closestPointOnSurface(ppid._1).point
      if (from.pointSet.findClosestPoint(clp).id.id == ppid._2.id) Option((ppid._2, clp)) else None
    }).toIndexedSeq
    (cor , from.pointSet.numberOfPoints - cor.length)
  }
}

/**
 * this strategy allows usage of the more general setup if correspondence is given due to synthetic data
 * @param mapping maps the pids from 'from' to 'to'
 */
case class Correspondence(mapping: (PointId) => Option[PointId]) extends SamplingStrategyUniform {
  override def establishCorrespondenceUniformBackground(from: TriangleMesh[_3D], to: TriangleMesh[_3D]): (IndexedSeq[(PointId, Point[_3D])], Int) = {
    val cor = from.triangulation.pointIds.flatMap(pid => {
      mapping(pid).map(cpid => (pid, to.pointSet.point(cpid)))
    })
    (cor, from.pointSet.numberOfPoints - cor.length)  // should have constant background
  }
}
