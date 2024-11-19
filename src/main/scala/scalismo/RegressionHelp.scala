package scalismo

import norms.L2norm
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.StatisticalMeshModel
import utility.MathHelp

object RegressionHelp {

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

  def getScalingEigf(mesh: TriangleMesh[_3D], center: Point[_3D]): (IndexedSeq[EuclideanVector[_3D]], Double) = {
    val vecs = mesh.pointSet.pointSequence.map(_ - center)
    val normz = 1.0 / L2norm[_3D]().norm2Vector(vecs)
    (vecs.map(_.*(normz)), 1.0 / normz) //also returns scalar needed for a single full scaling
  }

}
