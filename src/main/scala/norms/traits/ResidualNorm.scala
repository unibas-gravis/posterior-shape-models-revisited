package norms.traits

import scalismo.geometry.EuclideanVector

trait ResidualNorm[D] extends VectorNorm[D] with RealNorm[D] {
  override def norm2Vector(residual: IndexedSeq[EuclideanVector[D]]): Double
}
