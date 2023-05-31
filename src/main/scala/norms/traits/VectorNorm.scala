package norms.traits

import scalismo.geometry.EuclideanVector

trait VectorNorm[D] {
  def norm2Vector(residual: IndexedSeq[EuclideanVector[D]]): Double
  def normVector(residual: IndexedSeq[EuclideanVector[D]]): Double = math.sqrt(norm2Vector(residual))
}
