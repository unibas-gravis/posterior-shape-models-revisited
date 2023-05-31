package norms.traits

trait RealNorm[D] {
  def norm2(residual: IndexedSeq[Double]): Double
  def norm(residual: IndexedSeq[Double]): Double = math.sqrt(norm2(residual))
}
