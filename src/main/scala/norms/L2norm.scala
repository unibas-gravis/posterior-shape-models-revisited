package norms

import norms.traits.ResidualNorm
import scalismo.geometry.EuclideanVector

class L2norm[D] extends ResidualNorm[D] {
  override def norm2Vector(residual: IndexedSeq[EuclideanVector[D]]): Double = {
    var i = 0
    var s = 0.0
    while (i < residual.length){
      s+=residual(i).norm2
      i+=1
    }
    s
  }
  override def norm2(residual: IndexedSeq[Double]): Double = {
    var i = 0
    var s = 0.0
    while (i < residual.length){
      val r = residual(i)
      s+=r*r
      i+=1
    }
    s
  }
}

object L2norm {
  def apply[D](): L2norm[D] = new L2norm[D]()
}
