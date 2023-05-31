package utility

import utility.Control.TimeUnits.{Ms, Ns, S, TimeUnit}

import scala.util.control.NonFatal

object Control {
  def withResources[T <: AutoCloseable, V](r: => T)(f: T => V): V = {
    val resource: T = r
    require(resource != null, "resource is null")
    var exception: Throwable = null
    try {
      f(resource)
    } catch {
      case NonFatal(e) =>
        exception = e
        throw e
    } finally {
      closeAndAddSuppressed(exception, resource)
    }
  }

  private def closeAndAddSuppressed(e: Throwable, resource: AutoCloseable): Unit = {
    if (e != null) {
      try {
        resource.close()
      } catch {
        case NonFatal(suppressed) => e.addSuppressed(suppressed)
      }
    } else {
      resource.close()
    }
  }

  object TimeUnits extends Enumeration {
    type TimeUnit = Value

    val Ns, Ms, S = Value
  }

  def toNs(amount: Long, unit: TimeUnit): Long = {
    unit match {
      case Ns => amount
      case Ms => amount * 1000 * 1000
      case S => amount * 1000 * 1000 * 1000
    }
  }

  def toMs(amount: Long, unit: TimeUnit): Long = {
    unit match {
      case Ns => amount / (1000 * 1000)
      case Ms => amount
      case S => amount * 1000
    }
  }

  def toS(amount: Long, unit: TimeUnit): Double = {
    unit match {
      case Ns => amount / (1000.0 * 1000.0 * 1000.0)
      case Ms => amount / 1000.0
      case S =>amount
    }
  }

  /**
   * returns the value of f and the time in ns to execute f. use benchmarkit for benchmarking
   */
  def timeit[T](f: () => T): (T,Long) = {
    val time = System.nanoTime()
    val ret = f.apply()
    val measured = System.nanoTime()-time
    (ret, measured)
  }

  def timeitVerbose[T](f: () => T, prefix: String = ""): T = {
    val time = System.nanoTime()
    val ret = f.apply()
    val measured = System.nanoTime()-time
    println(s"measured $prefix : $measured ns ${measured/1e6.toInt} ms ${measured/1e9} s ${if (measured > 60e9) s"${measured/60e9} min" else ""}")
    ret
  }

  /**
   * benchmarks the function and returns the average time in nanoseconds
   */
  def benchmarkit[T](f: () => T, repeats: Int, warmup: Int = 0): Long = {
    for (_ <- 0 until warmup) f.apply()
    val time = System.nanoTime()
    for (_ <- 0 until repeats) f.apply()
    val measured = System.nanoTime()-time
    measured/repeats
  }
}