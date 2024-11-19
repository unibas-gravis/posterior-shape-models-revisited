package io

import breeze.linalg.{DenseMatrix, DenseVector}

import scala.{specialized => spec}
import io.RunWriter.{matrixToString, vectorToString}

import java.io.{BufferedWriter, File, FileWriter}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import scalismo.geometry.EuclideanVector

import scala.reflect.ClassTag

/**
 * contains various useful writing functions for data. basic writing is synchronized
 */
class RunWriter(writer: BufferedWriter, pythonStyle: Boolean){
  val bufferMat = scala.collection.mutable.HashMap.empty[String,Traversable[Traversable[Any]]]
  val bufferMatS = scala.collection.mutable.HashMap.empty[String,Traversable[(Traversable[(Any,Int)],Int)]]
  val bufferVec = scala.collection.mutable.HashMap.empty[String,Traversable[Any]]
  val bufferVecS = scala.collection.mutable.HashMap.empty[String,Traversable[(Any,Int)]]

  def writeSingle[A](data: A, name: String): Unit = {
    this.synchronized {
      writer.write(s"$name=${data.toString};\n")
      writer.flush()
    }
  }
  def writeVector[A](data: Traversable[A], name: String): Unit = {
    this.synchronized {
      writer.write(s"$name=[${vectorToString(data)}];\n")
      writer.flush()
    }
  }
  def writeVectorEV[D](data: Traversable[EuclideanVector[D]], name: String): Unit = {
    this.synchronized {
      writer.write(s"$name=[${vectorToString(data.flatMap(_.toArray))}];\n")
      writer.flush()
    }
  }
  def writeVectorDV[D](data: Traversable[DenseVector[D]], name: String): Unit = writeVector(data.flatMap(dv => dv.data),name)
  def writeMatrix[A](data: Traversable[Traversable[A]], name: String): Unit = {
    this.synchronized {
      writer.write(s"$name=[${matrixToString(data, pythonStyle)}];\n")
      writer.flush()
    }
  }
  def writeMatrix[@spec(Double, Int, Float, Long) A](data: DenseMatrix[A], name: String, rowFirst: Boolean=true): Unit = {
    this.synchronized {
      writer.write(s"$name=[${matrixToString(data, rowFirst, pythonStyle)}];\n")
      writer.flush()
    }
  }
  def writeMatrixEV[D](data: Traversable[Traversable[EuclideanVector[D]]], name: String): Unit = {
    this.synchronized {
      writer.write(s"$name=[${matrixToString(data.map(_.flatMap(_.toArray)), pythonStyle)}];\n")
      writer.flush()
    }
  }
  def writeCollected[A](data: A, name: String): Unit = {
    this.synchronized {
      bufferVec.get(name) match {
        case Some(cur) => bufferVec.update(name, cur++Traversable(data))
        case None => bufferVec.put(name, Traversable(data))
      }
    }
  }
  def writeCollectedSorted[A](data: A, index: Int, name: String): Unit = {
    this.synchronized {
      bufferVecS.get(name) match {
        case Some(cur) => bufferVecS.update(name, cur++Traversable((data,index)))
        case None => bufferVecS.put(name, Traversable((data,index)))
      }
    }
  }
  def writeCollected[A](data: Traversable[A], name: String): Unit = {
    this.synchronized {
      bufferMat.get(name) match {
        case Some(cur) => bufferMat.update(name, cur++Traversable(data))
        case None => bufferMat.put(name, Traversable(data))
      }
    }
  }
  def writeCollectedSorted[A](data: Traversable[(A,Int)], index: Int, name: String): Unit = {
    this.synchronized {
      bufferMatS.get(name) match {
        case Some(cur) => bufferMatS.update(name, cur++Traversable((data,index)))
        case None => bufferMatS.put(name, Traversable((data,index)))
      }
    }
  }
  def writeCollectedEV[D](data: Traversable[EuclideanVector[D]], name: String): Unit = writeCollected(data.flatMap(_.toArray), name)
  def writeCollectedEVsorted[D](data: Traversable[(EuclideanVector[D],Int)], index: Int, name: String): Unit = {
    writeCollectedSorted(data.flatMap(t => {
      t._1.toArray.zipWithIndex.map(d => (d._1,t._2*t._1.dimensionality+d._2))
    }), index, name)
  }
  def writeLiteral(data:String): Unit = {
    this.synchronized {
      writer.write(s"$data\n")
      writer.flush()
    }
  }
  def flushCollected(): Unit = {
    this.synchronized {
      bufferMat.foreach{case (name, data) => writeMatrix(data, name)}
      bufferMat.clear()
      bufferVec.foreach{case (name, data) => writeVector(data, name)}
      bufferVec.clear()
      bufferVecS.foreach{case (name, data) => writeVector(data.toIndexedSeq.sortBy(_._2).map(_._1), name)}
      bufferVecS.clear()
      bufferMatS.foreach{case (name, data) => writeMatrix(data.toIndexedSeq.sortBy(_._2).map(_._1.toIndexedSeq.sortBy(_._2).map(_._1)), name)}
      bufferMatS.clear()
      writer.flush()
    }
  }
  def flushCollected(keys: Set[String]): Unit = {
    this.synchronized {
      bufferMat.foreach { case (name, data) if keys.contains(name) => writeMatrix(data, name) case _ => }
      bufferMat.clear()
      bufferVec.foreach { case (name, data) if keys.contains(name) => writeVector(data, name) case _ => }
      bufferVec.clear()
      bufferVecS.foreach { case (name, data) if keys.contains(name) => writeVector(data.toIndexedSeq.sortBy(_._2).map(_._1), name) case _ => }
      bufferVecS.clear()
      bufferMatS.foreach { case (name, data) if keys.contains(name) => writeMatrix(data.toIndexedSeq.sortBy(_._2).map(_._1.toIndexedSeq.sortBy(_._2).map(_._1)), name) case _ => }
      bufferMatS.clear()
      writer.flush()
    }
  }
  def close(): Unit = {
    this.synchronized {
      flushCollected()
      writer.close()
    }
  }
}

object RunWriter {
  /**
   * replace the two characters $t with a timeStamp
   */
  def apply(file: String, pythonStyle: Boolean): (RunWriter, File) = {
    val f = new File(file.replace("$t",LocalDateTime.now.format(DateTimeFormatter.ofPattern("YYYYMMdd_HHmmss"))))
    (RunWriter(f, pythonStyle), f)
  }
  def apply(file: File, pythonStyle: Boolean): RunWriter = new RunWriter(new BufferedWriter(new FileWriter(file)), pythonStyle)

  def apply(pythonStyle: Boolean): RunWriter = {
    new RunWriter(new BufferedWriter(new FileWriter(new File(s"data/runData/${LocalDateTime.now.format(DateTimeFormatter.ofPattern("YYYYMMdd_HHmmss"))}"))), pythonStyle)
  }

  def vectorToString[A](data: Traversable[A]): String = {
    val sb = data.foldLeft(new StringBuilder())((sb, t) => sb.append(t.toString+","))
    sb.dropRight(1).toString()
  }

  def vectorToString(data: DenseVector[Double]): String = {
    vectorToString(data.toArray)
  }

  def matrixToString[A](data: Traversable[Traversable[A]], pythonStyle: Boolean): String = {
    if (pythonStyle){
      val sb = data.foldLeft(new StringBuilder())((sb, t) => {
        sb.append("[")
        val sb2 = t.foldLeft(sb)((sb2, t2) => sb2.append(t2.toString + ","))
        sb2.dropRight(1).append("],")
      })
      sb.dropRight(1).toString()
    } else {
      val sb = data.foldLeft(new StringBuilder())((sb, t) => {
        val sb2 = t.foldLeft(sb)((sb2, t2) => sb2.append(t2.toString + ","))
        sb2.dropRight(1).append(";")
      })
      sb.dropRight(1).toString()
    }
  }

  def matrixToString[@spec(Double, Int, Float, Long) A](data: DenseMatrix[A], rowFirst: Boolean, pythonStyle: Boolean): String = {
    val data2d = if (rowFirst) {
      (0 until data.cols).map(y => data(::,y).iterator.map(_._2).toIndexedSeq)
    } else {
      (0 until data.rows).map(x => data(x,::).t.iterator.map(_._2).toIndexedSeq)
    }
    matrixToString[A](data2d, pythonStyle)
  }
}
