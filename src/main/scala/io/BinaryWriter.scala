package io

import breeze.linalg.{DenseMatrix, DenseVector}

import java.io.{DataOutputStream, FileOutputStream}

class BinaryWriter(filePath: String) {
  val outputStream = new DataOutputStream(new FileOutputStream(filePath))

  def writedv[A](dv: DenseVector[A])(implicit writer: BinaryWritable[A]): Unit = {
    write1dSeq[A](dv.data.toIndexedSeq)
  }
  def writedm[A](dm: DenseMatrix[A])(implicit writer: BinaryWritable[A]): Unit = {
    this.synchronized {
      //dimensions
      outputStream.writeByte(2)
      outputStream.writeInt(dm.rows)
      outputStream.writeInt(dm.cols)
      outputStream.writeInt(writer.getBytes())
      //data
      dm.data.foreach(write[A](_))
      outputStream.flush()
    }
  }

  def write2dSeq[A](arr: IndexedSeq[IndexedSeq[A]])(implicit writer: BinaryWritable[A]): Unit = {
    this.synchronized {
      //dimensions
      outputStream.writeByte(2)
      outputStream.writeInt(arr.length)
      outputStream.writeInt(arr.head.length)
      outputStream.writeInt(writer.getBytes())
      //data
      arr.foreach(_.foreach(write[A](_)))
      outputStream.flush()
    }
  }

  def write1dSeq[A](arr: IndexedSeq[A])(implicit writer: BinaryWritable[A]): Unit = {
    this.synchronized {
      //dimension
      outputStream.writeByte(1)
      outputStream.writeInt(arr.length)
      outputStream.writeInt(writer.getBytes())
      //data
      arr.foreach(write[A](_))
      outputStream.flush()
    }
  }

  def write[A](item: A)(implicit writer: BinaryWritable[A]): Unit = {
    writer.write(outputStream, item)
  }

  def close(): Unit = {
    this.synchronized {
      outputStream.flush()
      outputStream.close()
    }
  }
}

object BinaryWriter {
  def apply(filePath: String): BinaryWriter = {
    new BinaryWriter(filePath)
  }
}

