package io

import java.io.DataOutputStream

trait BinaryWritable[A] {
  def write(outputStream: DataOutputStream, value: A): Unit
  def getBytes(): Int
}

object BinaryWritable {
  implicit val intWritable: BinaryWritable[Int] = new BinaryWritable[Int] {
    override def write(outputStream: DataOutputStream, value: Int): Unit = {
      outputStream.writeInt(value)
    }

    override def getBytes(): Int = 4
  }

  implicit val doubleWritable: BinaryWritable[Double] = new BinaryWritable[Double] {
    override def write(outputStream: DataOutputStream, value: Double): Unit = {
      outputStream.writeDouble(value)
    }
    override def getBytes(): Int = 8
  }

  implicit val booleanWritable: BinaryWritable[Boolean] = new BinaryWritable[Boolean] {
    override def write(outputStream: DataOutputStream, value: Boolean): Unit = {
      outputStream.writeBoolean(value)
    }

    override def getBytes(): Int = 1
  }
}