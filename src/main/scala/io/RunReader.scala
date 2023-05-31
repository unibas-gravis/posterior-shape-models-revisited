package io

import breeze.linalg.DenseVector
import scalismo.geometry.{EuclideanVector, _3D}
import utility.{Control, Tchange}

import java.io.File

class RunReader(file: File) {
  def getValue(name: String): Double = readValue(name).get
  def readValue(name: String): Option[Double] = {
    find(name).map(_.dropRight(1).toDouble)
  }
  def getVectorEV(name: String): IndexedSeq[EuclideanVector[_3D]] = Tchange.seqToVec3D(getVector(name))
  def getVectorDV(name: String, vectorLength: Int): IndexedSeq[DenseVector[Double]] = Tchange.seqToDenseVec(getVector(name),vectorLength)
  def getVector(name: String): IndexedSeq[Double] = readVector(name).get
  def readVector(name: String): Option[IndexedSeq[Double]] = {
    find(name).map(line => line.drop(1).dropRight(2).split(",").toIndexedSeq.map(_.toDouble))
  }
  def getMatrixEV(name: String): IndexedSeq[IndexedSeq[EuclideanVector[_3D]]] = getMatrix(name).map(Tchange.seqToVec3D)
  def getMatrix(name: String): IndexedSeq[IndexedSeq[Double]] = readMatrix(name).get
  def readMatrix(name: String): Option[IndexedSeq[IndexedSeq[Double]]] = {
    find(name).map(line => line.drop(1).dropRight(2).split(";").toIndexedSeq.map(_.split(",").toIndexedSeq.map(_.toDouble)))
  }
  def find(name: String): Option[String] = {
    Control.withResources(scala.io.Source.fromFile(file))(source => {
      source.getLines().find(line => line.startsWith(s"$name=")).map(_.drop(s"$name=".length))
    })
  }
  def getVariables(): Set[String] = {
    Control.withResources(scala.io.Source.fromFile(file))(source => {
      source.getLines().toSet
    })
  }
}

object RunReader {
  def apply(file: String): RunReader = RunReader(new File(file))
  def apply(file: File): RunReader = new RunReader(file)
}
