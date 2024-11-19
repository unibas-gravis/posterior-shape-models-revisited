package io

import io.jsonHandler.Instructions
import scalismo.common.PointId
import scalismo.faces.io.MoMoIO
import scalismo.faces.momo.MoMo
import scalismo.geometry.{Point, Point3D, _3D}
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.{ModelHelpers, StatisticalMeshModel}

import java.io.File
import scala.io.Source

object LoadData {
  def getModel(file:File=new File("./data/femur_gp_model_100-cholesky.h5")): StatisticalMeshModel = { //femur_gp_model_100-components.h5
    val loading = StatisticalModelIO.readStatisticalMeshModel(file)
    loading.getOrElse(throw new IllegalArgumentException(s"couldn't load ${file.getName}"))
  }
  def getMoMo(file:File=new File("./data/faceModel/model2019_bfm.h5")): MoMo = {
    val loading = MoMoIO.read(file)
    loading.getOrElse(throw new IllegalArgumentException(s"couldn't load ${file.getName}"))
  }
  def getFaceShapeModel(file:File=new File("./data/faceModel/model2019_bfm.h5")): StatisticalMeshModel = {
    val loading = MoMoIO.read(file)
    val model = loading.getOrElse(throw new IllegalArgumentException(s"couldn't load ${file.getName}"))
    val lrgp = ModelHelpers.pointToVectorDLRGP(model.expressionModel.get.shape.gpModel,model.referenceMesh)
    StatisticalMeshModel(model.referenceMesh, lrgp)
  }
  def getPathologicalTarget(file:File=new File("./data/syntheticPathologicTargets/0over.vtk")): TriangleMesh[_3D] = {
    MeshIO.readMesh(file).getOrElse(throw new IllegalArgumentException(s"couldn't load ${file.getName}"))
  }
  def getFaceSegmentation(file:File=new File("./data/faceModel/syntheticContours.json")): Map[String,Set[Int]] = {
    val inst = Instructions(file)
//    val perm = Seq(("reye",Seq(0,1,13)),("leye",Seq(2,3,14)),("nose",Seq(4,5,6)), ("mouth",Seq(7,8,15)), ("rear",Seq(9,10)), ("lear",Seq(11,12)))
    val perm = Seq(("eye",Seq(0,1,13,2,3,14)),("nose",Seq(4,5,6)), ("mouth",Seq(7,8,15)), ("ear",Seq(9,10,11,12)))
    perm.toMap.mapValues(sint => sint.map(i=>inst.idGroups(i)).reduce(_.union(_)))
  }
  def getFaceSegmentationSelection(selection:IndexedSeq[IndexedSeq[Int]], file:File=new File("./data/faceModel/syntheticContours.json")): IndexedSeq[Set[Int]] = {
    val inst = Instructions(file)
    selection.map(sint => sint.map(i=>inst.idGroups(i)).reduce(_.union(_)))
  }
  /**
   * read csv file of mesh. returns ids and point3d. fail fast
   */
  def readCsv(f: File): IndexedSeq[(PointId, Point[_3D])] = {
    val source = Source.fromFile(f)
    val ps = {
      val full = source.mkString
      val lines = full.split("\n").drop(1) //drop header information
      lines.map(l => l.split(",")).map(l => (PointId(l(0).toInt), Point3D(l(1).toDouble, l(2).toDouble,l(3).toDouble)))
    }
    source.close()
    ps.toIndexedSeq
  }
}
