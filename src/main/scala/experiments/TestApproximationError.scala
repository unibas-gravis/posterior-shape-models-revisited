package experiments

import hinge.Hinge
import hinge.Hinge.HingeShape
import io.RunWriter
import scalismo.geometry._3D
import scalismo.io.MeshIO
import scalismo.mesh.TriangleMesh
import scalismo.{ModelRecenter, ModelUtils}
import utility.Tchange
import utility.Tchange.MeanVar

import java.io.File

/**
 * this tests the approximation error for the hinge shape.
 */
object TestApproximationError {

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val meshPrinting = false  //because vtk cannot show such meshes needs to be saved to stl files for external viewing
    val ending = "stl" //for writing meshes
    val outputFolder = "./data/alignment/output/"
    implicit val rnd = scalismo.utils.Random(98237465L)

    //hinge parameters
    val (left, right) = (11,10)
    val n = 50 //number of model shapes
    val nt = 5000 //number of test shapes for error estimation

    //noise parameters
    val noisestdvs = (1 to 10).map(i => i.toDouble / 100 * 5)
    val noisemean = IndexedSeq.fill(noisestdvs.length)(math.Pi/2)

    val writer = RunWriter("./data/alignment/pythonScript/hingeRes.py", pythonStyle = true)._1
    val hinge = Hinge.getShape(left, right)
    val rotationAxis = {
      val cpoints = hinge.hingeCenter.map(hinge.mesh.pointSet.point)
      (cpoints.head - cpoints.last).normalize
    }

    noisemean.zip(noisestdvs).zipWithIndex.map(t => {
      val ((nmean, nstdv), i) = t
      val align = Hinge.getModel(nmean, nstdv, left, right, n, "left") //model calculated from shapes aligned on left points
      val apprx = {
        val orig = Hinge.getModel(nmean, nstdv, left, right, n, "all") //model calculated from shapes aligned on all points
        ModelRecenter.completeAlign(orig, hinge.leftIds, IndexedSeq(rotationAxis))
      }
      val trueShapes = {
        val unaligned = Hinge.variedShapes(hinge, nmean, nstdv, nt)
        ModelUtils.alignShapesGpa(unaligned, Option(hinge.leftIds))._1.map(_._1)
      }

      val ioshapes = 9
      val parAlign = {
        val alignshapes = (0 until nt).map(_=>align.sample())
        if (meshPrinting) alignshapes.zipWithIndex.take(ioshapes).foreach(t => MeshIO.writeMesh(t._1,new File(outputFolder, s"align${i}-${t._2}.${ending}")))
        recoverParameterDistribution(hinge, alignshapes)
      }
      val parApprx = {
        val apprxshapes = (0 until nt).map(_=>apprx.sample())
        if (meshPrinting) apprxshapes.zipWithIndex.take(ioshapes).foreach(t => MeshIO.writeMesh(t._1,new File(outputFolder, s"apprx${i}-${t._2}.${ending}")))
        recoverParameterDistribution(hinge, apprxshapes)
      }
      val parTrue = recoverParameterDistribution(hinge, trueShapes)
      if (meshPrinting) trueShapes.zipWithIndex.take(ioshapes).foreach(t => MeshIO.writeMesh(t._1,new File(outputFolder, s"true${i}-${t._2}.${ending}")))

      writer.writeCollectedSorted(parAlign.toIndexedSeq().zipWithIndex, i, "align")
      writer.writeCollectedSorted(parApprx.toIndexedSeq().zipWithIndex, i, "apprx")
      writer.writeCollectedSorted(parTrue.toIndexedSeq().zipWithIndex, i, "truep")
      println(s"tested for mean ${nmean} and stdv ${nstdv}")
      writer.writeCollectedSorted(nmean, i, "noisemean")
      writer.writeCollectedSorted(nstdv, i, "noisestdv")
      (parAlign, parApprx, parTrue)
    })
    println("done.\n")
    writer.close()
  }

  def recoverParameterDistribution(hinge: HingeShape, shapes: IndexedSeq[TriangleMesh[_3D]]): HingeParamaters = {
    //assuming that the characteristic of the hinge remains in the shape list.
    val hinges = shapes.map(s => HingeShape(s, hinge.left, hinge.right))
    val angles = hinges.map(_.angle)
    val leftl = hinges.map(_.leftLength)
    val rightl = hinges.map(_.rightLength)
    HingeParamaters(Tchange.getMeanVar(angles), Tchange.getMeanVar(leftl), Tchange.getMeanVar(rightl))
  }

  case class HingeParamaters(angle: MeanVar[Double], leftLength: MeanVar[Double], rightLength: MeanVar[Double]){
    def toIndexedSeq(): IndexedSeq[Double] = IndexedSeq(angle, leftLength, rightLength).flatMap(p => IndexedSeq(p.m, p.v))
  }
}
