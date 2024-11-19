package hinge

import scalismo.ModelUtils
import scalismo.common.PointId
import scalismo.geometry.{Point3D, SquareMatrix, _3D}
import scalismo.mesh.{TriangleCell, TriangleList, TriangleMesh, TriangleMesh3D}
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.transformations.Rotation3D
import scalismo.utils.Random
import utility.{MathHelp, Tchange}

object Hinge {

  /**
   * returns a model of a 2d 'road' in 3d space where there is a hinge after left number of points.
   *
   * that hinge is rotated by a desired angle stdv and n meshes are created from that.
   *
   * standard is a flat surface encompassing a half rotation 180°
   *
   * @param angleDev amount of angle gaussian standard deviation
   * @param left the amount of points moving on the hinge
   * @param right the amount of constant points on the right side of the hinge
   * @param n the number of samples used to create a model
   * @return empirical model calculated
   */
  def getModel(angleMean: Double, angleDev: Double, left: Int=5, right: Int=10, n: Int=20, align: String="all")(implicit rnd:Random): StatisticalMeshModel = {
    val hinge = getShape(left, right)
    val shapes = variedShapes(hinge, angleMean, angleDev, n)
    align match {
      case "all" => ModelUtils.pcaModelGpa(shapes, hinge.mesh.pointSet.pointIds.toIndexedSeq)
      case "left" => ModelUtils.pcaModelGpa(shapes, hinge.leftIds)
      case "right" => ModelUtils.pcaModelGpa(shapes, hinge.rightIds)
      case _ =>  throw new IllegalArgumentException(s"the align argument ${align} is not valid")
    }


  }

  def variedShapes(hinge: HingeShape, angleMean: Double, angleDev: Double, n: Int)(implicit rnd:Random): IndexedSeq[TriangleMesh[_3D]] = {
    val (pa, pb) = (hinge.mesh.pointSet.point(hinge.hingeCenter.head), hinge.mesh.pointSet.point(hinge.hingeCenter.last))
    val rotationVector = (pa - pb).normalize
    val rotpoints = hinge.leftIds.map(hinge.mesh.pointSet.point)
    val stablepoints = hinge.rightIds.map(hinge.mesh.pointSet.point)
    val meshes = (0 until n).map(_ => {
      val amp = angleMean + angleDev * rnd.scalaRandom.nextGaussian()
      val rotMat = MathHelp.rotMatrixFromAxisAngle(rotationVector, amp)
      val rot = Rotation3D.apply(SquareMatrix[_3D](rotMat.data), Tchange.getMean(IndexedSeq(pa, pb)))
      val rotated = rotpoints.map(rot.f)
      TriangleMesh3D(rotated ++ stablepoints, hinge.mesh.triangulation)
    })
    meshes
  }

  /**
   * returns mesh consisting of left (also has hinge rotation center) + right square fields -> twice the amount of triangles.
   * point ids:
   * 2468...
   * 1357...
   */
  def getShape(left: Int, right: Int): HingeShape = {
    val total = left + right

    //the y axis is used to add some deviation from perfect plane. needed for some optimizations
    val points = (-left+1 to right).flatMap(x => Seq(Point3D(x, ((x*100)%43.0)/100.0, -1), Point3D(x, ((x*100)%59)/100.0, 1)))
    val cells = (0 until total).flatMap(i => Seq(
      TriangleCell(PointId(i*2), PointId(i*2+2), PointId(i*2+1)),
      TriangleCell(PointId(i*2+3), PointId(i*2+1), PointId(i*2+2))
    )).dropRight(1)

    HingeShape(TriangleMesh3D(points, TriangleList(cells)), left, right)
  }

  case class HingeShape(mesh: TriangleMesh[_3D], left: Int, right: Int){
    val leftIds: IndexedSeq[PointId] = mesh.pointSet.pointIds.take(left*2).toIndexedSeq
    val rightIds: IndexedSeq[PointId] = mesh.pointSet.pointIds.toIndexedSeq.takeRight(right*2)
    val hingeCenter: IndexedSeq[PointId] = leftIds.takeRight(2)
    lazy val angle: Double = { //assuming a valid HingeShape mesh
      val (ida, idb, idc) = (leftIds.head, hingeCenter.head, rightIds(rightIds.length - 2))
      val vecHinge = (mesh.pointSet.point(ida) - mesh.pointSet.point(idb)).normalize
      val vecBase = (mesh.pointSet.point(idc) - mesh.pointSet.point(idb)).normalize
      MathHelp.getInAngle(vecHinge, vecBase)
    }
    lazy val (leftLength, rightLength): (Double, Double) = { //assuming a valid HingeShape mesh
      val (ida, idb, idc) = (leftIds.head, hingeCenter.head, rightIds(rightIds.length - 2))
      val ll = (mesh.pointSet.point(ida) - mesh.pointSet.point(idb)).norm
      val rl = (mesh.pointSet.point(idc) - mesh.pointSet.point(idb)).norm
      (ll, rl)
    }
  }

}
