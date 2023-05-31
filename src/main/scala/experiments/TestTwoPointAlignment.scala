package experiments

import io.RunWriter
import scalismo.common.{PointId, UnstructuredPointsDomain1D}
import scalismo.geometry.{Point1D, _3D}
import scalismo.io.MeshIO
import scalismo.mesh.TriangleMesh
import scalismo.numerics.PivotedCholesky.RelativeTolerance
import scalismo.statisticalmodel.DiscreteLowRankGaussianProcess
import scalismo.statisticalmodel.dataset.DataCollection
import scalismo.utils.Random

import java.io.File

object TestTwoPointAlignment {
  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    implicit val rnd = scalismo.utils.Random(98237465L)

    val femurFolder = "./data/femurData/registered/"

    val origShapes = new File(femurFolder).listFiles().map(MeshIO.readMesh).map(_.get).toIndexedSeq.zipWithIndex

    val writer = RunWriter(new File("./data/alignment/synthLengthTest/", "twopRes.txt"))

    IndexedSeq(false, true).foreach(aligned => {
      val label = if(aligned) "Align" else "Orig"
      IndexedSeq(1.0,10.0).foreach(sigma => {
        val res = origShapes.map{ case (target,i) =>
          val shapes = origShapes.filter(_._2 != i).map(_._1)
          check2Point1d(shapes,target, aligned, sigma)
        }
        val average = res.foldLeft((0.0,0.0))((t1,t2) => (t1._1+t2._1/res.length, t1._2+t2._2/res.length))
        writer.writeSingle(average._1,s"mean$label${sigma.toInt}")
        writer.writeSingle(average._2,s"var$label${sigma.toInt}")
      })
    })
    writer.close()
  }

  def check2Point1d(meshes: IndexedSeq[TriangleMesh[_3D]], fulltarget: TriangleMesh[_3D], aligned: Boolean, sigma: Double)(implicit rnd:Random): (Double, Double) = {
    //compare original model to translation aligned model. only looks at the variance of the prediction point
    val p1 = PointId(1840)
    val p2 = PointId(3130)
    val ps = meshes.map(mesh => (mesh.pointSet.point(p1) - mesh.pointSet.point(p2)).norm)
    val ts = (fulltarget.pointSet.point(p1) - fulltarget.pointSet.point(p2)).norm
    val shapes = if (aligned) {
      ps.map(d => UnstructuredPointsDomain1D(IndexedSeq(Point1D(0.0), Point1D(d))))
    } else ps.map(d => UnstructuredPointsDomain1D(IndexedSeq(Point1D(-d / 2), Point1D(d / 2))))
    val target = if (aligned) {
      UnstructuredPointsDomain1D(IndexedSeq(Point1D(0.0), Point1D(ts)))
    } else UnstructuredPointsDomain1D(IndexedSeq(Point1D(-ts / 2), Point1D(ts / 2)))

    val dgp = DiscreteLowRankGaussianProcess.createUsingPCA(DataCollection.fromUnstructuredPointsDomainSequence(shapes.head, shapes), RelativeTolerance(0.011))
    val obs = target.pointSet.point(PointId(0))
    val obsn = obs.map(d => d + rnd.scalaRandom.nextGaussian() * sigma)
    val pdgp = dgp.posterior(IndexedSeq((PointId(0), obsn - dgp.domain.pointSet.point(PointId(0)))), sigma * sigma)
    val pobsn = ((pdgp.domain.pointSet.point(PointId(0)) + pdgp.mean.data(0)) - (pdgp.domain.pointSet.point(PointId(1)) + pdgp.mean.data(1))).norm
    val pvar = breeze.linalg.sum(pdgp.variance)
    (pobsn, pvar)
  }

  def findMaxDistancePair(mesh: TriangleMesh[_3D]): (PointId, PointId) = {
    var pair = (PointId(0),PointId(1))
    var dist = 0.0
    (0 until mesh.pointSet.numberOfPoints).foreach(id => {
      val pid = PointId(id)
      val p = mesh.pointSet.point(pid)
      (id+1 until mesh.pointSet.numberOfPoints).foreach(aid => {
        val apid = PointId(aid)
        val d = (p - mesh.pointSet.point(apid)).norm
        if (d > dist){
          dist = d
          pair = (pid,apid)
        }
      })
    })
    pair
  }
}
