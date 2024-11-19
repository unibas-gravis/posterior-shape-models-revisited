package experiments

import experiments.TestPartialAlignment._
import io.RunWriter
import scalismo.ModelUtils
import scalismo.common.PointId
import scalismo.geometry._
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.ui.api.ScalismoUI
import utility.MeshUtils

import java.io.File

/**
 * performs loo reconstruction using the scaphoid dataset.
 */
object TestScaphoidPartialAlignment {

  var ui: ScalismoUI = null
  val useParallel = true
  val showRun = false
  val vizMeshFolder = "./data/alignment/scaphoids/"

  val runSettings = {
    val rotation = true
    val mcmc = true

    val rposet = 0.4
    val rposer = 0.02
    val iidNoise = 0.1
    val tildeLdmSigma = iidNoise
    val tScale = 1.0
    val rScale = 6.0
    val isoShapeScale = 0.2
    val shapeScale = 0.2
    val burnin = 500
    val numSamples = 7000
    val subsampling = 2
    val numRealignmentsMCMC = 5  //depending on what is observed, our naive estimation is exploring the space very broadly during burn-in. alternatively the initial pose guess should be worsened, as it provides the naive methods a good local minima to be stuck in
    val toUseLengthBasedRatio = true
    val icpIterations = 0 //this excludes icp methods. icp was never configured in this setup.
    val numOfLoo = 100 //there are 78 scaphoids.
    val sigma = 0.1
    val ratios = IndexedSeq(0.2, 0.3, 0.45, 0.6, 0.8)
    RecoSettings(rotation, mcmc, rposet, rposer, iidNoise, tildeLdmSigma, tScale, rScale, isoShapeScale, shapeScale, burnin, numSamples, subsampling, numRealignmentsMCMC, toUseLengthBasedRatio, icpIterations, numOfLoo, sigma, ratios)
  }


  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    //    ui = ScalismoUI()
    implicit val rnd = scalismo.utils.Random(98237465L)

    assert(!(useParallel && showRun), "should not use visualization when using parallel run")
    assert(!showRun || ui != null, "should activate ui to use visualization")

    val scaphoidFolder = "./data/scaphoidData/scaphoid_modelfiles/"

    assert((runSettings.rotation && runSettings.mcmc) || !runSettings.mcmc, "if mcmc is used rotation should be enabled")

    val origShapes = {
      println("loading data")
      val unaligned = new File(scaphoidFolder).listFiles().map(MeshIO.readMesh).map(_.get).toIndexedSeq
      val aligned = ModelUtils.alignShapesGpa(unaligned)._1.map(_._1).zipWithIndex
      aligned
    }
    println("loaded and aligned data")

    val methodDesc = if (runSettings.mcmc) "MCMC" else {if (runSettings.rotation) "Arot" else "At"}
    val (writer, fileForWriter) = RunWriter(s"./data/scaphoid/pythonScript/partialRes${methodDesc}.py", pythonStyle = true)
    val (writerLong, _) = RunWriter(s"./data/scaphoid/pythonScript/partialRes${methodDesc}Verbose.txt", pythonStyle = true)
    runSettings.print(writer, Option(writerLong), Option(fileForWriter), Option("Similar to TestPartialAlignment run to test analytical posterior, MCMC, Icp performance but for scaphoids"))

    val dirEncPids = (PointId(63), PointId(398)) //398 is observed

    val ldms = {
      val model = StatisticalModelIO.readStatisticalMeshModel(new File("./data/scaphoidData/empirical_model_R.h5")).get
      val ps = IndexedSeq( //these are not real landmarks. they are hard to identify and are unobservable for low ratios. they are only suitable for this abstracted setting
        Point3D(86.6594, 81.4422, 67.7105),
        Point3D(90.5881, 81.8316, 63.5196),
        Point3D(89.9648, 79.9865, 60.4572),
        Point3D(86.3238, 75.7147, 65.7089),
      )
      ps.map(model.mean.pointSet.findClosestPoint).map(_.id)
    }

    fieldNameSeq().foreach(s => {
      writerLong.writeLiteral(s"${s}${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${runSettings.ratios.length}))")
      writerLong.writeLiteral(s"${s}Obs${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${runSettings.ratios.length}))")
      writerLong.writeLiteral(s"${s}Pre${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${runSettings.ratios.length}))")
    })

    runSettings.ratios.zipWithIndex.foreach(ratio => {
      val toIterate = if (useParallel) origShapes.par else origShapes
      val res = toIterate.take(runSettings.numOfLoo).map { case (target, i) =>
        val (_, obs, _) = MeshUtils.getPartialMeshVector(target, target.pointSet.point(dirEncPids._1) - target.pointSet.point(dirEncPids._2), ratio._1, runSettings.toUseLengthBasedRatio)
        val shapes = origShapes.filter(_._2 != i).map(_._1)
        val res = checkPartial(shapes, target, ldms, obs, dirEncPids, runSettings)
        printOutDirect(writerLong, res, ratio._2, i, runSettings.sigma)
        res
      }.toIndexedSeq
      val average = res.reduce(_ + _).map(_ / res.length)
      printOutCollected(writer, average, ratio._2, runSettings.sigma)
      writer.writeCollectedSorted(ratio._1, ratio._2, s"obsRatio")
      println(s"fininshed for ratio ${(100 * ratio._1).toInt}% obs")
    })
    writer.close()
    writerLong.close()
  }

}
