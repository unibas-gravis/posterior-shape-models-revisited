package experiments

import fitting.FitScripts.Sample
import fitting.IcpFit.IcpSettings
import fitting._
import io.RunWriter
import norms.L2norm
import scalismo.{ModelRecenter, ModelUtils}
import scalismo.common.{PointId, ScalarMeshField}
import scalismo.geometry._
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.mesh.{TriangleMesh, TriangleMesh3D}
import scalismo.statisticalmodel.{PointDistributionModel, StatisticalMeshModel}
import scalismo.transformations.{Rotation3D, Translation3D, TranslationAfterRotation3D}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random
import utility.{MathHelp, MeshUtils, Tchange}

import java.io.File


/**
 * The main experiment of femur loo reconstruction. if you set mcmc=false you get the experiment in the appendix where
 * posterior shape models are used to perform reconstructions with fixed correspondence instead of the full fledged
 * optimization over pose and correspondence
 */
object TestPartialAlignment {

  var ui: ScalismoUI = null
  val useParallel = true
  val showRun = false
  val vizMeshFolder = "./data/alignment/scalarMeshes/"

  val runSettings = {
    val rotation = true
    val mcmc = true

    val rposet = 2.0
    val rposer = 0.02
    val iidNoise = 1.0
    val tildeLdmSigma = 1.0
    val tScale = 2.0
    val rScale = 6.0
    val isoShapeScale = 0.1
    val shapeScale = 0.3
    val burnin = 500
    val numSamples = 5000
    val subsampling = 4
    val numRealignmentsMCMC = 10
    val toUseLengthBasedRatio = true
    val icpIterations = 150
    val numOfLoo = 100 //100 there are 47 femurs. so any value > leads to all samples used.
    val sigma = 1.0
    val ratios = IndexedSeq(0.1,0.2,0.3,0.45,0.6,0.8)

    RecoSettings(rotation, mcmc, rposet, rposer, iidNoise, tildeLdmSigma, tScale, rScale, isoShapeScale, shapeScale, burnin, numSamples, subsampling, numRealignmentsMCMC, toUseLengthBasedRatio, icpIterations, numOfLoo, sigma, ratios)
  }


  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    //    ui = ScalismoUI()
    implicit val rnd = scalismo.utils.Random(98237465L)

    assert(!(useParallel&&showRun), "should not use visualization when using parallel run")
    assert(!showRun || ui!=null, "should activate ui to use visualization")

    val femurFolder = "./data/femurData/registered/"

    assert((runSettings.rotation && runSettings.mcmc) || !runSettings.mcmc, "if mcmc is used rotation should be enabled")

    val origShapes = {
      println("loading data")
      val unaligned = new File(femurFolder).listFiles().map(MeshIO.readMesh).map(_.get).toIndexedSeq
      val aligned = ModelUtils.alignShapesGpa(unaligned)._1.map(_._1).zipWithIndex
      aligned
    }
    println("loaded and aligned data")

    //val model = StatisticalMeshModel.createUsingPCA(DataCollection.fromTriangleMesh3DSequence(origShapes.head._1, origShapes.map(_._1))).get
    //    val means = (0 to 10).map(_ => Tchange.getMean(model.sample()))
    //    {
    //      val model = StatisticalMeshModel.createUsingPCA(DataCollection.fromTriangleMesh3DSequence(origShapes.head._1, origShapes.dropRight(1).map(_._1))).get
    //      val projLast = model.project(origShapes.last._1)
    //      val l2 = L2norm[_3D]()
    //      val diff = Tchange.getDef(projLast.pointSet, origShapes.last._1)
    //      val n = l2.norm2Vector(diff)
    //      println(s"norm2 ${n} leading to avgNorm2 of ${n/model.mean.pointSet.numberOfPoints}")
    //      println("test ended")
    //    }

    val methodDesc = if (runSettings.mcmc) "MCMC" else {if (runSettings.rotation) "Arot" else "At"}
    val (writer,fileForWriter) = RunWriter(s"./data/femur/pythonScript/partialRes${methodDesc}.py", pythonStyle = true)
    val (writerLong,_) = RunWriter(s"./data/femur/pythonScript/partialRes${methodDesc}Verbose.txt", pythonStyle = true)
    runSettings.print(writer, Option(writerLong), Option(fileForWriter), Option("TestPartialAlignment run to test analytical posterior, MCMC, Icp performance"))

    val dirEncPids = (PointId(1840), PointId(3130))

    val ldms = {
      val model = StatisticalModelIO.readStatisticalMeshModel(new File("./data/femurData/ssm.h5")).get
      val ps = IndexedSeq(
        Point3D(115.29815673828125,146.07376098632812,-101.59976959228516),
        Point3D(103.93435668945312,179.8994903564453,-100.22386932373047),
        Point3D(42.72134017944336,180.55906677246094,-117.84776306152344),
        Point3D(55.907447814941406,131.927978515625,-127.22222900390625)
      )
      ps.map(model.mean.pointSet.findClosestPoint).map(_.id)
      //      model.mean.pointSet.pointIds.toIndexedSeq
    }

    //    (1 to 30).map(i => { //ldms are contained with ratio >= 0.16
    //      val ratio = i/100.0
    //      val (_,obs,_) = MeshUtils.getPartialMeshVector(origShapes.head._1, origShapes.head._1.pointSet.point(dirEncPids._1) - origShapes.head._1.pointSet.point(dirEncPids._2), ratio)
    //      println(s"ldms in obs for ratio ${(ratio*100).toInt}%: ${ldms.forall(pid => obs.contains(pid))}")
    //    })
    //    println("done")
    //    {
    //      val ratios = IndexedSeq(0.1, 0.2, 0.3, 0.5, 0.65, 0.8)
    //      ratios.zipWithIndex /*.par*/ .foreach(ratio => { //TODO the landmark fields are aligned to landmarks -> not observed for 0.05 and 0.15 -> increased knowledge
    //        IndexedSeq(1.0 /*,10.0*/).foreach(sigma => {
    //            val (_, obs, _) = MeshUtils.getPartialMeshVector(origShapes.head._1, origShapes.head._1.pointSet.point(dirEncPids._1) - origShapes.head._1.pointSet.point(dirEncPids._2), ratio._1)
    //            ui.show(UnstructuredPointsDomain(obs.map(origShapes.head._1.pointSet.point)), s"ratio ${ratio._1}").opacity = 0.0
    //        })
    //      })
    //      Thread.sleep(1000000)
    //    }



    //    ratios.foreach(ratio => {
    //      val target = origShapes.head._1
    //      val (_, obs, _) = MeshUtils.getPartialMeshVector(target, target.pointSet.point(dirEncPids._1) - target.pointSet.point(dirEncPids._2), ratio, toUseLengthBasedRatio)
    //      println(s"ratio ${ratio} gives ${obs.length} out of ${target.pointSet.numberOfPoints}")
    //    })
    //    println("done")
    fieldNameSeq().foreach(s => {
      writerLong.writeLiteral(s"${s}${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${runSettings.ratios.length}))")
      writerLong.writeLiteral(s"${s}Obs${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${runSettings.ratios.length}))")
      writerLong.writeLiteral(s"${s}Pre${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${runSettings.ratios.length}))")
    })

    runSettings.ratios.zipWithIndex.par.foreach(ratio => {
      val toIterate = if(useParallel) origShapes.par else origShapes
      val res = toIterate.take(runSettings.numOfLoo).map{ case (target,i) =>
        val (_,obs,_) = MeshUtils.getPartialMeshVector(target, target.pointSet.point(dirEncPids._1) - target.pointSet.point(dirEncPids._2), ratio._1, runSettings.toUseLengthBasedRatio)
        val shapes = origShapes.filter(_._2 != i).map(_._1)
        val res = checkPartial(shapes, target, ldms, obs, dirEncPids, runSettings)
        printOutDirect(writerLong, res, ratio._2, i, runSettings.sigma)
        res
      }.toIndexedSeq
      val average = res.reduce(_+_).map(_/res.length)
      printOutCollected(writer, average, ratio._2, runSettings.sigma)
      writer.writeCollectedSorted(ratio._1,ratio._2,s"obsRatio")
      println(s"fininshed for ratio ${(100*ratio._1).toInt}% obs")
    })
    writer.close()
    writerLong.close()
  }

  def fieldNameSeq(): Seq[String] = {
    Seq(
      "varOrig", "varLdmAlign", "varLdmLdmAlign", "varLdmTildeLdmAlign",  //variance
      "varFullAlign", "varApprxAlign", "varAlignOnline",
      "lossOrig", "lossLdmAlign", "lossLdmLdmAlign", "lossLdmTildeLdmAlign", //loss clp
      "lossFullAlign", "lossApprxAlign", "lossAlignOnline",
      "lossIcpOrig", "lossIcpLdmAlign", "lossIcpLdmLdmAlign", "lossIcpLdmTildeLdmAlign", //icp loss clp
      "lossIcpFullAlign", "lossIcpApprxAlign", "lossIcpAlignOnline",
      "lossCorrOrig", "lossCorrLdmAlign", "lossCorrLdmLdmAlign", "lossCorrLdmTildeLdmAlign", //loss correspondence
      "lossCorrFullAlign", "lossCorrApprxAlign", "lossCorrAlignOnline",
      "lossCorrIcpOrig", "lossCorrIcpLdmAlign", "lossCorrIcpLdmLdmAlign", "lossCorrIcpLdmTildeLdmAlign", //loss icp correspondence
      "lossCorrIcpFullAlign", "lossCorrIcpApprxAlign", "lossCorrIcpAlignOnline",
    )
  }
  def toPrint(results: CompPostResults): Seq[(CompItem, String)] = {
    Seq(
      (results.original, "varOrig"), (results.ldmAlign, "varLdmAlign"), //variance
      (results.ldmAlignLdm, "varLdmLdmAlign"), (results.ldmAlignTildeLdm, "varLdmTildeLdmAlign"),
      (results.fullAlign, "varFullAlign"), (results.apprxAlign, "varApprxAlign"), (results.onlineAlign, "varAlignOnline"),
      (results.originalLoss, "lossOrig"), (results.ldmAlignLoss, "lossLdmAlign"), //loss clp
      (results.ldmAlignLossLdm, "lossLdmLdmAlign"), (results.ldmAlignLossTildeLdm, "lossLdmTildeLdmAlign"),
      (results.fullAlignLoss, "lossFullAlign"), (results.apprxAlignLoss, "lossApprxAlign"), (results.onlineAlignLoss, "lossAlignOnline"),
      (results.icpOrigLoss, "lossIcpOrig"), (results.icpLdmLoss, "lossIcpLdmAlign"), //icp loss clp
      (results.icpLdmLossLdm, "lossIcpLdmLdmAlign"), (results.icpLdmLossTildeLdm, "lossIcpLdmTildeLdmAlign"),
      (results.icpAlignLoss, "lossIcpFullAlign"), (results.icpApprxLoss, "lossIcpApprxAlign"), (results.icpOnlineLoss, "lossIcpAlignOnline"),
      (results.originalLossCorr, "lossCorrOrig"), (results.ldmAlignLossCorr, "lossCorrLdmAlign"), //loss correspondence
      (results.ldmAlignLossCorrLdm, "lossCorrLdmLdmAlign"), (results.ldmAlignLossCorrTildeLdm, "lossCorrLdmTildeLdmAlign"),
      (results.fullAlignLossCorr, "lossCorrFullAlign"), (results.apprxAlignLossCorr, "lossCorrApprxAlign"), (results.onlineAlignLossCorr, "lossCorrAlignOnline"),
      (results.icpOrigLossCorr, "lossCorrIcpOrig"), (results.icpLdmLossCorr, "lossCorrIcpLdmAlign"), //loss icp correspondence
      (results.icpLdmLossCorrLdm, "lossCorrIcpLdmLdmAlign"), (results.icpLdmLossCorrTildeLdm, "lossCorrIcpLdmTildeLdmAlign"),
      (results.icpAlignLossCorr, "lossCorrIcpFullAlign"), (results.icpApprxLossCorr, "lossCorrIcpApprxAlign"), (results.icpOnlineLossCorr, "lossCorrIcpAlignOnline"),
    )
  }
  def printOutDirect(writer: RunWriter, results: CompPostResults, index1: Int, index2: Int, sigma: Double, postFix:String=""): Unit = {
    def write(item: CompItem, name: String, sigma: Int): Unit = {
      writer.writeLiteral(s"${name}${sigma}[${index1}, ${index2}] = ${item.all}")
      writer.writeLiteral(s"${name}Obs${sigma}[${index1}, ${index2}] = ${item.obs}")
      writer.writeLiteral(s"${name}Pre${sigma}[${index1}, ${index2}] = ${item.pre}")
    }
    toPrint(results).foreach(t => write(t._1, t._2+postFix, sigma.toInt))
  }
  def printOutCollected(writer: RunWriter, results: CompPostResults, index: Int, sigma: Double, postFix:String=""): Unit = {
    def write(item: CompItem, index: Int, name: String, sigma: Int): Unit = {
      writer.writeCollectedSorted(item.all, index, s"${name}${sigma}")
      writer.writeCollectedSorted(item.obs, index, s"${name}Obs${sigma}")
      writer.writeCollectedSorted(item.pre, index, s"${name}Pre${sigma}")
    }
    toPrint(results).foreach(t => write(t._1, index, t._2+postFix, sigma.toInt))
  }

  def checkPartial(origShapes: IndexedSeq[TriangleMesh[_3D]], fulltarget: TriangleMesh[_3D], ldms: IndexedSeq[PointId], obs: IndexedSeq[PointId], pidsDirEncode: (PointId, PointId), runSettings: RecoSettings)(implicit rnd:Random): CompPostResults = {
    val (shapes, mean) = {
      val res = ModelUtils.alignShapesGpa(origShapes)
      (res._1.map(_._1), res._2)
    }
    val (modelLdm, modelAll, model) = if (runSettings.rotation) {
      val modelLdm = ModelUtils.pcaModelGpa(shapes, ldms)
      val modelAll = ModelUtils.pcaModelGpa(shapes, obs)
      val model = ModelUtils.pcaModelGpa(shapes, mean.pointSet.pointIds.toIndexedSeq)
      (modelLdm, modelAll, model)
    } else {
      val modelLdm = ModelUtils.pcaAlignIdModel(shapes.head, shapes, ldms, false)
      val modelAll = ModelUtils.pcaAlignIdModel(shapes.head, shapes, obs, false)
      val model = ModelUtils.pcaAlignIdModel(shapes.head, shapes, mean.pointSet.pointIds.toIndexedSeq, false)
      (modelLdm, modelAll, model)
    }

    val res = if (runSettings.mcmc) {
      val modelTildeLdm = ModelUtils.pcaModelGpaTildeLdm(shapes, ldms, runSettings.tildeLdmSigma)
      compareToAlignedMCMCicp(model, modelLdm, modelTildeLdm, modelAll, obs, ldms, fulltarget, pidsDirEncode, runSettings)
    } else {
      compareToAlignedAnalytical(model, modelLdm, modelAll, obs, ldms, fulltarget, pidsDirEncode, runSettings)
    }
    res
  }

  case class CompItem(all: Double, obs: Double, pre: Double){
    def +(that: CompItem): CompItem = CompItem(this.all+that.all, this.obs+that.obs, this.pre+that.pre)
    def map(f:Double=>Double): CompItem = CompItem(f(all), f(obs), f(pre))
    override def toString: String = s"(${this.all}, ${this.obs}, ${this.pre})"
  }
  case class CompPostResults(original: CompItem, ldmAlign: CompItem, ldmAlignLdm: CompItem, ldmAlignTildeLdm: CompItem, fullAlign: CompItem, apprxAlign: CompItem, onlineAlign: CompItem,
                             originalLoss: CompItem, ldmAlignLoss: CompItem, ldmAlignLossLdm: CompItem, ldmAlignLossTildeLdm: CompItem, fullAlignLoss: CompItem, apprxAlignLoss: CompItem, onlineAlignLoss: CompItem,
                             icpOrigLoss: CompItem, icpLdmLoss: CompItem, icpLdmLossLdm: CompItem, icpLdmLossTildeLdm: CompItem, icpAlignLoss: CompItem, icpApprxLoss: CompItem, icpOnlineLoss: CompItem,
                             originalLossCorr: CompItem, ldmAlignLossCorr: CompItem,  ldmAlignLossCorrLdm: CompItem, ldmAlignLossCorrTildeLdm: CompItem, fullAlignLossCorr: CompItem, apprxAlignLossCorr: CompItem, onlineAlignLossCorr: CompItem,
                             icpOrigLossCorr: CompItem, icpLdmLossCorr: CompItem, icpLdmLossCorrLdm: CompItem, icpLdmLossCorrTildeLdm: CompItem, icpAlignLossCorr: CompItem, icpApprxLossCorr: CompItem, icpOnlineLossCorr: CompItem
                            ){
    def +(that: CompPostResults): CompPostResults = {
      CompPostResults(
        this.original+that.original, this.ldmAlign+that.ldmAlign, this.ldmAlignLdm+that.ldmAlignLdm, this.ldmAlignTildeLdm+that.ldmAlignTildeLdm, this.fullAlign+that.fullAlign, this.apprxAlign+that.apprxAlign, this.onlineAlign+that.onlineAlign,
        this.originalLoss+that.originalLoss, this.ldmAlignLoss+that.ldmAlignLoss, this.ldmAlignLossLdm+that.ldmAlignLossLdm, this.ldmAlignLossTildeLdm+that.ldmAlignLossTildeLdm, this.fullAlignLoss+that.fullAlignLoss, this.apprxAlignLoss+that.apprxAlignLoss, this.onlineAlignLoss+that.onlineAlignLoss,
        this.icpOrigLoss+that.icpOrigLoss, this.icpLdmLoss+that.icpLdmLoss, this.icpLdmLossLdm+that.icpLdmLossLdm, this.icpLdmLossTildeLdm+that.icpLdmLossTildeLdm, this.icpAlignLoss+that.icpAlignLoss, this.icpApprxLoss+that.icpApprxLoss, this.icpOnlineLoss+that.icpOnlineLoss,
        this.originalLossCorr+that.originalLossCorr, this.ldmAlignLossCorr+that.ldmAlignLossCorr, this.ldmAlignLossCorrLdm+that.ldmAlignLossCorrLdm, this.ldmAlignLossCorrTildeLdm+that.ldmAlignLossCorrTildeLdm, this.fullAlignLossCorr+that.fullAlignLossCorr, this.apprxAlignLossCorr+that.apprxAlignLossCorr, this.onlineAlignLossCorr+that.onlineAlignLossCorr,
        this.icpOrigLossCorr+that.icpOrigLossCorr, this.icpLdmLossCorr+that.icpLdmLossCorr, this.icpLdmLossCorrLdm+that.icpLdmLossCorrLdm, this.icpLdmLossCorrTildeLdm+that.icpLdmLossCorrTildeLdm, this.icpAlignLossCorr+that.icpAlignLossCorr, this.icpApprxLossCorr+that.icpApprxLossCorr, this.icpOnlineLossCorr+that.icpOnlineLossCorr
      )
    }
    def map(f:Double=>Double): CompPostResults = {
      CompPostResults(
        original.map(f), ldmAlign.map(f), ldmAlignLdm.map(f), ldmAlignTildeLdm.map(f), fullAlign.map(f), apprxAlign.map(f), onlineAlign.map(f),
        originalLoss.map(f), ldmAlignLoss.map(f), ldmAlignLossLdm.map(f), ldmAlignLossTildeLdm.map(f), fullAlignLoss.map(f), apprxAlignLoss.map(f), onlineAlignLoss.map(f),
        icpOrigLoss.map(f), icpLdmLoss.map(f), icpLdmLossLdm.map(f), icpLdmLossTildeLdm.map(f), icpAlignLoss.map(f), icpApprxLoss.map(f), icpOnlineLoss.map(f),
        originalLossCorr.map(f), ldmAlignLossCorr.map(f), ldmAlignLossCorrLdm.map(f), ldmAlignLossCorrTildeLdm.map(f), fullAlignLossCorr.map(f), apprxAlignLossCorr.map(f), onlineAlignLossCorr.map(f),
        icpOrigLossCorr.map(f), icpLdmLossCorr.map(f), icpLdmLossCorrLdm.map(f), icpLdmLossCorrTildeLdm.map(f), icpAlignLossCorr.map(f), icpApprxLossCorr.map(f), icpOnlineLossCorr.map(f)
      )
    }
  }

  /**
   * returns case class containing fitting information. loss and variance.
   * * these values are distinguished for domains: all/observed/predicted
   */
  def compareToAlignedAnalytical(model: StatisticalMeshModel, modelLdm: StatisticalMeshModel, modelAlign: StatisticalMeshModel, obs: IndexedSeq[PointId], ldms: IndexedSeq[PointId], target: TriangleMesh[_3D], pidsDirEncode: (PointId, PointId), runSettings: RecoSettings): CompPostResults = {
    require(!runSettings.rotation || (obs.length>=3 && ldms.length>=3), "need at least three observations for rotation alignment")
    val nobs = {
      val sobs = obs.map(_.id).toSet
      model.mean.pointSet.pointIds.filter(pid => !sobs.contains(pid.id)).toIndexedSeq
    }
    def toVar(ssm: StatisticalMeshModel): CompItem = {
      val vars = ssm.mean.pointSet.pointIds.map(pid => breeze.linalg.sum(breeze.linalg.diag(ssm.cov(pid,pid)))).toIndexedSeq
      CompItem(vars.sum/vars.length, obs.map(pid => vars(pid.id)).sum/obs.length, nobs.map(pid => vars(pid.id)).sum/nobs.length)
    }
    //val tdList = obs.map(pid => (pid, target.pointSet.point(pid)))
    def getTd(mesh:TriangleMesh[_3D]): IndexedSeq[(PointId, Point[_3D])] = obs.map(pid => (pid, mesh.pointSet.point(pid)))
    val origTarget = ModelUtils.alignShape(target, model.mean, Option(obs), runSettings.rotation)._1
    val origTd = getTd(origTarget) //also used for apprx model
    val ldmTarget = ModelUtils.alignShape(target, modelLdm.mean, Option(obs), runSettings.rotation)._1
    val ldmTd = getTd(ldmTarget)
    val ldmTargetLdm = ModelUtils.alignShape(target, modelLdm.mean, Option(ldms), runSettings.rotation)._1
    val ldmldmTd = getTd(ldmTargetLdm)
    val alignTarget = ModelUtils.alignShape(target, modelAlign.mean, Option(obs), runSettings.rotation)._1
    val fullTd = getTd(alignTarget)

    val modelPost = model.posterior(origTd, runSettings.sigma)
    val orig = toVar(modelPost)
    val modelLdmPost = modelLdm.posterior(ldmTd, runSettings.sigma)
    val ldmalign = toVar(modelLdmPost)
    val modelLdmLdmPost = modelLdm.posterior(ldmldmTd, runSettings.sigma)
    val ldmldmalign = toVar(modelLdmLdmPost)

    val alignPost = modelAlign.posterior(fullTd,runSettings.sigma)
    val align = toVar(alignPost)
    val axisd = (model.mean.pointSet.point(pidsDirEncode._1) - model.mean.pointSet.point(pidsDirEncode._2)).normalize
    val axis = MathHelp.listOrthogonalVectors(axisd)
    val modelApprx = {
      val rec = ModelRecenter.recenterSsm(model, obs)
      if (runSettings.rotation) ModelRecenter.rerotate(rec, obs, axis) else rec
    }
    val apprxPost = modelApprx.posterior(origTd,runSettings.sigma)
    val apprx = toVar(apprxPost)

    //    ui.show(modelPost.mean,"postOrigmean").opacity = 0.0
    //    ui.show(origTarget, "origTarget").opacity = 0.0
    //    ui.show(alignPost.mean,"postAlignmean").opacity = 0.0
    //    ui.show(alignTarget, "alignTarget").opacity = 0.0
    //    Thread.sleep(1000000)

    //    ui.show(ui.createGroup("modelLdm"), modelLdm, "modelldm")
    //    ui.show(ui.createGroup("model"), model, "model")
    //    ui.show(ui.createGroup("modelFull"), modelAlign, "modelFull")
    //    ui.show(ui.createGroup("modelApprx"), modelApprx, "modelAprrx")
    //    Thread.sleep(1000000)


    //    ui.show(ui.createGroup("modelLdm"), modelLdmPost, "modelldm")
    //    ui.show(ui.createGroup("model"), modelPost, "model")
    //    ui.show(ui.createGroup("modelFull"), alignPost, "modelFull")
    //    ui.show(ui.createGroup("modelApprx"), apprxPost, "modelAprrx")
    //    Thread.sleep(10000000)

    //    ui.show(alignTarget, "aligntarget")
    //    ui.show(alignPost.mean, "alignPostMean")
    //    Thread.sleep(1000000000)

    //    { //sanity check
    //      val getTemp = (d: Double) => {
    //        val prealign = RegressionHelp.augmentPosteriorPreFix(model, getTd(model.mean), 1.0, d, withRotation = false)
    //        val postalign = ModelRecenter.recenter(prealign,obs)
    //        toVar(postalign)
    //      }
    //      val postalignvar = getTemp(100.0)
    //      val postalignvarHigh = getTemp(1000000000.0)
    //      println("test done")
    //    }

    val l2 = L2norm[_3D]()
    def loss(residual: IndexedSeq[EuclideanVector[_3D]]): CompItem =
      CompItem(
        l2.norm2Vector(residual)/residual.length,
        l2.norm2Vector(obs.map(pid => residual(pid.id)))/obs.length,
        l2.norm2Vector(nobs.map(pid => residual(pid.id)))/nobs.length
      )
    val zCompitem = CompItem(0.0,0.0,0.0)
    val res = CompPostResults(
      orig, ldmalign, ldmldmalign, zCompitem, align, apprx, zCompitem, //variance
      zCompitem,zCompitem,zCompitem,zCompitem,zCompitem,zCompitem,zCompitem, //mcmc losses based on clp
      zCompitem,zCompitem,zCompitem,zCompitem,zCompitem,zCompitem,zCompitem, //icp losses clp
      loss(Tchange.getDef(origTarget.pointSet, modelPost.mean)),//mean l2 loss with correct correspondence
      loss(Tchange.getDef(ldmTarget.pointSet, modelLdmPost.mean)),
      loss(Tchange.getDef(ldmTargetLdm.pointSet, modelLdmLdmPost.mean)),zCompitem, //empty entry for noisy ldm
      loss(Tchange.getDef(alignTarget.pointSet, alignPost.mean)),
      loss(Tchange.getDef(origTarget.pointSet, apprxPost.mean)), zCompitem,
      zCompitem,zCompitem,zCompitem,zCompitem,zCompitem,zCompitem,zCompitem,//again icp with corr
    )
    res
  }

  /**
   * returns case class containing fitting results. loss and variance for mcmc and just loss for icp.
   * these values are distinguished for domains: all/observed/predicted
   */
  def compareToAlignedMCMCicp(modelssm: StatisticalMeshModel, modelLdmssm: StatisticalMeshModel, modelTildeLdmssm: StatisticalMeshModel, modelAlignssm: StatisticalMeshModel, obs: IndexedSeq[PointId], ldms: IndexedSeq[PointId], target: TriangleMesh[_3D], pidsDirEncode: (PointId, PointId), runSettings: RecoSettings)(implicit rnd: Random): CompPostResults = {
    require(obs.length>=3 && ldms.length>=3, "need at least three observations for rotation alignment")
    val model = ModelUtils.ssmToPdm(modelssm)
    val modelLdm = ModelUtils.ssmToPdm(modelLdmssm)
    val modelTildeLdm = ModelUtils.ssmToPdm(modelTildeLdmssm)
    val modelAlign = ModelUtils.ssmToPdm(modelAlignssm)
    val axisd = (model.mean.pointSet.point(pidsDirEncode._1) - model.mean.pointSet.point(pidsDirEncode._2)).normalize
    val nobs = {
      val sobs = obs.map(_.id).toSet
      model.mean.pointSet.pointIds.filter(pid => !sobs.contains(pid.id)).toIndexedSeq
    }

    //val tdList = obs.map(pid => (pid, target.pointSet.point(pid)))
    val rndPose = (
      EuclideanVector3D(rnd.scalaRandom.nextGaussian()*runSettings.rposet,rnd.scalaRandom.nextGaussian()*runSettings.rposet,rnd.scalaRandom.nextGaussian()*runSettings.rposet),
      rnd.scalaRandom.nextGaussian()*runSettings.rposer,rnd.scalaRandom.nextGaussian()*runSettings.rposer,rnd.scalaRandom.nextGaussian()*runSettings.rposer
    )
    val origTarget = getNoisyTarget(ModelUtils.alignShape(target, model.mean, Option(obs), true)._1, obs, rndPose) //precise x target
    val ldmTarget = getNoisyTarget(ModelUtils.alignShape(target, modelLdm.mean, Option(obs), true)._1, obs, rndPose) //precise x target
    val ldmTargetLdm = getNoisyTarget(ModelUtils.alignShape(target, modelLdm.mean, Option(ldms), true)._1, obs, rndPose) //precise ldm target
    val ldmTargetTildeLdm = { //stdv tildeLdmSigma for landmark noise. still on surface                        //tilde ldm target
      val ldmsObs = ldms.map(pid => (target.pointSet.point(pid), modelLdm.mean.pointSet.point(pid)))
      val tildeLdms = ldmsObs.map(t => {
        val rndv = EuclideanVector3D(rnd.scalaRandom.nextGaussian(),rnd.scalaRandom.nextGaussian(),rnd.scalaRandom.nextGaussian())
        (target.operations.closestPointOnSurface(t._1 + rndv * runSettings.tildeLdmSigma).point, t._2)
      })
      val transform = ModelUtils.alignLdmsToLdms(tildeLdms, true)._2
      getNoisyTarget(target.transform(transform), obs, rndPose) //tilde ldm target
    }
    val alignTarget = getNoisyTarget(ModelUtils.alignShape(target, modelAlign.mean, Option(obs), true)._1, obs, rndPose) //precise x target

    val obsset = obs.map(_.id).toSet

    //creation of the apprxmodel with gt correspondence
    val axis = MathHelp.listOrthogonalVectors(axisd)
    val modelApprxssm = {
      val rec = ModelRecenter.recenterSsm(modelssm, obs)
      ModelRecenter.rerotate(rec, obs, axis)
    }
    val modelApprx = ModelUtils.ssmToPdm(modelApprxssm)

    //starting the fitting processes
    val zCompitem = (CompItem(0.0,0.0,0.0), CompItem(0.0,0.0,0.0), CompItem(0.0,0.0,0.0)) //no run object
    //precise omega model, precise x target
    val (orig, origl, origlc) = if (runSettings.numSamples > 0) {
      val (partial, rotpoint) = prepareTarget(origTarget, obsset)
      val modelChain = getChain(model, partial, rotpoint, axisd, runSettings)
      handleChain(model, modelChain, origTarget, obs, nobs, "original")
    } else zCompitem
    //precise ldm model, precise x target
    val (ldmalign, ldmalignl, ldmalignlc) = if (runSettings.numSamples > 0) {
      val modelLdmChain = {
        val (partial, _) = prepareTarget(ldmTarget, obsset)
        val rotpoint = Tchange.getMean(ldms.map(modelLdm.mean.pointSet.point))
        getChain(modelLdm, partial, rotpoint, axisd, runSettings)
      }
      handleChain(modelLdm, modelLdmChain, ldmTarget, obs, nobs, "ldm")
    } else zCompitem
    //precise ldm model, precise ldm target
    val (ldmalignldm, ldmalignldml, ldmalignldmlc) = if (runSettings.numSamples > 0) {
      val modelLdmChain = {
        val (partial, _) = prepareTarget(ldmTargetLdm, obsset)
        val rotpoint = Tchange.getMean(ldms.map(modelLdm.mean.pointSet.point))
        getChain(modelLdm, partial, rotpoint, axisd, runSettings)
      }
      handleChain(modelLdm, modelLdmChain, ldmTarget, obs, nobs, "ldmldm")
    } else zCompitem
    //tilde ldm model, tilde ldm target
    val (ldmaligntildeldm, ldmaligntildeldml, ldmaligntildeldmlc) = if (runSettings.numSamples > 0) {
      val modelLdmChain = {
        val (partial, _) = prepareTarget(ldmTargetTildeLdm, obsset)
        val rotpoint = Tchange.getMean(ldms.map(modelTildeLdm.mean.pointSet.point))
        getChain(modelTildeLdm, partial, rotpoint, axisd, runSettings)
      }
      handleChain(modelTildeLdm, modelLdmChain, ldmTarget, obs, nobs, "ldmtildeldm")
    } else zCompitem
    //precise x model, precise x target
    val (align, alignl, alignlc) = if (runSettings.numSamples > 0) {
      val (partial, rotpoint) = prepareTarget(alignTarget, obsset)
      val modelAlignChain = getChain(modelAlign, partial, rotpoint, axisd, runSettings)
      handleChain(modelAlign, modelAlignChain, alignTarget, obs, nobs, "aligned")
    } else zCompitem
    //my x model, precise x target
    val (apprx, apprxl, apprxlc) = if (runSettings.numSamples > 0) {
      val (partial, rotpoint) = prepareTarget(origTarget, obsset)
      val modelApprxChain = getChain(modelApprx, partial, rotpoint, axisd, runSettings)
      handleChain(modelApprx, modelApprxChain, origTarget, obs, nobs, "apprx")
    } else zCompitem
    //my tilde x model, precise x target
    val (online, onlinel, onlinelc) = if (runSettings.numSamples > 0) {
      val (partial, rotpoint) = prepareTarget(origTarget, obsset)
      val onlineChain = FitScripts.fitPartialWithAlignment(model, partial, axisd, rotpoint, iidNoiseSigma = runSettings.iidNoise, burnin = runSettings.burnin, numSamples = runSettings.numSamples, subsampling = runSettings.subsampling, tScale = runSettings.tScale, rScale = runSettings.rScale, shapeScale = runSettings.shapeScale, realignments = runSettings.numRealignmentsMCMC, samplingStrategy = TargetSampling(), correspondenceStrategy = TargetSamplingUnique())
      handleChainMeshes(model, onlineChain, origTarget, obs, nobs, "iterApprx")
    } else zCompitem

    if (showRun) {
      println("stopping for visualization. 'go' for continue")
      while (scala.io.StdIn.readLine() != "go") {
        println("sleeping 1 minute")
        Thread.sleep(60000)
      }
    }

    val icpSettings = IcpSettings[SamplingStrategy](BidirectionalSamplingFromTarget(), true, runSettings.icpIterations, sigma2 = runSettings.sigma*runSettings.sigma)
    val (icp,icpit) = if (runSettings.icpIterations > 0) IcpFit.apply(model, prepareTarget(origTarget, obsset)._1, icpSettings) else (model.mean, 0)
    val (icpAlign,icpAlignit) = if (runSettings.icpIterations > 0) IcpFit.apply(modelAlign, prepareTarget(alignTarget, obsset)._1, icpSettings) else (modelAlign.mean, 0)
    val (icpLdm,icpLdmit) = if (runSettings.icpIterations > 0) IcpFit.apply(modelLdm, prepareTarget(ldmTarget, obsset)._1, icpSettings) else (modelLdm.mean, 0)
    val (icpLdmLdm,icpLdmLdmit) = if (runSettings.icpIterations > 0) IcpFit.apply(modelLdm, prepareTarget(ldmTargetLdm, obsset)._1, icpSettings) else (modelLdm.mean, 0)
    val (icpLdmTildeLdm,icpLdmTildeLdmit) = if (runSettings.icpIterations > 0) IcpFit.apply(modelTildeLdm, prepareTarget(ldmTargetTildeLdm, obsset)._1, icpSettings) else (modelTildeLdm.mean, 0)
    val (icpApprx,icpApprxit) = if (runSettings.icpIterations > 0) IcpFit.apply(modelApprx, prepareTarget(origTarget, obsset)._1, icpSettings) else (modelApprx.mean, 0)
    val (icpOnline,icpOnlineit) = if (runSettings.icpIterations > 0) IcpFit.fitWithOnlineAlignment(model, prepareTarget(origTarget, obsset)._1, icpSettings, true, axisd, Option(TargetSamplingUnique())) else (model.mean, 0)

    def loss(residual: IndexedSeq[EuclideanVector[_3D]]): CompItem = lossFunction(residual, obs, nobs)

    val res = CompPostResults(
      original=orig, ldmAlign=ldmalign, ldmAlignLdm=ldmalignldm, ldmAlignTildeLdm=ldmaligntildeldm, fullAlign=align, apprxAlign=apprx, onlineAlign=online, //variance
      originalLoss=origl, ldmAlignLoss=ldmalignl, ldmAlignLossLdm=ldmalignldml, ldmAlignLossTildeLdm=ldmaligntildeldml, fullAlignLoss=alignl, apprxAlignLoss=apprxl, onlineAlignLoss=onlinel, //mean l2 loss with clp
      icpOrigLoss=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(origTarget, icp)),
      icpLdmLoss=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(ldmTarget, icpLdm)),
      icpLdmLossLdm=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(ldmTargetLdm, icpLdmLdm)),
      icpLdmLossTildeLdm=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(ldmTargetTildeLdm, icpLdmTildeLdm)),
      icpAlignLoss=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(alignTarget, icpAlign)),
      icpApprxLoss=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(origTarget, icpApprx)),
      icpOnlineLoss=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(origTarget, icpOnline)),
      originalLossCorr=origlc, ldmAlignLossCorr=ldmalignlc, ldmAlignLossCorrLdm=ldmalignldmlc, ldmAlignLossCorrTildeLdm=ldmaligntildeldmlc, fullAlignLossCorr=alignlc, apprxAlignLossCorr=apprxlc, onlineAlignLossCorr=onlinelc, //mean l2 loss with correspondence
      icpOrigLossCorr=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(origTarget.pointSet, icp)),
      icpLdmLossCorr=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(ldmTarget.pointSet, icpLdm)),
      icpLdmLossCorrLdm=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(ldmTargetLdm.pointSet, icpLdmLdm)),
      icpLdmLossCorrTildeLdm=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(ldmTargetTildeLdm.pointSet, icpLdmTildeLdm)),
      icpAlignLossCorr=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(alignTarget.pointSet, icpAlign)),
      icpApprxLossCorr=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(origTarget.pointSet, icpApprx)),
      icpOnlineLossCorr=if (runSettings.icpIterations == 0) zCompitem._1 else loss(Tchange.getDef(origTarget.pointSet, icpOnline))
    )
    res
  }

  def getNoisyTarget(mesh: TriangleMesh[_3D], obs: IndexedSeq[PointId], rndPose: (EuclideanVector[_3D], Double, Double, Double)) = {
    val rotp = Tchange.getMean(mesh, Option(obs))
    val noisyPoseT = TranslationAfterRotation3D.apply(Translation3D(rndPose._1), Rotation3D(rndPose._2, rndPose._3, rndPose._4, rotp))
    val noisyAlignedTarget = mesh.transform(noisyPoseT)
    noisyAlignedTarget
  }

  def handleChain(model: PointDistributionModel[_3D, TriangleMesh], chain: IndexedSeq[(Sample, Double)], target: TriangleMesh[_3D], obs: IndexedSeq[PointId], nobs: IndexedSeq[PointId], name: String): (CompItem, CompItem, CompItem) =
    handleChainMeshes(model, chain.map(t => (t._1.calcInstance(model), t._2)), target, obs, nobs, name)

  def makePartial(mesh: TriangleMesh[_3D], obsset: Set[Int]) = mesh.operations.maskPoints(pid => obsset.contains(pid.id)).transformedMesh

  def prepareTarget(mesh: TriangleMesh[_3D], obsset: Set[Int]) = {
    val partial = makePartial(mesh, obsset)
    (partial, Tchange.getMean(partial.pointSet.pointSequence))
  }

  def getChain(model: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], rotPoint: Point[_3D], axisd:EuclideanVector[_3D], runSettings: RecoSettings)(implicit rnd: Random) = {
    FitScripts.fitPartial(model, target, axisd, rotPoint, iidNoiseSigma = runSettings.iidNoise, burnin = runSettings.burnin, numSamples = runSettings.numSamples, subsampling = runSettings.subsampling, tScale = runSettings.tScale, rScale = runSettings.rScale, shapeScale = runSettings.shapeScale, isoShapeScale = runSettings.isoShapeScale)
  }

  def lossFunction(residual: IndexedSeq[EuclideanVector[_3D]], obs: IndexedSeq[PointId], nobs: IndexedSeq[PointId]): CompItem = {
    val l2 = L2norm[_3D]()
    CompItem(
      l2.norm2Vector(residual) / residual.length,
      l2.norm2Vector(obs.map(pid => residual(pid.id))) / obs.length,
      l2.norm2Vector(nobs.map(pid => residual(pid.id))) / nobs.length
    )
  }

  def handleChainMeshes(model: PointDistributionModel[_3D, TriangleMesh], chain: IndexedSeq[(TriangleMesh[_3D], Double)], target: TriangleMesh[_3D], obs: IndexedSeq[PointId], nobs: IndexedSeq[PointId], name: String): (CompItem, CompItem, CompItem) = {
    //returns the MAP and the average point variance
    val best = chain.maxBy(_._2)._1
    val meshpoints = chain.map(c => c._1.pointSet.pointSequence)
    val mean = {
      val scaledPoints = meshpoints.map(points => points.map(p => p.toVector.map(_ * (1.0 / chain.length))))
      TriangleMesh3D(scaledPoints.transpose.map(_.reduce(_ + _).toPoint), model.mean.triangulation)
    }
    val variance = {
      meshpoints.foldLeft(IndexedSeq.fill(model.mean.pointSet.numberOfPoints)(0.0)) { case (vs, points) => {
        points.zip(mean.pointSet.pointSequence).map(t => {
          val v = (t._1 - t._2).toBreezeVector
          //breeze.linalg.sum(breeze.linalg.diag(v*v.t)) //extended notation
          v.t * v
        }).zip(vs).map(t => t._1 + t._2)
      }
      }
    }.map(d => d / meshpoints.length)
    val residual = Tchange.getDef(target, best) //TODO discuss if mean should be used here.
    val residualCorr = Tchange.getDef(target.pointSet, best)

    if (showRun) {
      val dv = EuclideanVector3D(0.0, 1.0, 0.4).map(_ * 100)
      val tf = name match {
        case "original" => Translation3D.apply(EuclideanVector3D.zero)
        case "ldm" => Translation3D.apply(dv)
        case "ldmldm" => Translation3D.apply(dv.map(_ * 2.0))
        case "ldmtildeLdm" => Translation3D.apply(dv.map(_ * 3.0))
        case "aligned" => Translation3D.apply(dv.map(_ * 4.0))
        case "apprx" => Translation3D.apply(dv.map(_ * 5.0))
        case "iterApprx" => Translation3D.apply(dv.map(_ * 6.0))

        case _ => Translation3D.apply(EuclideanVector3D.zero)
      }
      val varianceMesh = ScalarMeshField(mean.transform(tf), variance)
      val lossMesh = ScalarMeshField(mean.transform(tf), residual.map(_.norm))
      val lossCorrMesh = ScalarMeshField(mean.transform(tf), residualCorr.map(_.norm))
      Seq((varianceMesh, "Var"), (lossMesh, "Loss"), (lossCorrMesh, "LossCorr")).map { case (scalarMesh, postfix) => {
        MeshIO.writeScalarMeshField(scalarMesh, new File(vizMeshFolder, s"${name}${postfix}.vtk"))
        ui.show(scalarMesh, name + postfix)
      }
      }
    }

    (
      CompItem( //variance
        variance.sum / model.mean.pointSet.numberOfPoints,
        obs.map(pid => variance(pid.id)).sum / obs.length,
        nobs.map(pid => variance(pid.id)).sum / nobs.length
      ),
      lossFunction(residual, obs, nobs), //loss clp
      lossFunction(residualCorr, obs, nobs) //loss correspondence
    )
  }

  case class RecoSettings(rotation: Boolean, mcmc: Boolean, rposet: Double, rposer: Double, iidNoise: Double, tildeLdmSigma: Double, tScale: Double, rScale: Double, isoShapeScale: Double, shapeScale: Double, burnin: Int, numSamples: Int, subsampling: Int, numRealignmentsMCMC: Int, toUseLengthBasedRatio: Boolean, icpIterations: Int, numOfLoo: Int, sigma: Double, ratios: IndexedSeq[Double]) {
    def print(writer: RunWriter, writerLong: Option[RunWriter]=None, fileForWriter: Option[File]=None, msg: Option[String]=None): Unit = {
      if (writerLong.isDefined) {
        if (msg.isDefined) writerLong.get.writeLiteral(s"# ${msg.get}")
        writerLong.get.writeLiteral(s"# for infos on the hyperparameters look at the associated result file ${fileForWriter.getOrElse(new File("")).getName}")
      }

      writer.writeLiteral(s"rotation = ${if(rotation) 1 else 0};")
      writer.writeLiteral(s"mcmc = ${if(mcmc) 1 else 0};")
      writer.writeLiteral(s"iidNoise = $iidNoise;")
      writer.writeLiteral(s"tildeLdmSigma = $tildeLdmSigma;")
      writer.writeLiteral(s"tScale = $tScale;")
      writer.writeLiteral(s"rScale = $rScale;")
      writer.writeLiteral(s"shapeScale = $shapeScale;")
      writer.writeLiteral(s"isoShapeScale = $isoShapeScale;")
      writer.writeLiteral(s"rposet = $rposet;")
      writer.writeLiteral(s"rposer = $rposer;")
      writer.writeLiteral(s"burnin = $burnin;")
      writer.writeLiteral(s"numSamples = $numSamples;")
      writer.writeLiteral(s"subsampling = $subsampling;")
      writer.writeLiteral(s"numRealignmentsMCMC = $numRealignmentsMCMC;")
      writer.writeLiteral(s"toUseLengthBasedRatio = ${if(toUseLengthBasedRatio) 1 else 0};")
      writer.writeLiteral(s"icpIterations = $icpIterations;")
      writer.writeLiteral(s"numOfLoo = $numOfLoo;")
      writer.writeLiteral(s"sigma = ${sigma};")
    }
  }

}
