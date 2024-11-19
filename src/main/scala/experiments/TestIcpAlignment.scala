package experiments

import experiments.TestPartialAlignment._
import fitting.IcpFit.IcpSettings
import fitting.{BidirectionalSamplingFromTarget, IcpFit, SamplingStrategy, SamplingStrategyNormals}
import io.RunWriter
import scalismo.common.PointId
import scalismo.geometry.{EuclideanVector, EuclideanVector3D, Point3D, _3D}
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.{PointDistributionModel, StatisticalMeshModel}
import scalismo.utils.Random
import scalismo.{ModelRecenter, ModelUtils}
import utility.{MathHelp, MeshUtils, Tchange}

import java.io.File


/**
 * this test shows the icp performance with respect to iterations.
 */
object TestIcpAlignment {

  val settings = TestPartialAlignment.runSettings
  val icpIterationsToConv = 1 * runSettings.icpIterations
  val icpSamples = 2 * icpIterationsToConv
  val useParallel = TestPartialAlignment.useParallel
  val ratios = IndexedSeq(0.2)//,0.3)

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    implicit val rnd = scalismo.utils.Random(98237465L)

    assert(!(useParallel && showRun), "should not use visualization when using parallel run")
    assert(!showRun || ui != null, "should activate ui to use visualization")

    val femurFolder = "./data/femurData/registered/"

    val rotation = true
    val mcmc = true

    assert((rotation && mcmc) || !mcmc, "if mcmc is used rotation should be enabled")

    val origShapes = {
      println("loading data")
      val unaligned = new File(femurFolder).listFiles().map(MeshIO.readMesh).map(_.get).toIndexedSeq
      val aligned = ModelUtils.alignShapesGpa(unaligned)._1.map(_._1).zipWithIndex
      aligned
    }

    println("loaded and aligned data")

    val methodDesc = "ICPpath"
    val (writer,fileForWriter) = RunWriter(s"./data/femur/pythonScript/partialRes${methodDesc}.py", pythonStyle = true)
    val (writerLong,_) = RunWriter(s"./data/femur/pythonScript/partialRes${methodDesc}Verbose.txt", pythonStyle = true)
    runSettings.print(writer, Option(writerLong), Option(fileForWriter), Option("TestIcpAlignment run to test icp iteration performance"))

    val dirEncPids = (PointId(1840), PointId(3130))

    val ldms = {
      val model = StatisticalModelIO.readStatisticalMeshModel(new File("./data/femurData/ssm.h5")).get
      val ps = IndexedSeq( //TODO change to 'real used' ldms as these also add some additional information for some targets with under ratio 0.2
        //        Point3D(66.61833190917969,201.5722198486328,-98.29879760742188),
        //        Point3D(78.68072509765625,121.51073455810547,-120.26091766357422),
        Point3D(115.29815673828125, 146.07376098632812, -101.59976959228516),
        Point3D(103.93435668945312, 179.8994903564453, -100.22386932373047),
        Point3D(42.72134017944336, 180.55906677246094, -117.84776306152344),
        Point3D(55.907447814941406, 131.927978515625, -127.22222900390625)
      )
      ps.map(model.mean.pointSet.findClosestPoint).map(_.id)
      //      model.mean.pointSet.pointIds.toIndexedSeq
    }

    fieldNameSeq().foreach(s => {
      writerLong.writeLiteral(s"${s}${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${ratios.length}, ${icpIterationsToConv}))")
      writerLong.writeLiteral(s"${s}Obs${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${ratios.length}, ${icpIterationsToConv}))")
      writerLong.writeLiteral(s"${s}Pre${runSettings.sigma.toInt} = np.zeros((${math.min(runSettings.numOfLoo, origShapes.length)}, ${ratios.length}, ${icpIterationsToConv}))")
    })

    ratios.zipWithIndex.foreach(ratio => {
      val toIterate = if (useParallel) origShapes.par else origShapes
      val res = toIterate.take(runSettings.numOfLoo).map { case (target, i) =>
        val (_, obs, _) = MeshUtils.getPartialMeshVector(target, target.pointSet.point(dirEncPids._1) - target.pointSet.point(dirEncPids._2), ratio._1, runSettings.toUseLengthBasedRatio)
        val shapes = origShapes.filter(_._2 != i).map(_._1)
        val res = checkPartialWithIcp(shapes, target, ldms, obs, dirEncPids, runSettings)
        printOutDirect(writerLong, res, i, ratio._2, runSettings.sigma)
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

  def toPrint(results: CompPostIcpResults): Seq[(IndexedSeq[CompItem], String)] = {
    Seq(
      (results.icpOrigLoss, "lossIcpOrig"), (results.icpLdmLoss, "lossIcpLdmAlign"), //icp loss clp
      (results.icpLdmLossLdm, "lossIcpLdmLdmAlign"), (results.icpLdmLossTildeLdm, "lossIcpLdmTildeLdmAlign"),
      (results.icpAlignLoss, "lossIcpFullAlign"), (results.icpApprxLoss, "lossIcpApprxAlign"), (results.icpOnlineLoss, "lossIcpAlignOnline"),
      (results.icpOrigLinScaleLoss, "lossOrigLinScaleAlign"), (results.icpLinScaleLoss, "lossLinScaleAlign"),
      (results.icpOrigLossCorr, "lossCorrIcpOrig"), (results.icpLdmLossCorr, "lossCorrIcpLdmAlign"), //loss icp correspondence
      (results.icpLdmLossCorrLdm, "lossCorrIcpLdmLdmAlign"), (results.icpLdmLossCorrTildeLdm, "lossCorrIcpLdmTildeLdmAlign"),
      (results.icpAlignLossCorr, "lossCorrIcpFullAlign"), (results.icpApprxLossCorr, "lossCorrIcpApprxAlign"), (results.icpOnlineLossCorr, "lossCorrIcpAlignOnline"),
      (results.icpOrigLinScaleLossCorr, "lossCorrOrigLinScaleAlign"), (results.icpLinScaleLossCorr, "lossCorrLinScaleAlign"),
    )
  }

  def printOutDirect(writer: RunWriter, results: CompPostIcpResults, index1: Int, index2: Int, sigma: Double, postFix: String = ""): Unit = {
    def write(item: IndexedSeq[CompItem], name: String, sigma: Int): Unit = {
      item.zipWithIndex.foreach(it => {
        writer.writeLiteral(s"${name}${sigma}[${index1}, ${index2}, ${it._2}] = ${it._1.all}")
        writer.writeLiteral(s"${name}Obs${sigma}[${index1}, ${index2}, ${it._2}] = ${it._1.obs}")
        writer.writeLiteral(s"${name}Pre${sigma}[${index1}, ${index2}, ${it._2}] = ${it._1.pre}")
      })
    }
    toPrint(results).foreach(t => write(t._1, t._2 + postFix, sigma.toInt))
  }

  def printOutCollected(writer: RunWriter, results: CompPostIcpResults, index: Int, sigma: Double, postFix: String = ""): Unit = {
    def write(item: IndexedSeq[CompItem], index: Int, name: String, sigma: Int): Unit = {
      writer.writeCollectedSorted(item.map(_.all).zipWithIndex, index, s"${name}${sigma}")
      writer.writeCollectedSorted(item.map(_.obs).zipWithIndex, index, s"${name}Obs${sigma}")
      writer.writeCollectedSorted(item.map(_.pre).zipWithIndex, index, s"${name}Pre${sigma}")
    }
    toPrint(results).foreach(t => write(t._1, index, t._2 + postFix, sigma.toInt))
  }


  def checkPartialWithIcp(origShapes: IndexedSeq[TriangleMesh[_3D]], fulltarget: TriangleMesh[_3D], ldms: IndexedSeq[PointId], obs: IndexedSeq[PointId], pidsDirEncode: (PointId, PointId), runSettings: RecoSettings)(implicit rnd:Random): CompPostIcpResults = {
    val (shapes, mean) = {
      val res = ModelUtils.alignShapesGpa(origShapes)
      (res._1.map(_._1), res._2)
    }
    val (modelLdm, modelAll, model) = {
      val modelLdm = ModelUtils.pcaModelGpa(shapes, ldms)
      val modelAll = ModelUtils.pcaModelGpa(shapes, obs)
      val model = ModelUtils.pcaModelGpa(shapes, mean.pointSet.pointIds.toIndexedSeq)
      (modelLdm, modelAll, model)
    }
    val modelTildeLdm = ModelUtils.pcaModelGpaTildeLdm(shapes, ldms, runSettings.tildeLdmSigma)

    val res = compareToAligned(model, modelLdm, modelTildeLdm, modelAll, obs, ldms, fulltarget, pidsDirEncode, runSettings)
    res
  }

  def compareToAligned(modelssm: StatisticalMeshModel, modelLdmssm: StatisticalMeshModel, modelTildeLdmssm: StatisticalMeshModel, modelAlignssm: StatisticalMeshModel, obs: IndexedSeq[PointId], ldms: IndexedSeq[PointId], target: TriangleMesh[_3D], pidsDirEncode: (PointId, PointId), runSettings: RecoSettings)(implicit rnd:Random): CompPostIcpResults = {
    require(obs.length >= 3 && ldms.length >= 3, "need at least three observations for rotation alignment")
    val model = ModelUtils.ssmToPdm(modelssm)
    val modelLdm = ModelUtils.ssmToPdm(modelLdmssm)
    val modelTildeLdm = ModelUtils.ssmToPdm(modelTildeLdmssm)
    val modelAlign = ModelUtils.ssmToPdm(modelAlignssm)
    val nobs = {
      val sobs = obs.map(_.id).toSet
      model.mean.pointSet.pointIds.filter(pid => !sobs.contains(pid.id)).toIndexedSeq
    }
    val obsset = obs.map(_.id).toSet

    val rndPose = (
      EuclideanVector3D(rnd.scalaRandom.nextGaussian() * runSettings.rposet, rnd.scalaRandom.nextGaussian() * runSettings.rposet, rnd.scalaRandom.nextGaussian() * runSettings.rposet),
      rnd.scalaRandom.nextGaussian() * runSettings.rposer, rnd.scalaRandom.nextGaussian() * runSettings.rposer, rnd.scalaRandom.nextGaussian() * runSettings.rposer
    )
    val origTarget = getNoisyTarget(ModelUtils.alignShape(target, model.mean, Option(obs), true)._1, obs, rndPose) //precise x target
    val ldmTarget = getNoisyTarget(ModelUtils.alignShape(target, modelLdm.mean, Option(obs), true)._1, obs, rndPose) //precise x target
    val ldmTargetLdm = getNoisyTarget(ModelUtils.alignShape(target, modelLdm.mean, Option(ldms), true)._1, obs, rndPose) //precise ldm target
    val ldmTargetTildeLdm = { //stdv tildeLdmSigma for landmark noise. still on surface                        //tilde ldm target
      val ldmsObs = ldms.map(pid => (target.pointSet.point(pid), modelLdm.mean.pointSet.point(pid)))
      val tildeLdms = ldmsObs.map(t => {
        val rndv = EuclideanVector3D(rnd.scalaRandom.nextGaussian(), rnd.scalaRandom.nextGaussian(), rnd.scalaRandom.nextGaussian())
        (target.operations.closestPointOnSurface(t._1 + rndv * runSettings.tildeLdmSigma).point, t._2)
      })
      val transform = ModelUtils.alignLdmsToLdms(tildeLdms, true)._2
      getNoisyTarget(target.transform(transform), obs, rndPose) //tilde ldm target
    }
    val alignTarget = getNoisyTarget(ModelUtils.alignShape(target, modelAlign.mean, Option(obs), true)._1, obs, rndPose) //precise x target

    val axisd = (model.mean.pointSet.point(pidsDirEncode._1) - model.mean.pointSet.point(pidsDirEncode._2)).normalize
    val axis = MathHelp.listOrthogonalVectors(axisd)
    val modelApprxssm = {
      val rec = ModelRecenter.recenterSsm(modelssm, obs)
      ModelRecenter.rerotate(rec, obs, axis)
    }
    val modelApprx = ModelUtils.ssmToPdm(modelApprxssm)

    def loss(residual: IndexedSeq[EuclideanVector[_3D]]): CompItem = lossFunction(residual, obs, nobs)

    def evalIt(it: Iterator[(TriangleMesh[_3D], Boolean)]) = {
      val states = it.take(icpSamples).toIndexedSeq
      (states.map(t => loss(Tchange.getDef(origTarget, t._1))), states.map(t => loss(Tchange.getDef(origTarget.pointSet, t._1))))
    }

    //val icpSettings: IcpSettings[SamplingStrategyNormals] = IcpSettings(NormalSamplingSimpleExtension(TargetSamplingUnique()), true, icpIterationsToConv, sigma2 = sigma * sigma)
    val icpSettings: IcpSettings[SamplingStrategy] = IcpSettings(BidirectionalSamplingFromTarget(), true, icpIterationsToConv, sigma2 = runSettings.sigma * runSettings.sigma)
    val (icpOrigLoss, icpOrigLossCorr) = {
      val it = getIt(model, prepareTarget(origTarget, obsset)._1, icpSettings, useRealingment = false, lookAheadRegression = true, mainAxis = axisd)
      evalIt(it.map(t => (t._1, t._2)))
    }
    val (icpAlignLoss, icpAlignLossCorr) = {
      val it = getIt(modelAlign, prepareTarget(alignTarget, obsset)._1, icpSettings, useRealingment = false, lookAheadRegression = true, mainAxis = axisd)
      evalIt(it.map(t => (t._1, t._2)))
    }
//    val (icpLdmLoss, icpLdmLossCorr) = {
//      val it = getIt(modelLdm, prepareTarget(ldmTarget, obsset)._1, icpSettings, useRealingment = false, lookAheadRegression = true, mainAxis = axisd)
//      evalIt(it.map(t => (t._1, t._2)))
//    }
//    val (icpLdmLossLdm, icpLdmLossCorrLdm) = {
//      val it = getIt(modelLdm, prepareTarget(ldmTargetLdm, obsset)._1, icpSettings, useRealingment = false, lookAheadRegression = true, mainAxis = axisd)
//      evalIt(it.map(t => (t._1, t._2)))
//    }
//    val (icpLdmLossTildeLdm, icpLdmLossCorrTildeLdm) = {
//      val it = getIt(modelTildeLdm, prepareTarget(ldmTargetTildeLdm, obsset)._1, icpSettings, useRealingment = false, lookAheadRegression = true, mainAxis = axisd)
//      evalIt(it.map(t => (t._1, t._2)))
//    }
//    val (icpApprxLoss, icpApprxLossCorr) = {
//      val it = getIt(modelApprx, prepareTarget(origTarget, obsset)._1, icpSettings, useRealingment = false, lookAheadRegression = true, mainAxis = axisd)
//      evalIt(it.map(t => (t._1, t._2)))
//    }
    val (icpOnlineLoss, icpOnlineLossCorr) = {
      val it = getIt(model, prepareTarget(origTarget, obsset)._1, icpSettings, useRealingment = true, lookAheadRegression = true, mainAxis = axisd)
      evalIt(it.map(t => (t._1, t._2)))
    }
//    val (icpOrigLinScalingLoss, icpOrigLinScalingLossCorr) = { //linear scaling with unaligned basis
//      val it = getIt(ModelUtils.getLinearScalingModel(model.mean), prepareTarget(origTarget, obsset)._1, icpSettings, useRealingment = false, lookAheadRegression = true, mainAxis = axisd)
//      evalIt(it.map(t => (t._1, t._2)))
//    }
    val (icpLinScalingLoss, icpLinScalingLossCorr) = { //linear scaling with aligned basis
      val it = getIt(ModelUtils.getLinearScalingModel(model.mean, Option(obs)), prepareTarget(origTarget, obsset)._1, icpSettings, useRealingment = false, lookAheadRegression = true, mainAxis = axisd)
      evalIt(it.map(t => (t._1, t._2)))
    }

    val zres = {
      val zc = CompItem(0.0,0.0,0.0)
      IndexedSeq(zc,zc,zc)
    }
    val res = CompPostIcpResults(
      icpOrigLoss = icpOrigLoss,
      icpLdmLoss = zres,//icpLdmLoss,
      icpLdmLossLdm = zres,//icpLdmLossLdm,
      icpLdmLossTildeLdm = zres,//icpLdmLossTildeLdm,
      icpAlignLoss = icpAlignLoss,
      icpApprxLoss = zres,//icpApprxLoss,
      icpOnlineLoss = icpOnlineLoss,
      icpOrigLinScaleLoss = zres,//icpOrigLinScalingLoss,
      icpLinScaleLoss = icpLinScalingLoss,
      icpOrigLossCorr = icpOrigLossCorr,
      icpLdmLossCorr = zres,//icpLdmLossCorr,
      icpLdmLossCorrLdm = zres,//icpLdmLossCorrLdm,
      icpLdmLossCorrTildeLdm = zres,//icpLdmLossCorrTildeLdm,
      icpAlignLossCorr = icpAlignLossCorr,
      icpApprxLossCorr = zres,//icpApprxLossCorr,
      icpOnlineLossCorr = icpOnlineLossCorr,
      icpOrigLinScaleLossCorr = zres,//icpOrigLinScalingLossCorr,
      icpLinScaleLossCorr = icpLinScalingLossCorr
    )
    res
  }

  def getIt[A <: SamplingStrategy : Manifest](model: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], settings: IcpSettings[A], useRealingment: Boolean, lookAheadRegression: Boolean, mainAxis: EuclideanVector[_3D]): Iterator[(TriangleMesh[_3D], Boolean, Int)] = {
    settings.sampling match {
      case _: SamplingStrategyNormals => IcpFit.fitSurfaceAwareIteratorWithIndex(model, target, settings.asInstanceOf[IcpSettings[SamplingStrategyNormals]], useRealingment = useRealingment, lookAheadRegression = lookAheadRegression)
      case _: SamplingStrategy =>
        if (useRealingment)
          IcpFit.fitWithOnlineAlignmentIteratorWithIndex(model, target, settings.asInstanceOf[IcpSettings[SamplingStrategy]], withRotation = true, mainAxis = mainAxis)
        else
          IcpFit.iteratorWithIndex(model, target, settings.asInstanceOf[IcpSettings[SamplingStrategy]])
    }
  }

  case class CompPostIcpResults(
                             icpOrigLoss: IndexedSeq[CompItem], icpLdmLoss: IndexedSeq[CompItem], icpLdmLossLdm: IndexedSeq[CompItem], icpLdmLossTildeLdm: IndexedSeq[CompItem], icpAlignLoss: IndexedSeq[CompItem], icpApprxLoss: IndexedSeq[CompItem], icpOnlineLoss: IndexedSeq[CompItem], icpOrigLinScaleLoss: IndexedSeq[CompItem], icpLinScaleLoss: IndexedSeq[CompItem],
                             icpOrigLossCorr: IndexedSeq[CompItem], icpLdmLossCorr: IndexedSeq[CompItem], icpLdmLossCorrLdm: IndexedSeq[CompItem], icpLdmLossCorrTildeLdm: IndexedSeq[CompItem], icpAlignLossCorr: IndexedSeq[CompItem], icpApprxLossCorr: IndexedSeq[CompItem], icpOnlineLossCorr: IndexedSeq[CompItem], icpOrigLinScaleLossCorr: IndexedSeq[CompItem], icpLinScaleLossCorr: IndexedSeq[CompItem]
                            ) {
    def +(that: CompPostIcpResults): CompPostIcpResults = {
      CompPostIcpResults(
        this.icpOrigLoss.zip(that.icpOrigLoss).map(t => t._1+t._2), this.icpLdmLoss.zip(that.icpLdmLoss).map(t => t._1+t._2), this.icpLdmLossLdm.zip(that.icpLdmLossLdm).map(t => t._1+t._2), this.icpLdmLossTildeLdm.zip(that.icpLdmLossTildeLdm).map(t => t._1+t._2), this.icpAlignLoss.zip(that.icpAlignLoss).map(t => t._1+t._2), this.icpApprxLoss.zip(that.icpApprxLoss).map(t => t._1+t._2), this.icpOnlineLoss.zip(that.icpOnlineLoss).map(t => t._1+t._2), this.icpOrigLinScaleLoss.zip(that.icpOrigLinScaleLoss).map(t => t._1+t._2), this.icpLinScaleLoss.zip(that.icpLinScaleLoss).map(t => t._1+t._2),
        this.icpOrigLossCorr.zip(that.icpOrigLossCorr).map(t => t._1+t._2), this.icpLdmLossCorr.zip(that.icpLdmLossCorr).map(t => t._1+t._2), this.icpLdmLossCorrLdm.zip(that.icpLdmLossCorrLdm).map(t => t._1+t._2), this.icpLdmLossCorrTildeLdm.zip(that.icpLdmLossCorrTildeLdm).map(t => t._1+t._2), this.icpAlignLossCorr.zip(that.icpAlignLossCorr).map(t => t._1+t._2), this.icpApprxLossCorr.zip(that.icpApprxLossCorr).map(t => t._1+t._2), this.icpOnlineLossCorr.zip(that.icpOnlineLossCorr).map(t => t._1+t._2), this.icpOrigLinScaleLossCorr.zip(that.icpOrigLinScaleLossCorr).map(t => t._1+t._2), this.icpLinScaleLossCorr.zip(that.icpLinScaleLossCorr).map(t => t._1+t._2)
      )
    }

    def map(f: Double => Double): CompPostIcpResults = {
      CompPostIcpResults(
        icpOrigLoss.map(_.map(f)), icpLdmLoss.map(_.map(f)), icpLdmLossLdm.map(_.map(f)), icpLdmLossTildeLdm.map(_.map(f)), icpAlignLoss.map(_.map(f)), icpApprxLoss.map(_.map(f)), icpOnlineLoss.map(_.map(f)), icpOrigLinScaleLoss.map(_.map(f)), icpLinScaleLoss.map(_.map(f)),
        icpOrigLossCorr.map(_.map(f)), icpLdmLossCorr.map(_.map(f)), icpLdmLossCorrLdm.map(_.map(f)), icpLdmLossCorrTildeLdm.map(_.map(f)), icpAlignLossCorr.map(_.map(f)), icpApprxLossCorr.map(_.map(f)), icpOnlineLossCorr.map(_.map(f)), icpOrigLinScaleLossCorr.map(_.map(f)), icpLinScaleLossCorr.map(_.map(f))
      )
    }
  }

}
