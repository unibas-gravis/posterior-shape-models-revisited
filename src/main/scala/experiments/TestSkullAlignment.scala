package experiments

import breeze.linalg.DenseVector
import fitting.FitScripts.Sample
import fitting.IcpFit.IcpSettings
import fitting._
import io.RunWriter
import norms.L2norm
import scalismo.common.interpolation.NearestNeighborInterpolator3D
import scalismo.common.{PointId, ScalarMeshField, UnstructuredPoints}
import scalismo.geometry._
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.mesh.{SurfacePointProperty, TriangleList, TriangleMesh}
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.utils.Random
import scalismo.{ModelRecenter, ModelUtils}
import utility.{MeshUtils, Tchange}

import java.io.File

/**
 * fitting code for the skull targets. tests a number of methods, all results are difficult to interpret without gt.
 * we do not provide the data for this experiment. Nonetheless, the code can be inspected for transparency.
 */
object TestSkullAlignment {

  //changing some running parameters as the skulls have a lot of points with a somewhat significant amount of noise in some areas
  val runSettings = {
    val iidNoise = 6.0*TestPartialAlignment.runSettings.iidNoise
    val isoShapeScale = 0.5*TestPartialAlignment.runSettings.isoShapeScale
    val shapeScale = 1.0*TestPartialAlignment.runSettings.shapeScale
    val rScale = 1.0*TestPartialAlignment.runSettings.rScale
    val tScale = 1.0*TestPartialAlignment.runSettings.tScale
    val sigma = iidNoise
    val icpIterations = 50
    val numSamples = 10*TestPartialAlignment.runSettings.numSamples
    val burnin = TestPartialAlignment.runSettings.burnin

    TestPartialAlignment.runSettings.copy(iidNoise=iidNoise, isoShapeScale=isoShapeScale, shapeScale=shapeScale, rScale=rScale, tScale=tScale, sigma=sigma, icpIterations = icpIterations, numSamples=numSamples, burnin=burnin)
  }
  val perpendicularFactor = 10
  val backgroundll = -0.06  //is applied for every point not in correspondence -> depending on correspondence heuristic will not influence anything
  val ldmSigmaFactor = 0.1  //in case landmarks need to have different associated certainty
  val origModelRank = 40
  val augModelRank = origModelRank * 2
  val cardAxis = IndexedSeq(EuclideanVector3D.unitX, EuclideanVector3D.unitY, EuclideanVector3D.unitZ)

  val dataFolder = "./data/skulls/skull-osteotomy/"
  val outputFolder = new File(dataFolder, "output/")

  //some landmarks to help initialization
  val psldms = {
    val ssmldms = IndexedSeq(
      Point3D(67.7,408.8,261.2),
      Point3D(66.2,412.4,302.1),
      Point3D(-41.3,407.0,257.1),
      Point3D(-44.7,411.6,293.9),
      Point3D(15.2,390.0,222.0),
      Point3D(-42.4, 432.8, 321.9),
      Point3D(61.5, 434.7, 327.4),
      Point3D(8.3, 431.0, 348.8),
//      Point3D(-39.0,409.5,307.9), //no great correspondence
//      Point3D(59.8,412.4,313.4), //no great correspondence
//      Point3D(-40.7,419.5,310.3), //no great correspondence
//      Point3D(61.6,420.8,315.4), //no great correspondence
    )
    val targetldms = IndexedSeq(
      Point3D(69.4,409.5,264.0),
      Point3D(69.1,410.7,304.3),
      Point3D(-43.6,406.9,259.9),
      Point3D(-45.7,408.6,301.0),
      Point3D(14.4,387.8,226.7),
      Point3D(-40.1, 432.6, 324.3),
      Point3D(63.2, 432.2, 326.4),
      Point3D(10.3, 426.2, 354.4),
//      Point3D(-41.6,406.0,312.8), //no great correspondence
//      Point3D(67.7,409.1,319.1), //no great correspondence
//      Point3D(-36.3,415.1,314.8), //no great correspondence
//      Point3D(59.2,416.8,316.4), //no great correspondence
    )
    ssmldms.zip(targetldms)
  }


  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    implicit val rnd = scalismo.utils.Random(98237465L)

    val useParallelModel = true
    val useParallel = true
    val useOrigModel = false
    val useAugModel = true
    val useLandmarks = true
    if (!outputFolder.exists()) outputFolder.mkdir()

    val (writer, _) = RunWriter("./data/skulls/skull-osteotomy/skullRecoRes$t.txt", pythonStyle = true)
    runSettings.print(writer, None, None, None)
    writer.writeLiteral(s"ldmSigmaFactor = ${ldmSigmaFactor};")
    writer.writeLiteral(s"origRank = ${origModelRank};")
    writer.writeLiteral(s"augRank = ${augModelRank};")

    val (lscModel, fscModel) = {
      val file = new File(dataFolder, "skullmodel.h5")
      val modelarbitrary = StatisticalModelIO.readStatisticalMeshModel(file).get
      // this is the omega model aligned over the full domain. this is a good approximation of the true model
      val model = ModelRecenter.completeAlign(modelarbitrary, modelarbitrary.mean.pointSet.pointIds.toIndexedSeq, cardAxis).truncate(origModelRank)
      (model.decimate(1000), model)
    }
    val (lscModelAug, fscModelAug) = {
      val file = new File(dataFolder, "skullmodelGaussAugmented.h5")
      val modelarbitrary = StatisticalModelIO.readStatisticalMeshModel(file).get
      // this is the omega model aligned over the full domain. this is a good approximation of the true model.
      // this inclues some augmentation synthetic kernels which can also add pose to the shape. These issues are apprx removed.
      val model = ModelRecenter.completeAlign(modelarbitrary, modelarbitrary.mean.pointSet.pointIds.toIndexedSeq, cardAxis).truncate(augModelRank)
      (model.decimate(1000), model)
    }

    println(s"rank of origmodel ${lscModel.rank}, rank of augmented model ${lscModelAug.rank}")

    val ttargetNoLdms = {
      val file = new File(dataFolder, "skull-osteotomy-2-aligned-cleaned.stl")//"skull-osteotomy-2-aligned-simplified-cleaned.stl")
      MeshIO.readMesh(file).get.operations.decimate(1000)
    }

    //check for issues with mesh
    val issues = Seq(lscModel.mean, ttargetNoLdms).map(mesh => mesh.vertexNormals.pointData.exists(v => v.x.isNaN))
    val (fixedlscModel, fixedlscModelAug) = if (issues.head) { //cut the points with issue out of mesh
      val newMesh = MeshUtils.cutPointsWithUndefinedNormals(lscModel.mean, 10)
      val newmodel = ModelUtils.ssmToPdm(lscModel).newReference(newMesh, NearestNeighborInterpolator3D())
      val newmodelpreaug = ModelUtils.ssmToPdm(lscModelAug).newReference(newMesh, NearestNeighborInterpolator3D())
      (ModelUtils.pdmToSsm(newmodel), ModelUtils.pdmToSsm(newmodelpreaug))
    } else (lscModel, lscModelAug)
    assert(!fixedlscModel.mean.vertexNormals.pointData.exists(v => v.x.isNaN)) //the model should be clean
    val ttarget = if (issues.last){ //cut the points with issue out of mesh
      MeshUtils.cutPointsWithUndefinedNormals(ttargetNoLdms, 10)
    } else ttargetNoLdms
    assert(!ttarget.vertexNormals.pointData.exists(v => v.x.isNaN))

    val lldmsUnaligned = if (useLandmarks) { //to be clear: low scale ldms -> high scale probably not needed.
      val lldms = psldms.map(ldm => (lscModel.mean.pointSet.findClosestPoint(ldm._1).id, ldm._2))
      val lldmsAug = psldms.map(ldm => (lscModelAug.mean.pointSet.findClosestPoint(ldm._1).id, ldm._2))
      //non perfect safety check for correspondence between models
      assert(lldms.map(_._1.id).zip(lldmsAug.map(_._1.id)).forall(t => t._1==t._2))
      lldms
    } else IndexedSeq.empty

    val (target, lldms) = if (useLandmarks) {
      val l2norm = L2norm[_3D]()
      val ldmCorrp = lldmsUnaligned.map(t => (t._2, lscModel.mean.pointSet.point(t._1)))
      val tr = ModelUtils.alignLdmsToLdms(ldmCorrp)
      val transformation = tr._2
      val targetMesh = ttarget.transform(transformation)
      MeshIO.writeMesh(targetMesh, new File(outputFolder, "ldmAlignedTarget.stl"))
      val ldmCorrpAfter = lldmsUnaligned.map(t => (lscModel.mean.pointSet.point(t._1), transformation.f(t._2)))
      println(s"avg ldm distance before: ${l2norm.normVector(ldmCorrp.map(t=>t._1-t._2))} and after: ${l2norm.normVector(ldmCorrpAfter.map(t=>t._1-t._2))} ldm alignment")
      (targetMesh, lldmsUnaligned.map(t => (t._1, transformation.f(t._2))))
    } else (ttarget, lldmsUnaligned)

    //val corSampling = NormalSamplingSimpleExtension(BidirectionalSamplingFromOrigin())//NormalSamplingSimpleExtension(TargetSamplingUnique())
    //writer.writeLiteral(s"sampler =  '${"NormalExtension-BidirectionalSamplingFromOrigin"}';")//'${"NormalExtension-TargetSamplingUnique"}';")
    val corSampling = BidirectionalSamplingFromOrigin() //TargetSampling()
    writer.writeLiteral(s"sampler =  '${"BidirectionalSamplingFromOrigin"}';")//'${"NormalExtension-TargetSamplingUnique"}';")
    val icpSettings = IcpSettings[SamplingStrategy](
      corSampling, optAlignment = true, runSettings.icpIterations, sigma2 = runSettings.sigma * runSettings.sigma, perpendicularFactor = perpendicularFactor)

    val corSamplingInitialX = TargetSamplingUnique()
    writer.writeLiteral("samplerInitialX = 'TargetSamplingUnique';")

    val modelsSeq = {
      val ll1: IndexedSeq[(StatisticalMeshModel, StatisticalMeshModel, String)] = if (useOrigModel) IndexedSeq((fixedlscModel, fscModel, "orig")) else IndexedSeq.empty
      val ll2: IndexedSeq[(StatisticalMeshModel, StatisticalMeshModel, String)] = if (useAugModel) IndexedSeq((fixedlscModelAug, fscModelAug, "aug")) else IndexedSeq.empty
      if (useParallelModel) (ll1++ll2).par else ll1++ll2
    }

    modelsSeq.foreach{ case (model, umodel, modelLabel) =>
      //the target and model are adequately alignment -> estimate an initial correspondence
      val initialX = corSamplingInitialX.establishCorrespondenceUniform(model.mean, target).map(_._1)

      val axis = cardAxis
      val setOff = FitInfo(icpSettings, axis, writer, online = false)
      val setOn = FitInfo(icpSettings, axis, writer, online = true)

      //in case some other initialization is tested
      val initMesh = umodel.mean
      println(s"starting data preparation for ${modelLabel}...")
      val reconstructionfs = IndexedSeq(
        (icpReconstruction(model, target, setOff, umodel, lldms, s"${modelLabel}ICP", corSamplingInitialX), s"${modelLabel}ICP"),
        (icpReconstruction(ModelRecenter.completeAlign(model, initialX, axis), target, setOff, ModelRecenter.completeAlign(umodel, initialX, axis), lldms, s"${modelLabel}ApprxICP", corSamplingInitialX), s"${modelLabel}ApprxICP"),
        (icpReconstruction(model, target, setOn, umodel, lldms, s"${modelLabel}OnlineICP", corSamplingInitialX), s"${modelLabel}OnlineICP"),
        (mcmcReconstruction(model, initMesh, target, (setOff, s"${modelLabel}MCMC"), umodel, lldms, corSamplingInitialX), s"${modelLabel}MCMC"),
        (mcmcReconstruction(ModelRecenter.completeAlign(model, initialX, axis), initMesh, target, (setOff, s"${modelLabel}ApprxMCMC"), ModelRecenter.completeAlign(umodel, initialX, axis), lldms, corSamplingInitialX), s"${modelLabel}ApprxMCMC"),
        (mcmcReconstruction(model, initMesh, target, (setOn, s"${modelLabel}OnlineMCMC"), umodel, lldms, corSamplingInitialX), s"${modelLabel}OnlineMCMC"),
      )

      println(s"\npreparation complete for ${modelLabel}, starting reconstruction procedures...\n\n")
      (if (useParallel) reconstructionfs.par else reconstructionfs).foreach(t => {
        t._1.apply().zipWithIndex.foreach(meshi => {
          MeshIO.writeMesh(meshi._1, new File(outputFolder, t._2 + s"-${meshi._2+1}.vtk"))
        })
        println(s"done for ${t._2}")
      })
    }

    writer.writeLiteral("# Completed run and saved output meshes")
  }

  def mcmcReconstruction(model: StatisticalMeshModel, initm: TriangleMesh[_3D], target: TriangleMesh[_3D], ssettings: (FitInfo[SamplingStrategy],String), fullscModel: StatisticalMeshModel, ldms: IndexedSeq[(PointId, Point[_3D])], corrHeuristic: SamplingStrategyUniform)(implicit rnd: Random): () => IndexedSeq[TriangleMesh[_3D]] = {
    val (settings, sname) = ssettings
    val rp = Tchange.getMean(model.mean, Option(corrHeuristic.establishCorrespondenceUniform(model.mean, target).map(_._1)))
    val transf = ModelUtils.alignShape(fullscModel.mean, initm)
    val alignedFullscModel = fullscModel.transform(transf._2)
    val coeff = alignedFullscModel.coefficients(initm)
    val pmodel = ModelUtils.ssmToPdm(model.transform(transf._2))
    val pmodelf = ModelUtils.ssmToPdm(alignedFullscModel)
    val initSample = Sample(coeff, EuclideanVector3D.zero, EuclideanVector3D.zero, rp, settings.axis)
    if (settings.online)
      () => {
        MeshIO.writeMesh(initSample.instance(pmodelf), new File(outputFolder, ssettings._2+"-init.vtk"))
        val chain = FitScripts.fitPartialWithAlignment(pmodel, target, settings.axis.head, rp, iidNoiseSigma = runSettings.iidNoise, burnin = runSettings.burnin, initial = Option(initSample), numSamples = runSettings.numSamples, subsampling = runSettings.subsampling, tScale = runSettings.tScale, rScale = runSettings.rScale, shapeScale = runSettings.shapeScale, isoShapeScale = runSettings.isoShapeScale, ldms = ldms, ldmSigmaFactor = ldmSigmaFactor, backgroundll = backgroundll, perpendicularFactor = perpendicularFactor, correspondenceStrategy = corrHeuristic, samplingStrategy = settings.icpSettings.sampling)
        val variance = {
          val filtered = chain.zipWithIndex.filter(t => t._2 % (chain.length / 100 + 1) == 0)
          val vectorized = filtered.map(t => DenseVector(Tchange.getDef(t._1._1.pointSet, pmodel.mean).map(_.norm).toArray))
          Tchange.getMeanVarDV(vectorized).variance
        }
        settings.writer.writeSingle[Double](Tchange.getMeanVar(variance.data).mean, sname+"var")
        settings.writer.writeVector[Double](variance.data, sname+"varVerbose")
        val lcbest = chain.maxBy(_._2)._1
        val (best, pids) = {
          val pids = Option(corrHeuristic.establishCorrespondenceUniform(lcbest, target).map(_._1))
          (ModelUtils.upscaleFreeMesh(lcbest, pids, pmodel, pmodelf, alignModel = true, settings.axis), pids)
        }
        val (fsmesh, fsscmesh) = {
          val fcmesh = best
          val sprop = SurfacePointProperty.apply[Double](lcbest.triangulation, variance.data)
          val fcvalues = Tchange.handleUcorSurface(sprop, fcmesh.pointSet.pointSequence.map(p => lcbest.operations.closestPointOnSurface(p)))
          (fcmesh, ScalarMeshField(fcmesh, fcvalues))
        }
        MeshIO.writeMesh(fsmesh, new File(outputFolder, ssettings._2+"-best.vtk"))
        MeshIO.writeMesh(best.copy(pointSet = UnstructuredPoints(corrHeuristic.establishCorrespondenceUniform(best, target).map(_._1).map(best.pointSet.point)), triangulation = TriangleList.empty), new File(outputFolder, ssettings._2+"-best-pcld.vtk"))
        MeshIO.writeScalarMeshField(fsscmesh, new File(outputFolder, ssettings._2+"-best-var.vtk"))
        val ms = ((0 until 10).map(_ => rnd.scalaRandom.nextInt(chain.length-2))++IndexedSeq(chain.length-1)).map(chain(_)._1)
        ms ++ ms.map(m =>
          ModelUtils.upscaleFreeMesh(m, pids, pmodel, pmodelf, alignModel = true, settings.axis)) //recalculate a model for each sample
      }
    else
      () => {
        MeshIO.writeMesh(initSample.instance(pmodelf), new File(outputFolder, ssettings._2+"-init.vtk"))
        val chain = FitScripts.fitPartialStrat(pmodel, target, settings.axis.head, rp, iidNoiseSigma = runSettings.iidNoise, burnin = runSettings.burnin, initial = Option(initSample), numSamples = runSettings.numSamples, subsampling = runSettings.subsampling, tScale = runSettings.tScale, rScale = runSettings.rScale, shapeScale = runSettings.shapeScale, isoShapeScale = runSettings.isoShapeScale, ldms = ldms, ldmSigmaFactor = ldmSigmaFactor, backgroundll = backgroundll, perpendicularFactor = perpendicularFactor, samplingStrategy = settings.icpSettings.sampling)
        val variance = {
          val filtered = chain.zipWithIndex.filter(t => t._2 % (chain.length / 100 + 1) == 0).map(_._1._1.instance(pmodel))
          val vectorized = filtered.map(t => DenseVector(Tchange.getDef(t.pointSet, pmodel.mean).map(_.norm).toArray))
          Tchange.getMeanVarDV(vectorized).variance
        }
        settings.writer.writeSingle[Double](Tchange.getMeanVar(variance.data).mean, sname+"var")
        settings.writer.writeVector[Double](variance.data, sname+"varVerbose")
        val best = chain.maxBy(_._2)
        val bestm = best._1.instance(pmodel)
        val (fsmesh, fsscmesh) = {
          val fcmesh = best._1.instance(pmodelf)
          val lcmesh = bestm
          val sprop = SurfacePointProperty.apply[Double](best._1.instance(pmodel).triangulation, variance.data)
          val fcvalues = Tchange.handleUcorSurface(sprop, fcmesh.pointSet.pointSequence.map(p => lcmesh.operations.closestPointOnSurface(p)))
          (fcmesh, ScalarMeshField(fcmesh, fcvalues))
        }
        MeshIO.writeMesh(fsmesh, new File(outputFolder, ssettings._2+"-best.vtk"))
        MeshIO.writeMesh(bestm.copy(pointSet = UnstructuredPoints(corrHeuristic.establishCorrespondenceUniform(bestm, target).map(_._1).map(bestm.pointSet.point)), triangulation = TriangleList.empty), new File(outputFolder, ssettings._2+"-best-pcld.vtk"))
        MeshIO.writeScalarMeshField(fsscmesh, new File(outputFolder, ssettings._2+"-best-var.vtk"))
        val ms = ((0 until 10).map(_ => rnd.scalaRandom.nextInt(chain.length-2))++IndexedSeq(chain.length-1)).map(chain(_)._1)
        ms.map(_.instance(pmodel)) ++ ms.map(_.instance(pmodelf))
      }
  }

  def icpReconstruction(model: StatisticalMeshModel, target: TriangleMesh[_3D], settings: FitInfo[SamplingStrategy], fullscModel: StatisticalMeshModel, ldms: IndexedSeq[(PointId, Point[_3D])], outputName: String, corrHeuristic: SamplingStrategyUniform): () => IndexedSeq[TriangleMesh[_3D]] = {
    val pdm = ModelUtils.ssmToPdm(model)
    val pdmf = ModelUtils.ssmToPdm(fullscModel)
    if (settings.online)
      () => {
        val converged = {
          val it = settings.icpSettings match {
            case _:SamplingStrategyNormals => IcpFit.fitSurfaceAwareIteratorWithIndex(pdm, target, settings.icpSettings.asInstanceOf[IcpSettings[SamplingStrategyNormals]], settings.axis.head, ldms = ldms, alignmentSettings = Option(corrHeuristic)).map(t => (t._1,t._2))
            case _:SamplingStrategy => IcpFit.fitWithOnlineAlignmentIterator(pdm, target, settings.icpSettings, withRotation = true, mainAxis = settings.axis.head, ldms = ldms, alignmentSettings = Option(corrHeuristic))
          }
          it.drop(settings.icpSettings.maxIterations-1).next()._1
        }
        val obsPids = corrHeuristic.establishCorrespondenceUniform(converged, target).map(_._1)
        val upscale = ModelUtils.upscaleFreeMesh(converged, Option(obsPids), pdm, pdmf, alignModel = true, axis = settings.axis)
        //val obsPidsUpscale = obsPids.map(pid => upscale.pointSet.findClosestPoint(converged.pointSet.point(pid)).id)
        val analyticalVariance = {
          val model = pdmf
          //val allObsPidsUpscale = settings.icpSettings.sampling.establishCorrespondence(upscale, target).map(_._1)
          //val posterior = model.posterior(allObsPidsUpscale.map(pid => (pid, upscale.pointSet.point(pid))), sigma * sigma)
          val posterior = FitScripts.getSurfaceAwarePosterior(model, target, settings.icpSettings.sampling, Option(upscale), varianceAlongNormal = runSettings.sigma * runSettings.sigma)
          pdmf.mean.pointSet.pointIds.map(pid => {
            val dat = breeze.linalg.diag(posterior.cov(pid, pid)).toArray
            dat.sum / dat.length
          }).toIndexedSeq
        }
        val mesh = ScalarMeshField.apply(upscale, analyticalVariance)
        MeshIO.writeScalarMeshField(mesh, new File(outputFolder, outputName+"-best-var.vtk"))
        IndexedSeq(converged, upscale)
      }
    else {
      () => {
        val converged = IcpFit.apply(ModelUtils.ssmToPdm(model), target, settings.icpSettings, ldms)._1
        val obsPids = corrHeuristic.establishCorrespondenceUniform(converged, target).map(_._1)
        val upscale = ModelUtils.upscaleFreeMesh(converged, Option(obsPids), pdm, pdmf, alignModel = false)
        val obsPidsUpscale = obsPids.map(pid => upscale.pointSet.findClosestPoint(converged.pointSet.point(pid)).id)
        val analyticalVariance = {
          val model = ModelRecenter.completeAlignPdm(pdmf, obsPidsUpscale, cardAxis)
          //val allObsPidsUpscale = settings.icpSettings.sampling.establishCorrespondence(upscale, target).map(_._1)
          //val posterior = model.posterior(allObsPidsUpscale.map(pid => (pid, upscale.pointSet.point(pid))), sigma * sigma)
          val posterior = FitScripts.getSurfaceAwarePosterior(model, target, settings.icpSettings.sampling, Option(upscale), varianceAlongNormal = runSettings.sigma * runSettings.sigma)
          pdmf.mean.pointSet.pointIds.map(pid => {
            val dat = breeze.linalg.diag(posterior.cov(pid, pid)).toArray
            dat.sum / dat.length
          }).toIndexedSeq
        }
        val mesh = ScalarMeshField.apply(upscale, analyticalVariance)
        MeshIO.writeScalarMeshField(mesh, new File(outputFolder, outputName+"-best-var.vtk"))
        IndexedSeq(converged, upscale)
      }
    }
  }

  case class FitInfo[A <: SamplingStrategy](icpSettings: IcpSettings[A], axis: IndexedSeq[EuclideanVector[_3D]], writer: RunWriter, online: Boolean)
}

