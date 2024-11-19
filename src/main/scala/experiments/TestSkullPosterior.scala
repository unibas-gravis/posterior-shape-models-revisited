package experiments

import fitting.BidirectionalSamplingFromOrigin
import io.RunReader
import scalismo.common.ScalarMeshField
import scalismo.common.interpolation.NearestNeighborInterpolator3D
import scalismo.geometry.{EuclideanVector3D, Point3D, _3D}
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.mesh.{SurfacePointProperty, TriangleMesh}
import scalismo.numerics.UniformMeshSampler3D
import scalismo.statisticalmodel.{PointDistributionModel, StatisticalMeshModel}
import scalismo.{ModelRecenter, ModelUtils}
import utility.Tchange

import java.io.File

/**
 * This is the qualitative skull reconstruction experiment.
 * we do not provide the data for this experiment. Nonetheless, the code can be inspected for transparency.
 */
object TestSkullPosterior {
  val iidNoise = 1.0*TestPartialAlignment.runSettings.iidNoise
  val perpendicularFactor = 1
  val origModelRank = 50
  val augModelRank = origModelRank * 3
  val cardAxis = IndexedSeq(EuclideanVector3D.unitX, EuclideanVector3D.unitY, EuclideanVector3D.unitZ)

  val dataFolder = "./data/skulls/skull-osteotomy/"
  val outputFolder = new File(dataFolder, "output/")

  val withTeeth = true
  val augmentedModel = false

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    implicit val rnd = scalismo.utils.Random(98237465L)

    val modelarbitraryF = if (augmentedModel) {
      val file = new File(dataFolder, "skullmodelGaussAugmented.h5")
      StatisticalModelIO.readStatisticalMeshModel(file).get
    } else {
      val file = new File(dataFolder, "skullmodel.h5")
      StatisticalModelIO.readStatisticalMeshModel(file).get
    }

    val (ulscmodel, ufscmodel, decimatedomegadomain, omegadomain) = {
      MeshIO.writeMesh(modelarbitraryF.referenceMesh, new File(dataFolder, "referenceMesh.stl"))
      println("parsing text files")
      val outerDomain = {
        val txtdomains = IndexedSeq("skullrightOuter.txt", "skullleftOuter.txt", "skullcap.txt", "skullteeth.txt")
        val files = (if (withTeeth) txtdomains else txtdomains.dropRight(1)).map(s => new File(dataFolder, s))
        val pointclouds = files.par.map(f => {
          val vars = RunReader(f).getVariables()
          vars.drop(1).map(s => s.split(",").drop(1).map(_.toDouble)).map(t => Point3D(t(0),t(1),t(2)))
        })
        val pids = pointclouds.flatMap(pcl => pcl.map(p => modelarbitraryF.mean.pointSet.findClosestPoint(p).id)) //or reference mesh
        pids.map(_.id).toIndexedSeq.toSet
      }
      println("done")
      val reducedMesh = modelarbitraryF.mean.operations.maskPoints(pid => outerDomain.contains(pid.id)).transformedMesh
      MeshIO.writeMesh(reducedMesh, new File(dataFolder, "reducedMesh.stl"))
      println("interpolating and discretizing to cleaner mesh...")
      val modelarbitrary = {
        val lgp = modelarbitraryF.gp.interpolate(NearestNeighborInterpolator3D()).discretize(reducedMesh)
        StatisticalMeshModel(reducedMesh, lgp)
      }
      println("done")
      val almostUniformSampledPids = {
        val uniformSampler = UniformMeshSampler3D(reducedMesh, 500)
        val pointsOnReduced = uniformSampler.sample()
        pointsOnReduced.map(t => modelarbitrary.mean.pointSet.findClosestPoint(t._1).id).distinct
      }
      val model = ModelRecenter.completeAlign(modelarbitrary, modelarbitrary.mean.pointSet.pointIds.toIndexedSeq, cardAxis).truncate(if (augmentedModel) augModelRank else origModelRank)
      //val model = ModelRecenter.recenterSsm(modelarbitrary, almostUniformSampledPids).truncate(if (augmentedModel) augModelRank else origModelRank)
      val decimated = model.decimate(1000)
      val decimatedAlmostUniformSampledPids = almostUniformSampledPids.map(pid => decimated.mean.pointSet.findClosestPoint(modelarbitrary.mean.pointSet.point(pid)).id).distinct
      val lsc = ModelRecenter.recenterSsm(decimated, decimatedAlmostUniformSampledPids)
      (lsc, model, decimatedAlmostUniformSampledPids, almostUniformSampledPids)
    }
    println("the models have updated references and are nicely aligned")
    //println(s"all pids aligned? : lsc ${ulscmodel.gp.klBasis.map(t => Tchange.getMeanCovDV(t.eigenfunction.data.map(_.toBreezeVector))).forall(t => t.mean.data.map(math.abs).sum < 1e-6)} fsc ${ufscmodel.gp.klBasis.map(t => Tchange.getMeanCovDV(t.eigenfunction.data.map(_.toBreezeVector))).forall(t => t.mean.data.map(math.abs).sum < 1e-6)}")

    val target = MeshIO.readMesh(new File(outputFolder, "ldmAlignedTarget.stl")).get

    val ostate = MeshIO.readMesh(new File(outputFolder, "origICP-2.vtk")).get
    val (lscmodel, modelFull) = {
      val transformation = ModelUtils.alignLdmsToLdms(modelarbitraryF.mean.pointSet.pointSequence.zip(ostate.pointSet.pointSequence))._2
      (ModelUtils.ssmToPdm(ulscmodel.transform(transformation)), ModelUtils.ssmToPdm(modelarbitraryF.transform(transformation)))
    }

    val stateFull = modelFull.project(ostate)
    val pidCorr = {
      val clst = lscmodel.mean.pointSet.pointsWithId.map(t => (t._2, modelFull.mean.operations.closestPointOnSurface(t._1))).toIndexedSeq
      val cps = Tchange.applyClpsOnMesh(stateFull, clst.map(_._2))
      clst.zip(cps).map(t => (t._1._1, t._2))
    }
    val state = lscmodel.posterior(pidCorr, 0.01).mean
    MeshIO.writeMesh(state, new File(outputFolder, "restructuredMesh.stl")).get
    println("finished state transformations")

    val corrSampling = BidirectionalSamplingFromOrigin()
    val obs = corrSampling.establishCorrespondenceUniform(state, target)
    println(s"domain size ${state.pointSet.numberOfPoints} and observations size ${obs.length}")

    //val origPost = FitScripts.getSurfaceAwarePosterior(lscmodel, target, corrSampling, Option(state), varianceAlongNormal = iidNoise*iidNoise, perpendicularFactor = perpendicularFactor)
    val origPost = lscmodel.posterior(obs, sigma2 = iidNoise*iidNoise)
    val apprxlscmodel = ModelRecenter.completeAlignPdm(lscmodel, obs.map(_._1), axis = cardAxis)//, rp = Option(Tchange.getMean(lscmodel.mean.pointSet.pointSequence)))
    //val apprxlscmodel = ModelRecenter.recenterPdm(lscmodel, obs.map(_._1))
    //val apprxPost = FitScripts.getSurfaceAwarePosterior(apprxlscmodel, target, corrSampling, Option(state), varianceAlongNormal = iidNoise*iidNoise, perpendicularFactor = perpendicularFactor)
    val apprxPost = apprxlscmodel.posterior(obs, sigma2 = iidNoise*iidNoise)

    def getVar(model: PointDistributionModel[_3D,TriangleMesh]): Double = breeze.linalg.sum(model.gp.variance)
    println(s"variance: origmodel ${getVar(lscmodel)} apprxModel ${getVar(apprxlscmodel)}\norigPost ${getVar(origPost)} apprxPost ${getVar(apprxPost)}")
//    def getEigfNorms(model: PointDistributionModel[_3D,TriangleMesh]): IndexedSeq[Double] = {
//      val norms = breeze.linalg.norm(model.gp.basisMatrix, Axis._0)
//      norms.inner.data.toIndexedSeq
//    }
//    println(s"eignorms: origmodel ${getEigfNorms(lscmodel)}\napprxModel ${getEigfNorms(apprxlscmodel)}\norigPost ${getEigfNorms(origPost)}\napprxPost ${getEigfNorms(apprxPost)}")

    val origvariance = state.pointSet.pointIds.map(pid => Tchange.getMeanVar(breeze.linalg.diag(origPost.cov(pid,pid)).data).mean).toIndexedSeq
    val apprxvariance = state.pointSet.pointIds.map(pid => Tchange.getMeanVar(breeze.linalg.diag(apprxPost.cov(pid,pid)).data).mean).toIndexedSeq
    IndexedSeq((origvariance,"orig"), (apprxvariance,"apprx")).foreach(t => {
      val sprop = SurfacePointProperty.apply[Double](state.triangulation, t._1)
      val fcvalues = Tchange.handleUcorSurface(sprop, stateFull.pointSet.pointSequence.map(p => state.operations.closestPointOnSurface(p)))
      val smesh = ScalarMeshField(stateFull, fcvalues)
      MeshIO.writeScalarMeshField(smesh, new File(outputFolder, s"fairSimplePosterior${t._2}.vtk"))
    })

//    val uniformSampler = UniformMeshSampler3D(lscmodel.mean, 1000)

//    val ui = ScalismoUI()
//    ui.show(ui.createGroup("origModel"), lscmodel, "origModel")
//    ui.show(ui.createGroup("apprxModel"), apprxlscmodel, "apprxModel")
//    ui.show(ui.createGroup("estimated domain"), UnstructuredPoints3D(obs.map(t => state.pointSet.point(t._1))), "estimated domain")
//    ui.show(ui.createGroup("uniform surface domain"), UnstructuredPoints3D(uniformSampler.sample().map(_._1)), "uniform surface domain")
//    ui.show(ui.createGroup("model domain"), UnstructuredPoints3D(lscmodel.mean.pointSet.pointSequence), "model domain")
//    ui.show(ui.createGroup("origPost"), origPost, "origPost")
//    ui.show(ui.createGroup("apprxPost"), apprxPost, "apprxPost")
  }
}
