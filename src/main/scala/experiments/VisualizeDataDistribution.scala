package experiments

import breeze.linalg.{DenseMatrix, DenseVector}
import fitting.{Correspondence, FitScripts}
import io.{LoadData, RunWriter}
import scalismo.{ModelRecenter, ModelUtils}
import scalismo.common.{PointId, ScalarMeshField, UnstructuredPoints3D}
import scalismo.geometry.{EuclideanVector3D, Point3D}
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.numerics.PivotedCholesky
import scalismo.registration.LandmarkRegistration
import scalismo.transformations.Rotation3D
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random
import utility.{Control, MathHelp, MeshUtils, Tchange}

import java.io.File

/**
 * used for a few visualization.
 */
object VisualizeDataDistribution {

  var ui: ScalismoUI = null

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    ui = ScalismoUI()
    implicit val rnd = scalismo.utils.Random(437342612L)

    val femurFolder = "./data/femurData/registered/"
    new File("./data/alignment/figures/randPoseMeshes/").mkdir()
    new File("./data/alignment/figures/partialTargets/").mkdir()

    val createRandPoses = true  //not used in the paper -> simply creates meshes that follow our pose assumption
    val createPartialTarget = true
    val createPriorPosteriorPcs = true
    val createSkullPoseAdj = false  //data is not available and this will fail
    val createFemurPoseAdj = true
    val createFaceViz = false  //needs bfm to work. produces - among other things - a projected eigenfunction.


    val origShapes = {
      println("loading data")
      val unalignedUncentered = new File(femurFolder).listFiles().map(MeshIO.readMesh).map(_.get).toIndexedSeq
      val unaligned = unalignedUncentered.map(mesh => {
        val center = Tchange.getMean(mesh).toVector
        mesh.transform(p => p - center)
      })
      val aligned = ModelUtils.alignShapesGpa(unaligned)._1.map(_._1)
      aligned
    }
    println("loaded and aligned data")

    val rpt = (300.0, 100.0, 5.0)
    val tf = 4.0
    val rpr = 2 * math.Pi

    if (createRandPoses) {
      origShapes.zipWithIndex.foreach { case (mesh, i) => {
        val t = EuclideanVector3D(tf * rnd.scalaRandom.nextDouble() * rpt._1, tf * rnd.scalaRandom.nextDouble() * rpt._2, tf * rnd.scalaRandom.nextDouble() * rpt._3)
        val r = (rnd.scalaRandom.nextDouble() * rpr, rnd.scalaRandom.nextDouble() * rpr, rnd.scalaRandom.nextDouble() * rpr)
        val R = Rotation3D.apply(r._1, r._2, r._3, Point3D.origin)
        val randpose = mesh.transform(p => R.apply(p) + t)
        ui.show(randpose, "mesh")
        MeshIO.writeMesh(randpose, new File(s"./data/alignment/figures/randPoseMeshes/mesh${i}.vtk"))
        MeshIO.writeMesh(randpose, new File(s"./data/alignment/figures/randPoseMeshes/mesh${i}.stl"))
      }
      }

      println("visualized uniform pose prior done")
    }


    if (createPartialTarget) {
      val pids = IndexedSeq(PointId(4054), PointId(4871), PointId(144))
      val pidsP = IndexedSeq(Point3D(1500, 0, -3000), Point3D(-2000, 1000, -3000), Point3D(0, 0, 6000))
      val mesh = origShapes.head
      val transf = LandmarkRegistration.rigid3DLandmarkRegistration(pids.map(id => mesh.pointSet.point(id)).zip(pidsP), Tchange.getMean(mesh))
      val m = mesh.transform(transf)
      val ratios = IndexedSeq(0.1, 0.2, 0.3, 0.45, 0.6, 0.8)
      val dirEncPids = (PointId(1840), PointId(3130))
      ratios.zipWithIndex.foreach { case (r, i) => {
        val (pm, _, _) = MeshUtils.getPartialMeshVector(m, m.pointSet.point(dirEncPids._1) - m.pointSet.point(dirEncPids._2), r, true)
        val t = EuclideanVector3D(i * 105.0, 0.0, 0.0)
        val pmm = pm.transform(p => p + t)
        val mm = m.transform(p => p + t)
        ui.show(pmm, s"ratio ${r}")
        ui.show(mm, s"ratio ${r} full").opacity = 0.2
        MeshIO.writeMesh(pmm, new File(s"./data/alignment/figures/partialTargets/meshRat${(r * 100).toInt}.vtk"))
        MeshIO.writeMesh(pmm, new File(s"./data/alignment/figures/partialTargets/meshRat${(r * 100).toInt}.stl"))
        MeshIO.writeMesh(mm, new File(s"./data/alignment/figures/partialTargets/meshRat${(r * 100).toInt}full.vtk"))
        MeshIO.writeMesh(mm, new File(s"./data/alignment/figures/partialTargets/meshRat${(r * 100).toInt}full.stl"))
      }
      }

      println("visualized partial targets done")
    }

    if (createPriorPosteriorPcs) {
      val priorFile = new File("./data/alignment/figures/pcPrior/")
      val posteriorFile = new File("./data/alignment/figures/pcPosterior25/")
      priorFile.mkdir()
      posteriorFile.mkdir()

      // the posterior calculation can sometimes struggle with numerical precision.
      val sigma2 = 10.0 // high value to make it easier to see. also when calculating the posterior with lower rrank this helps
      val perpendicularFactor = 100.0
      val rrank = 40
      val rratio = 0.05 // this is approximately the amount used in the figure
      //think about using surfaceaware as this helps target aware generally more
      val useSurfaceAware = false // if normals should be used during posterior calc

      /**
       * val decimationLevel = 500
       * val (shapes, dirEncPids) = {
       * val mmean = Tchange.getMean(origShapes)
       * val model = ModelUtils.pcaModel(mmean, origShapes).decimate(decimationLevel)
       * println(s"orig variance: ${breeze.linalg.sum(model.gp.variance)}")
       * val dirEncPids = IndexedSeq(PointId(1840), PointId(3130))
       * val newPids = dirEncPids.map(pid => model.mean.pointSet.findClosestPoint(mmean.pointSet.point(pid)).id)
       * (
       * (0 until model.rank).flatMap(i => IndexedSeq(-1, 1).map(d => model.instance(DenseVector.tabulate[Double](model.rank)(j => if (i == j) d else 0.0)))),
       * (newPids.head, newPids.last)
       * )
       * }
       * */
      val shapes = origShapes
      val dirEncPids = (PointId(1840), PointId(3130))
      val dirEncPidsSeq = Seq(dirEncPids._1, dirEncPids._2)
      val meanshape = Tchange.getMean(shapes)
      val naiveModel = ModelUtils.pcaModelGpa(shapes, shapes.head.triangulation.pointIds).truncate(rrank)

      val target = {
        val sample = { //a large bone is taken as example
          val samples = (0 to 10).map(_ => naiveModel.sample()).map(mesh => (mesh, dirEncPidsSeq.map(mesh.pointSet.point).map(_.toVector).reduce(_ - _).norm))
          samples.maxBy(_._2)._1
        }
        val targetInfo = MeshUtils.getPartialMeshVector(meanshape, meanshape.pointSet.point(dirEncPids._1) - meanshape.pointSet.point(dirEncPids._2), rratio, true)
        ModelUtils.alignShape(sample, meanshape, Option(targetInfo._2))._1
      }
      val targetInfo = MeshUtils.getPartialMeshVector(target, meanshape.pointSet.point(dirEncPids._1) - meanshape.pointSet.point(dirEncPids._2), rratio, true)
      val targetToModel = targetInfo._3 //Tchange.reverseMapping(meanshape.triangulation.pointIds, targetInfo._3)
      val alignedModel = ModelUtils.pcaModelGpa(shapes, targetInfo._2).truncate(rrank)
      assert(naiveModel.rank == alignedModel.rank, "easier that way")
      MeshIO.writeMesh(targetInfo._1, new File(s"./data/alignment/figures/target_25.stl"))

      /**
       * //to check if the mapping is correct
       * {
       * val trueCorrsFullToPartial = meanshape.pointSet.pointsWithId.map(t => (t._2, t._1, targetInfo._1.pointSet.findClosestPoint(t._1))).filter(t => (t._3.point - t._2).norm < 1e-4).map(t => (t._1, t._3.id)).toIndexedSeq
       * val trueCorrsPartialToFull = trueCorrsFullToPartial.map(t => (t._2, t._1)).toIndexedSeq
       * val test = trueCorrsFullToPartial.map(t => (t, targetInfo._3(t._1)))
       * val test2 = trueCorrsPartialToFull.map(t => (t, targetToModel(t._1)))
       * println("done")
       * }
       * */

      println(s"true length: ${dirEncPidsSeq.map(target.pointSet.point).map(_.toVector).reduce(_ - _).norm}")
      val v = DenseVector.tabulate[Double](rrank)(i => if (i == 0) 1.0 else 0.0)
      IndexedSeq((naiveModel, "naive"), (alignedModel, "aligned")).foreach { case (model, name) => {
        val posterior = if (useSurfaceAware) {
          val posteriorPdm = FitScripts.getSurfaceAwarePosterior(ModelUtils.ssmToPdm(model), target, Correspondence(targetToModel), varianceAlongNormal = sigma2, asSymmetricProposal = true, perpendicularFactor = perpendicularFactor)
          ModelUtils.pdmToSsm(posteriorPdm)
        } else model.posterior(targetInfo._2.map(pid => (pid, target.pointSet.point(pid))), sigma2)
        println(s"predicted ${name} map length: ${dirEncPidsSeq.map(posterior.mean.pointSet.point).map(_.toVector).reduce(_ - _).norm}")
        val dir = dirEncPidsSeq.map(target.pointSet.point).map(_.toVector).reduce(_ - _).normalize.toBreezeVector
        val dotLengthUncertainty = dir.t * posterior.cov(dirEncPids._1, dirEncPids._1) * dir
        println(s"${name} length uncertainty ${dotLengthUncertainty}")
        ui.show(meanshape, s"meanshape-${name}")
        ui.show(UnstructuredPoints3D(targetInfo._2.map(meanshape.pointSet.point)), s"pcl-${name}")

        println(s"prior and posterior variance left for ${name}: ${breeze.linalg.sum(model.gp.variance)} | ${breeze.linalg.sum(posterior.gp.variance)}")
        ui.show(ui.createGroup(s"prior model ${name}"), model, "prior")
        ui.show(ui.createGroup(s"posterior model ${name}"), posterior, "posterior")
        IndexedSeq((model, priorFile), (posterior, posteriorFile)).foreach { case (model, folder) => {
          IndexedSeq(-3.0, 0.0, 3.0).foreach(alpha => {
            MeshIO.writeMesh(model.instance(alpha * v), new File(folder, s"${name}_${if (alpha < 0) "m" else if (alpha > 0) "p" else ""}${math.abs(alpha.toInt)}.stl"))
          })
        }
        }
      }
      }

      println("visualizing prior and posterior pcs done")
    }

    if (createSkullPoseAdj){
      shapePCviz(true)
    }
    if (createFemurPoseAdj){
      shapePCviz(false)
    }
    if (createFaceViz){
      faceRealignmentTest()
    }

  }

  def shapePCviz(skull: Boolean)(implicit rng: Random): Unit = {
    val modelpath = if (skull) {
      "./data/skulls/skull-osteotomy/skullmodel.h5"
    } else "./data/femurData/ssm.h5"
    val outPath = if (skull) {
      "./data/skulls/skullViz/"
    } else "./data/femurData/femurViz/"
    val model = StatisticalModelIO.readStatisticalMeshModel(new File(modelpath)).get
    val (xdomain, infSet) = if (skull) {
      val meshInf = LoadData.readCsv(new File("./data/skulls/skullViz/left_cheek.csv"))
      //should be +-1e-2 to 0
      val checkCorr = meshInf.map(t => (model.mean.pointSet.point(t._1), t._2)).map(t => t._1-t._2).map(_.norm)
      val xdomain = meshInf.map(_._1)
      val infSet = meshInf.map(_._1.id).toSet
      (xdomain, infSet)
    } else {
      val dirEncPids = (PointId(1840), PointId(3130))  //(PointId(3130), PointId(1840))  // observe femur head
      val (_,obs,_) = MeshUtils.getPartialMeshVector(model.mean, model.mean.pointSet.point(dirEncPids._1) - model.mean.pointSet.point(dirEncPids._2), 0.2, true)
      (obs, obs.map(_.id).toSet)
    }
    println("read all data")

    val untilRank = 5
    val std = 2.0
    val nsamples = 100
    val createMeshForViz = false

    //create samples along pcs
    val zeroCoeff = DenseVector.zeros[Double](model.rank)
    val axis = IndexedSeq(EuclideanVector3D.unitX, EuclideanVector3D.unitY, EuclideanVector3D.unitZ)
    Seq((model, "naive"), (ModelRecenter.completeAlign(model, xdomain, axis), "align")).foreach(mmodel => {
      val outd = s"${outPath}/"
      new File(outd).mkdirs()
      val model = mmodel._1
      val smallAid = if (mmodel._2.equals("naive")) "n" else "a"
      if (createMeshForViz) {
        val meshesLabeled = (0 until untilRank).map(i => {
          val vec = zeroCoeff
          vec(i) = 1.0
          IndexedSeq(std*vec, -std*vec).zip(IndexedSeq(s"${std.toInt}p${i}d", s"${std.toInt}m${i}d"))
        })
        def writemeshes(t: (DenseVector[Double], String)): Unit = {
          val mesh = model.instance(t._1)
          //save standard mesh
          MeshIO.writeMesh(mesh, new File(s"${outd}${smallAid}${t._2}.stl")).get
          //MeshIO.writeMesh(mesh, new File(s"${outd}${smallAid}${t._2}.vtk")).get
          //cheek mesh slightly moved along normal -> makes viz nicer. check if mesh changed (could have changing signs in regions)
          val cheek = {
            val cut = mesh.operations.maskPoints(pid => infSet.contains(pid.id)).transformedMesh
            cut.transform(p => {
              val ppid = cut.pointSet.findClosestPoint(p)
              p + cut.vertexNormals(ppid.id).*(1e-1)
            })
          }
          MeshIO.writeMesh(cheek, new File(s"${outd}/${smallAid}${t._2}.stl")).get
          //MeshIO.writeMesh(cheek, new File(s"${outd}/${smallAid}${t._2}.vtk")).get
          println(s"done with ${mmodel._2}: ${t._2}")
        }
        //meshesLabeled.foreach(_.foreach(writemeshes))
        //writemeshes((zeroCoeff, "mean"))
      }
      //create visualization of translation and rotation space
      val poseParas = (0 until nsamples).map(i => {
        val sample = model.sample()
        val transform = ModelUtils.alignShape(sample, model.mean, Option(xdomain))._2
        if (i % 5 == 0) println(s"creating pose data for ${mmodel._2}: ${i} out of ${nsamples}")
        transform.parameters
      })
      val poseMat = DenseMatrix.apply(poseParas:_*)
      val writer = RunWriter.apply(new File(s"${outd}/poseSamples_${mmodel._2}.py"), pythonStyle = true)
      writer.writeMatrix[Double](poseMat, "pmat")
      writer.close()
    })
    if (createMeshForViz) println("created visualization meshes")
    println("done")
  }


  def faceRealignmentTest()(implicit rng: Random): Unit = {
    //a test to check the difference in samples from realigned models compared to original meshes
    //val ui = ScalismoUI()

    val model = LoadData.getFaceShapeModel()

    val namedSegmWithoutDef = LoadData.getFaceSegmentationSelection(IndexedSeq(IndexedSeq(4,5,6),IndexedSeq(9,10)))
    val segm = namedSegmWithoutDef.map(_.toIndexedSeq)++IndexedSeq(IndexedSeq(model.mean.pointSet.findClosestPoint(Point3D(0.0,-75.0,110.0)).id.id))

    //    println(s"average: ${Control.toS(Control.benchmarkit(()=>ModelRecenter.recenter(model,segm.head.map(PointId)),10,2),Control.TimeUnits.Ns)}")

    //    val sampleGroup = ui.createGroup("ExtremeSamples")
    val sampleCoef = DenseVector.tabulate[Double](model.rank){_ => rng.scalaRandom.nextGaussian()*2.0}//to visualize potential differences better
    val selSim = rng.scalaRandom.shuffle(model.mean.pointSet.pointIds.toIndexedSeq).take(400)
    val selSimRotCenter = selSim.map(model.mean.pointSet.point).map(_.toVector).reduce(_+_).map(_./(selSim.length.toDouble))
    val axis = IndexedSeq(EuclideanVector3D.unitX, EuclideanVector3D.unitY, EuclideanVector3D.unitZ)
    Seq(("default",IndexedSeq.empty[PointId]),("nose", segm(0).map(PointId)),("rightear", segm(1).map(PointId))/*,("rightearSample", segm(1).map(PointId))*//*,("chin", segm(2).map(PointId))*/).foreach(t => {
      //      val group = ui.createGroup(t._1)
      val pselect = t._2//.zipWithIndex.filter(_._2%(t._2.length/50)==0).map(_._1)
      val ssm = if (pselect.nonEmpty) {
        //        ui.show(group, UnstructuredPoints(pselect.map(model.referenceMesh.pointSet.point)), "centeringPoints").opacity=0.0
        println(s"transforming ${t._1} ...")
        val rerotateCenter = Option(Tchange.getMean(model.mean, Option(pselect)))
        if (t._1.contains("Sample")){
          val sr = 800  //the implementations can be optimized quite a bit more. this just as a starting point
          //has higher rank (potentially up to discretization/sample limit) than original model -> only calculate 'model.rank' number of eigenfunctions.
          Control.timeitVerbose(()=>ModelUtils.pdmToSsm(ModelUtils.buildSampleModelId(ModelUtils.ssmToPdm(model), Option(pselect), sr, PivotedCholesky.NumberOfEigenfunctions(math.min(model.rank, sr)-1))), s"${t._1} realignmentSamples")
        } else {
          //Control.timeitVerbose(()=>ModelRecenter.recenterSsm(model, pselect, rediag=false), s"${t._1} recenter")
          Control.timeitVerbose(()=>ModelRecenter.completeAlign(model, pselect, axis, rerotateCenter, rediag=false), s"${t._1} realignment")
        }
      } else model
      //      ui.show(group, ssm, t._1).meshView.opacity=0.0

      //saves a visualization of the alignment area
      if (pselect.nonEmpty){
        val usedIds = pselect.map(_.id).toSet
        val mesh = ScalarMeshField.apply(ssm.mean, ssm.mean.pointSet.pointIds.map(pid => if (usedIds.contains(pid.id)) 1.0 else 0.0).toIndexedSeq)
        MeshIO.writeScalarMeshField(mesh, new File("./data/alignment/faces/", s"alignDomain-${t._1}.vtk"))
      } else {
        //printing out the default rotation detected on other domains
        val domains = IndexedSeq(("nose", segm(0).map(PointId)), ("rightear", segm(1).map(PointId)))
        val angleInfos = (0 until 20).flatMap(_ => {
          val mesh = ssm.sample()
          domains.map(d => {
            val mps = d._2.map(ssm.mean.pointSet.point(_))
            val rp = Tchange.getMean(mps)
            val transf = LandmarkRegistration.rigid3DLandmarkRegistration(d._2.map(pid => mesh.pointSet.point(pid)).zip(mps), rp)
            val tmesh = mesh.transform(transf)
            val angles = d._2.map(pid => tmesh.pointSet.point(pid)).zip(mps).map(t => {
              val (v1, v2) = (t._1-rp, t._2-rp)
              MathHelp.getInAngle(v1,v2)
            })
            val meanangle = Tchange.getMeanVar(angles).mean
            (d._1, meanangle)
          })
        })
        val grouped = angleInfos.groupBy(_._1)
        grouped.values.foreach(domangle => println(s"${domangle.head._1} with an average angle of ${Tchange.getMeanVar(domangle.map(_._2)).mean}"))
      }
      //saves model variance
      val mvariance = ssm.mean.pointSet.pointIds.map(pid => {
        val dat = breeze.linalg.diag(ssm.cov(pid, pid)).toArray
        dat.sum / dat.length
      })
      val mesh = ScalarMeshField.apply(ssm.mean, mvariance.toIndexedSeq)
      MeshIO.writeScalarMeshField(mesh, new File("./data/alignment/faces/", s"modelVar-${t._1}.vtk"))

      //saving the first pc for visualization
      (0 until 5).foreach(pc => {
        val ppc = ssm.instance(DenseVector.tabulate(ssm.rank)(i => if (i == pc) 3.0 else 0.0))
        val mpc = ssm.instance(DenseVector.tabulate(ssm.rank)(i => if (i == pc) -3.0 else 0.0))
        MeshIO.writeMesh(ppc, new File("./data/alignment/faces/", s"pcViz${t._1}${pc}Ppc3.vtk"))
        MeshIO.writeMesh(mpc, new File("./data/alignment/faces/", s"pcViz${t._1}${pc}Mpc3.vtk"))
      })

      if (pselect.nonEmpty){ //if there is an alignment happening
        println("measuring avg mesh distances")
        val meshDistances = if (t._1.contains("Sample")) {
          (0 until 20).map(_ => {
            //average both directions
            val inst = ssm.sample()
            val align = ModelUtils.alignShape(inst, model.mean)
            val instAlign = model.project(align._1)
            val (m1,m2) = (align._1, instAlign)
            val apprxerr1 = Tchange.getDef(m2.pointSet, m1).map(_.norm)
            val mean1 = Tchange.getMeanVar(apprxerr1).mean //avg point error
            val inst2 = model.sample()
            val align2 = ModelUtils.alignShape(inst2, ssm.mean, Option(pselect))
            val instAlign2 = ssm.project(align2._1)
            val (m21,m22) = (align2._1, instAlign2)
            val apprxerr2 = Tchange.getDef(m22.pointSet, m21).map(_.norm)
            val mean2 = Tchange.getMeanVar(apprxerr2).mean //avg point error
            println(s"mean1 $mean1 mean2 $mean2")
            ((mean1 + mean2) / 2.0, apprxerr1.zip(apprxerr2).map(t => (t._1+t._2)/2))
          })
        } else {
          (0 until 20).map(_ => {
            val sampleCoef = DenseVector.tabulate[Double](model.rank){_ => rng.scalaRandom.nextGaussian()} //normal coefficient sampling
            val inst = ssm.instance(sampleCoef)
            val instAlign = model.instance(sampleCoef)
            val simTransform = LandmarkRegistration.rigid3DLandmarkRegistration(selSim.map(pid => (inst.pointSet.point(pid),instAlign.pointSet.point(pid))), selSimRotCenter.toPoint)
            val (m1,m2) = (inst.transform(simTransform), instAlign)
            val apprxerr = Tchange.getDef(m2.pointSet, m1).map(_.norm)
            (Tchange.getMeanVar(apprxerr).mean, apprxerr)
          })
        }
        val sameMeshDistances = meshDistances.map(_._1)
        println(s"${t._1}: mean overall distance between technically same mesh: ${Tchange.getMeanVar(sameMeshDistances).mean}")
        val allDistances = meshDistances.map(_._2)
        val averaged = allDistances.transpose.map(t => Tchange.getMeanVar(t).mean)
        val mesh = ScalarMeshField.apply(ssm.mean, averaged)
        MeshIO.writeScalarMeshField(mesh, new File("./data/alignment/faces/", s"apprxErr-${t._1}.vtk"))
      }

      val inst = ssm.instance(sampleCoef)
      val instAlign = model.instance(sampleCoef)
      val simTransform = LandmarkRegistration.rigid3DLandmarkRegistration(selSim.map(pid => (inst.pointSet.point(pid),instAlign.pointSet.point(pid))), selSimRotCenter.toPoint)
      //      ui.show(sampleGroup, inst.transform(simTransform), t._1).opacity=0.0
      val predVar = ssm.posterior(segm(1).map(PointId).map(pid => (pid, ssm.mean.pointSet.point(pid))),1.0).gp.variance.data.sum / ssm.mean.pointSet.numberOfPoints
      println(s"posterior variance for ${t._1} aligned model given ear points: ${predVar}")
    })
  }
}
