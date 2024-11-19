package experiments

import breeze.linalg.{DenseMatrix, DenseVector}
import io.BinaryWriter
import scalismo.{ModelRecenter, ModelUtils}
import scalismo.common.PointId
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.io.MeshIO
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.utils.Random
import utility.MeshUtils

import java.io.File

/**
 * this showcases the bias in posterior shape models.
 */
object TestPosteriorReconstructions {

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    implicit val rng = scalismo.utils.Random(1234983L)
    val useFemur = true

    val input = if (useFemur) "./data/femurData/registered/" else "../data/scaphoid/scaphoid_registered/"
    val output = {
      val n = "./data/alignment/posteriorReconstructionExperiment/"
      val f = new File(n)
      if (!f.exists()) f.mkdirs()
      n
    }

    posteriorRecoTest(input, output)
  }

  def posteriorRecoTest(inputdir: String, outputdir: String)(implicit rng: Random): Unit = {
    // writes a list of posterior mean and covs out into a file for later processing
    val decimate = 200
    val iterations = 100000
    val ratioBounds = IndexedSeq(
      (1.5 / decimate, 4.5 / decimate), // leads to very few observations
      (1.5 / decimate, 0.25), // small parts observed
      // the more peaked a posterior is the more samples are required to adequately cover the divergence
      //(1.5 / decimate, 1.0), // if all observed, all models will perform the same -> not very interesting
    )
    // TODO think about using surface aware noise. that would favor the aligned methods even more
    val sigma2 = 48.0 // significant noise such that it does not requires an extreme number of posteriors for a reasonably good approximation of the true distribution in higher dimensions
    println("preparing data")
    // create a first model to decimate -> then align them again to form a proper decimated model
    val (naiveModel, alignedData) = {
      val unaligned = new File(inputdir).listFiles().map(MeshIO.readMesh).map(_.get).toIndexedSeq
      val roughAlign = ModelUtils.alignShapesGpa(unaligned)._1.map(_._1)
      val fullmodel = ModelUtils.pcaModel(roughAlign.head, roughAlign)
      val decimateModel = fullmodel.decimate(decimate)
      val decimatedData = roughAlign.map(mesh => decimateModel.instance(fullmodel.coefficients(mesh)))
      val alignedData = ModelUtils.alignShapesGpa(decimatedData)._1.map(_._1)
      (ModelUtils.pcaModel(alignedData.head, alignedData), alignedData)
    }
    // same ratio and obsDir for each posterior
    val ratioObs = (1 to iterations).map(_ => {
      val v = scalismo.geometry.EuclideanVector3D(rng.scalaRandom.nextGaussian(), rng.scalaRandom.nextGaussian(), rng.scalaRandom.nextGaussian())
      val rb = ratioBounds(rng.scalaRandom.nextInt(ratioBounds.length))
      (rng.scalaRandom.nextDouble() * (rb._2 - rb._1) + rb._1, v.normalize, PointId(rng.scalaRandom.nextInt(naiveModel.mean.pointSet.numberOfPoints)))
    })
    val ratioPoints = (1 to iterations).map(_ => rng.scalaRandom.nextBoolean())
    val rngLongs = (1 to iterations).map(_ => rng.scalaRandom.nextLong())
    val configs = {
      val dims = IndexedSeq(2, 3, 5, 10, 45) //the posteriors start to have some issues for very low ranks. -> divergence doesnt fully converge for target specific models
      // the (true, false, d), (true, true, d) will produce the same output because the target specific model space is orthogonal to x-observed translation
      dims.flatMap(d => IndexedSeq((false, false, d), (false, true, d), (true, false, d), (true, true, d)))
    }
    println("beginning simulation")
    configs.par.foreach { case (align, tdAlign, rank) =>
      val naiveModelTruncated = naiveModel.truncate(rank)
      val writer = BinaryWriter(s"${outputdir}dat_${if (align) "align" else "naive"}_${if (tdAlign) "tdA_" else ""}${rank}.txt")
      writer.outputStream.writeInt(iterations)
      writer.writedv(naiveModelTruncated.gp.variance)
      ratioObs.zip(rngLongs.zip(ratioPoints)).zipWithIndex.foreach(t => {
        val obsTuple = t._1._1.copy(_3 = naiveModelTruncated.mean.pointSet.point(t._1._1._3))
        val ids = getObsIds(naiveModelTruncated.mean, t._1._2._2, obsTuple)
        val model = if (align) getAlignModel(alignedData, ids).truncate(rank) else naiveModelTruncated
        val result = handlePosterior(model, naiveModelTruncated, ids, align, sigma2, tdAlign)(scalismo.utils.Random(t._1._2._1))
        writer.writedv(result._1)
        writer.writedm(result._2)
      })
      writer.close()
      println(s"done with config (${align}, ${tdAlign}, ${rank})")
    }
  }

  def getObsIds(mesh: TriangleMesh[_3D], pointSource: Boolean, obs: (Double, EuclideanVector[_3D], Point[_3D])): IndexedSeq[PointId] = {
    if (pointSource) {
      MeshUtils.getPartialMeshPoint(mesh, obs._3, obs._1)._2
    } else {
      MeshUtils.getPartialMeshVector(mesh, obs._2, obs._1)._2
    }
  }

  def getAlignModel(meshes: IndexedSeq[TriangleMesh[_3D]], ids: IndexedSeq[PointId]): StatisticalMeshModel = {
    ModelUtils.pcaModelGpa(meshes, ids, false)
  }

  def handlePosterior(model: StatisticalMeshModel, naiveModel: StatisticalMeshModel, ids: IndexedSeq[PointId], align: Boolean, sigma2: Double, tdAlign: Boolean)(implicit rng: Random): (DenseVector[Double], DenseMatrix[Double]) = {
    //calculate the analytical posterior. be careful to use the correct observations -> including noise
    val td = getTrainingData(naiveModel, ids, math.sqrt(sigma2), tdAlign)
    val posteriorModel = model.posterior(td, sigma2)

    //map the noise back to the naiveModel - no projection required for naive model
    val projectPosteriorModel = if (align) ModelRecenter.recenterSsm(posteriorModel, model.mean.pointSet.pointIds.toIndexedSeq) else posteriorModel
    val posteriorCov = ModelUtils.getEmbeddedPosteriorCoeffCovariance(ModelUtils.ssmToPdm(naiveModel), ModelUtils.ssmToPdm(projectPosteriorModel))
    val posteriorMean = naiveModel.coefficients(projectPosteriorModel.mean)

    //return mean and covariance of the coefficients in the naiveModel space
    (posteriorMean, posteriorCov)
  }

  def getTrainingData(naiveModel: StatisticalMeshModel, ids: IndexedSeq[PointId], sigma: Double, tdAlign: Boolean)(implicit rng: Random): IndexedSeq[(PointId, Point[_3D])] = {
    val shape = naiveModel.sample()
    val tdshape = if (tdAlign) ModelUtils.alignShape(shape, naiveModel.mean, Option(ids), rotation = false)._1 else shape
    val td = ids.map(id => (id, tdshape.pointSet.point(id) + EuclideanVector(rng.scalaRandom.nextGaussian(), rng.scalaRandom.nextGaussian(), rng.scalaRandom.nextGaussian()).map(_ * sigma)))
    td
  }
}
