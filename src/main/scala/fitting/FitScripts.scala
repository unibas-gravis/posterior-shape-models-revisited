package fitting

import breeze.linalg.{DenseVector, diag, norm}
import norms.L2norm
import scalismo.common.PointId
import scalismo.faces.render.Translation3D
import scalismo.geometry._
import scalismo.mesh.TriangleMesh
import scalismo.sampling.algorithms.MetropolisHastings
import scalismo.sampling.loggers.AcceptRejectLogger
import scalismo.sampling.proposals.MixtureProposal
import scalismo.sampling.{DistributionEvaluator, ProposalGenerator, SymmetricTransitionRatio}
import scalismo.statisticalmodel.{MultivariateNormalDistribution, PointDistributionModel}
import scalismo.transformations.Rotation3D
import scalismo.utils.{Memoize, Random}
import scalismo.{ModelRecenter, ModelUtils}
import utility.{MathHelp, Tchange}

object FitScripts {

  /**
   * optimizing over shape parameters with free correspondence. pose optimization uses estimated pose posterior.
   */
  def fitPartial(model: PointDistributionModel[_3D, TriangleMesh],
                 target: TriangleMesh[_3D],
                 mainAxis: EuclideanVector[_3D],
                 rotPoint: Point[_3D],
                 numSamples: Int = 2000, //number of final samples -> mh samples numSamples * subsampling
                 iidNoiseSigma: Double = 1.0,
                 shapeScale: Double = 1.0,
                 shapeUpdateNormalScale: Double = 1.0,
                 isoShapeScale: Double = 1.0,
                 tScale: Double = 1.0,
                 rScale: Double = 1.0,
                 subsampling: Int = 5,
                 burnin: Int = 100,
                 initial: Option[Sample] = None
                )(implicit rnd: Random): IndexedSeq[(Sample, Double)] = {
    //TODO apply pose of initial to these initial calculations
    //TODO correspondence in this posterior calc and the isoGenerator below can be unified
    val postModel = getSurfaceAwarePosterior(model, target, iidNoiseSigma * shapeUpdateNormalScale)
    val axis = MathHelp.listOrthogonalVectors(mainAxis)
    val stscale = IndexedSeq(0.2,0.6,0.2).zip(IndexedSeq(0.5, 1.0, 2.0))
    val alignIds = {
      val mapped = target.pointSet.pointsWithId.map(t => (t._2, model.mean.pointSet.findClosestPoint(t._1).id)).toIndexedSeq
      mapped.map(_._2.id).distinct.map(PointId)
    }

    val shapeProposal = getShapeProposalGenerator(model, postModel, stscale.map(t => (t._1,t._2*shapeScale)))
    val isoProposal = getScaledIsoGenerator(model, iidNoiseSigma, alignIds, stscale.map(t => (t._1,t._2*isoShapeScale)))
    val poseProposal = {
      val (tScales, rScales) = stscale.map(t => ((t._1,t._2*tScale),(t._1,t._2*rScale))).unzip
      getPoseProposalGenerator(model, alignIds, tScales, rScales, mainAxis, rotPoint)
    }
    val proposals = MixtureProposal.fromSymmetricProposals(IndexedSeq((0.5, shapeProposal), (0.2, isoProposal), (0.3, poseProposal)): _*)
    val eval = getPostEval(iidNoiseSigma, model, target)

    val mh = MetropolisHastings.apply(proposals, eval)
    val trackAcc = true
    val logger = new CustomLogger[Sample](false, trackAcc, model)
    val init = initial.getOrElse(Sample(DenseVector.zeros[Double](model.rank), EuclideanVector3D.zero, EuclideanVector3D.zero, rotPoint, axis))

    val chain = mh.iterator(init, logger).drop(burnin).zipWithIndex.filter(_._2 % subsampling == 0).map(t => (t._1, eval.logValue(t._1))).take(numSamples)


    val resChain = chain.toIndexedSeq
    if (trackAcc) println(s"logger ${logger.toString}")
    resChain
  }

  /**
   * optimizing over shape parameters with free correspondence. pose optimization uses estimated pose posterior.
   * the model is aligned every number of iterations
   */
  def fitPartialWithAlignment(model: PointDistributionModel[_3D, TriangleMesh],
                                 target: TriangleMesh[_3D],
                                 mainAxis: EuclideanVector[_3D],
                                 rotPoint: Point[_3D],
                                 numSamples: Int = 2000, //number of final samples -> mh samples numSamples * subsampling
                                 iidNoiseSigma: Double = 1.0,
                                 shapeScale: Double = 1.0,
                                 shapeUpdateNormalScale: Double = 1.0,
                                 isoShapeScale: Double = 1.0,
                                 tScale: Double = 1.0,
                                 rScale: Double = 1.0,
                                 subsampling: Int = 5,
                                 burnin: Int = 100,
                                 initial: Option[Sample] = None,
                                 realignments: Int = 5,
                                 realignmentRotation: Boolean = true,
                                 samplingStrategy: SamplingStrategy = BidirectionalSamplingFromOrigin() //this sampling strategy is quite expensive but safe, TargetSamplingUnique is cheaper option
                                )(implicit rnd: Random): IndexedSeq[(TriangleMesh[_3D], Double)] = {
    val axis = MathHelp.listOrthogonalVectors(mainAxis)
    val zsample = Sample(DenseVector.zeros[Double](model.rank), EuclideanVector3D.zero, EuclideanVector3D.zero, rotPoint, axis)
    val init = initial.getOrElse(zsample)
    (0 until realignments).scanLeft((IndexedSeq((init.instance(model),init, Double.MinValue)),model)){case ((pchain,cmodel),_) => {
      //create aligned model
      val clps = samplingStrategy.establishCorrespondence(pchain.last._1, target)
      val alignids = clps.map(_._1.id).distinct.map(PointId)
      val recentered = ModelRecenter.recenter(ModelUtils.pdmToSsm(model), alignids)
      val ssm = if (realignmentRotation) ModelRecenter.rerotate(recentered, alignids, axis, Option(rotPoint)) else recentered
      val amodel = ModelUtils.ssmToPdm(ssm)

      //recreate current state in new model
      val noPose = zsample.setShape(pchain.last._2.shape)
      val shapecoeff = amodel.coefficients(noPose.instance(cmodel))
      val ainit = pchain.last._2.setShape(shapecoeff)

      //getChain. see model, numSamples, burnin and init for differences
      val chain = fitPartial(amodel, target, mainAxis, rotPoint,
        (numSamples+burnin/subsampling)/realignments, iidNoiseSigma, shapeScale, shapeUpdateNormalScale, isoShapeScale, tScale, rScale, subsampling,
        0, Option(ainit))
      (chain.map(t => (t._1.instance(amodel), t._1, t._2)), amodel)
    }}.drop(1).flatMap(_._1.map(t => (t._1,t._3))).drop(burnin/subsampling)
  }

  /**
   * symmetric version of the informed proposal structure in
   * "A closest point proposal for MCMC-based probabilistic surface registration"
   * by Madsen et al. (2020)
   */
  def getSurfaceAwarePosterior(model: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], uncertaintyAlongNormal: Double, perpendicularFactor: Double = 1e4): PointDistributionModel[_3D, TriangleMesh] = {
    val alignIds = {
      val mapped = target.pointSet.pointsWithId.map(t => (t._2, model.mean.pointSet.findClosestPoint(t._1).id)).toIndexedSeq
      mapped.map(_._2.id).distinct.map(PointId)
    }
    val td = {
      val clps = alignIds.map(pid => {
        target.operations.closestPointOnSurface(model.mean.pointSet.point(pid))
      })
      val normals = Tchange.handleUcorSurface(target.vertexNormals, clps)
      val mvns = normals.map(n => {
        val nv = n.toBreezeVector
        val nvsvd = breeze.linalg.svd(nv * nv.t)
        val cov = nvsvd.leftVectors * diag(DenseVector.tabulate(3)(i => uncertaintyAlongNormal * (if (i > 0) perpendicularFactor else 1.0))) * nvsvd.rightVectors
        MultivariateNormalDistribution(DenseVector.zeros[Double](3), cov)
      })
      alignIds.zip(clps.zip(mvns)).map(t => {
        (t._1, t._2._1.point, t._2._2)
      })
    }
    model.posterior(td) //if cheaper required use the RegressionHelp.regressMean
  }

  /**
   * returns a coefficient proposal function that creates samples from the posterior in the model coefficient space.
   * useful to build more complex proposal functions
   */
  def getShapeProposal(model: PointDistributionModel[_3D, TriangleMesh], posterior: PointDistributionModel[_3D, TriangleMesh])(implicit rnd: Random): () => DenseVector[Double] = {
    //scaling the posterior up so the space is correctly represented and dividing the model space to calculate the coefficients in the original space
    val projMatrix = diag(model.gp.variance.map(d => 1.0 / math.sqrt(d))) * model.gp.basisMatrix.t * posterior.gp.basisMatrix * diag(posterior.gp.variance.map(math.sqrt))
    val f = () => projMatrix * DenseVector.tabulate[Double](model.rank)(_ => rnd.scalaRandom.nextGaussian())
    f
  }

  /**
   * uses the posterior to get adjusted proposals that better represent the target distribution.
   * additionally applies a function to the proposal scale that favors the dominant pcs to improve sampling.
   */
  def getShapeProposalGenerator(model: PointDistributionModel[_3D, TriangleMesh], posterior: PointDistributionModel[_3D, TriangleMesh], scale: IndexedSeq[(Double, Double)], decayFactor: Option[Int => Double] = None)(implicit rnd: Random): ProposalGenerator[Sample] with SymmetricTransitionRatio[Sample] = {
    val prop = getShapeProposal(model, posterior)
    val decay = decayFactor.getOrElse((_:Int)=>1.0)
    MixtureProposal.fromSymmetricProposals(
      scale.map(t => (t._1,
        new ProposalGenerator[Sample] with SymmetricTransitionRatio[Sample] {
          override def propose(current: Sample): Sample =
            current.changeShape(t._2 * (if (decayFactor.isDefined) prop.apply().mapPairs((i, d) => d * decay(i)) else prop.apply()))
          override def toString: String = s"${t._1} shapeProp sc ${t._2}"
        }
      )): _*
    )
  }

  /**
   * returns an isotropic shape proposal that is scaled to a certain domain.
   */
  def getScaledIsoGenerator(model: PointDistributionModel[_3D, TriangleMesh], noiseSigma: Double, alignIds: IndexedSeq[PointId], scale: IndexedSeq[(Double, Double)])(implicit rnd: Random): ProposalGenerator[Sample] with SymmetricTransitionRatio[Sample] = {
    val l2 = L2norm[_3D]()
    val totalv = model.gp.klBasis.map(kl => {
      val funcValues = alignIds.map(pid => kl.eigenfunction(pid))
      l2.norm2Vector(funcValues)*kl.eigenvalue
    }).sum
    val ratio =  math.sqrt(alignIds.length*noiseSigma*noiseSigma / totalv)
    MixtureProposal.fromSymmetricProposals(
      scale.map(t => (t._1,
        new ProposalGenerator[Sample] with SymmetricTransitionRatio[Sample] {
          override def propose(current: Sample): Sample =
            current.changeShape(t._2 * ratio * DenseVector.tabulate(current.shape.length)(_ => rnd.scalaRandom.nextGaussian()))
          override def toString: String = s"${t._1} shapeIsoProp sc ${t._2}"
        }
      )): _*
    )
  }

  /**
   * returns a function returning a translation vector. It is scaled with the number of points.
   * The general assumption is that there is one degree of freedom that increases the loss for each point.
   */
  def getTranslationProposal(numberPoints: Int)(implicit rnd: Random): () => EuclideanVector[_3D] = {
    val nsqrt = 1.0 / math.sqrt(numberPoints)
    val f = () => EuclideanVector3D(nsqrt * rnd.scalaRandom.nextGaussian(), nsqrt * rnd.scalaRandom.nextGaussian(), nsqrt * rnd.scalaRandom.nextGaussian())
    f
  }

  /**
   * returns a function returning rotation parameters. It is scaled depending on the eigenfunction.
   * The general assumption is that there is one degree of freedom that increases the loss for each point.
   */
  def getRotationProposal(model: PointDistributionModel[_3D, TriangleMesh], pids: IndexedSeq[PointId], mainAxis: EuclideanVector[_3D], rotPoint: Point[_3D], ignoreGeom: Boolean = false)(implicit rnd: Random): () => EuclideanVector[_3D] = {
    val ssm = ModelUtils.pdmToSsm(model)
    val orthv = MathHelp.getOrthogonalVectors(mainAxis)
    val axis = IndexedSeq(mainAxis, orthv._1, orthv._2)

    val std = axis.map(ax => {
      val eigf = ModelRecenter.getRf(ssm, ax, rotPoint)
      //find alignment to vertex normals
      val weights = if (!ignoreGeom) {
        pids.map(pid => {
          val efv = eigf(pid.id).normalize
          val vn = model.mean.vertexNormals(pid)
          val vnefvdot = vn.dot(efv)
          val res = math.abs(vnefvdot) //efv.norm - math.abs(vnefvdot)
          res
        })
      } else IndexedSeq.fill(pids.length)(1.0)
      assert(ignoreGeom || weights.sum > 1e-4, s"the geometry allows for a strong rotation along ${ax}")

      val angw = pids.zip(weights).map(t => {
        val pid = t._1
        val weight = t._2
        val p = model.mean.pointSet.point(pid)
        //distance to rotation axis
        val dist = (MathHelp.projectToPlane(p, rotPoint, ax) - rotPoint).norm
        //weighted angular speed
        eigf(pid.id).norm * weight / dist
      })
      val angwbar = angw.sum / angw.length
      //how much angular speed is required to reach sigma amount of deviation on average
      val ratio = math.sqrt(eigf.map(_.norm2).sum / pids.map(pid => eigf(pid.id).norm2).sum)

      angwbar*ratio
    })

    val f = () => EuclideanVector3D(std(0) * rnd.scalaRandom.nextGaussian(), std(1) * rnd.scalaRandom.nextGaussian(), std(2) * rnd.scalaRandom.nextGaussian())
    f
  }

  /**
   * uses the adjusted pose proposals to build a potentially mixture proposal distribution.
   * tscale = indexedseq(weight, scale), rscale = indexedseq(weight, scale)
   */
  def getPoseProposalGenerator(model: PointDistributionModel[_3D, TriangleMesh], alignIds: IndexedSeq[PointId], tscale: IndexedSeq[(Double, Double)], rscale: IndexedSeq[(Double, Double)], mainAxis: EuclideanVector[_3D], rotPoint: Point[_3D])(implicit rnd: Random): ProposalGenerator[Sample] with SymmetricTransitionRatio[Sample] = {
    val tprop = getTranslationProposal(alignIds.length)
    val rprop = getRotationProposal(model, alignIds, mainAxis, rotPoint)
    MixtureProposal.fromSymmetricProposals[Sample](
      (tscale.map(t => (t._1,
        new ProposalGenerator[Sample] with SymmetricTransitionRatio[Sample] {
          override def propose(current: Sample): Sample = current.changeT(tprop.apply().map(_ * t._2))
          override def toString: String = s"${t._1} tProp sc ${t._2}"
        }
      )) ++
        rscale.map(r => (r._1,
          new ProposalGenerator[Sample] with SymmetricTransitionRatio[Sample] {
            override def propose(current: Sample): Sample = current.changeR(rprop.apply().map(_ * r._2))
            override def toString: String = s"${r._1} rProp sc ${r._2}"
          }
        ))
        ): _*
    )
  }

  case class Sample(shape: DenseVector[Double], t: EuclideanVector[_3D], r: EuclideanVector[_3D], rp: Point[_3D], axis: IndexedSeq[EuclideanVector[_3D]]) {
    def calcInstance(model: PointDistributionModel[_3D, TriangleMesh]): TriangleMesh[_3D] = {
      val noPose = model.instance(shape)
      val rotated = if (r.norm > 1e-13) {
        val m = r.toArray.zip(axis).map(t => {
          SquareMatrix.apply[_3D](MathHelp.rotMatrixFromAxisAngle(t._2, t._1).data)
        }).reverse.reduce(_ * _)
        noPose.transform(Rotation3D(m, rp).apply)
      } else noPose

      val mesh = if (t.norm > 1e-13) rotated.transform(Translation3D.apply(t).apply) else rotated
      mesh
    }

    def instance(model: PointDistributionModel[_3D, TriangleMesh]): TriangleMesh[_3D] = calcInstance(model)

    def setShape(shape: DenseVector[Double]): Sample = this.copy(shape = shape)

    def changeShape(dshape: DenseVector[Double]): Sample = setShape(this.shape + dshape)

    def setT(t: EuclideanVector[_3D]): Sample = this.copy(t = t)

    def changeT(dt: EuclideanVector[_3D]): Sample = setT(this.t + dt)

    def setR(r: EuclideanVector[_3D]): Sample = this.copy(r = r)

    def changeR(dr: EuclideanVector[_3D]): Sample = setR(this.r + dr)
  }

  /**
   * targetsampling. cached
   */
  def getPostEval(sigma: Double, model: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D]): DistributionEvaluator[Sample] = {
    new DistributionEvaluator[Sample] {
      val l2Norm = L2norm[_3D]()
      val mem = Memoize(calcLogValue, 3)

      def calcLogValue(sample: Sample): Double = {
        val current = sample.instance(model)
        val ll = l2Norm.norm2Vector(Tchange.getDef(current, target)) / (2 * sigma * sigma)
        val lp = 0.5 * math.pow(norm(sample.shape), 2)
        -ll - lp //likelihood + prior
      }

      override def logValue(sample: Sample): Double = mem(sample)
    }
  }

  class CustomLogger[A](verbose:Boolean, trackAcc: Boolean, model: PointDistributionModel[_3D, TriangleMesh]) extends AcceptRejectLogger[A] {
    val accdata = scala.collection.mutable.HashMap[String, Int]()
    val rejdata = scala.collection.mutable.HashMap[String, Int]()

    override def accept(current: A, sample: A, generator: ProposalGenerator[A], evaluator: DistributionEvaluator[A]): Unit = {
      if (trackAcc) {
        val gs = generator.toString
        accdata.update(gs, accdata.getOrElse(gs, 0) + 1)
      }
      if (verbose) println(s"acc: ${evaluator.logValue(current) - evaluator.logValue(sample)}")
    }

    override def reject(current: A, sample: A, generator: ProposalGenerator[A], evaluator: DistributionEvaluator[A]): Unit = {
      if (trackAcc) {
        val gs = generator.toString
        rejdata.update(gs, rejdata.getOrElse(gs, 0) + 1)
      }
      if (verbose) println(s"rej: ${evaluator.logValue(current) - evaluator.logValue(sample)}")
    }

    override def toString: String = if (trackAcc) {
      accdata.map(t => {
        s"${t._1}: acc ${t._2} / ${t._2 + rejdata.getOrElse(t._1, 0)}\n"
      }).reduce(_ + _)
    } else "custom logger - activate trackAcc for more detailed info"
  }

}
