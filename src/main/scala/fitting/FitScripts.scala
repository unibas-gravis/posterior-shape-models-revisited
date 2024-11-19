package fitting

import breeze.linalg.{DenseMatrix, DenseVector, diag, norm}
import norms.L2norm
import scalismo.common.PointId
import scalismo.geometry._
import scalismo.mesh.TriangleMesh
import scalismo.sampling.algorithms.MetropolisHastings
import scalismo.sampling.loggers.AcceptRejectLogger
import scalismo.sampling.proposals.MixtureProposal
import scalismo.sampling.{DistributionEvaluator, ProposalGenerator, SymmetricTransitionRatio}
import scalismo.statisticalmodel.{MultivariateNormalDistribution, PointDistributionModel}
import scalismo.transformations._
import scalismo.utils.{Memoize, Random}
import scalismo.{ModelRecenter, ModelUtils, RegressionHelp}
import utility.{MathHelp, Tchange}

object FitScripts {

  /**
   * optimizing over shape parameters with free correspondence. pose optimization uses estimated pose posterior.
   * uses fixed targetSampling
   */
  def fitPartial(model: PointDistributionModel[_3D, TriangleMesh],
                 target: TriangleMesh[_3D],
                 mainAxis: EuclideanVector[_3D],
                 rotPoint: Point[_3D],
                 numSamples: Int = 2000, //number of final samples -> mh samples numSamples * subsampling
                 iidNoiseSigma: Double = 1.0,
                 shapeScale: Double = 1.0,
                 perpendicularFactor: Double = 1.0,
                 isoShapeScale: Double = 1.0,
                 tScale: Double = 1.0,
                 rScale: Double = 1.0,
                 subsampling: Int = 5,
                 burnin: Int = 100,
                 initial: Option[Sample] = None
                )(implicit rnd: Random): IndexedSeq[(Sample, Double)] = {
    //TODO apply pose of initial to these initial calculations
    //TODO correspondence in this posterior calc and the isoGenerator below can be unified
    //val temp = TestPartialAlignment.ui.show(target, "target")
    val postModel = getSurfaceAwarePosterior(model, target, TargetSamplingUnique(), varianceAlongNormal = iidNoiseSigma * iidNoiseSigma, perpendicularFactor = perpendicularFactor, asSymmetricProposal = true)
    val axis = MathHelp.listOrthogonalVectors(mainAxis)
    val stscale = IndexedSeq(0.2, 0.6, 0.2).zip(IndexedSeq(0.5, 1.0, 2.0))
    val alignIds = {
      val mapped = target.pointSet.pointsWithId.map(t => (t._2, model.mean.pointSet.findClosestPoint(t._1).id)).toIndexedSeq
      mapped.map(_._2.id).distinct.map(PointId)
    }

    val shapeProposal = getShapeProposalGenerator(model, postModel, stscale.map(t => (t._1, t._2 * shapeScale)))
    val isoProposal = getScaledIsoGenerator(model, iidNoiseSigma, alignIds, stscale.map(t => (t._1, t._2 * isoShapeScale)))
    val poseProposal = {
      val (tScales, rScales) = stscale.map(t => ((t._1, t._2 * tScale * iidNoiseSigma), (t._1, t._2 * rScale * iidNoiseSigma))).unzip
      getPoseProposalGenerator(model, alignIds, tScales, rScales, mainAxis, rotPoint)
    }
    val proposals = MixtureProposal.fromSymmetricProposals(IndexedSeq((0.5, shapeProposal), (0.2, isoProposal), (0.3, poseProposal)): _*)
    val eval = getPostEval(iidNoiseSigma, model, target)

    val mh = MetropolisHastings.apply(proposals, eval)
    val trackAcc = true
    val logger = new CustomLogger[Sample](false, trackAcc, model)
    val init = initial.getOrElse(Sample(DenseVector.zeros[Double](model.rank), EuclideanVector3D.zero, EuclideanVector3D.zero, rotPoint, axis))
    //    val temp0 = (0 to 10).map(_ => poseProposal.asInstanceOf[MixtureProposal[Sample]].generators(4).propose(init))
    //    val temp = temp0.map(t0 => Tchange.getDef(t0.instance(model).pointSet, model.mean))
    //    val temp2 = temp.map(t => t.map(_.norm2).sum)
    val chain = mh.iterator(init, logger).drop(burnin).zipWithIndex.filter(_._2 % subsampling == 0).map(t => (t._1, eval.logValue(t._1))).take(numSamples)

    val resChain = chain.toIndexedSeq
    //    temp.remove()
    //    logger.last.get.remove()
    if (trackAcc) println(s"logger ${logger.toString}")
    resChain
  }

  /**
   * optimizing over shape parameters with free correspondence. pose optimization uses estimated pose posterior.
   * allows to have arbitrary sampling strategies
   */
  def fitPartialStrat(model: PointDistributionModel[_3D, TriangleMesh],
                      target: TriangleMesh[_3D],
                      mainAxis: EuclideanVector[_3D],
                      rotPoint: Point[_3D],
                      numSamples: Int = 2000, //number of final samples -> mh samples numSamples * subsampling
                      iidNoiseSigma: Double = 1.0,
                      shapeScale: Double = 1.0,
                      perpendicularFactor: Double = 1.0,
                      isoShapeScale: Double = 1.0,
                      tScale: Double = 1.0,
                      rScale: Double = 1.0,
                      subsampling: Int = 5,
                      burnin: Int = 100,
                      mixtures: (Double,Double,Double) = (0.5, 0.2, 0.3),
                      initial: Option[Sample] = None,
                      ldms: IndexedSeq[(PointId, Point[_3D])] = IndexedSeq.empty,
                      ldmSigmaFactor: Double = 1.0,
                      backgroundll: Double = 0.0,
                      samplingStrategy: SamplingStrategy = BidirectionalSamplingFromTarget() //expensive but stable sampling strategy
                     )(implicit rnd: Random): IndexedSeq[(Sample, Double)] = {
    assert(backgroundll <= 0.0, "backgroundll is in log likelihood format and shoulb be <=0")
    //TODO apply pose of initial to these initial calculations
    val postModel = getSurfaceAwarePosterior(model, target, samplingStrategy, varianceAlongNormal = iidNoiseSigma * iidNoiseSigma, perpendicularFactor = perpendicularFactor, asSymmetricProposal = true)
    val axis = MathHelp.listOrthogonalVectors(mainAxis)
    val stscale = IndexedSeq(0.2, 0.6, 0.2).zip(IndexedSeq(0.5, 1.0, 2.0))
    val alignIds = {
      val mapped = target.pointSet.pointsWithId.map(t => (t._2, model.mean.pointSet.findClosestPoint(t._1).id)).toIndexedSeq
      mapped.map(_._2.id).distinct.map(PointId)
    }

    val shapeProposal = getShapeProposalGenerator(model, postModel, stscale.map(t => (t._1, t._2 * shapeScale)))
    val isoProposal = getScaledIsoGenerator(model, iidNoiseSigma, alignIds, stscale.map(t => (t._1, t._2 * isoShapeScale)))
    val poseProposal = {
      val (tScales, rScales) = stscale.map(t => ((t._1, t._2 * tScale * iidNoiseSigma), (t._1, t._2 * rScale * iidNoiseSigma))).unzip
      getPoseProposalGenerator(model, alignIds, tScales, rScales, mainAxis, rotPoint)
    }
    val proposals = MixtureProposal.fromSymmetricProposals(IndexedSeq((mixtures._1, shapeProposal), (mixtures._2, isoProposal), (mixtures._3, poseProposal)): _*)
    val eval = getPostEvalStrat(iidNoiseSigma, model, target, ldms, ldmSigmaFactor * iidNoiseSigma, backgroundll, samplingStrategy)

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
                              perpendicularFactor: Double = 1.0,
                              isoShapeScale: Double = 1.0,
                              tScale: Double = 1.0,
                              rScale: Double = 1.0,
                              subsampling: Int = 5,
                              burnin: Int = 100,
                              mixtures: (Double,Double,Double) = (0.5, 0.2, 0.3),
                              initial: Option[Sample] = None,
                              ldms: IndexedSeq[(PointId, Point[_3D])] = IndexedSeq.empty,
                              ldmSigmaFactor: Double = 1.0,
                              realignments: Int = 5,
                              realignmentRotation: Boolean = true,
                              backgroundll: Double = 0.0,
                              samplingStrategy: SamplingStrategy = BidirectionalSamplingFromTarget(), //this sampling strategy is quite expensive but safe, TargetSamplingUnique is cheaper option
                              correspondenceStrategy: SamplingStrategyUniform = BidirectionalSamplingFromOrigin(),
                             )(implicit rnd: Random): IndexedSeq[(TriangleMesh[_3D], Double)] = {
    val axis = MathHelp.listOrthogonalVectors(mainAxis)
    val zsample = Sample(DenseVector.zeros[Double](model.rank), EuclideanVector3D.zero, EuclideanVector3D.zero, rotPoint, axis)
    val init = initial.getOrElse(zsample)
    (0 until realignments).scanLeft((IndexedSeq((init.instance(model), init, Double.MinValue)), model)) { case ((pchain, cmodel), _) => {
      //create aligned model
      val clps = correspondenceStrategy.establishCorrespondenceUniform(pchain.last._1, target)
      val alignids = clps.map(_._1.id).distinct.map(PointId)
      val recentered = ModelRecenter.recenterSsm(ModelUtils.pdmToSsm(model), alignids)
      val ssm = if (realignmentRotation && alignids.size >= 3) ModelRecenter.rerotate(recentered, alignids, axis, Option(rotPoint)) else recentered
      val amodel = ModelUtils.ssmToPdm(ssm)

      //recreate current state in new model
      val noPose = zsample.setShape(pchain.last._2.shape)
      val shapecoeff = amodel.coefficients(noPose.instance(cmodel))
      val ainit = pchain.last._2.setShape(shapecoeff)

      //getChain. see model, numSamples, burnin and init for differences
      val chain = fitPartialStrat(amodel, target, mainAxis, rotPoint,
        (numSamples + burnin / subsampling) / realignments, iidNoiseSigma, shapeScale, perpendicularFactor, isoShapeScale, tScale, rScale, subsampling,
        0, mixtures, Option(ainit), ldms, ldmSigmaFactor, backgroundll, samplingStrategy)
      (chain.map(t => (t._1.instance(amodel), t._1, t._2)), amodel)
    }
    }.drop(1).flatMap(_._1.map(t => (t._1, t._3))).drop(burnin / subsampling)
  }

  /**
   * an approximation of Madsen, Dennis, et al. "A closest point proposal for MCMC-based probabilistic surface registration." Computer Vision–ECCV 2020
   * used for posteriors better suited for surfaces which allow shifting correspondence. this is intended to provide
   * the posterior implied by meshes. if asSymmetricProposal is true then it will only observe
   */
  def getSurfaceAwarePosterior(model: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], samplingStrategy: SamplingStrategy, startingMesh: Option[TriangleMesh[_3D]]=None, varianceAlongNormal: Double, perpendicularFactor: Double = 1e4, asSymmetricProposal:Boolean = false): PointDistributionModel[_3D, TriangleMesh] = {
    val mesh = startingMesh.getOrElse(model.mean)
    samplingStrategy match {
      case samplingStrategy: SamplingStrategyNormals =>
        val corrPoints = samplingStrategy.establishCorrespondenceNormal(mesh, target)
        val td = getSurfaceAwareTd(model, corrPoints, varianceAlongNormal, perpendicularFactor, asSymmetricProposal)
        model.posterior(td)
      case samplingStrategy: SamplingStrategy =>
        val normExtended = NormalSamplingSimpleExtension(samplingStrategy)
        val corrPoints = normExtended.establishCorrespondenceNormal(mesh, target)
        val td = getSurfaceAwareTd(model, corrPoints, varianceAlongNormal, perpendicularFactor, asSymmetricProposal)
        model.posterior(td)
    }
  }

  /**
   * similar to getSurfaceAwarePosterior but uses given correspondence and each observation uncertainty can be modified
   * by additional scalar given.
   * if meanObs is used -> then only the mean is observed. useful for symmetric proposals
   */
  def getSurfaceAwareTd(model: PointDistributionModel[_3D, TriangleMesh], obs: IndexedSeq[(PointId, Point[_3D], EuclideanVector[_3D], Double)], uncertaintyAlongNormal: Double, perpendicularFactor: Double = 1e4, meanObs: Boolean=false): IndexedSeq[(PointId, Point[_3D], MultivariateNormalDistribution)] = {
    //TODO there is sometimes a numerical issue here -> compare normal posterior with surface aware posterior
    val td = {
      val mvns = obs.map(t => {
        val n = t._3
        val nv = n.toBreezeVector
        val cov = try {
          val nvsvd = breeze.linalg.svd(nv * nv.t) //+ ddiag)
          nvsvd.leftVectors * diag(DenseVector.tabulate(3)(i => t._4 * uncertaintyAlongNormal * (if (i > 0) perpendicularFactor else 1.0))) * nvsvd.rightVectors
        } catch {
          //nan values in nv replace with normal iid
          case _: breeze.linalg.NotConvergedException => t._4 * uncertaintyAlongNormal * perpendicularFactor * DenseMatrix.eye[Double](3)
        }
        MultivariateNormalDistribution(DenseVector.zeros[Double](3), cov)
      })
      obs.zip(mvns).map(t => {
        (t._1._1, if (meanObs) model.mean.pointSet.point(t._1._1) else t._1._2, t._2)
      })
    }
    td
  }

  /**
   * returns a coefficient proposal function that creates samples from the posterior in the model coefficient space.
   * useful to build more complex proposal functions
   */
  def getShapeProposal(model: PointDistributionModel[_3D, TriangleMesh], posterior: PointDistributionModel[_3D, TriangleMesh])(implicit rnd: Random): () => DenseVector[Double] = {
    //scaling the posterior up so the space is correctly represented and dividing the model space to calculate the coefficients in the original space
    val projMatrix = ModelUtils.getRootMatrixMappingSpaces(model, posterior)
    val f = () => projMatrix * DenseVector.tabulate[Double](model.rank)(_ => rnd.scalaRandom.nextGaussian())
    f
  }

  /**
   * uses the posterior to get adjusted proposals that better represent the target distribution.
   * additionally applies a function to the proposal scale that favors the dominant pcs to improve sampling.
   */
  def getShapeProposalGenerator(model: PointDistributionModel[_3D, TriangleMesh], posterior: PointDistributionModel[_3D, TriangleMesh], scale: IndexedSeq[(Double, Double)], decayFactor: Option[Int => Double] = None)(implicit rnd: Random): ProposalGenerator[Sample] with SymmetricTransitionRatio[Sample] = {
    val prop = getShapeProposal(model, posterior)
    val decay = decayFactor.getOrElse((_: Int) => 1.0)
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
      l2.norm2Vector(funcValues) * kl.eigenvalue
    }).sum
    val ratio = math.sqrt(alignIds.length * noiseSigma * noiseSigma / totalv)
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
      val eigf = RegressionHelp.getRf(ssm, ax, rotPoint)
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

      angwbar * ratio
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
      val transformation = getTransform()
      noPose.transform(transformation)
    }

    def calcLdmInstance(model: PointDistributionModel[_3D, TriangleMesh], ldm: IndexedSeq[PointId]): IndexedSeq[Point[_3D]] = {
      val points = ldm.map(pid => model.mean.pointSet.point(pid) + model.gp.instanceAtPoint(shape, pid))
      val transformation = getTransform()
      points.map(transformation.f)
    }

    def getTransform(): RigidTransformation[_3D] = {
      val rotation: Rotation[_3D] = if (r.norm > 1e-13) {
        //        r.toArray.zip(axis).map(t => { //not so efficient, but clearer
        //          Rotation3D(SquareMatrix.apply[_3D](MathHelp.rotMatrixFromAxisAngle(t._2, t._1).data), rp)
        //        }).foldLeft(noPose)((mesh, r) => mesh.transform(r.apply))
        val m = r.toArray.zip(axis).map(t => {
          SquareMatrix.apply[_3D](MathHelp.rotMatrixFromAxisAngle(t._2, t._1).data)
        }).reverse.reduce(_ * _)
        Rotation3D(m, rp)
      } else RotationSpace3D(Point3D.origin).identityTransformation //Identity transformation
      val translation: Translation[_3D] = if (t.norm > 1e-13) scalismo.transformations.Translation(t) else TranslationSpace3D.identityTransformation //Identity transformation
      TranslationAfterRotation3D(translation, rotation)
    }

    def instance(model: PointDistributionModel[_3D, TriangleMesh]): TriangleMesh[_3D] = calcInstance(model)

    def ldmInstance(model: PointDistributionModel[_3D, TriangleMesh], ldm: IndexedSeq[PointId]): IndexedSeq[Point[_3D]] = calcLdmInstance(model, ldm)

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

  /**
   * arbitrary correspondence strategy. cached
   */
  def getPostEvalStrat(sigma: Double, model: PointDistributionModel[_3D, TriangleMesh], target: TriangleMesh[_3D], ldms: IndexedSeq[(PointId, Point[_3D])], sigmaLdms: Double, backgroundll: Double, strat: SamplingStrategy): DistributionEvaluator[Sample] = {
    new DistributionEvaluator[Sample] {
      val l2Norm = L2norm[_3D]()
      val mem = Memoize(calcLogValue, 3)

      def calcLogValue(sample: Sample): Double = {
        val current = sample.instance(model)
        val (dvec, bkgn) = {
          val (corr, bkgn) = strat.establishCorrespondenceBackground(current, target)
          (corr.map(t => (current.pointSet.point(t._1) - t._2).*(t._3)), bkgn)
        }

        val ll = l2Norm.norm2Vector(dvec) / (2.0 * sigma * sigma)
        val lldm = if (ldms.nonEmpty) {
          val dvecldm = ldms.map(t => current.pointSet.point(t._1) - t._2)
          l2Norm.norm2Vector(dvecldm) / (2.0 * sigmaLdms * sigmaLdms)
        } else 0.0
        val lbkg = bkgn * backgroundll
        val lp = 0.5 * math.pow(norm(sample.shape), 2)
        val total = -ll - lp - lldm + lbkg    //likelihood + prior + landmarks + background
        total
      }

      override def logValue(sample: Sample): Double = mem(sample)
    }
  }

  class CustomLogger[A](verbose: Boolean, trackAcc: Boolean, model: PointDistributionModel[_3D, TriangleMesh]) extends AcceptRejectLogger[A] {
    //var last: Option[ShowInScene.ShowInSceneMesh.View] = None
    val accdata = scala.collection.mutable.HashMap[String, Int]()
    val rejdata = scala.collection.mutable.HashMap[String, Int]()

    override def accept(current: A, sample: A, generator: ProposalGenerator[A], evaluator: DistributionEvaluator[A]): Unit = {
      /*
      sample match {
        case s: Sample =>
          val mesh = s.calcInstance(model)
          val temp = last
          last = Option(TestPartialAlignment.ui.show(mesh, "newState"))
          if (temp.isDefined) temp.get.remove()
      }
      */
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
