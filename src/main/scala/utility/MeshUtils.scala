package utility

import scalismo.common.{DiscreteDomain, DiscreteField, PointId, UnstructuredPoints, UnstructuredPoints3D, UnstructuredPointsDomain, UnstructuredPointsDomain3D}
import scalismo.geometry.{EuclideanVector, EuclideanVector3D, NDSpace, Point, _3D}
import scalismo.mesh.{TriangleCell, TriangleList, TriangleMesh, TriangleMesh3D}

object MeshUtils {

  /**
   * sorts points based on distance to growP
   */
  def getPartialMeshPoint(mesh: TriangleMesh[_3D], growP: Point[_3D], ratio: Double, ratioForInterval: Boolean = true): (TriangleMesh[_3D], IndexedSeq[PointId], PointId => Option[PointId]) = {
    cutTo(mesh, (p,_) => (p-growP).norm2, ratio, ratioForInterval)
  }

  /**
   * sorts points based on dot product with growV
   */
  def getPartialMeshVector(mesh: TriangleMesh[_3D], growV: EuclideanVector[_3D], ratio: Double, ratioForInterval: Boolean = true): (TriangleMesh[_3D], IndexedSeq[PointId], PointId => Option[PointId]) = {
    cutTo(mesh, (p,_) => p.toVector.dot(growV), ratio, ratioForInterval)
  }

  /**
   * cuts given mesh down to provided point ids
   */
  def makePartial(mesh: TriangleMesh[_3D], obsset: Set[Int]): TriangleMesh[_3D] =
    mesh.operations.maskPoints(pid => obsset.contains(pid.id)).transformedMesh

  /**
   * returns (new mesh, old ids that are kept, a lazy mapping from old to new ids)
   * ratioForInterval if true ratio is relative to length; if false number of points
   */
  def cutTo(mesh: TriangleMesh[_3D], order: (Point[_3D],PointId)=>Double, ratio: Double, ratioForInterval: Boolean): (TriangleMesh[_3D], IndexedSeq[PointId], PointId => Option[PointId]) = {
    val ids = if (ratioForInterval) {
      val ordered = mesh.pointSet.pointsWithId.toIndexedSeq.map(t => (t, order(t._1,t._2))).sortBy(_._2)
      val h = (ordered.last._2-ordered.head._2)*ratio+ordered.head._2
      ordered.takeWhile(t => t._2 < h).map(_._1._2)
    } else {
      mesh.pointSet.pointsWithId.toIndexedSeq.sortBy(t => order(t._1,t._2))
        .take((mesh.pointSet.numberOfPoints*ratio).toInt).map(_._2)
    }
    val idset = ids.map(_.id).toSet
    val tm = mesh.operations.maskPoints(pid => idset.contains(pid.id)).transformedMesh
    (tm, ids, pid => {
      if (idset.contains(pid.id)) {
        Option(tm.pointSet.findClosestPoint(mesh.pointSet.point(pid)).id)
      } else None
    })
  }

  /**
   * transforms the points to a mesh. either as cubes or spheres
   */
  def pointsToMesh(points: IndexedSeq[Point[_3D]], radius: Double, cubes: Boolean = true): TriangleMesh[_3D] = {
    if (cubes){
      val (ps, triangles) = points.zipWithIndex.map(t => {
        val cubePoints = IndexedSeq(
          EuclideanVector3D(1,1,1),EuclideanVector3D(-1,1,1),EuclideanVector3D(1,-1,1),EuclideanVector3D(1,1,-1),
          EuclideanVector3D(-1,-1,1),EuclideanVector3D(-1,1,-1),EuclideanVector3D(1,-1,-1),EuclideanVector3D(-1,-1,-1)
        ).map(e => t._1 + e*radius)
        val faces = IndexedSeq(
          IndexedSeq(0,1,4,2),IndexedSeq(0,1,5,3),IndexedSeq(0,2,6,3), //faces connected to 111
          IndexedSeq(7,6,3,5),IndexedSeq(7,6,2,4),IndexedSeq(7,5,1,4) //faces connected to -1-1-1
        )
        val indices = faces.flatMap(list => {
          val oneSided = IndexedSeq(IndexedSeq(list(0),list(1),list(2)), IndexedSeq(list(0),list(3),list(2))) //one side
          oneSided.flatMap(one => IndexedSeq(one,IndexedSeq(one(0),one(2),one(1)))) //double sided
        })
        val trianglelist = indices.map(_.map(i => PointId(t._2*8+i))).map(tr => TriangleCell(tr(0),tr(1),tr(2)))
        (cubePoints, trianglelist)
      }).unzip
      TriangleMesh3D(ps.flatten, TriangleList(triangles.flatten))
    } else {
      throw new NotImplementedError("not implemented for spheres")
    }
  }

  /**
   * transforms a discreteVectorField to a mesh. TODO add length modifier
   */
  def fieldToMesh(vecField: DiscreteField[_3D, UnstructuredPointsDomain, EuclideanVector[_3D]], thicknessMod: Double = 0.05, headMod: Double = 2.0): TriangleMesh[_3D] = {
    val plusVec: (IndexedSeq[EuclideanVector[_3D]],IndexedSeq[IndexedSeq[Int]]) = {
      val hs = 0.62 //apprx 1/goldenRatio
      val cubePoints = IndexedSeq(
        EuclideanVector3D(1,1,hs),EuclideanVector3D(-1,1,hs),EuclideanVector3D(1,-1,hs),EuclideanVector3D(1,1,0),
        EuclideanVector3D(-1,-1,hs),EuclideanVector3D(-1,1,0),EuclideanVector3D(1,-1,0),EuclideanVector3D(-1,-1,0)
      )
      val arrowPoints = IndexedSeq(
        EuclideanVector3D(headMod,headMod,hs),EuclideanVector3D(-headMod,headMod,hs),EuclideanVector3D(headMod,-headMod,hs),
        EuclideanVector3D(-headMod,-headMod,hs), EuclideanVector3D(0,0,1)
      )
      val allPoints = (cubePoints++arrowPoints).map(v => EuclideanVector3D(v.x*thicknessMod,v.y*thicknessMod,v.z))
      val faces = IndexedSeq( //cube with one side open and then a pyramid
        IndexedSeq(0,1,4,2),IndexedSeq(0,1,5,3),IndexedSeq(0,2,6,3), //faces connected to 111
        IndexedSeq(7,6,3,5),IndexedSeq(7,6,2,4),IndexedSeq(7,5,1,4), //faces connected to -1-1-1
        IndexedSeq(8,9,11,10), //base of pyramid
        IndexedSeq(8,9,12),IndexedSeq(8,10,12),IndexedSeq(9,11,12),IndexedSeq(10,11,12) //sides of pyramid
      )
      val indices = faces.flatMap(list => {
        val oneSided = if (list.length == 3) IndexedSeq(list) else {
          IndexedSeq(IndexedSeq(list(0), list(1), list(2)), IndexedSeq(list(0), list(3), list(2))) //one side
        }
        oneSided.flatMap(one => IndexedSeq(one, IndexedSeq(one(0), one(2), one(1)))) //double sided
      })
      (allPoints,indices)
    }

    val (arrows, triangleLists) = vecField.domain.pointSet.pointIds.map(pid => {
      val point = vecField.domain.pointSet.point(pid)
      val vector = vecField(pid)
      require(math.abs(vector.x) + math.abs(vector.y) < 1e-10, "only implemented for z vectors")
      val goodPoints = plusVec._1.map(v => point+v*vector.z)
      val goodTrianglelist = plusVec._2.map(_.map(i => PointId(pid.id*13+i))).map(tr => TriangleCell(tr(0),tr(1),tr(2)))
      (goodPoints, goodTrianglelist)
    }).toIndexedSeq.unzip
    TriangleMesh3D(arrows.flatten, TriangleList(triangleLists.flatten))
  }

  /**
   * can be used to do basic deformations (like scaling) to meshes
   * if no point is provided the mean point of the mesh is used
   * pnew = f(pold-point) + point
   * internally:
   * pnew = pold + f(pold - point) - (pold - point) -> to only scale z by factor of 2 => f(x,y,z)=(x,y,2z)
   * TODO should scaling be changed to f(x,y,z) = (0,0,z) -> remove implicit -1?
   */
  def deformMesh(mesh: TriangleMesh[_3D], f: EuclideanVector[_3D]=>EuclideanVector[_3D], point: Option[Point[_3D]] = None): TriangleMesh[_3D] = {
    val p = point.getOrElse(Tchange.getMean(mesh))
    Tchange.meshPlusEv(mesh, mesh.pointSet.pointSequence.map(ap => ap-p).map(v => f(v)-v))
  }

}
