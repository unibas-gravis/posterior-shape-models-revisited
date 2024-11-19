package io.jsonHandler

import spray.json.DefaultJsonProtocol._
import spray.json._

import java.io.File
import scala.io.Source

class Instructions(val idGroups: Seq[Set[Int]], val groups: Seq[ContourGroup]) {

}

object Instructions {
  def apply(f: File): Instructions = {
    val source = Source.fromFile(f)
    val instructions = source.mkString.parseJson.convertTo[Map[String,JsValue]]
    source.close()
    val files = instructions("idxFiles").convertTo[Seq[String]].map(idxf => new File(f.getParent+"/"+idxf))
    val idGroups = files.map(getIds).map(_.toSet)
    val groups = instructions("contours").convertTo[Seq[JsValue]].map(ContourGroup(_))
    new Instructions(idGroups,groups)
  }

  private def getIds(f: File): Seq[Int] = {
    val source = Source.fromFile(f)
    val result = source.mkString.parseJson.convertTo[Map[String,Seq[Int]]].apply("pointIds")
    source.close()
    result
  }
}

