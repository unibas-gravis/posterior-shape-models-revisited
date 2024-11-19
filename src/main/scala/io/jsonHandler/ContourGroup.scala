package io.jsonHandler

import spray.json.DefaultJsonProtocol._
import spray.json._

class ContourGroup(js: JsValue) {
  private val info = js.convertTo[Map[String,JsValue]]
  val always = info("alwaysContour").convertTo[Boolean]
  val groups = info("groups").convertTo[Seq[Int]]
  val onto = info("onto").convertTo[Seq[Int]]
}

object ContourGroup {
  def apply(js: JsValue): ContourGroup = new ContourGroup(js)
}
