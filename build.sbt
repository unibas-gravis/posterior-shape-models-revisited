organization  := "ch.unibas.cs.gravis"

name := "eigenfunctionInspection"

version       := "0.1"

scalaVersion  := "2.12.8"

scalacOptions := Seq("-unchecked", "-deprecation", "-encoding", "utf8")

resolvers += Resolver.bintrayRepo("unibas-gravis", "maven")


libraryDependencies ++= Seq(
  "ch.unibas.cs.gravis" %% "scalismo-ui" % "0.90.0",
  "ch.unibas.cs.gravis" %% "scalismo-faces" % "0.90.0",
  "ch.qos.logback" % "logback-classic" % "1.2.3",
  "com.typesafe.scala-logging" %% "scala-logging" % "3.9.0"
).map(_.force())

libraryDependencies ~= { _.map(_.exclude("org.slf4j", "slf4j-nop")) }

fork := true

javaOptions += "-Djava.awt.headless=true"
javaOptions += "-Xmx72g"
//javaOptions += "-Xms2g"
//javaOptions += "-Xmn2g"
//javaOptions += "-Xss4g"
javaOptions += "-XX:+UseParallelGC"
//javaOptions += "-XX:+HeapDumpOnOutOfMemoryError"


