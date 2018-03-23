package fir

import chisel3._
import chisel3.util._
import chisel3.iotesters._
import chisel3.iotesters.{ChiselFlatSpec, PeekPokeTester, SteppedHWIOTester}
import org.scalatest.{FreeSpec, Matchers}
import java.io.File
import java.io.PrintWriter
import scala.io.Source

class IRCSfilterTester(c: IRCSfilter) extends PeekPokeTester(c) {
	val isig_file = "./src/test/scala/fir/I_out.csv"
	val qsig_file = "./src/test/scala/fir/Q_out.csv"
	val isig = Source.fromFile(isig_file).getLines().next.split(",").map.(_.trim).map(_.toInt).toSeq
	val qsig = Source.fromFile(qsig_file).getLines().next.split(",").map.(_.trim).map(_.toInt).toSeq
	assert(isig.length == qsig.length)
	// Error message: Something's not right about your input vectors
	isig_out = Seq[Int]()


	reset(5)
	// +47 because I'm expecting 47 delays
	for(v <- 1 until isig.length+47) {
		poke(c.io.isig, isig(v))
		poke(c.io.qsig, qsit(v))
		step(1)
		isig_out = isig_out :+ peek(c.io.isig_out)
	}
	
	val writer = new PrintWriter(new File("isig_out.csv"))

	for(v <- isig_out) {
		writer.write(v.toString + ", ")
	}
	writer.close()


}


class IRCSfilterTesterSpec extends FreeSpec with Matchers {
	val rcoeffs_file = "./src/test/scala/fir/rcoeffs.csv"
	val icoeffs_file = "./src/test/scala/fir/icoeffs.csv"
	val rcoeffs:Seq[SInt] = Source.fromFile(rcoeffs_file).getLines().next.split(",").map(_.trim).map(_.toInt).map(_.asSint).toSeq
	val icoeffs:Seq[SInt] = Source.fromFile(icoeffs_file).getLines().next.split(",").map(_.trim).map(_.toInt).map(_.asSint).toSeq
	"tester should filter things" in {
		iotesters.Driver.execute(Array("--backend-name", "firrtl", "--target-dir", "test_run_dir", "--fint-write-vcd"), () => new IRCSfilter(length = rcoeffs.length, inputwidth = 5, filterwidth = 8, import_rcoeffs = rcoeffs, import_icoeffs = icoeffs)) { c =>
			new IRCSfilterTester(c)
		} should be (true)
	}
}