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
	// val isig_file = "./src/test/scala/fir/I_side.csv"
	// val qsig_file = "./src/test/scala/fir/Q_side.csv"
	// val writer = new PrintWriter(new File("isig_side.csv"))
	// val isig_file = "./src/test/scala/fir/I_image.csv"
	// val qsig_file = "./src/test/scala/fir/Q_image.csv"
	// val writer = new PrintWriter(new File("isig_image.csv"))
	val isig_file = "./src/test/scala/fir/I_both.csv"
	val qsig_file = "./src/test/scala/fir/Q_both.csv"
	val writer = new PrintWriter(new File("isig_both.csv"))
	val isig = Source.fromFile(isig_file).getLines().next.split(",").map(_.trim).map(_.toInt).toSeq
	val qsig = Source.fromFile(qsig_file).getLines().next.split(",").map(_.trim).map(_.toInt).toSeq
	assert(isig.length == qsig.length)
	// Error message: Something's not right about your input vectors
	var isig_out = Seq[Int]()


	reset(5)
	for(v <- 1 until isig.length) {
		poke(c.io.isig, isig(v))
		poke(c.io.qsig, qsig(v))
		step(1)
		isig_out = isig_out :+ peek(c.io.isig_out).toInt
	}
	// 49-tap means 49 delays, methinks.
	for(v <- 1 until 49) {
		step(1)
		isig_out = isig_out :+ peek(c.io.isig_out).toInt		
	}
	

	for(v <- isig_out) {
		writer.write(v.toString + ",")
	}
	writer.close()


}


class IRCSfilterTesterSpec extends FreeSpec with Matchers {
	val rcoeffs_file = "./src/test/scala/fir/rcoeffs.csv"
	val icoeffs_file = "./src/test/scala/fir/icoeffs.csv"
	val rcoeffs:Seq[SInt] = Source.fromFile(rcoeffs_file).getLines().next.split(",").map(_.trim).map(_.toInt).map(_.asSInt).toSeq
	val icoeffs:Seq[SInt] = Source.fromFile(icoeffs_file).getLines().next.split(",").map(_.trim).map(_.toInt).map(_.asSInt).toSeq
	"tester should filter things" in {
		iotesters.Driver.execute(Array("--backend-name", "firrtl", "--target-dir", "test_run_dir", "--fint-write-vcd"), () => new IRCSfilter(length = rcoeffs.length, inputwidth = 5, filterwidth = 8, import_rcoeffs = rcoeffs, import_icoeffs = icoeffs)) { c =>
			new IRCSfilterTester(c)
		} should be (true)
	}
}