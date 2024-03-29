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
	val isig_file = "./src/test/scala/fir/test/I_out_withImage.csv"
	val qsig_file = "./src/test/scala/fir/test/Q_out_withImage.csv"
	val writer = new PrintWriter(new File("./src/test/scala/fir/test/I_out_noImage.csv"))

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
	val rcoeffs_file = "./src/test/scala/fir/test/rcoeffs.csv"
	val icoeffs_file = "./src/test/scala/fir/test/icoeffs.csv"
	val rcoeffs:Seq[SInt] = Source.fromFile(rcoeffs_file).getLines().next.split(",").map(_.trim).map(_.toInt).map(_.asSInt).toSeq
	val icoeffs:Seq[SInt] = Source.fromFile(icoeffs_file).getLines().next.split(",").map(_.trim).map(_.toInt).map(_.asSInt).toSeq
	"tester should filter things" in {
		iotesters.Driver.execute(Array("--backend-name", "firrtl", "--target-dir", "test_run_dir", "--fint-write-vcd"), () => new IRCSfilter(length = rcoeffs.length, inputwidth = 5, filterwidth = 5, scalingfactor = 8, import_rcoeffs = rcoeffs, import_icoeffs = icoeffs)) { c =>
			new IRCSfilterTester(c)
		} should be (true)
	}
}