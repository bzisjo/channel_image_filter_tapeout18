package fir

import chisel3._
import chisel3.util._
import chisel3.iotesters._
import chisel3.iotesters.{ChiselFlatSpec, PeekPokeTester, SteppedHWIOTester}
import org.scalatest.{FreeSpec, Matchers}


class FIRfilterTesterSpec extends FreeSpec with Matchers {
	val rcoeffs_file = "./src/test/scala/fir/rcoeffs.csv"
	val icoeffs_file = "./src/test/scala/fir/icoeffs.csv"
	val rcoeffs:Seq[SInt] = Source.fromFile(filename).getLines().next.split(",").map(_.trim).map(_.toInt).map(_.asSint).toSeq
	"tester should filter things" in {
		iotesters.Driver.execute(Array("--backend-name", "firrtl", "--target-dir", "test_run_dir", "--fint-write-vcd"), () => new FIRfilter()) { c =>
			new FIRfilterTester(c)
		} should be (true)
	}
}