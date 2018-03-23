import chisel3._
import chisel3._
import chisel3.util._
import chisel3.iotesters.{ChiselFlatSpec, Driver, PeekPokeTester}


class FIRfilter(length: Int, bitwidth: Int, import_coeffs: Seq[Int]) extends Module{
	val io = IO(new Bundle {
		val in = Input(SInt(bitwidth.W))
		val out = Output(SInt(bitwidth*2+length-1).W))	// figure out bit growth later
	})


	val coeffs = import_coeffs.map(_.S)

	val delays = Seq.fill(length)(Wire(UInt(bitwidth.W))).scan(io.in)(prev: SInt, next: SInt) => {
		next := RegNext(prev)
		next
	}

	val mults = delays.zip(coeffs).map { case(delay: UInt, coeff: UInt) => delay * coeff }

	val result = mults.reduce(_+&_)

	io.out := result

}