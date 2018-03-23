package fir

import chisel3._
import chisel3._
import chisel3.util._
import chisel3.iotesters.{ChiselFlatSpec, Driver, PeekPokeTester}


class FIRfilter(length: Int, inputwidth: Int, filterwidth: Int, import_coeffs: Seq[SInt]) extends Module{
	val io = IO(new Bundle {
		val in = Input(SInt(inputwidth.W))
		val out = Output(SInt((filterwidth*2+length-1).W))	// figure out bit growth later
	})

	// map coeffs to SInts
	// val coeffs = import_coeffs.map(_.S)
	val coeffs = import_coeffs

	// create an array holding the output of the delays
	val delays = Seq.fill(length)(Wire(SInt(inputwidth.W))).scan(io.in)( (prev: SInt, next: SInt) => {
		next := RegNext(prev)
		next
	})

	// multiply, storing results in 'mults'
	val mults = delays.zip(coeffs).map { case(delay: SInt, coeff: SInt) => delay * coeff }

	// add up multiplier outputs with bit growth
	val result = mults.reduce(_+&_)

	io.out := result

}