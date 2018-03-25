package fir

import chisel3._
import chisel3.util._
import chisel3.iotesters.{ChiselFlatSpec, Driver, PeekPokeTester}

class IRCSfilter(length: Int, inputwidth: Int, filterwidth: Int, import_rcoeffs: Seq[SInt], import_icoeffs: Seq[SInt]) extends Module {

	val io = IO(new Bundle {
		val isig = Input(SInt(inputwidth.W))
		val qsig = Input(SInt(inputwidth.W))
		val isig_out = Output(SInt(inputwidth.W))
	})

	val isig_filter = Module(new FIRfilter(length = length, inputwidth = inputwidth, filterwidth = filterwidth, import_coeffs = import_rcoeffs))
	val x1 = Wire(SInt((filterwidth*2+length-1).W))
	isig_filter.io.in := io.isig
	x1 := isig_filter.io.out


	val qsig_filter = Module(new FIRfilter(length = length, inputwidth = inputwidth, filterwidth = filterwidth, import_coeffs = import_icoeffs))
	val x2 = Wire(SInt((filterwidth*2+length-1).W))
	qsig_filter.io.in := io.qsig
	x2 := qsig_filter.io.out

	io.isig_out := (x1 +& x2) >> 11.U

}