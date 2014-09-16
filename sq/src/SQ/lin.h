/*
	Siam Quantum [SQ] is a program pacakage that performs quantum
	mechanical calculations on atomic systems.

	Copyright (C) 2008  Computational Physics @ KKU Group 

	This file is a part of SQ.
	                                                       
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2, or (at your option) 
	any later version.                                                  
*/

#ifndef LIN_H
#define LIN_H
int gen_sym_eigen(
	int nDim,                           // array dimension
	const double *A, const double *S,   // input A and S
	double *e, double *C);              // output value and vector

int linear_solver(int nRow,         // number of row
                  int nCol,         // number column
                  const double *A,  // input matrix A
                  double *B);       // input matrix B, output matrix X
#endif
