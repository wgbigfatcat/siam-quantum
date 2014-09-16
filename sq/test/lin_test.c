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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "lin.h"

// print_matrix : print matrix of dimension nRow x nCol 
//
// Dec 13, 2008 - Teepanis Chachiyo
// 		Initial implementation
//
void print_matrix(double *A, int nRow, int nCol){
	int i,j;

	for(i=0; i < nRow; i++){
		for(j=0; j < nCol; j++){
			printf(" %18.5e ", A[i*nRow+j]);
		}
		printf("\n");
	}
}

// lin_test : testing subroutine lin.c which solves
// linear algebra system of equation by calling
// LAPACK library.
//
// Dec 13, 2008 - Teepanis Chachiyo
// 		Initial implementation
//
int main(int argc, char *argv[]){

	// define matrix
	double A[16] = {1.0, 2.0, 0.0, 0.0,
	                2.0, 1.0, 0.0, 0.0,
	                0.0, 0.0, 1.0, 0.0,
	                0.0, 0.0, 0.0, 1.0};
	double S[16] = {1.0, 0.0, 0.0, 0.0,
	                0.0, 1.0, 0.0, 0.0,
	                0.0, 0.0, 1.0, 0.0,
	                0.0, 0.0, 0.0, 2.0};
	double e[4], X[16];

	// define answer using MathCad14
	double MCad_e[4]  = {-1, 3, 1, 0.5};
	double MCad_X[16] = { 1.0, 1.0, 0.0, 0.0,
	                     -1.0, 1.0, 0.0, 0.0,
	                      0.0, 0.0, 1.0, 0.0,
	                      0.0, 0.0, 0.0, 1.0};

	printf("Solving generalized eigen value problem Ax=eSx\n");
	printf("Matrix A = \n");
	print_matrix(A, 4, 4);
	printf("Matrix S = \n");
	print_matrix(S, 4, 4);

	gen_sym_eigen(4, A, S, e, X);

	printf("Lapack Answer -- Eigen Row Vector -------\n");
	printf("Eigen Value = \n");
	print_matrix(e, 1, 4);
	printf("Matrix X = \n");
	print_matrix(X, 4, 4);
	
	printf("MathCAD14 Answer -- Eigen Column Vector ------\n");
	printf("Eigen Value = \n");
	print_matrix(MCad_e, 1, 4);
	printf("Matrix X = \n");
	print_matrix(MCad_X, 4, 4);

}
