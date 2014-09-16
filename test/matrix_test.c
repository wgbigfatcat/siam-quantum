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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basis.h"
#include "mol.h"
#include "int.h"
#include "matrix.h"

// matrix_test : test subroutines related to matrix.c
// which calculates kinetic and potential matrix element.
// Note. Electron-Electron matrix elements are tested
// elsewhere.
//
// The test is done by solving 1-electron hamiltonian
// system like H2+ , and He+. Then, compare the results
// with other Abinitio software package.
//
// Dec  13, 2008 - Teepanis Chachiyo
// 	Initial Implementation
//
int main(int argc, char *argv[]){

	struct Molecule_t *mol;      // pointer to molecule
	struct GTOBasis_t *bas;      // pointer to set of basis functions
	int    nBasis;               // number of basis function
	struct GTOBasisSet_t *dbSet; // pointer to basis-set database
	int    dbSetItem;            // number of item in basis-set database
	FILE   *dbSetFd;             // pointer to basis-set file
	FILE   *xyzFd;               // pointer to XYZ file

	double *H;                   // hamiltonian matrix
	double *S;                   // overlap matrix element
	double *E;                   // energy eigen values
	double *C;                   // eigen vectors

	int i,j;                     // indexes

	// check argument
	if(argc != 3){
		printf("usage: %s <Basis Set File>  <XYZ File>\n"
		       "<Basis Set File> - Basis set in GAMESS format\n"
		       "<XYZ File>       - Cartesian coordinate of atoms\n", argv[0]);
		exit(-1);
	}

	// open file
	dbSetFd = fopen(argv[1],"r");
	if(dbSetFd == NULL){
		printf("Cannot open file %s\n", argv[1]);
		exit(-1);
	}
	xyzFd = fopen(argv[2],"r");
	if(xyzFd == NULL){
		printf("Cannot open file %s\n", argv[2]);
		exit(-1);
	}

	// parse xyz file
	mol = readMolecule_XYZ(xyzFd);

	// parse basis set database
	dbSet = read_GAMESS_BasisSet(dbSetFd, argv[2], &dbSetItem);

	// generate basis set function
	bas  = genBasis(mol, &nBasis, dbSetItem, dbSet);

	// print basis set function
	printGTOBasis(stdout, nBasis, bas);

	// allocate matrix
	H = calloc(nBasis*nBasis, sizeof(double));
	S = calloc(nBasis*nBasis, sizeof(double));
	C = calloc(nBasis*nBasis, sizeof(double));
	E = calloc(nBasis,        sizeof(double));
	if(H==NULL || S==NULL || C==NULL || E==NULL){
		printf("Cannot allocate matrix\n");
		exit(-1);
	}

	// compute symmetric matrix element
	for(i=0; i < nBasis; i++)
	for(j=i; j < nBasis; j++){
		H[i*nBasis+j] = GTO_kinetic(i, j, bas) + GTO_nuclei(i, j, bas, mol);
		S[i*nBasis+j] = GTO_overlap(i, j, bas);

		H[j*nBasis+i] = H[i*nBasis+j];
		S[j*nBasis+i] = S[i*nBasis+j];
	}

	// solve generalize eigen value problem
	gen_sym_eigen(nBasis, H, S, E, C);

	// print out result
	printf("Ground state energy = %15.8E [Hartrees]\n", 
	                                            E[0] + nuclei_coulomb(mol));	

	// clean memory
	cleanGTOBasis( bas, nBasis);
	cleanMolecule( mol );	

	free(H);
	free(S);
	free(C);
	free(E);

	// close file
	fclose(dbSetFd);
	fclose(xyzFd);
}
