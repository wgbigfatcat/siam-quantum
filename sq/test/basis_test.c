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

// basis_test : test subroutines related to basis function and
// molecule managment. 
//
// Dec  12, 2008 - Teepanis Chachiyo
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

	// print basis set information
	print_BasisSet(dbSetItem, dbSet);

	// generate basis set function
	bas  = genBasis(mol, &nBasis, dbSetItem, dbSet);

	// print basis set function
	printGTOBasis(stdout, nBasis, bas);

	// clean memory
	cleanGTOBasis( bas, nBasis);
	cleanMolecule( mol );	

	// close file
	fclose(dbSetFd);
	fclose(xyzFd);
}
