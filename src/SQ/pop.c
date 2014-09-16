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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "int.h"
#include "basis.h"
#include "matrix.h"
#include "mol.h"
#include "option.h"
#include "util.h"
#include "uhf.h"
#include "pop.h"

// mulliken: compute Mulliken population analysis
//
// Aug 11, 2010 - Theerapon Khamla
//    Initial implementation
// 
// Oct 22, 2010 - Teepanis Chachiyo
//    Adjust for unrestricted case
//
void mulliken(
	int nBasis,               // number of basis function
	struct GTOBasis_t * gto,  // pointer to basis set information
	struct Molecule_t * mol,  // pointer to molecular structure
	int nEA,                  // number of spin up electron
	int nEB,                  // number of spin down electron
	double *CA,               // calculated molecular orbital coeff.
	double *CB,               // save but for beta spin
	double *eA,               // spin up eigen value
	double *eB,               // spin down eigen value
	struct option_t *option){ // print according to options

	double *PA;      // alpha electron density matrix
	double *PB;      // beta electron density matrix 
	double *PT;      // total density matrix
	double *PS;      // spin density matrix
	int *basis2Atom; // basis 2 atom mapping
	double *pop;     // mulliken population
	int i,j;         // generic loop indexes
	char str[256];   // atom name

	// memory allocation
	PA         = calloc(nBasis*nBasis, sizeof(double));
	PB         = calloc(nBasis*nBasis, sizeof(double));
	PT         = calloc(nBasis*nBasis, sizeof(double));
	PS         = calloc(nBasis*nBasis, sizeof(double));
	pop        = calloc(mol->nAtom,    sizeof(double));
	basis2Atom = calloc(nBasis,        sizeof(int));
	if(PA==NULL || PB==NULL || PT==NULL || PS==NULL || 
	   basis2Atom==NULL || pop==NULL){
		printf("mulliken: error - cannot allocate memory\n");
		exit(-1);
	}

	// determine atom center mapping
	for(i=0; i < nBasis; i++)
	for(j=0; j < mol->nAtom; j++)
		if(gto[i].x0==mol->x[j] && 
		   gto[i].y0==mol->y[j] && 
		   gto[i].z0==mol->z[j]){
			basis2Atom[i] = j;
			break;
		}

	// compute various types density matrix
	uhf_getDMatrix(nBasis, nEA, CA, PA);
	uhf_getDMatrix(nBasis, nEB, CB, PB);
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		PT[i*nBasis+j] = PA[i*nBasis+j] + PB[i*nBasis+j];
		PS[i*nBasis+j] = PA[i*nBasis+j] - PB[i*nBasis+j];
	}

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----               MULLIKEN POPULATION                 -----\n"
	"-------------------------------------------------------------\n"
	);

	// compute total mulliken charge due to electron contribution
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++)
		pop[basis2Atom[i]] -= PT[i*nBasis+j] * GTO_overlap(i, j, gto);

	// add nuclei contribution
	for(i=0; i < mol->nAtom; i++)
		pop[i] += mol->Z[i];

	// print value
#define CAPVALUE(a) (fabs(a)<1.0E-6)?0.0:a
	printf("       Atom       Mulliken Charge\n");
	for(i=0; i < mol->nAtom; i++){
		Z2SymShort(mol->Z[i],str);
		printf("%4d %5s    %15.6f\n",i+1,str,CAPVALUE(pop[i]));
	}

	// unrestricted case
	if(option->UHF){

		// reset zero
		for(i=0; i < mol->nAtom; i++) pop[i] = 0.0;

		// compute total mulliken spin population due to electron contribution
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++)
			pop[basis2Atom[i]] += PS[i*nBasis+j] * GTO_overlap(i, j, gto);

		// report value
		printf("       Atom       Mulliken Spin\n");
		for(i=0; i < mol->nAtom; i++){
			Z2SymShort(mol->Z[i],str);
			printf("%4d %5s    %15.6f\n",i+1,str,CAPVALUE(pop[i]));
		}

	}
#undef CAPVALUE

	// free memory
	free(PA);
	free(PB);
	free(PT);
	free(PS);
	free(basis2Atom);
	free(pop);
}


// electrostatic: print eletric potential and electric field at the nuclei
// positions in AU units.
//
// Nov 22, 2010 - Nanta Sophonrat
//         Initial Implementation     
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//         Added to Siam Quantum source tree
//
void electrostatic(
	int nBasis,               // number of basis function
	struct GTOBasis_t * gto,  // pointer to basis set information
	struct Molecule_t * mol,  // pointer to molecular structure
	int nEA,                  // number of spin up electron
	int nEB,                  // number of spin down electron
	double *CA,               // calculated molecular orbital coeff.
	double *CB,               // save but for beta spin
	double *eA,               // spin up eigen value
	double *eB,               // spin down eigen value
	struct option_t *option){ // print according to options

	int    i;             // loop index
	double *PA, *PB;      // density matrix for spin up and down
	double E[3];          // electric field
	double pot;           // potential

	// memory allocation
	PA = calloc(nBasis*nBasis, sizeof(double));
	PB = calloc(nBasis*nBasis, sizeof(double));
	if(PA==NULL || PB==NULL){
		printf("electrostatic : errror cannot allocate memroy\n");
		exit(-1);
	}

	// compute density matrix
	uhf_getDMatrix(nBasis, nEA, CA, PA);
	uhf_getDMatrix(nBasis, nEB, CB, PB);

	printf("\n\n"
       "-------------------------------------------------------------\n"
       "-----          ELECTROSTATIC PROPERTIES (AU)            -----\n"
       "-------------------------------------------------------------\n"
       "  Atom    Potential                 Electric Field           \n"
       "                             X             Y             Z   \n"
       "-------------------------------------------------------------\n");

	for(i=0; i<mol->nAtom; i++){
		// compute electric potential
		pot = uhf_potential(nBasis, gto, mol, PA, PB,
		                    mol->x[i], mol->y[i], mol->z[i]);

		// comput electric field
		uhf_efield(nBasis, gto, mol, PA, PB,
		           mol->x[i], mol->y[i], mol->z[i], E+0, E+1, E+2);

		printf("%4d %14.6f%14.6f%14.6f%14.6f\n", i+1, pot, E[0], E[1], E[2]);
	}
	printf(
       "-------------------------------------------------------------\n");
	free(PA);
	free(PB);
}

// electric_multipole: print eletric multipole moments of the molecule 
// The moments are calculated relative to the molecule coordinate. 
//
// As of Oct 13, 2012, there are some problems with quadrupole moments.
// If the molecular charge is zero, the Siam Quantum's results are identical
// to Gaussian's results. But if the net molecular charge is not zero, the
// results are different. I still need to sort this out.
//
// 2011 - Aniwat Kesorn
//      Initial implementation and testing
//
// Oct 13, 2012 - Teepanis Chachiyo
//      Added to Siam Quantum 
//
void electric_multipole(
	int nBasis,               // number of basis function
	struct GTOBasis_t * gto,  // pointer to basis set information
	struct Molecule_t * mol,  // pointer to molecular structure
	int nEA,                  // number of spin up electron
	int nEB,                  // number of spin down electron
	double *CA,               // calculated molecular orbital coeff.
	double *CB,               // save but for beta spin
	double *eA,               // spin up eigen value
	double *eB,               // spin down eigen value
	struct option_t *option){ // print according to options

	int    i,j;               // loop index
	double *PA, *PB;          // density matrix for spin up and down
	double px,py,pz;          // dipole moment

	//double qxx,qyy,qzz;       // diagonoal quadrupole moment
	//double qxy,qxz,qyz;       // off-diagonal quadrupole moment

	// memory allocation
	PA = calloc(nBasis*nBasis, sizeof(double));
	PB = calloc(nBasis*nBasis, sizeof(double));
	if(PA==NULL || PB==NULL){
		printf("electric_multipole : errror cannot allocate memory\n");
		exit(-1);
	}

	// compute density matrix
	uhf_getDMatrix(nBasis, nEA, CA, PA);
	uhf_getDMatrix(nBasis, nEB, CB, PB);

	printf("\n\n"
       "-------------------------------------------------------------\n"
       "-----           ELECTRIC MULTIPOLE MOMENT               -----\n"
       "-------------------------------------------------------------\n");

	//////////////////////////
	// compute dipole moment
	////////////////////////// 
	printf("\n"
       "                   Dipole Moment (Debye)                     \n"
       "                   ---------------------                     \n");

	px = 0.0; py = 0.0; pz = 0.0;

	// electronic contribution
	for(i=0; i<nBasis; i++)
	for(j=0; j<nBasis; j++){
		px -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,1,0,0,0.0,0.0,0.0);
		py -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,0,1,0,0.0,0.0,0.0);
		pz -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,0,0,1,0.0,0.0,0.0);
	}

	// nuclei contribution
	for(i=0; i<mol->nAtom; i++){
		px += mol->Z[i] * mol->x[i];
		py += mol->Z[i] * mol->y[i];
		pz += mol->Z[i] * mol->z[i];
	}

#define CAPVALUE(value,max) ((fabs(value)>max)?value:0.0)
	// convert to debye unit
	px *= AU2DEBYE; py *= AU2DEBYE; pz *= AU2DEBYE;
	printf("   X = %12.4f     Y = %12.4f     Z = %12.4f\n",
	       CAPVALUE(px,1E-4), CAPVALUE(py,1E-4), CAPVALUE(pz,1E-4));

	///////////////////////////////
	// compute quadrupole moment
	///////////////////////////////
	/*
	// EXPERIMENTAL
	printf("\n\n"
       "              Quadrupole Moment (Debye-Angstrom)             \n"
       "              ----------------------------------             \n");

	qxx = qyy = qzz = qxy = qxz = qyz = 0.0;

	// electronic contribution
	for(i=0; i<nBasis; i++)
	for(j=0; j<nBasis; j++){
		qxx -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,2,0,0,0.0,0.0,0.0);
		qyy -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,0,2,0,0.0,0.0,0.0);
		qzz -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,0,0,2,0.0,0.0,0.0);
		qxy -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,1,1,0,0.0,0.0,0.0);
		qxz -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,1,0,1,0.0,0.0,0.0);
		qyz -= (PA[i*nBasis+j]+PB[i*nBasis+j])*GTO_moment(i,j,gto,0,1,1,0.0,0.0,0.0);
	}

	// nuclei contribution
	for(i=0; i<mol->nAtom; i++){
		qxx += mol->Z[i] * mol->x[i] * mol->x[i];
		qyy += mol->Z[i] * mol->y[i] * mol->y[i];
		qzz += mol->Z[i] * mol->z[i] * mol->z[i];
		qxy += mol->Z[i] * mol->x[i] * mol->y[i];
		qxz += mol->Z[i] * mol->x[i] * mol->z[i];
		qyz += mol->Z[i] * mol->y[i] * mol->z[i];
	}

#define CAPVALUE(value,max) ((fabs(value)>max)?value:0.0)
	// convert to debye unit
	qxx *= AU2DEBYE*BOHR2ANGSTROM;
	qyy *= AU2DEBYE*BOHR2ANGSTROM;
	qzz *= AU2DEBYE*BOHR2ANGSTROM;
	qxy *= AU2DEBYE*BOHR2ANGSTROM;
	qxz *= AU2DEBYE*BOHR2ANGSTROM;
	qyz *= AU2DEBYE*BOHR2ANGSTROM;

	printf("  XX = %12.4f    YY = %12.4f    ZZ = %12.4f\n",
	       CAPVALUE(qxx,1E-4),
	       CAPVALUE(qyy,1E-4),
	       CAPVALUE(qzz,1E-4));
	printf("  XY = %12.4f    XZ = %12.4f    YZ = %12.4f\n",
	       CAPVALUE(qxy,1E-4),
	       CAPVALUE(qxz,1E-4),
	       CAPVALUE(qyz,1E-4));
	*/

	printf(
       "-------------------------------------------------------------\n");

	free(PA);
	free(PB);
}
