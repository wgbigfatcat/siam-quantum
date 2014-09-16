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
#ifndef MECP_H
#define MECP_H

#include "basis.h"
#include "mol.h"
#include "option.h"

// convergence criterion
#define CONV_FORCEMAX   0.000450
#define CONV_FORCERMS   0.000300
#define CONV_DISPMAX    0.001800
#define CONV_DISPRMS    0.001200
#define MAXSTEPSIZE     0.5
#define MAXOPTSTEP      30

void Hessian_BFGS(
	int nDim,            // matrix dimension
	const double *dR,    // displacement vector
	const double *dGrad, // change of gradient vector,
	double *H);          // original Hessian and the returned update values

void mecp(
	int dbSize,                         // number of record in basis set db
	const struct GTOBasisSet_t *basisDB,// pointer to basis set database
	struct Molecule_t *mol,             // returned molecular coordinate
	struct option_t *opt);              // options
#endif
