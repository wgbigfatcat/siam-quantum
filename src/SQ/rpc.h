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
#ifndef RPC_H
#define RPC_H

#include "basis.h"
#include "option.h"

int rpc_GTO_JK_Matrix_Quartet_Parallel(
	int status,                    // current rpc status
	int childID,                   // childID to call
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *schwarz_basis,   // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB,                    // return G for spin down
	struct option_t *opt);         // global option

int rpc_mp2_rhf_aqij_Parallel_Contribute(
	int status,                    // current rpc status
	int childID,                   // childID to call
	double *mp2,                   // returned contribution
	int a,                         // starting orbital index
	int maxCorr,                   // number of correlated orbitals
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const double *Schwarz,         // schwarz inequality
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

int rpc_GradEri_ShellSet_Parallel(
	int status,                     // current rpc status
	int nAtom,                      // number of atoms
	int childID,                    // this child id
	int nCPU,                       // number of CPUs
	int nBasis,                     // number of basis function
	const struct GTOBasis_t *gto,   // basis function database
	const int *basis2Atom,          // basis to atom mapping
	const double *PT,               // total density matrix
	const double *PS,               // spin density matrix
	double *Gx,                     // returned gradient in x direction
	double *Gy,                     // returned gradient in y direction
	double *Gz,                     // returned gradient in z direction
	const struct option_t *opt);    // global options

void trapChild(int argc, char *argv[]);
void rpcSpawnChildren(int argc, char *argv[], struct option_t *opt);
void rpcExitChildren(struct option_t *opt);

// rpc message header
struct rpcHdr_t{
	int id;     // message id
	int type;   // type of this message
	int len;    // length in bytes of the following message
};

// rpc message types
#define RPC_CHILDEXIT            0
#define RPC_GTOJKMATRIXQUARTET   1
#define RPC_MP2RHFAQIJCONTRIBUTE 2
#define RPC_GRADERISHELLSET      3

// rpc status
#define RPC_DONE -1
#define RPC_IDLE  0

#define RPC_DELAY 200 // delay time in milli-seconds

#endif // RPC_H
