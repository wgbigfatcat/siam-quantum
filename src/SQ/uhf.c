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
#include "basis.h"
#include "matrix.h"
#include "lin.h"
#include "int.h"
#include "fgamma.h"
#include "util.h"
#include "option.h"
#include "rhf.h"
#include "uhf.h"
#include "conv.h"
#include "check.h"
#include "rpc.h"

// normalizeC : normalized eigen vector 
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static void normalizeC(int nBasis, double *S, double *C){
	double sum;
	int i,j,k;

	for(k=0; k < nBasis; k++){
		sum = 0.0;
		for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++)
				sum += C[k*nBasis+i]
				      *C[k*nBasis+j]
				      *S[i*nBasis+j];
		sum = 1.0/sqrt(sum);
		for(i=0; i < nBasis; i++)
			C[k*nBasis+i] = sum * C[k*nBasis+i];
	}
	return;
}

// uhf_getDMatrix : compute density matrix
//
// July 10, 2010 - Teepanis Chachiyo
//    Migrate to uhf scheme
//
// 2008 - Teepanis Chachiyo
//    Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//    Not using static anymore
//
void uhf_getDMatrix(int nBasis, int nOcc, double *C, double *P){
	int i,j,k;
	double sum;

	for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			sum = 0.0;
			for(k=0; k < nOcc; k++){
				sum += C[k*nBasis+i] * C[k*nBasis+j];
			}
			P[i*nBasis+j] = sum;
		}
	return;
}

// uhf_getGMatrix : compute G matrix element for both alpha and beta spin
//
// Mar 6, 2013 - Teepanis Chachiyo
//  Use parallel version
//
// Nov 19, 2012 - Teepanis Chachiyo
//  Passing cutoff value as an argument
//
// July 10, 2010 - Teepanis Chachiyo
//  Migrate to unrestricted calculations
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static void uhf_getGMatrix(
	int nBasis,
	struct GTOBasis_t *gto,
	double *Schwarz,
	double cutoff,
	double *PA, double *PB,
	double *GA, double *GB,
	struct option_t *opt){

	int i,j;
	
	// reset to zero
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		GA[i*nBasis+j] = 0.0;
		GB[i*nBasis+j] = 0.0;
	}

	/////////////////////
	// parallel version
	/////////////////////
	double *GAset, *GBset;   // set of GA,GB for each cpu
	int *status;             // status for each cpu
	int alldone;             // all idle flag
	int n;                   // cpu counter

	// allocate memory
	status=calloc(opt->nCPU,sizeof(int));
	GAset =calloc(nBasis*nBasis*opt->nCPU,sizeof(double));
	GBset =calloc(nBasis*nBasis*opt->nCPU,sizeof(double));
	if(status==NULL || GAset==NULL || GBset==NULL){
		printf("uhf_getGMatrix - error cannot allocate memory\n");
		exit(-1);
	}

	// reset status to idle
	for(n=(opt->nCPU-1);n>=0;n--) status[n] = RPC_IDLE;

	// loop thru all cpu and compute GA and GB
	do{
		// remote call to all cpu except childID=0
		for(n=(opt->nCPU-1);n>0;n--)
			if(status[n] != RPC_DONE)
			status[n] = rpc_GTO_JK_Matrix_Quartet_Parallel(status[n],n,nBasis,
			                                               PA, PB, gto, Schwarz,
			                                               cutoff,
			                                               GAset+nBasis*nBasis*n,
			                                               GBset+nBasis*nBasis*n,
			                                               opt);

		// local call for childID=0
		if(status[0]==RPC_IDLE){
			GTO_JK_Matrix_Quartet_Parallel(0, nBasis, PA, PB, gto, Schwarz, cutoff,
			                               GAset, GBset, opt);
			status[0]=RPC_DONE;
		}

		// check if all done
		alldone=1;
		for(n=(opt->nCPU-1);n>=0;n--) if(status[n] != RPC_DONE) alldone=0;

	}while(!alldone);

	// accumulate GA and GB
	for(n=(opt->nCPU-1);n>=0;n--){
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			GA[i*nBasis+j] += GAset[nBasis*nBasis*n +i*nBasis+j];
			GB[i*nBasis+j] += GBset[nBasis*nBasis*n +i*nBasis+j];
		}
	}

	// clean memory
	free(status);
	free(GAset);
	free(GBset);

	///////////////////////////////////
	// perform G matrix computation  //
	///////////////////////////////////
	//GTO_JK_Matrix_Quartet(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB, opt);
	//GTO_JK_Matrix_ShellSet(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB, opt);
	//GTO_JK_Matrix_NoSymm(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB);
	//GTO_JK_Matrix_PrimeSpace(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB);
	//GTO_JK_Matrix_PrimeSpaceShellSet(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB);
}

// getEtotal : compute total energy this is equal to
// electronic part + nuclei part
//
// July 12, 2010 - Teepanis Chachiyo
//  Migrate to unrestricted calculations
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static double uhf_getEtotal(
	int nBasis, struct Molecule_t *mol,
	double *PA, double *FA,
	double *PB, double *FB, 
	double *H){

	double E=0.0;
	int i,j;

	// compute nuclei repulsion energy
	E += nuclei_coulomb(mol);

	// include electron energy given in (Szabo and Ostlund)
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		E += 0.5*((H[i*nBasis+j]+FA[i*nBasis+j])*PA[i*nBasis+j] +
			      (H[i*nBasis+j]+FB[i*nBasis+j])*PB[i*nBasis+j]);
	}

	return E;
}

// uhf_rho : computes and return electron density at
// specified Cartesian point (x,y,z)
//
// August 15, 2010 Nanta Sophonrat
//  Fix bug when computing electron density. Forgot to multiply by 2
//
// July 12, 2010 - Teepanis Chachiyo
//  Migrate to unrestricted calculation
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//  The subroutine now takes density matrix as an argument
//  as supposed to molecular orbital coefficient. This helps
//  increase speed.
//
// March 18, 2010 - Teepanis Chachiyo
//  Store eval_chi() in memory to reduce the number of calculations
//  from nBasis*nBasis to nBasis only. Also exploit P[] symmetrical
//  nature to reduce calculation time by a factor of two.
//
double uhf_rho(int nBasis,              // number of basis function
               struct GTOBasis_t *gto,  // function structure
               double *PA, double *PB,  // density matrix
               double x, double y, double z){

	int i,j;
	double sum=0.0;
	double *chi;

/* //////// very inefficient //////////

	// evaluate density
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		sum += P[i*nBasis+j]*eval_chi(i,gto,x,y,z)
		                    *eval_chi(j,gto,x,y,z);
	}
///////////////////////////////////// */

	// allocate memory
	chi=calloc(nBasis, sizeof(double));
	if(chi==NULL){
		printf("uhf_rho - error cannot allocate memory\n");
		exit(-1);
	}

	// evaluate chi at a point
	for(i=0; i < nBasis; i++) chi[i] = eval_chi(i,gto,x,y,z);

	//
	// evaluate density
	//

	// off diagonal
	for(i=0; i < nBasis; i++)
	for(j=0; j < i; j++)
		sum += (PA[i*nBasis+j]+PB[i*nBasis+j])*chi[i]*chi[j];
	sum = 2.0*sum;

	// diagonal
	for(i=0; i < nBasis; i++)
		sum += (PA[i*nBasis+i]+PB[i*nBasis+i])*chi[i]*chi[i];

	// clean up memory
	free(chi);

	return sum;
}


//
// uhf_rho_XYPlane : compute electron density for all points on a grid
// on xy-plane. This should be faster than compute each point one-by-one.
//
// Oct 8, 2012 - Teepanis Chachiyo
//     Initial implementationa from uhf_rho(...) and testing
//
void uhf_rho_XYPlane(
	int nBasis,                         // number of basis function
    const struct GTOBasis_t *gto,       // basis function structure
    const double *PA, const double *PB, // density matrix
	double cutoff,                      // ignore prefactor
	double x0, double y0, double z0,    // origin of the grid
	double dx, double dy,               // step size
	int nx, int ny,                     // number of points
	double *rhoXY){                     // returned matrix of size ny*nx

	int p,q;                  // basis function index
	int cP, cQ;               // contracted function index
	int i,j;                  // grid point index

	double *X, *Y, Z;         // buffer in x, y, and z directions
	double x,y;               // current coordinates
	double prefactor;         // prefactor
	double rab2;              // distance between two basis function
	double x1,x2,y1,y2,z1,z2; // center of each basis functions

	// allocation
	if((X=calloc(nx,sizeof(double)))==NULL ||
	   (Y=calloc(ny,sizeof(double)))==NULL){
		printf("uhf_rho_XYPlane : error cannot allocate memory\n");
		exit(-1);
	}

	// reset values
	for(j=0; j < ny; j++)
	for(i=0; i < nx; i++)
		rhoXY[j*nx+i] = 0.0;

	// loop through all basis functions
	for(p=0; p < nBasis; p++)
	for(q=0; q <= p; q++){

		// set center variables
		x1 = gto[p].x0; y1 = gto[p].y0; z1 = gto[p].z0;
		x2 = gto[q].x0; y2 = gto[q].y0; z2 = gto[q].z0;
#define DIST2(x1,y1,z1,x2,y2,z2) ((x1-x2)*(x1-x2)+\
                                  (y1-y2)*(y1-y2)+\
                                  (z1-z2)*(z1-z2))
		rab2 = DIST2(x1,y1,z1,x2,y2,z2);

		// loop through contractions
		for(cP=0; cP < gto[p].nContract; cP++)
		for(cQ=0; cQ < gto[q].nContract; cQ++){

			// screening
			prefactor = exp(-gto[p].exp[cP]*gto[q].exp[cQ]*rab2/
			                (gto[p].exp[cP]+gto[q].exp[cQ]));
			if(fabs(prefactor) < cutoff) continue;

			// compute prefactor which is indenpendent of points
			prefactor = (PA[p*nBasis+q] + PB[p*nBasis+q]) *
			            gto[p].norm[cP] * gto[q].norm[cQ] *
			            gto[p].coef[cP] * gto[q].coef[cQ];

			// apply symmetry
			if(p!=q) prefactor = prefactor+prefactor;

			// compute Z specific value
			Z = pow_int(z0-gto[p].z0,gto[p].n)*pow_int(z0-gto[q].z0,gto[q].n)*
			    exp( - gto[p].exp[cP] * (z0-gto[p].z0) * (z0-gto[p].z0)
			         - gto[q].exp[cQ] * (z0-gto[q].z0) * (z0-gto[q].z0));

			// absorb Z into prefactor
			prefactor = prefactor * Z;

			// compute Y specfic values
			for(j=0; j < ny; j++){
				y    = y0+j*dy;
				Y[j] = pow_int(y-gto[p].y0,gto[p].m)*pow_int(y-gto[q].y0,gto[q].m)*
				       exp( - gto[p].exp[cP] * (y-gto[p].y0) * (y-gto[p].y0)
				            - gto[q].exp[cQ] * (y-gto[q].y0) * (y-gto[q].y0));
			}

			// compute X specific values
			for(i=0; i < nx; i++){
				x    = x0+i*dx;
				X[i] = pow_int(x-gto[p].x0,gto[p].l)*pow_int(x-gto[q].x0,gto[q].l)*
				       exp( - gto[p].exp[cP] * (x-gto[p].x0) * (x-gto[p].x0)
				            - gto[q].exp[cQ] * (x-gto[q].x0) * (x-gto[q].x0));
			}

			// loop through all points on xy-plane
			for(j=0; j < ny; j++)
			for(i=0; i < nx; i++)
				rhoXY[j*nx+i] += prefactor * X[i] * Y[j];
		}
	}

	// free memory
	free(X);
	free(Y);
}


//
// uhf_potential_XYPlane : compute electric potential for all points on a grid
// on xy-plane. This should be faster than compute each point one-by-one.
//
// Oct 9, 2012 - Teepanis Chachiyo
//     Initial implementationa from uhf_rho_XYPlane(...) and testing
//
void uhf_potential_XYPlane(
	int nBasis,                         // number of basis function
    const struct GTOBasis_t *gto,       // basis function structure
	const struct Molecule_t *mol,       // molecule structure 
    const double *PA, const double *PB, // density matrix
	double cutoff,                      // ignore prefactor
	double x0, double y0, double z0,    // origin of the grid
	double dx, double dy,               // step size
	int nx, int ny,                     // number of points
	double *phiXY){                     // returned matrix of size ny*nx

	int p,q;             // basis function index
	int cP, cQ;          // contracted function index
	int i,j;             // grid point index
	int A;               // atom index
	double rr;           // distance from point to nuclei
	double r[3];         // position of nucleus

	double prefactor;    // prefactor

#define MAXL 8
	// nai integral sections
	double eta,xp,yp,zp,sum,rab2,rcp2;
	double *Ax,*Ay,Az[2*MAXL],*tAx,*tAy;
	double F[2*MAXL];
	int kX, kY, kZ;
	double x1,y1,z1,x2,y2,z2;
	int l1,m1,n1,l2,m2,n2;

	// allocation
	if((Ax=calloc(nx*2*MAXL,sizeof(double)))==NULL ||
	   (Ay=calloc(ny*2*MAXL,sizeof(double)))==NULL){
		printf("uhf_potential_XYPlane : error cannot allocate memory\n");
		exit(-1);
	}

	// reset values
	for(j=0; j < ny; j++)
	for(i=0; i < nx; i++)
		phiXY[j*nx+i] = 0.0;

	// loop through all basis functions
	for(p=0; p < nBasis; p++)
	for(q=0; q <= p; q++){

		// set angular index
		l1 = gto[p].l; m1 = gto[p].m; n1 = gto[p].n;
		l2 = gto[q].l; m2 = gto[q].m; n2 = gto[q].n;

		// set center variables
		x1 = gto[p].x0; y1 = gto[p].y0; z1 = gto[p].z0;
		x2 = gto[q].x0; y2 = gto[q].y0; z2 = gto[q].z0;

		rab2 = DIST2(x1,y1,z1,x2,y2,z2);

		// loop through contractions
		for(cP=0; cP < gto[p].nContract; cP++)
		for(cQ=0; cQ < gto[q].nContract; cQ++){

			// prepare for nai1d
			eta = 1.0/(gto[p].exp[cP]+gto[q].exp[cQ]);

			// prefactor and screening
			prefactor = exp(-gto[p].exp[cP]*gto[q].exp[cQ]*rab2*eta);

			// screening
			if(fabs(prefactor) < cutoff) continue;

			// compute prefactor which is indenpendent of points
			prefactor = prefactor                         * 
			            (PA[p*nBasis+q] + PB[p*nBasis+q]) *
			            gto[p].norm[cP] * gto[q].norm[cQ] *
			            gto[p].coef[cP] * gto[q].coef[cQ] *
	                    2.0 * M_PI * eta;

			// apply symmetry
			if(p!=q) prefactor = prefactor+prefactor;

			// compute center of two-gaussian function
			xp   = (gto[p].exp[cP]*x1+gto[q].exp[cQ]*x2)*eta;
			yp   = (gto[p].exp[cP]*y1+gto[q].exp[cQ]*y2)*eta;
			zp   = (gto[p].exp[cP]*z1+gto[q].exp[cQ]*z2)*eta;

			// compute Z specific value
			naiA_1d(Az,n1,n2,zp-z1,zp-z2,zp-z0,eta);

			// compute Y specfic values
			for(j=0; j < ny; j++)
				naiA_1d(Ay+j*2*MAXL,m1,m2,yp-y1,yp-y2,yp-y0-j*dy,eta);

			// compute X specific values
			for(i=0; i < nx; i++)
				naiA_1d(Ax+i*2*MAXL,l1,l2,xp-x1,xp-x2,xp-x0-i*dx,eta);

			// loop through all points on xy-plane
			for(j=0; j < ny; j++)
			for(i=0; i < nx; i++){

				rcp2 = DIST2(x0+i*dx,y0+j*dy,z0,xp,yp,zp);
				fgamma_set(l1+l2+m1+m2+n1+n2,rcp2/eta,F);

				sum = 0.0;
				tAx = Ax+i*2*MAXL;
				tAy = Ay+j*2*MAXL;
				for(kX=0; kX < l1+l2+1; kX++)
				for(kY=0; kY < m1+m2+1; kY++)
				for(kZ=0; kZ < n1+n2+1; kZ++)
					sum += tAx[kX]*tAy[kY]*Az[kZ]*F[kX+kY+kZ];

				phiXY[j*nx+i] -= prefactor * sum;
			}
		}
	}

	// include potential from nuclei
	for(j=0; j < ny; j++)
	for(i=0; i < nx; i++)
	for(A=0; A<mol->nAtom; A++){
		r[0]  = x0+i*dx - mol->x[A];
		r[1]  = y0+j*dy - mol->y[A];
		r[2]  = z0      - mol->z[A];

		rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

		if(rr == 0.0) continue;
		else phiXY[j*nx+i] += mol->Z[A]/rr;
	}

	// free memory
	free(Ax);
	free(Ay);
}


// uhf_mo : computes and return molecular orbital at
// specified Cartesian point (x,y,z)
//
// July 12, 2010 - Teepanis Chachiyo
//    Migrate to unrestricted calculation
//
// Dec 2009 - Theerapon Khamla
//    Initial implementation
//
// Mar 05, 2010 - Teepanis Chachiyo
//    Bug fix, invalid range of molecular orbital should be
//    "n >= nBasis" not "n > nBasis" as before.
//
double uhf_mo(int nBasis,              // number of basis function
              struct GTOBasis_t *gto,  // function structure
              double *C,               // molecular orbital
              int n,                   // orbital index (1st orbital is zero)
              double x, double y, double z){
	int i;
	double sum;

	// sanitize
	if(n >= nBasis || n < 0){
		printf("uhf_mo: error invalid molecular orbital index\n");
		exit(EXIT_FAILURE);
	}

	// evaluate molecular orbital
	sum = 0.0;
	for(i=0; i < nBasis; i++){
		sum += C[n*nBasis+i]*eval_chi(i,gto,x,y,z);
	}

	return sum;

}

// uhf_potential : computes and return electric potential at the specified 
// Cartesian point (x,y,z). It also includes the contribution due to the
// nuclei in the system, except when the point (x,y,z) is right on top of
// a nucleus, the contribution from that particular nucleus is ignored.
//
// Oct 26, 2010 - Nanta Sophonrat
//         Initial Implementation     
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//         Added to Siam Quantum source tree
//
double uhf_potential(int nBasis,                    // number of basis function
                     const struct GTOBasis_t *gto,  // function structure
                     const struct Molecule_t *mol,  // molecular structure info
                     const double *PA,              // alpha density matrix
                     const double *PB,              // beta density matrix
                     double x, double y, double z){

	int    i=0,j=0;     // loop index
	double rr;          // distance between nucleus and point(x1, x2, x3)
	double pot = 0.0;   // electric potential
	double r[3];        // distances between atoms

	// multiply by 2 due to the symmetry
	// off-diagonal elements
	for(i=1; i<nBasis; i++)
	for(j=0; j<i; j++){
		pot += 2.0*(PA[i*nBasis+j]+PB[i*nBasis+j])
		          * GTO_nai(i, j, gto, x, y, z);
	}

	// diagonal elements
	for(i=0; i<nBasis; i++)
		pot += (PA[i*nBasis+i]+PB[i*nBasis+i])*GTO_nai(i, i, gto, x, y, z);

	// include potential from nuclei
	for(i=0; i<mol->nAtom; i++){
		r[0]  = x - mol->x[i];
		r[1]  = y - mol->y[i];
		r[2]  = z - mol->z[i];

		rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

		if(rr == 0.0) continue;
		else pot += mol->Z[i]/rr;
	}
	return pot;
}

// uhf_efield: compute electric field at a specified point (x,y,z). It also 
// includes the nuclei contribution, except when the point (x,y,z) is right
// on top of the nucles, the contribution from that particular nucleus is
// ignored.
//
// Nov 22, 2010 - Nanta Sophonrat
//         Initial Implementation     
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//         Added to Siam Quantum source tree
//
void uhf_efield(int nBasis,                          // number of basis function
                const struct GTOBasis_t *gto,        // function structure
                const struct Molecule_t *mol,        // molecular structure info
                const double *PA,                    // alpha density matrix
                const double *PB,                    // beta density matrix
                double x, double y, double z,        // specified point
                double *ex, double *ey, double *ez){ // return field

	int    i=0,j=0;       // basis function index
	double rr;            // distance between the nucleus and (x,y,z)
	double r[3];          // displacement between the nucleus and (x,y,z)
	double tex, tey, tez; // temporary field

	// initial value
	*ex = 0.0; *ey = 0.0; *ez = 0.0;

	// multiply by 2 due to the symmetry
	// off-diagonal elements
	for(i=1; i<nBasis; i++)
	for(j=0; j<i; j++){
			GTO_efi(i,j,gto,x,y,z,&tex,&tey,&tez);
			*ex -= 2.0*(PA[i*nBasis+j]+PB[i*nBasis+j])*tex;
			*ey -= 2.0*(PA[i*nBasis+j]+PB[i*nBasis+j])*tey;			
			*ez -= 2.0*(PA[i*nBasis+j]+PB[i*nBasis+j])*tez;
	}

	// diagonal elements
	for(i=0; i<nBasis; i++){
		GTO_efi(i,i,gto,x,y,z,&tex,&tey,&tez);
		*ex -= (PA[i*nBasis+i]+PB[i*nBasis+i])*tex;
		*ey -= (PA[i*nBasis+i]+PB[i*nBasis+i])*tey;
		*ez -= (PA[i*nBasis+i]+PB[i*nBasis+i])*tez;
	}

	// include nuclei contribution
	for(i=0; i<mol->nAtom; i++){
		r[0]  = x - mol->x[i];
		r[1]  = y - mol->y[i];
		r[2]  = z - mol->z[i];

		rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
		rr = rr*rr*rr;

		if(rr == 0.0) continue;
		else{
			*ex += r[0]*mol->Z[i]/rr;
			*ey += r[1]*mol->Z[i]/rr;
			*ez += r[2]*mol->Z[i]/rr;
		}
	}
}

// uhf : carries out Unrestricted Hartree-Fock calculation until
// convergence is reached. It returns total electronic energy of the
// systems.
//
// If the scf convergence is not reached, it returns zero instead.
//
// In this version, the e-e integral is calculated using direct
// integration every time it is needed. They are not stored in
// the memory, which is not a good idea!
//
// Feb 2013 - Teepanis Chachiyo
//    - No longer need to rebuild fock matrix when switching accuracy
//      This is only for the case when the initial accuracy is <= 10^-3
//
// Jan 26, 2013 - Teepanis Chachiyo
//    - Rewrite the scf convergence structure
//
// Jan 9, 2013 - Teepanis Chachiyo
//    - Forgot to check avgdP criteria when switching accuracy
//
// Nov 26, 2012 - Teepanis Chachiyo
//    - Using 3 steps integral accuracy
//    - Check point is saved even if the run is not converged
//    - Add option to save checkpoint every scf cycle
//
// Nov 19, 2012 - Teepanis Chachiyo
//    - Schwarz cutoff is calculated according to convergence requirement
//      and the density matrix
//    - do not report initial non-interacting electronic energy
// 
// Oct 21, 2010 - Teepanis Chachiyo
//    Use option SCFConv to set convergence threshold
//    The function now returns zero if the convergence failed.
//
// Sep 6, 2010 - Teepanis Chachiyo
//    Use density matrix difference (insted of absolute density matrix value)
//    to compute the change in G matrix
//    Use idensity matrix as initial density matrix
//
// Sep 6, 2010 - Teepanis Chachiyo
//    Take the convergence part outside "uhf" function
//
// July 12, 2010 - Teepanis Chachiyo
//    Migrate from RHF to UHF
//
// 2008 - Teepanis Chachiyo
// 	  Initial implementation
//
// Dec 14, 2008 - Teepanis Chachiyo
// 	  Simplify and use Direct integral for the first release
//	  version.
//
// March 11, 2010 - Teepanis Chachiyo
//    Add reading and saving density matrix
//
// May 21, 2011 - Teepanis Chachiyo
//    Add diis option
//
// Mar 6, 2013 - Teepanis Chachiyo
//     Use parallel version
//
// Mar 19, 2013 - Teepanis Chachiyo
//     Revamp convergence checking
//
double uhf(
	int nBasis,              // number of basis functions
	struct GTOBasis_t * gto, // pointer to function structure
	struct Molecule_t * mol, // pointer to molecule structure
	int nEA,                 // total number of spin up electrons
	int nEB,                 // total number of spin down electrons
	double *CA,              // returned molecular alpha spin orbital
	double *CB,              // returned molecular beta spin orbital 
	double *eA,              // returned eigen values
	double *eB,              // returned eigen values
	struct option_t *opt){   // global option

	double Etot=0.0;      // total electronic energy
	double dE=0.0;        // energy change
	double avgdP=0.0;     // average change in density matrix
	double gamma=1.0;     // update coefficient
	double sum=0.0;       // generic summation variable
	double realCutoff;    // schwarz cutoff
	double thisCutoff;    // current cutoff
	double cutoffA=0.0;   // first stage cutoff 
	double cutoffB=0.0;   // intermediate cutoff
	double cutoffC=0.0;   // final stage cutoff
	double sumA=0.0;      // sum of alpha density matrix
	double sumB=0.0;      // sum of beta density matrix

	// matrix elements
	double *PA=NULL;     // density matrix
	double *GA=NULL;     // EE matrix
	double *dPA=NULL;    // change in density matrix
	double *dGA=NULL;    // change in G matrix
	double *FA=NULL;     // fock matrix
	double *PB=NULL;     // density matrix
	double *GB=NULL;     // EE matrix
	double *dPB=NULL;    // change in density matrix
	double *dGB=NULL;    // change in G matrix
	double *FB=NULL;     // fock matrix
	double *H=NULL;      // h core matrix
	double *T=NULL;      // kinetic matrix
	double *V=NULL;      // nuclei potential matrix
	double *S=NULL;      // overlap matrix
	double *Schwarz=NULL;// Schwarz upper bound matrix

	int i,j,iter=0;
	int notConverged=1;

	FILE *fd;            // file descriptor for density matrix

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----            SOLVING HARTREE-FOCK EQUATION          -----\n"
	"-------------------------------------------------------------\n"
	);
	fflush(stdout);

	// report
	if(opt->RHF && opt->UHF){
		printf("uhf - error detect both RHF and UHF activated\n");
		exit(-1);
	}

	if(opt->UHF){
		printf("Requested Unrestricted Hartree-Fock calculations\n");
		printf("There are %d alpha spin and %d beta spin electrons\n", nEA, nEB);
	}

	if(opt->RHF){
		printf("Requested Restricted Hartree-Fock calculations\n");
		printf("There are %d electrons in the density matrix\n", nEA+nEB);
		if(nEA!=nEB){
			printf("uhf - error number of electron in each spin is not the same\n");
			exit(-1);
		}
	}
	printMolecule_XYZ(mol, stdout);

#define ALLOCATE(P)                                      \
P = calloc(nBasis*nBasis, sizeof(double));               \
if(P==NULL){                                             \
	printf("rhf: Error - Cannot allocate memory\n"); \
	exit(EXIT_FAILURE);                              \
}

	// memory allocation
	ALLOCATE(PA);
	ALLOCATE(GA);
	ALLOCATE(dPA);
	ALLOCATE(dGA);
	ALLOCATE(FA);
	ALLOCATE(PB);
	ALLOCATE(GB);
	ALLOCATE(dPB);
	ALLOCATE(dGB);
	ALLOCATE(FB);
	ALLOCATE(H);
	ALLOCATE(T);
	ALLOCATE(V);
	ALLOCATE(S);

#undef  ALLOCATE

	/////////////////////////////////////////////
	// Building necessary matrix elements     ///
	/////////////////////////////////////////////

	// report
	printf(
	"Computing 1-electron matrix elements ... \n"
	"    H - Core Hamiltonian                 \n"
	"    S - Overlap Matrix                   \n");

	// get kinetic matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		T[i*nBasis+j] = GTO_kinetic(i,j,gto);
		// symmetrize matrix
		T[j*nBasis+i] = T[i*nBasis+j];
	}

	// get nuclei matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		V[i*nBasis+j] = GTO_nuclei(i,j,gto,mol);
		// symmetrize matrix
		V[j*nBasis+i] = V[i*nBasis+j];
	}

	// get overlap matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		S[i*nBasis+j] = GTO_overlap(i,j,gto);
		// symmetrize matrix
		S[j*nBasis+i] = S[i*nBasis+j];
	}

	// build Hcore 
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		H[i*nBasis+j] = T[i*nBasis+j] + V[i*nBasis+j];
	}

	//
	// Diagonalize core hamiltonian to guess density matrix
	//
	if(opt->SCFGuess == SCFGUESS_CORE){
		fflush(stdout);
		printf("Diagonalizing H for initial density matrix ...\n");
		gen_sym_eigen(nBasis, H, S, eA, CA);
		gen_sym_eigen(nBasis, H, S, eB, CB);
		normalizeC(nBasis, S, CA);
		normalizeC(nBasis, S, CB);
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);
	}

	//
	// set density matrix to diagonal as a guess
	//
	if(opt->SCFGuess == SCFGUESS_DIAG){
		fflush(stdout);
		printf("Use identity matrix as initial density matrix ...\n");
		sum=0.0;
		for(i=0; i<nBasis; i++) sum+=S[i*nBasis+i];
		for(i=0; i<nBasis; i++)
		for(j=0; j<nBasis; j++){
			if(i==j){
				PA[i*nBasis+j] = nEA/sum;
				PB[i*nBasis+j] = nEB/sum;
			}else{
				PA[i*nBasis+j] = 0.0;
				PB[i*nBasis+j] = 0.0;
			}
		}
	}

	//
	// read orbital from checkpoint
	//
	if(opt->SCFGuess == SCFGUESS_CHECK){
		fflush(stdout);
		printf("Read orbital from checkpoint for initial density matrix ...\n");
		guess_checkpoint_orbital(nBasis, gto, mol, opt, S, CA, CB);
		normalizeC(nBasis, S, CA);
		normalizeC(nBasis, S, CB);
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);
	}

	// load density matrix if requested
	if(opt->loadDMatrix){

		// openfile
		fd=fopen(opt->DMatrixFile,"r");
		if(fd==NULL){
			printf("uhf - error cannot open file %s for reading\n", 
			       opt->DMatrixFile);
			exit(-1);
		}

		//////////////////
		// loop and read
		//////////////////
		// handle RHF case
		if(opt->RHF){
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				if(fscanf(fd, "%lf", PA+i*nBasis+j)!=1){
					printf("uhf - error reading density matrix\n");
					exit(-1);
				}
				PB[i*nBasis+j] = PA[i*nBasis+j];
			}
		}
		// handle UHF case - read 2 times
		if(opt->UHF){
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				if(fscanf(fd, "%lf", PA+i*nBasis+j)!=1){
					printf("uhf - error reading alpha density matrix\n");
					exit(-1);
				}
			}
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				if(fscanf(fd, "%lf", PB+i*nBasis+j)!=1){
					printf("uhf - error reading beta spin density matrix\n");
					exit(-1);
				}
			}
		}

		// close file
		fclose(fd);

		// report status
		printf("Density matrix loaded from %s\n", opt->DMatrixFile);
	}

	// compute cutoff
#define GETCUTOFF(accuracy) ((accuracy)/(sumA*sumA+sumA*sumB+sumB*sumB)/50.0)
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		sumA += PA[i*nBasis+j];
		sumB += PB[i*nBasis+j];
	}
	realCutoff = GETCUTOFF(opt->SCFConv);
	opt->SCFCutoff = realCutoff;

	// 2-electron integral information
	printf("Processing 2E integrals ...\n");
	printf("Schwarz inequality screening cut-off %.8E\n", opt->SCFCutoff);
	printf("Primitive prefactor cut-off %0.8E\n", PRIMITIVE_CUTOFF);
	Schwarz = create_Schwarz(nBasis, gto);

	// scf loop
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   printf("Use 4-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DIIS3:   printf("Use 3-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DIIS2:   printf("Use 2-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DAMPING: printf("Use simple weighting convergence method\n"); break;
	default:
		printf("uhf - error no specific request for convergence method\n");
		exit(-1);
	break;
	}

	printf("Drag coefficient %f\n", opt->SCFDrag);
	printf("SCFConv %.2E\n", opt->SCFConv);
	printf("SCFMax %d iterations\n", opt->SCFMax);
	printf("Enter SCF loop ... \n");
	printf(
	"Iteration  Total Energy [Hartrees]  RMSD Density\n"
	"------------------------------------------------\n");
	fflush(stdout);
	
	// oscillation drag coefficient
	gamma = opt->SCFDrag;

	// call convergence function for the first time to initialize it
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   conv_diis4(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DIIS3:   conv_diis3(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DIIS2:   conv_diis2(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DAMPING: conv_damping(nBasis, gamma, PA, PB); break;
	default:
		printf("uhf - error unknown opt->convMethod\n");
		exit(-1);
	break;
	}

	// preparations
	for(i=0; i<nBasis; i++)
	for(j=0; j<nBasis; j++){

		// set G matrix to zero
		GA[i*nBasis+j] = 0.0;
		GB[i*nBasis+j] = 0.0;

		// set delta density to the initial density
		dPA[i*nBasis+j] = PA[i*nBasis+j];
		dPB[i*nBasis+j] = PB[i*nBasis+j];
	}

	// manage integral accuracy
	switch(opt->SCFAccuracy){
	case SCFACCURACY_1STEP:
		thisCutoff = realCutoff;
	break;
	case SCFACCURACY_3STEP:
		cutoffA    = GETCUTOFF(SCFACCURACY_3STEP_A);
		cutoffB    = GETCUTOFF(SCFACCURACY_3STEP_B);
		cutoffC    = GETCUTOFF(SCFACCURACY_3STEP_C);
		thisCutoff = cutoffA;
	break;
	default:
		printf("uhf - error unknown SCFAccuracy\n");
		exit(-1);
	break;
	}

	do{
		iter++;

		// compute delta G matrix
		if(opt->UHF)
			uhf_getGMatrix(nBasis, gto, Schwarz, thisCutoff, dPA, dPB, dGA, dGB, opt);
		if(opt->RHF)
			uhf_getGMatrix(nBasis, gto, Schwarz, thisCutoff, dPA, dPA, dGA, dGB, opt);

		// updates
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			// update G matrix
			GA[i*nBasis+j] = GA[i*nBasis+j] + dGA[i*nBasis+j];
			GB[i*nBasis+j] = GB[i*nBasis+j] + dGB[i*nBasis+j];

			// update fock matrix
			FA[i*nBasis+j] = H[i*nBasis+j] + GA[i*nBasis+j];
			FB[i*nBasis+j] = H[i*nBasis+j] + GB[i*nBasis+j];

			// saving current density matrix to dPA and dPB
			dPA[i*nBasis+j] = PA[i*nBasis+j];
			dPB[i*nBasis+j] = PB[i*nBasis+j];
		}

		// solve generalized eigen value problem and normalize orbital
		gen_sym_eigen(nBasis, FA, S, eA, CA);
		gen_sym_eigen(nBasis, FB, S, eB, CB);
		normalizeC(nBasis, S, CA);
		normalizeC(nBasis, S, CB);

		// get new P matrix
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);

		// compute energy and energ difference
		dE     = Etot;
		Etot   = uhf_getEtotal(nBasis, mol, PA, FA, PB, FB, H);
		dE     = Etot - dE;

		// update P matrix using convergence method
		switch(opt->convMethod){
		case CONVMETHOD_DIIS4:   avgdP = conv_diis4(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DIIS3:   avgdP = conv_diis3(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DIIS2:   avgdP = conv_diis2(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DAMPING: avgdP = conv_damping(nBasis, gamma, PA, PB); break;
		}

		// check convergence
		notConverged = fabs(dE) > opt->SCFConv || avgdP    > opt->SCFConv;

		// compute delta density matrix for the next step
		for(i=0; i<nBasis; i++)
		for(j=0; j<nBasis; j++){
			dPA[i*nBasis+j] = PA[i*nBasis+j] - dPA[i*nBasis+j];
			dPB[i*nBasis+j] = PB[i*nBasis+j] - dPB[i*nBasis+j];
		}

		printf(" %5d %20.8f %20.4E\n", iter, Etot, avgdP);
		
		// check if we have reached scfmax limit
		if(iter >= opt->SCFMax) break;

		// flush output
		fflush(stdout);

		// manage integral accuracy
		if(opt->SCFAccuracy == SCFACCURACY_3STEP){
			
			// check convergence
			if(!notConverged) break;

			// check if we need to switch accuracy
			if((thisCutoff == cutoffA && fabs(dE) <= SCFACCURACY_3STEP_A && avgdP <= SCFACCURACY_3STEP_A) ||
			   (thisCutoff == cutoffB && fabs(dE) <= SCFACCURACY_3STEP_B && avgdP <= SCFACCURACY_3STEP_B) ||
			   (thisCutoff == cutoffC && fabs(dE) <= SCFACCURACY_3STEP_C && avgdP <= SCFACCURACY_3STEP_C)){

				// switch accuracy
				     if(thisCutoff == cutoffA) thisCutoff = cutoffB;
				else if(thisCutoff == cutoffB) thisCutoff = cutoffC;
				else if(thisCutoff == cutoffC) thisCutoff = realCutoff;

				// rebuilding fock matrix
				//for(i=0; i < nBasis; i++)
				//for(j=0; j < nBasis; j++){
				//	 GA[i*nBasis+j] = 0.0;
				//	 GB[i*nBasis+j] = 0.0;
				//	dPA[i*nBasis+j] = PA[i*nBasis+j];
				//	dPB[i*nBasis+j] = PB[i*nBasis+j];  
				//}
				printf("................. switch accuracy ..............\n");
				//printf("........ switch accuracy and rebuild fock matrix\n");
				fflush(stdout);

			}

		}

		// save checkpoint if requested
		if(opt->saveCheckAll){
			//printf("\n[Begin] saving checkpoint file to %s\n", opt->CheckFile);
			save_checkpoint(nBasis, gto, mol, Etot, dE, avgdP, CA, CB, eA, eB, opt);
			//printf("[Done]  saving checkpoint file to %s\n\n", opt->CheckFile);
			//fflush(stdout);
		}

	}while(notConverged);

	// report
	if(notConverged){
		printf("SCF have not converged because iter >= SCFMax\n");
	}else{
		printf("Done SCF Loop Total Energy is %20.8f Hartrees\n", Etot);
	}
	fflush(stdout);

	// save density matrix if requested
	if(opt->saveDMatrix){

		// openfile
		fd=fopen(opt->DMatrixFile,"w");
		if(fd==NULL){
			printf("uhf - error cannot open file %s for writing\n", 
			       opt->DMatrixFile);
		}

		///////////////////////////////////////////////////
		// loop and save each element, always write twice
		///////////////////////////////////////////////////
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			if((i*nBasis+j)%5==0) fprintf(fd,"\n");
			fprintf(fd," %15.8E ", PA[i*nBasis+j]);
		}
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			if((i*nBasis+j)%5==0) fprintf(fd,"\n");
			fprintf(fd," %15.8E ", PB[i*nBasis+j]);
		}
		
		// close file
		fclose(fd);

		// report status
		printf("Density matrix saved to %s\n", opt->DMatrixFile);
	}

	// clear convergence routine
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   conv_diis4(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DIIS3:   conv_diis3(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DIIS2:   conv_diis2(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DAMPING: conv_damping(0, 0.0, NULL, NULL); break;
	}

	// clean memory
	free(Schwarz);
	free(PA);
	free(GA);
	free(dPA);
	free(dGA);
	free(FA);
	free(PB);
	free(GB);
	free(dPB);
	free(dGB);
	free(FB);
	free(H);
	free(T);
	free(V);
	free(S);

	// clean electron storage in JK integral subroutine
	//GTO_JK_Matrix_ShellSet(0, NULL, NULL, NULL, NULL, 0.0, NULL, NULL, NULL);
	//GTO_JK_Matrix_Quartet(0, NULL, NULL, NULL, NULL, 0.0, NULL, NULL, NULL);

	//
	// parallel version of cleaning 
	//
	int *status;             // status for each cpu
	int alldone;             // all idle flag

	// allocate memory
	status=calloc(opt->nCPU,sizeof(int));
	if(status==NULL){
		printf("uhf - error cannot allocate memory\n");
		exit(-1);
	}

	// reset status to idle
	for(i=(opt->nCPU-1);i>=0;i--) status[i] = RPC_IDLE;

	// clean all child processes
	do{
		// remote call to all cpu except childID=0
		for(i=(opt->nCPU-1);i>0;i--)
			status[i] = rpc_GTO_JK_Matrix_Quartet_Parallel(status[i],i,0,
			                                               NULL, NULL, NULL, NULL,
			                                               0.0, NULL, NULL, opt);
		// local call for childID=0
		if(status[0]==RPC_IDLE){
			GTO_JK_Matrix_Quartet_Parallel(0, 0, NULL, NULL, NULL, NULL,
			                               0.0, NULL, NULL, opt);
			status[0] = RPC_DONE;
		}

		// check if all done
		alldone=1;
		for(i=(opt->nCPU-1);i>=0;i--) if(status[i] != RPC_DONE) alldone=0;

	}while(!alldone);

	// clean memory
	free(status);

	// save checkpoint if requested
	if(opt->saveCheck){
		printf("\n[Begin] saving checkpoint file to %s\n", opt->CheckFile);
		save_checkpoint(nBasis, gto, mol, Etot, dE, avgdP, CA, CB, eA, eB, opt);
		printf("[Done]  saving checkpoint file to %s\n\n", opt->CheckFile);
		fflush(stdout);
	}

	if(notConverged) return 0.0;
	else             return Etot;
}
