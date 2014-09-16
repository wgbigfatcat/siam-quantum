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
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "fgamma.h"
#include "util.h" 

/**********************************************************************
 * cints.c  C implementation of simple math functions in pyutil.
 *
 * The equations herein are based upon
 * 'Gaussian Expansion Methods for Molecular Orbitals.' H. Taketa,
 * S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
 * [THO paper].
 *
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
 **********************************************************************/
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define SMALL 0.00000001

static void gser(double *gamser, double a, double x, double *gln){
  int n;
  double sum,del,ap;

  *gln=lgamma(a);
  if (x <= 0.0) {
    assert(x>=0.);
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    printf("a too large, ITMAX too small in routine gser");
    return;
  }
}
 
static void gcf(double *gammcf, double a, double x, double *gln){
  int i;
  double an,b,c,d,del,h;
  
  *gln=lgamma(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  assert(i<=ITMAX);
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}


static double gamm_inc(double a, double x){ /* Taken from NR routine gammap */
  double gamser,gammcf,gln;
  
  assert (x >= 0.);
  assert (a > 0.);
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return exp(gln)*gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return exp(gln)*(1.0-gammcf);
  }
}

static double fgamma_pyquante(double m, double x){
  double val;
  if (fabs(x) < SMALL) x = SMALL;
  val = gamm_inc(m+0.5,x);
  /* if (val < SMALL) return 0.; */ /* Gives a bug for D orbitals. */
  return 0.5*pow(x,-m-0.5)*val; 
}

// fgamma_test : testing fgamma implementation
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
int main(int argc, char *argv[]){
	int m,i,j,k,l;
	double maxerr, err;
	double fset[20];
	char str[256];

	printf("Testing FGamma using PyQuante by Richard Muller as references\n");
	printf("Case 1: double fgamma(int m, double x);\n");
	for(m=0; m <= 15; m++){
		maxerr = 0.0;
		for(i=0; i < 100000; i++){
			err = fgamma_pyquante(m,i/1000.)-fgamma(m,i/1000.);
			err = fabs(err);
			if(err > maxerr) maxerr = err;
		}
		printf("m = %2d maxerr = %20.8E\n", m, maxerr);
	}

	printf("Case 2: void fgamma_set(int max, double x);\n");
	for(m=0; m <= 15; m++){
		maxerr = 0.0;
		for(i=0; i < 100000; i++){
			fgamma_set(19, i/1000., fset);
			err = fgamma_pyquante(m,i/1000.)-fset[m];
			err = fabs(err);
			if(err > maxerr) maxerr = err;
		}
		printf("m = %2d maxerr = %20.8E\n", m, maxerr);
	}

	printf("Evaluate fgamma(0,...) 100,000,000 times\n");
	time_str(256,str);
	printf("Start: %s\n", str);
	m   = 100;
	err = 0.0;
	for(i=0; i < m; i++)
	for(j=0; j < m; j++)
	for(k=0; k < m; k++)
	for(l=0; l < m; l++){
		err += fgamma(0,err);
	}
	time_str(256,str);
	printf("Done:  %s\n", str);
}
