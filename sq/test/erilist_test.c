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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/int.h"

struct gto_t{
double x,y,z; // center position
double exp;   // exponent
int    l,m,n; // angular index
};

// initialize the Gaussian orbial using 
// random function
void random_gto(struct gto_t *gto){
#define MAX_R   5.0 
#define MIN_EXP 0.5
#define MAX_EXP 2.0
#define MAX_L   3.0
	gto->x   = MAX_R*rand()/RAND_MAX;
	gto->y   = MAX_R*rand()/RAND_MAX;
	gto->z   = MAX_R*rand()/RAND_MAX;
	gto->exp = MAX_EXP*rand()/RAND_MAX + MIN_EXP;
	do{
		gto->l   = (int)((MAX_L+1)*rand()/RAND_MAX);
		gto->m   = (int)((MAX_L+1)*rand()/RAND_MAX);
		gto->n   = (int)((MAX_L+1)*rand()/RAND_MAX);
	}while(((gto->l)+(gto->m)+(gto->n))>MAX_L);
	return;
}

// int_test : using PyQuante by Richard Muller int package
// as reference integrals.
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
#define MAXITEM 4096
int main(int argc, char *argv[]){
	double max_err     = 0.0;
	double max_percent = 0.0;
	unsigned long long n,i;
	double err;
	double percent;
	struct gto_t a,b,c,d;     // using 4 gto
	int la[MAXITEM], ma[MAXITEM], na[MAXITEM];
	int lb[MAXITEM], mb[MAXITEM], nb[MAXITEM];
	int lc[MAXITEM], mc[MAXITEM], nc[MAXITEM];
	int ld[MAXITEM], md[MAXITEM], nd[MAXITEM];
	double results[MAXITEM], eriresult;
	char str[256];

	printf("----------------------------------------------------\n"
	       "<================  Accuracy  Testing  =============>\n"
	       "                                                    \n");

	
	// electron repulsion integrals testing
	printf("Testing electron repulsion integral (ab|1/r12|cd)\n"
	       "1,000,000 integrals generated randomly.\n");
	max_err     = 0.0;
	max_percent = 0.0;
	for(n=0; n < 1000; n++){
		random_gto(&a);
		random_gto(&b);
		random_gto(&c);
		random_gto(&d);

		// load eri_list
		la[0] = a.l; ma[0] = a.m; na[0] = a.n;
		lb[0] = b.l; mb[0] = b.m; nb[0] = b.n;
		lc[0] = c.l; mc[0] = c.m; nc[0] = c.n;
		ld[0] = d.l; md[0] = d.m; nd[0] = d.n;
		eri_list(a.x,a.y,a.z,1.0,la,ma,na,a.exp,
		         b.x,b.y,b.z,1.0,lb,mb,nb,b.exp,
		         c.x,c.y,c.z,1.0,lc,mc,nc,c.exp,
		         d.x,d.y,d.z,1.0,ld,md,nd,d.exp,
		         1,results);
		eriresult = eri(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,c.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,d.m,d.n,d.exp);

		err = fabs(results[0] - eriresult);
		percent = err*100/eriresult;
		percent = fabs(percent);

		if(err > max_err){
			max_err = err;
			max_percent = percent;
			printf("   max error is %18.6E  (%18.6E Percent )"
			       " @ (lmn)(lmn)(lmn)(lmn) = (%d%d%d)(%d%d%d)(%d%d%d)(%d%d%d)\n",
			       max_err,max_percent,a.l,a.m,a.n,
			                           b.l,b.m,b.n,
			                           c.l,c.m,c.n,
			                           d.l,d.m,d.n);
		}
	}
	printf("Done testing: max error is %18.6E (%18.6E Percent )\n"
	       "-------------------\n\n",
	                                 max_err,max_percent);
	
	
	printf("----------------------------------------------------\n"
	       "<================ Performance Testing =============>\n"
	       "                                                    \n");

	time_str(256,str);
	printf("%s\n",str);

	for(i=0; i < 40000; i++){
	random_gto(&a);
	random_gto(&b);
	random_gto(&c);
	random_gto(&d);

	n = 0;
	for(a.l = 0; a.l <= 1; a.l++)
	for(a.m = 0; a.m <= 1; a.m++)
	for(a.n = 0; a.n <= 1; a.n++)
	for(b.l = 0; b.l <= 1; b.l++)
	for(b.m = 0; b.m <= 1; b.m++)
	for(b.n = 0; b.n <= 1; b.n++)
	for(c.l = 0; c.l <= 1; c.l++)
	for(c.m = 0; c.m <= 1; c.m++)
	for(c.n = 0; c.n <= 1; c.n++)
	for(d.l = 0; d.l <= 1; d.l++)
	for(d.m = 0; d.m <= 1; d.m++)
	for(d.n = 0; d.n <= 1; d.n++){
		if(a.l+a.m+a.n <= 1 && 
		   b.l+b.m+b.n <= 1 && 
		   c.l+c.m+c.n <= 1 &&
		   d.l+d.m+d.n <= 1){

			// store integral list
			la[n] = a.l; ma[n] = a.m; na[n] = a.n;
			lb[n] = b.l; mb[n] = b.m; nb[n] = b.n;
			lc[n] = c.l; mc[n] = c.m; nc[n] = c.n;
			ld[n] = d.l; md[n] = d.m; nd[n] = d.n;

// NORMAL ERI 
		eriresult = eri(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,c.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,d.m,d.n,d.exp);

			n++;
		}
	}
//  ERI_LIST TEST
//	eri_list(a.x,a.y,a.z,1.0,la,ma,na,a.exp,
//	         b.x,b.y,b.z,1.0,lb,mb,nb,b.exp,
//		     c.x,c.y,c.z,1.0,lc,mc,nc,c.exp,
//		     d.x,d.y,d.z,1.0,ld,md,nd,d.exp,
//		     n,results);

	}

	time_str(256,str);
	printf("%s\n",str);

	/*
	printf("Computing 1,000,000 (SS|SS)\n");
	for(n=0; n < 1000000; n++){
		eri(
		    a.x,a.y,a.z,1.0,0,0,0,a.exp,
		    b.x,b.y,b.z,1.0,0,0,0,b.exp,
		    c.x,c.y,c.z,1.0,0,0,0,c.exp,
		    d.x,d.y,d.z,1.0,0,0,0,d.exp);
	}
	printf("Computing 1,000,000 (PP|PP)\n");
	for(n=0; n < 1000000; n++){
		eri(
		    a.x,a.y,a.z,1.0,1,0,0,a.exp,
		    b.x,b.y,b.z,1.0,1,0,0,b.exp,
		    c.x,c.y,c.z,1.0,1,0,0,c.exp,
		    d.x,d.y,d.z,1.0,1,0,0,d.exp);
	}
	*/
}
