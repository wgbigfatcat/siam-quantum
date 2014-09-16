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
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/fgamma.h"
#include "../src/int.h"
#include "../src/basis.h"
#include "../src/intmd.h"

/*************************************************************************
 *
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
 **************************************************************************/
/* My routines */
static double fB(int i, int l1, int l2, double px, double ax, double bx, 
	  int r, double g);
static double Bfunc(int i, int r, double g);
double contr_coulomb(int ia, double *aexps, double *acoefs, double *anorms,
			    double xa, double ya, double za, int la, int ma, int na, 
			    int ib, double *bexps, double *bcoefs, double *bnorms,
			    double xb, double yb, double zb, int lb, int mb, int nb, 
			    int ic, double *cexps, double *ccoefs, double *cnorms,
			    double xc, double yc, double zc, int lc, int mc, int nc, 
			    int id, double *dexps, double *dcoefs, double *dnorms,
			    double xd, double yd, double zd, int ld, int md, int nd);

double coulomb_repulsion(double xa, double ya, double za, double norma,
				int la, int ma, int na, double alphaa,
				double xb, double yb, double zb, double normb,
				int lb, int mb, int nb, double alphab,
				double xc, double yc, double zc, double normc,
				int lc, int mc, int nc, double alphac,
				double xd, double yd, double zd, double normd,
				int ld, int md, int nd, double alphad);

static double *B_array(int l1, int l2, int l3, int l4, double p, double a,
		double b, double q, double c, double d,
		double g1, double g2, double delta);

static double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
		     int l3, int l4, double Px, double Ax, double Bx,
		     double Qx, double Cx, double Dx, double gamma1,
		     double gamma2, double delta);
double kinetic_py(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb);
double overlap_py(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb);
static double overlap_1D(int l1, int l2, double PAx,
			 double PBx, double gamma);
double nuclear_attraction(double x1, double y1, double z1, double norm1,
				 int l1, int m1, int n1, double alpha1,
				 double x2, double y2, double z2, double norm2,
				 int l2, int m2, int n2, double alpha2,
				 double x3, double y3, double z3);
static double A_term(int i, int r, int u, int l1, int l2,
		     double PAx, double PBx, double CPx, double gamma);
static double *A_array(int l1, int l2, double PA, double PB,
		       double CP, double g);

static int fact(int n);
static int fact2(int n);
static double dist2(double x1, double y1, double z1, 
		    double x2, double y2, double z2);
static double dist(double x1, double y1, double z1, 
		   double x2, double y2, double z2);
static double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb);
static int binomial(int a, int b);

double Fgamma(double m, double x);
static double gamm_inc(double a, double x);

static int ijkl2intindex(int i, int j, int k, int l);

static int fact_ratio2(int a, int b);

static double product_center_1D(double alphaa, double xa, 
			 double alphab, double xb);

static double three_center_1D(double xi, int ai, double alphai,
			      double xj, int aj, double alphaj,
			      double xk, int ak, double alphak);

/* Routines from Numerical Recipes */
static void gser(double *gamser, double a, double x, double *gln);
static void gcf(double *gammcf, double a, double x, double *gln);


#if defined(_WIN32)
double lgamma(double x);
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define SMALL 0.00000001

static double fB(int i, int l1, int l2, double px, double ax, double bx, 
		 int r, double g){
  return binomial_prefactor(i,l1,l2,px-ax,px-bx)*Bfunc(i,r,g);
}

static double Bfunc(int i, int r, double g){
  return fact_ratio2(i,r)*pow(4*g,r-i);
}

double contr_coulomb(int lena, double *aexps, double *acoefs,
			    double *anorms, double xa, double ya, double za,
			    int la, int ma, int na, 
			    int lenb, double *bexps, double *bcoefs,
			    double *bnorms, double xb, double yb, double zb,
			    int lb, int mb, int nb, 
			    int lenc, double *cexps, double *ccoefs,
			    double *cnorms, double xc, double yc, double zc,
			    int lc, int mc, int nc, 
			    int lend, double *dexps, double *dcoefs,
			    double *dnorms, double xd, double yd, double zd,
			    int ld, int md, int nd){

  int i,j,k,l;
  double Jij = 0.,incr=0.;

  for (i=0; i<lena; i++)
    for (j=0; j<lenb; j++)
      for (k=0; k<lenc; k++)
	for (l=0; l<lend; l++){
	  incr = coulomb_repulsion(xa,ya,za,anorms[i],la,ma,na,aexps[i],
			      xb,yb,zb,bnorms[j],lb,mb,nb,bexps[j],
			      xc,yc,zc,cnorms[k],lc,mc,nc,cexps[k],
			      xd,yd,zd,dnorms[l],ld,md,nd,dexps[l]);
	  
	  Jij += acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]*incr;
	}
  return Jij;
}

double coulomb_repulsion(double xa, double ya, double za, double norma,
				int la, int ma, int na, double alphaa,
				double xb, double yb, double zb, double normb,
				int lb, int mb, int nb, double alphab,
				double xc, double yc, double zc, double normc,
				int lc, int mc, int nc, double alphac,
				double xd, double yd, double zd, double normd,
				int ld, int md, int nd, double alphad){

  double rab2, rcd2,rpq2,xp,yp,zp,xq,yq,zq,gamma1,gamma2,delta,sum;
  double *Bx, *By, *Bz;
  int I,J,K;

  rab2 = dist2(xa,ya,za,xb,yb,zb);
  rcd2 = dist2(xc,yc,zc,xd,yd,zd);
  xp = product_center_1D(alphaa,xa,alphab,xb);
  yp = product_center_1D(alphaa,ya,alphab,yb);
  zp = product_center_1D(alphaa,za,alphab,zb);
  xq = product_center_1D(alphac,xc,alphad,xd);
  yq = product_center_1D(alphac,yc,alphad,yd);
  zq = product_center_1D(alphac,zc,alphad,zd);
  rpq2 = dist2(xp,yp,zp,xq,yq,zq);
  gamma1 = alphaa+alphab;
  gamma2 = alphac+alphad;
  delta = (1./gamma1+1./gamma2)/4.;

  Bx = B_array(la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta);
  By = B_array(ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta);
  Bz = B_array(na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta);

  sum = 0.;
  for (I=0; I<la+lb+lc+ld+1;I++)
    for (J=0; J<ma+mb+mc+md+1;J++)
      for (K=0; K<na+nb+nc+nd+1;K++)
	sum += Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,0.25*rpq2/delta);

  free(Bx);
  free(By);
  free(Bz);  
  
  return 2.*pow(M_PI,2.5)/(gamma1*gamma2*sqrt(gamma1+gamma2))
    *exp(-alphaa*alphab*rab2/gamma1) 
    *exp(-alphac*alphad*rcd2/gamma2)*sum*norma*normb*normc*normd;
}

static double *B_array(int l1, int l2, int l3, int l4, double p, double a,
		double b, double q, double c, double d,
		double g1, double g2, double delta){
  int Imax,i1,i2,r1,r2,u,I,i;
  double *B;
  Imax = l1+l2+l3+l4+1;
  B = (double *)malloc(Imax*sizeof(double));
  for (i=0; i<Imax; i++) B[i] = 0.;

  for (i1=0; i1<l1+l2+1; i1++)
    for (i2=0; i2<l3+l4+1; i2++)
      for (r1=0; r1<i1/2+1; r1++)
	for (r2=0; r2<i2/2+1; r2++)
	  for (u=0; u<(i1+i2)/2-r1-r2+1; u++){
	    I = i1+i2-2*(r1+r2)-u;
	    B[I] = B[I] + B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
				 p,a,b,q,c,d,g1,g2,delta);
	  }

  return B;
}

static double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
	      int l3, int l4, double Px, double Ax, double Bx,
	      double Qx, double Cx, double Dx, double gamma1,
	      double gamma2, double delta){
  /* THO eq. 2.22 */
  return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)
    *pow(-1,i2)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)
    *pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u)
    *pow(Qx-Px,i1+i2-2*(r1+r2)-2*u)
    /pow(delta,i1+i2-2*(r1+r2)-u);
}


double kinetic_py(double alpha1, int l1, int m1, int n1,
	       double xa, double ya, double za,
	       double alpha2, int l2, int m2, int n2,
	       double xb, double yb, double zb){

  double term0,term1,term2;
  term0 = alpha2*(2*(l2+m2+n2)+3)*
    overlap_py(alpha1,l1,m1,n1,xa,ya,za,
		   alpha2,l2,m2,n2,xb,yb,zb);
  term1 = -2*pow(alpha2,2)*
    (overlap_py(alpha1,l1,m1,n1,xa,ya,za,
		    alpha2,l2+2,m2,n2,xb,yb,zb)
     + overlap_py(alpha1,l1,m1,n1,xa,ya,za,
		      alpha2,l2,m2+2,n2,xb,yb,zb)
     + overlap_py(alpha1,l1,m1,n1,xa,ya,za,
		      alpha2,l2,m2,n2+2,xb,yb,zb));
  term2 = -0.5*(l2*(l2-1)*overlap_py(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2-2,m2,n2,xb,yb,zb) +
		m2*(m2-1)*overlap_py(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2,m2-2,n2,xb,yb,zb) +
		n2*(n2-1)*overlap_py(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2,m2,n2-2,xb,yb,zb));
  return term0+term1+term2;
}

double overlap_py(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb){
  /*Taken from THO eq. 2.12*/
  double rab2,gamma,xp,yp,zp,pre,wx,wy,wz;

  rab2 = dist2(xa,ya,za,xb,yb,zb);
  gamma = alpha1+alpha2;
  xp = product_center_1D(alpha1,xa,alpha2,xb);
  yp = product_center_1D(alpha1,ya,alpha2,yb);
  zp = product_center_1D(alpha1,za,alpha2,zb);

  pre = pow(M_PI/gamma,1.5)*exp(-alpha1*alpha2*rab2/gamma);

  wx = overlap_1D(l1,l2,xp-xa,xp-xb,gamma);
  wy = overlap_1D(m1,m2,yp-ya,yp-yb,gamma);
  wz = overlap_1D(n1,n2,zp-za,zp-zb,gamma);
  return pre*wx*wy*wz;
}

static double overlap_1D(int l1, int l2, double PAx,
			 double PBx, double gamma){
  /*Taken from THO eq. 2.12*/
  int i;
  double sum;
  sum = 0.;
  for (i=0; i<(1+floor(0.5*(l1+l2))); i++)
    sum += binomial_prefactor(2*i,l1,l2,PAx,PBx)* 
      fact2(2*i-1)/pow(2*gamma,i);
  return sum;
}
    
double nuclear_attraction(double x1, double y1, double z1, double norm1,
				 int l1, int m1, int n1, double alpha1,
				 double x2, double y2, double z2, double norm2,
				 int l2, int m2, int n2, double alpha2,
				 double x3, double y3, double z3){
  int I,J,K;
  double gamma,xp,yp,zp,sum,rab2,rcp2;
  double *Ax,*Ay,*Az;

  gamma = alpha1+alpha2;

  xp = product_center_1D(alpha1,x1,alpha2,x2);
  yp = product_center_1D(alpha1,y1,alpha2,y2);
  zp = product_center_1D(alpha1,z1,alpha2,z2);

  rab2 = dist2(x1,y1,z1,x2,y2,z2);
  rcp2 = dist2(x3,y3,z3,xp,yp,zp);

  Ax = A_array(l1,l2,xp-x1,xp-x2,xp-x3,gamma);
  Ay = A_array(m1,m2,yp-y1,yp-y2,yp-y3,gamma);
  Az = A_array(n1,n2,zp-z1,zp-z2,zp-z3,gamma);

  sum = 0.;
  for (I=0; I<l1+l2+1; I++)
    for (J=0; J<m1+m2+1; J++)
      for (K=0; K<n1+n2+1; K++)
	sum += Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K,rcp2*gamma);

  free(Ax);
  free(Ay);
  free(Az);
  return -norm1*norm2*
    2*M_PI/gamma*exp(-alpha1*alpha2*rab2/gamma)*sum;
}
    
static double A_term(int i, int r, int u, int l1, int l2,
		     double PAx, double PBx, double CPx, double gamma){
  /* THO eq. 2.18 */
  return pow(-1,i)*binomial_prefactor(i,l1,l2,PAx,PBx)*
    pow(-1,u)*fact(i)*pow(CPx,i-2*r-2*u)*
    pow(0.25/gamma,r+u)/fact(r)/fact(u)/fact(i-2*r-2*u);
}

static double *A_array(int l1, int l2, double PA, double PB,
		double CP, double g){
  /* THO eq. 2.18 and 3.1 */
  int Imax,i,r,u,I;
  double *A;

  Imax = l1+l2+1;
  A = (double *)malloc(Imax*sizeof(double));
  for (i=0; i<Imax; i++) A[i] = 0.;
  for (i=0; i<Imax; i++)
    for (r=0; r<floor(i/2)+1;r++)
      for (u=0; u<floor((i-2*r)/2.)+1; u++){
	I = i-2*r-u;
	A[I] += A_term(i,r,u,l1,l2,PA,PB,CP,g);
      }
  return A;
}


static int fact(int n){
  if (n <= 1) return 1;
  return n*fact(n-1);
}

static int fact2(int n){ /* double factorial function = 1*3*5*...*n */
  if (n <= 1) return 1;
  return n*fact2(n-2);
}

static double dist2(double x1, double y1, double z1,
		    double x2, double y2, double z2){
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}
static double dist(double x1, double y1, double z1,
		   double x2, double y2, double z2){
  return sqrt(dist2(x1,y1,z1,x2,y2,z2));
}

static double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb){
  int t;
  double sum=0.;
  for (t=0; t<s+1; t++)
    if ((s-ia <= t) && (t <= ib)) 
      sum += binomial(ia,s-t)*binomial(ib,t)*pow(xpa,ia-s+t)*pow(xpb,ib-t);
  return sum;
} 

static int binomial(int a, int b){return fact(a)/(fact(b)*fact(a-b));}

double Fgamma(double m, double x){
  double val;
  if (fabs(x) < SMALL) x = SMALL;
  val = gamm_inc(m+0.5,x);
  /* if (val < SMALL) return 0.; */ /* Gives a bug for D orbitals. */
  return 0.5*pow(x,-m-0.5)*val; 
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

static int ijkl2intindex(int i, int j, int k, int l){
  int tmp,ij,kl;
  if (i<j){
    tmp = i;
    i = j;
    j = tmp;
  }
  if (k<l){
    tmp = k;
    k = l;
    l = tmp;
  }
  ij = i*(i+1)/2+j;
  kl = k*(k+1)/2+l;
  if (ij<kl){
    tmp = ij;
    ij = kl;
    kl = tmp;
  }
  return ij*(ij+1)/2+kl;
}

static int fact_ratio2(int a, int b){ return fact(a)/fact(b)/fact(a-2*b); }

static double product_center_1D(double alphaa, double xa, 
			 double alphab, double xb){
  return (alphaa*xa+alphab*xb)/(alphaa+alphab);
}

static double three_center_1D(double xi, int ai, double alphai,
			      double xj, int aj, double alphaj,
			      double xk, int ak, double alphak){

  double gamma, dx, px, xpi,xpj,xpk,intgl;
  int q,r,s,n;
  
  gamma = alphai+alphaj+alphak;
  dx = exp(-alphai*alphaj*pow(xi-xj,2)/gamma) *
    exp(-alphai*alphak*pow(xi-xk,2)/gamma) *
    exp(-alphaj*alphak*pow(xj-xk,2)/gamma);
  px = (alphai*xi+alphaj*xj+alphak*xk)/gamma;
    
  xpi = px-xi;
  xpj = px-xj;
  xpk = px-xk;
  intgl = 0;
  for (q=0; q<ai+1; q++){
    for (r=0; r<aj+1; r++){
      for (s=0; s<ak+1; s++){
	n = (q+r+s)/2;
	if ((q+r+s)%2 == 0) {
	  intgl += binomial(ai,q)*binomial(aj,r)*binomial(ak,s)*
	    pow(xpi,ai-q)*pow(xpj,aj-r)*pow(xpk,ak-s)*
	    fact2(2*n-1)/pow(2*gamma,n)*sqrt(M_PI/gamma);
	}
      }
    }
  }
  return dx*intgl;
}

#undef ITMAX
#undef EPS
#undef FPMIN

#define MAXL 5
struct eriOSData_t{
	double alphaa,alphab; // exponent for a and b
	double alphac,alphad; // exponent for c and d
	double zeta;          // alphaa+alphab
	double eta;           // alphac+alphad
	double xa,ya,za;      // center for a
	double xb,yb,zb;      // center for b
	double xc,yc,zc;      // center for c
	double xd,yd,zd;      // center for d
	double xp,yp,zp;      // P = (A*alphaa+B*alphab)/zeta
	double xq,yq,zq;      // Q = (C*alphac+D*alphad)/eta
	double xw,yw,zw;      // W = (P*zeta  +   Q*eta)/(zeta+eta) 
	double ssss[4*MAXL+1];// ssss integrals
};

double eriOS_R(int la, int ma, int na,
               int lb, int mb, int nb,
               int lc, int mc, int nc,
               int ld, int md, int nd,
               int m, struct eriOSData_t *info){

	if(la<0 || ma<0 || na<0 || 
	   lb<0 || mb<0 || nb<0 ||
	   lc<0 || mc<0 || nc<0 ||
	   ld<0 || md<0 || nd<0) return 0.0;

	else if(   la==0 && ma==0 && na==0 
	        && lb==0 && mb==0 && nb==0
	        && lc==0 && mc==0 && nc==0
	        && ld==0 && md==0 && nd==0) return info->ssss[m];

	else if(la>0) return 
	   eriOS_R(la-1,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*(info->xp-info->xa)
	  +eriOS_R(la-1,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*(info->xw-info->xp)
	  +eriOS_R(la-2,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*(la-1)*0.5/info->zeta
	  -eriOS_R(la-2,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*(la-1)*0.5/info->zeta
	  +eriOS_R(la-1,ma,na,lb-1,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*lb*0.5/info->zeta
	  -eriOS_R(la-1,ma,na,lb-1,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*lb*0.5/info->zeta
	  +eriOS_R(la-1,ma,na,lb,mb,nb,lc-1,mc,nc,ld,md,nd,m+1,info)*lc*0.5/(info->zeta+info->eta)
	  +eriOS_R(la-1,ma,na,lb,mb,nb,lc,mc,nc,ld-1,md,nd,m+1,info)*ld*0.5/(info->zeta+info->eta);
	else if(ma>0) return
	   eriOS_R(la,ma-1,na,lb,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*(info->yp-info->ya)
	  +eriOS_R(la,ma-1,na,lb,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*(info->yw-info->yp)
	  +eriOS_R(la,ma-2,na,lb,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*(ma-1)*0.5/info->zeta
	  -eriOS_R(la,ma-2,na,lb,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*(ma-1)*0.5/info->zeta
	  +eriOS_R(la,ma-1,na,lb,mb-1,nb,lc,mc,nc,ld,md,nd,m  ,info)*mb*0.5/info->zeta
	  -eriOS_R(la,ma-1,na,lb,mb-1,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*mb*0.5/info->zeta
	  +eriOS_R(la,ma-1,na,lb,mb,nb,lc,mc-1,nc,ld,md,nd,m+1,info)*mc*0.5/(info->zeta+info->eta)
	  +eriOS_R(la,ma-1,na,lb,mb,nb,lc,mc,nc,ld,md-1,nd,m+1,info)*md*0.5/(info->zeta+info->eta);
	else if(na>0) return
	   eriOS_R(la,ma,na-1,lb,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*(info->zp-info->za)
	  +eriOS_R(la,ma,na-1,lb,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*(info->zw-info->zp)
	  +eriOS_R(la,ma,na-2,lb,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*(na-1)*0.5/info->zeta
	  -eriOS_R(la,ma,na-2,lb,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*(na-1)*0.5/info->zeta
	  +eriOS_R(la,ma,na-1,lb,mb,nb-1,lc,mc,nc,ld,md,nd,m  ,info)*nb*0.5/info->zeta
	  -eriOS_R(la,ma,na-1,lb,mb,nb-1,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*nb*0.5/info->zeta
	  +eriOS_R(la,ma,na-1,lb,mb,nb,lc,mc,nc-1,ld,md,nd,m+1,info)*nc*0.5/(info->zeta+info->eta)
	  +eriOS_R(la,ma,na-1,lb,mb,nb,lc,mc,nc,ld,md,nd-1,m+1,info)*nd*0.5/(info->zeta+info->eta);
	else if(lb>0) return
	   eriOS_R(la,ma,na,lb-1,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*(info->xp-info->xb)
	  +eriOS_R(la,ma,na,lb-1,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*(info->xw-info->xp)
	  +eriOS_R(la,ma,na,lb-2,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*(lb-1)*0.5/info->zeta
	  -eriOS_R(la,ma,na,lb-2,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*(lb-1)*0.5/info->zeta
	  +eriOS_R(la-1,ma,na,lb-1,mb,nb,lc,mc,nc,ld,md,nd,m  ,info)*la*0.5/info->zeta
	  -eriOS_R(la-1,ma,na,lb-1,mb,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*la*0.5/info->zeta
	  +eriOS_R(la,ma,na,lb-1,mb,nb,lc-1,mc,nc,ld,md,nd,m+1,info)*lc*0.5/(info->zeta+info->eta)
	  +eriOS_R(la,ma,na,lb-1,mb,nb,lc,mc,nc,ld-1,md,nd,m+1,info)*ld*0.5/(info->zeta+info->eta);
	else if(mb>0) return
	   eriOS_R(la,ma,na,lb,mb-1,nb,lc,mc,nc,ld,md,nd,m  ,info)*(info->yp-info->yb)
	  +eriOS_R(la,ma,na,lb,mb-1,nb,lc,mc,nc,ld,md,nd,m+1,info)*(info->yw-info->yp)
	  +eriOS_R(la,ma,na,lb,mb-2,nb,lc,mc,nc,ld,md,nd,m  ,info)*(mb-1)*0.5/info->zeta
	  -eriOS_R(la,ma,na,lb,mb-2,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*(mb-1)*0.5/info->zeta
	  +eriOS_R(la,ma-1,na,lb,mb-1,nb,lc,mc,nc,ld,md,nd,m  ,info)*ma*0.5/info->zeta
	  -eriOS_R(la,ma-1,na,lb,mb-1,nb,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*ma*0.5/info->zeta
	  +eriOS_R(la,ma,na,lb,mb-1,nb,lc,mc-1,nc,ld,md,nd,m+1,info)*mc*0.5/(info->zeta+info->eta)
	  +eriOS_R(la,ma,na,lb,mb-1,nb,lc,mc,nc,ld,md-1,nd,m+1,info)*md*0.5/(info->zeta+info->eta);
	else if(nb>0) return
	   eriOS_R(la,ma,na,lb,mb,nb-1,lc,mc,nc,ld,md,nd,m  ,info)*(info->zp-info->zb)
	  +eriOS_R(la,ma,na,lb,mb,nb-1,lc,mc,nc,ld,md,nd,m+1,info)*(info->zw-info->zp)
	  +eriOS_R(la,ma,na,lb,mb,nb-2,lc,mc,nc,ld,md,nd,m  ,info)*(nb-1)*0.5/info->zeta
	  -eriOS_R(la,ma,na,lb,mb,nb-2,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*(nb-1)*0.5/info->zeta
	  +eriOS_R(la,ma,na-1,lb,mb,nb-1,lc,mc,nc,ld,md,nd,m  ,info)*na*0.5/info->zeta
	  -eriOS_R(la,ma,na-1,lb,mb,nb-1,lc,mc,nc,ld,md,nd,m+1,info)*info->eta/(info->zeta+info->eta)*na*0.5/info->zeta
	  +eriOS_R(la,ma,na,lb,mb,nb-1,lc,mc,nc-1,ld,md,nd,m+1,info)*nc*0.5/(info->zeta+info->eta)
	  +eriOS_R(la,ma,na,lb,mb,nb-1,lc,mc,nc,ld,md,nd-1,m+1,info)*nd*0.5/(info->zeta+info->eta);
	else if(lc>0) return
	   eriOS_R(la,ma,na,lb,mb,nb,lc-1,mc,nc,ld,md,nd,m  ,info)*(info->xq-info->xc)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc-1,mc,nc,ld,md,nd,m+1,info)*(info->xw-info->xq)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc-2,mc,nc,ld,md,nd,m  ,info)*(lc-1)*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc-2,mc,nc,ld,md,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*(lc-1)*0.5/info->eta
	  +eriOS_R(la,ma,na,lb,mb,nb,lc-1,mc,nc,ld-1,md,nd,m  ,info)*ld*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc-1,mc,nc,ld-1,md,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*ld*0.5/info->eta
	  +eriOS_R(la-1,ma,na,lb,mb,nb,lc-1,mc,nc,ld,md,nd,m+1,info)*la*0.5/(info->eta+info->zeta)
	  +eriOS_R(la,ma,na,lb-1,mb,nb,lc-1,mc,nc,ld,md,nd,m+1,info)*lb*0.5/(info->eta+info->zeta);
	else if(mc>0) return
	   eriOS_R(la,ma,na,lb,mb,nb,lc,mc-1,nc,ld,md,nd,m  ,info)*(info->yq-info->yc)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc-1,nc,ld,md,nd,m+1,info)*(info->yw-info->yq)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc-2,nc,ld,md,nd,m  ,info)*(mc-1)*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc-2,nc,ld,md,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*(mc-1)*0.5/info->eta
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc-1,nc,ld,md-1,nd,m  ,info)*md*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc-1,nc,ld,md-1,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*md*0.5/info->eta
	  +eriOS_R(la,ma-1,na,lb,mb,nb,lc,mc-1,nc,ld,md,nd,m+1,info)*ma*0.5/(info->eta+info->zeta)
	  +eriOS_R(la,ma,na,lb,mb-1,nb,lc,mc-1,nc,ld,md,nd,m+1,info)*mb*0.5/(info->eta+info->zeta);
	else if(nc>0) return
	   eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc-1,ld,md,nd,m  ,info)*(info->zq-info->zc)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc-1,ld,md,nd,m+1,info)*(info->zw-info->zq)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc-2,ld,md,nd,m  ,info)*(nc-1)*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc-2,ld,md,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*(nc-1)*0.5/info->eta
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc-1,ld,md,nd-1,m  ,info)*nd*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc-1,ld,md,nd-1,m+1,info)*info->zeta/(info->eta+info->zeta)*nd*0.5/info->eta
	  +eriOS_R(la,ma,na-1,lb,mb,nb,lc,mc,nc-1,ld,md,nd,m+1,info)*na*0.5/(info->eta+info->zeta)
	  +eriOS_R(la,ma,na,lb,mb,nb-1,lc,mc,nc-1,ld,md,nd,m+1,info)*nb*0.5/(info->eta+info->zeta);
	else if(ld>0) return
	   eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld-1,md,nd,m  ,info)*(info->xq-info->xd)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld-1,md,nd,m+1,info)*(info->xw-info->xq)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld-2,md,nd,m  ,info)*(ld-1)*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld-2,md,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*(ld-1)*0.5/info->eta
	  +eriOS_R(la,ma,na,lb,mb,nb,lc-1,mc,nc,ld-1,md,nd,m  ,info)*lc*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc-1,mc,nc,ld-1,md,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*lc*0.5/info->eta
	  +eriOS_R(la-1,ma,na,lb,mb,nb,lc,mc,nc,ld-1,md,nd,m+1,info)*la*0.5/(info->eta+info->zeta)
	  +eriOS_R(la,ma,na,lb-1,mb,nb,lc,mc,nc,ld-1,md,nd,m+1,info)*lb*0.5/(info->eta+info->zeta);
	else if(md>0) return
	   eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md-1,nd,m  ,info)*(info->yq-info->yd)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md-1,nd,m+1,info)*(info->yw-info->yq)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md-2,nd,m  ,info)*(md-1)*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md-2,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*(md-1)*0.5/info->eta
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc-1,nc,ld,md-1,nd,m  ,info)*mc*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc-1,nc,ld,md-1,nd,m+1,info)*info->zeta/(info->eta+info->zeta)*mc*0.5/info->eta
	  +eriOS_R(la,ma-1,na,lb,mb,nb,lc,mc,nc,ld,md-1,nd,m+1,info)*ma*0.5/(info->eta+info->zeta)
	  +eriOS_R(la,ma,na,lb,mb-1,nb,lc,mc,nc,ld,md-1,nd,m+1,info)*mb*0.5/(info->eta+info->zeta);
	else if(nd>0) return
	   eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd-1,m  ,info)*(info->zq-info->zd)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd-1,m+1,info)*(info->zw-info->zq)
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd-2,m  ,info)*(nd-1)*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd-2,m+1,info)*info->zeta/(info->eta+info->zeta)*(nd-1)*0.5/info->eta
	  +eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc-1,ld,md,nd-1,m  ,info)*nc*0.5/info->eta
	  -eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc-1,ld,md,nd-1,m+1,info)*info->zeta/(info->eta+info->zeta)*nc*0.5/info->eta
	  +eriOS_R(la,ma,na-1,lb,mb,nb,lc,mc,nc,ld,md,nd-1,m+1,info)*na*0.5/(info->eta+info->zeta)
	  +eriOS_R(la,ma,na,lb,mb,nb-1,lc,mc,nc,ld,md,nd-1,m+1,info)*nb*0.5/(info->eta+info->zeta);
	else{
		printf("eriOS_R : Error unexpected event has occurred\n");
		exit(-1);
	}
}

double eriOS(double xa, double ya, double za, double norma,
             int la, int ma, int na, double alphaa,
             double xb, double yb, double zb, double normb,
             int lb, int mb, int nb, double alphab,
             double xc, double yc, double zc, double normc,
             int lc, int mc, int nc, double alphac,
             double xd, double yd, double zd, double normd,
             int ld, int md, int nd, double alphad){

	struct eriOSData_t info;
	int m;
	double T;

	info.alphaa=alphaa; info.alphab=alphab;
	info.alphac=alphac; info.alphad=alphad;

	info.zeta = alphaa+alphab;
	info.eta  = alphac+alphad;

	info.xa=xa; info.ya=ya; info.za=za;
	info.xb=xb; info.yb=yb; info.zb=zb;
	info.xc=xc; info.yc=yc; info.zc=zc;
	info.xd=xd; info.yd=yd; info.zd=zd;

	info.xp = (xa*alphaa+xb*alphab)/(alphaa+alphab);
	info.yp = (ya*alphaa+yb*alphab)/(alphaa+alphab);
	info.zp = (za*alphaa+zb*alphab)/(alphaa+alphab);

	info.xq = (xc*alphac+xd*alphad)/(alphac+alphad);
	info.yq = (yc*alphac+yd*alphad)/(alphac+alphad);
	info.zq = (zc*alphac+zd*alphad)/(alphac+alphad);

	fgamma_set(la+ma+na+lb+mb+nb+lc+mc+nc+ld+md+nd
	           ,info.zeta*info.eta/(info.zeta+info.eta)
	                              *( (info.xp-info.xq)*(info.xp-info.xq)
	                                +(info.yp-info.yq)*(info.yp-info.yq)
	                                +(info.zp-info.zq)*(info.zp-info.zq))
	           ,info.ssss);

#define TWO_PI_POW2_5 34.986836655249725
	T = TWO_PI_POW2_5 / (info.zeta*info.eta*sqrt(info.zeta+info.eta))
	    * exp( -info.alphaa*info.alphab/info.zeta
	                       *( (info.xa-info.xb)*(info.xa-info.xb)
	                         +(info.ya-info.yb)*(info.ya-info.yb)
	                         +(info.za-info.zb)*(info.za-info.zb))
	           -info.alphac*info.alphad/info.eta
	                       *( (info.xc-info.xd)*(info.xc-info.xd)
	                         +(info.yc-info.yd)*(info.yc-info.yd)
	                         +(info.zc-info.zd)*(info.zc-info.zd))  )
	    * norma * normb * normc * normd;
#undef  TWO_PI_POW2_5

	for(m=(la+ma+na+lb+mb+nb+lc+mc+nc+ld+md+nd); m>=0; m--)
		info.ssss[m] = T*info.ssss[m];

	info.xw = (info.xp*info.zeta+info.xq*info.eta)/(info.zeta+info.eta);
	info.yw = (info.yp*info.zeta+info.yq*info.eta)/(info.zeta+info.eta);
	info.zw = (info.zp*info.zeta+info.zq*info.eta)/(info.zeta+info.eta);

	return eriOS_R(la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd,0,&info);
}
#undef MAXL

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
int main(int argc, char *argv[]){
	double max_err     = 0.0;
	double max_percent = 0.0;
	unsigned long long n;
	double err;
	double percent;
	struct gto_t a,b,c,d;     // using 4 gto

	printf("----------------------------------------------------\n"
	       "<================  Accuracy  Testing  =============>\n"
	       "                                                    \n");

	// overlap integrals testing
	printf("Testing overlap integral (a|b)\n"
	       "1,000,000 integrals generated randomly.\n");
	max_err     = 0.0;
	max_percent = 0.0;
	for(n=0; n < 1000000; n++){
		random_gto(&a);
		random_gto(&b);

		err = fabs(
		    overlap(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z) -
			overlap_py(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z)
		          );
		percent = err*100/overlap(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z);
		percent = fabs(percent);

		if(err > max_err){
			max_err = err;
			max_percent = percent;
			printf("   max error is %20.8E  (%18.6E Percent )\n",
			                         max_err,max_percent);
		}
	}
	printf("Done testing: max error is %18.6E (%18.6E Percent )\n"
	       "-------------------\n\n",
	                                 max_err,max_percent);


	// overlap integrals MD Scheme testing
	printf("Testing overlap integral (a|b) using McMurchie-Davidson Scheme\n"
	       "1,000,000 integrals generated randomly.\n");
	max_err     = 0.0;
	max_percent = 0.0;
	for(n=0; n < 1000000; n++){
		random_gto(&a);
		random_gto(&b);

		err = fabs(
		    overlap_MD(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z) -
			overlap_py(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z)
		          );
		percent = err*100/overlap(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z);
		percent = fabs(percent);

		if(err > max_err){
			max_err = err;
			max_percent = percent;
			printf("   max error is %20.8E  (%18.6E Percent )\n",
			                         max_err,max_percent);
		}
	}
	printf("Done testing: max error is %18.6E (%18.6E Percent )\n"
	       "-------------------\n\n",
	                                 max_err,max_percent);


	// kinetic integrals testing
	printf("Testing kinetic integral (a|-1/2 Laplacian|b)\n"
	       "1,000,000 integrals generated randomly.\n");
	max_err     = 0.0;
	max_percent = 0.0;
	for(n=0; n < 1000000; n++){
		random_gto(&a);
		random_gto(&b);

		err = fabs(
		    kinetic(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z) -
			kinetic_py(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z)
		          );
		percent = err*100/kinetic(
		              a.exp,a.l,a.m,a.n,a.x,a.y,a.z,
		              b.exp,b.l,b.m,b.n,b.x,b.y,b.z);
		percent = fabs(percent);

		if(err > max_err){
			max_err = err;
			max_percent = percent;
			printf("   max error is %20.6E  (%18.6E Percent )\n",
			                         max_err,max_percent);
		}
	}
	printf("Done testing: max error is %18.6E (%18.6E Percent )\n"
	       "-------------------\n\n",
	                                 max_err,max_percent);
	
	// nuclear attraction integrals testing
	printf("Testing nuclear attraction integral (a|-1/Rc|b)\n"
	       "1,000,000 integrals generated randomly.\n");
	max_err     = 0.0;
	max_percent = 0.0;
	for(n=0; n < 1000000; n++){
		random_gto(&a);
		random_gto(&b);
		random_gto(&c);

		err = fabs(
		    nuclear_attraction(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z) -
		    nai(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z)
		          );
		percent = err*100/
		    nuclear_attraction(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z);

		percent = fabs(percent);

		if(err > max_err){
			max_err = err;
			max_percent = percent;
			printf("   max error is %18.6E  (%18.6E Percent ) "
			       "@ (lmn)-(lmn) = (%d%d%d)-(%d%d%d)\n",
			       max_err,max_percent,a.l,a.m,a.n,b.l,b.m,b.n);
		}
	}
	printf("Done testing: max error is %18.6E (%18.6E Percent )\n"
	       "-------------------\n\n",
	                                 max_err,max_percent);
    
	// electron repulsion integrals testing		
	printf("Testing electron repulsion integral (ab|1/r12|cd)\n"
	       "1,000,000 integrals generated randomly.\n");
	fflush(stdout);
	max_err     = 0.0;
	max_percent = 0.0;
	for(n=0; n < 1000000; n++){
		random_gto(&a);
		random_gto(&b);
		random_gto(&c);
		random_gto(&d);

		err = fabs(
		    coulomb_repulsion(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,d.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,c.m,d.n,d.exp)
		    -
		              eri(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,d.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,c.m,d.n,d.exp)
		          );
		percent = err*100/
		    coulomb_repulsion(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,d.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,c.m,d.n,d.exp);

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
	fflush(stdout);
/*
MDScheme:

	// electron repulsion integrals testing	using McMurchie-Davidson Scheme	
	printf("Testing electron repulsion integral (ab|1/r12|cd) using McMurchie-Davidson\n"
	       "1,000,000 integrals generated randomly.\n");
	fflush(stdout);
	max_err     = 0.0;
	max_percent = 0.0;
	for(n=0; n < 1000000; n++){
		random_gto(&a);
		random_gto(&b);
		random_gto(&c);
		random_gto(&d);

		err = fabs(
		    //coulomb_repulsion(
		             eri(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,d.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,c.m,d.n,d.exp)
		    -
		              eri_MD(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,d.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,c.m,d.n,d.exp)
		          );
		percent = err*100/
		    //coulomb_repulsion(
		            eri(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,d.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,c.m,d.n,d.exp);

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
	fflush(stdout);
*/

	// electron repulsion integrals testing	using Obara-Saika Scheme	
	printf("Testing electron repulsion integral (ab|1/r12|cd) using Obara-Saika recursive\n"
	       "1,000,000 integrals generated randomly.\n");
	fflush(stdout);
	max_err     = 0.0;
	max_percent = 0.0;
	for(n=0; n < 1000000; n++){
		random_gto(&a);
		random_gto(&b);
		random_gto(&c);
		random_gto(&d);

		err = fabs(
		    //coulomb_repulsion(
		         eri(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,c.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,d.m,d.n,d.exp)
		    -
		              eriOS(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,b.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,c.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,d.m,d.n,d.exp)
		          );
		percent = err*100/
		    //coulomb_repulsion(
		         eri(
		              a.x,a.y,a.z,1.0,a.l,a.m,a.n,a.exp,
		              b.x,b.y,b.z,1.0,b.l,b.m,a.n,b.exp,
		              c.x,c.y,c.z,1.0,c.l,c.m,c.n,c.exp,
			          d.x,d.y,d.z,1.0,d.l,d.m,d.n,d.exp);

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
	fflush(stdout);

	printf("----------------------------------------------------\n"
	       "<================ Performance Testing =============>\n"
	       "                                                    \n");
	printf("Computing 1,000,000 (SS|SS)\n");
	for(n=0; n < 1000000; n++){
		eriOS(
		    a.x,a.y,a.z,1.0,0,0,0,a.exp,
		    b.x,b.y,b.z,1.0,0,0,0,b.exp,
		    c.x,c.y,c.z,1.0,0,0,0,c.exp,
		    d.x,d.y,d.z,1.0,0,0,0,d.exp);
	}
	printf("Computing 1,000,000 (PP|PP)\n");
	for(n=0; n < 1000000; n++){
		eriOS(
		    a.x,a.y,a.z,1.0,1,0,0,a.exp,
		    b.x,b.y,b.z,1.0,1,0,0,b.exp,
		    c.x,c.y,c.z,1.0,1,0,0,c.exp,
		    d.x,d.y,d.z,1.0,1,0,0,d.exp);
	}

	return 0;
}
