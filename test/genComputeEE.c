#include <stdio.h>
#include <stdlib.h>

#define MAXBASIS    16
#define MAXCONTRACT 16
struct GTOShell_t{
	double x,y,z;                             // shell center
	int nBasis;                               // number of basis in the shell
	int nContract;                            // number contract functions
	int maxL;                                 // maximum angular index
	int min,max;                              // basis index
	int l[MAXBASIS],m[MAXBASIS],n[MAXBASIS];  // angular index
	double exps[MAXCONTRACT];                 // exponents
	double coef[MAXCONTRACT][MAXBASIS];       // coefficients
	int type;                                 // predefined type 
#define TYPES 0
#define TYPEL 1
#define TYPEP 2
#define TYPED 3
#define TYPEF 4
};


#define MAXL          4
#define MAXSHELLINT  50625
int iBx[MAXSHELLINT];   // index for Bx storage
int iBy[MAXSHELLINT];   // index for By storage
int iBz[MAXSHELLINT];   // index for Bz storage
int  mX[MAXSHELLINT];   // maximum in x direction
int  mY[MAXSHELLINT];   // maximum in y direction
int  mZ[MAXSHELLINT];   // maximum in z direction
int  mT[MAXSHELLINT];   // type of mX,mY,mZ for loop optimization

double Bx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double By[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double Bz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double Sx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double Sy[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double Sz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double  F[4*(MAXL+1)];         // Boys functions


// genSetB : generate the "unrolled" version of genSetBxyzSxyz
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void genSetB(
	char *str,                 // name string
	struct GTOShell_t *P,      // pointer to P shell
	struct GTOShell_t *Q,      // pointer to Q shell
	struct GTOShell_t *I,      // pointer to I shell
	struct GTOShell_t *J){     // pointer to J shell

	int maxa,maxb,maxc,maxd;
	maxa = P->maxL;
	maxb = Q->maxL;
	maxc = I->maxL;
	maxd = J->maxL;

	// header
	printf(
	"double  genSet%sBxyzSxyzF(\n"
	"    double xa,  double ya,  double za, double alphaa,\n"
	"    double xb,  double yb,  double zb, double alphab,\n"
	"    double xc,  double yc,  double zc, double alphac,\n"
	"    double xd,  double yd,  double zd, double alphad,\n"
	"    double *Bx, double *By, double *Bz,\n"
	"    double *Sx, double *Sy, double *Sz,\n"
	"    double *F){\n"

	"    double rab2,rcd2,rpq2,xp,yp,zp,xq,yq,zq;   // auxilary\n"
	"    double eta1, eta2, zeta;                   // auxilary distance\n"
	"    double t;                                  // argument of fgamma\n"
	"    register double sum;\n"

	"    eta1 = 1.0/(alphaa+alphab);\n"
	"    eta2 = 1.0/(alphac+alphad);\n"
	"    zeta = 4.0*(alphaa+alphab)*(alphac+alphad)\n"
	"              /(alphaa+alphab+alphac+alphad);\n"

	"    xp   = (alphaa*xa+alphab*xb)*eta1;\n"
	"    yp   = (alphaa*ya+alphab*yb)*eta1;\n"
	"    zp   = (alphaa*za+alphab*zb)*eta1;\n"
	"    xq   = (alphac*xc+alphad*xd)*eta2;\n"
	"    yq   = (alphac*yc+alphad*yd)*eta2;\n"
	"    zq   = (alphac*zc+alphad*zd)*eta2;\n"

	"#define DIST2(x1,y1,z1,x2,y2,z2) ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))\n"

	"    rab2 = DIST2(xa,ya,za,xb,yb,zb);\n"
	"    rcd2 = DIST2(xc,yc,zc,xd,yd,zd);\n"
	"    rpq2 = DIST2(xp,yp,zp,xq,yq,zq);\n"

	"    t     = 0.25*rpq2*zeta;\n"
	"    fgamma_set(%d, t, F);\n"

	"#define TWO_PI_POW2_5 34.986836655249725\n"
	"    t = TWO_PI_POW2_5\n"
	"        *eta1*eta2\n"
	"        /sqrt(alphaa+alphab+alphac+alphad)\n"
	"        *exp(-alphaa*alphab*rab2*eta1\n"
	"             -alphac*alphad*rcd2*eta2);\n"
	"#undef  TWO_PI_POW2_5\n"
	"    if(t < PRIMITIVE_CUTOFF) return 0.0;\n"

	,str,maxa+maxb+maxc+maxd);

	int a,b,c,d;                               // angular index
	int i,j,k,n;                               // generic index
	int na,nb,nc;                              // index locator

	nc = (maxd+1); nb = nc*(maxc+1); na = nb*(maxb+1);
	for(a=0; a < maxa+1; a++)
	for(b=0; b < maxb+1; b++)
		if(a!=0  && b!=maxb){
			for(c=0; c < maxc+1; c++)
			for(d=0; d < maxd+1; d++){
				i =4*(MAXL+1)*(     a*na +     b*nb + c*nc + d);
				j =4*(MAXL+1)*( (a-1)*na + (b+1)*nb + c*nc + d);
				k =4*(MAXL+1)*( (a-1)*na +     b*nb + c*nc + d);
				n = a+b+c+d;

				printf(
				"    Bx[%d] = Bx[%d];" "By[%d] = By[%d];" "Bz[%d] = Bz[%d];\n"
				,i+n,j+n,  i+n,j+n,  i+n,j+n
				);

				                              printf("    sum=xa-xb;\n");
				for(n=0; n<=(a+b+c+d-1); n++) printf("    Bx[%d] = Bx[%d] - sum*Bx[%d];\n",i+n,j+n,k+n);

				                              printf("    sum=ya-yb;\n");
				for(n=0; n<=(a+b+c+d-1); n++) printf("    By[%d] = By[%d] - sum*By[%d];\n",i+n,j+n,k+n);

				                              printf("    sum=za-zb;\n");
				for(n=0; n<=(a+b+c+d-1); n++) printf("    Bz[%d] = Bz[%d] - sum*Bz[%d];\n",i+n,j+n,k+n);

				n = a+b+c+d;
				for(k=0; k<=maxa+maxb+maxc+maxd-n; k++){

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=Bx[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sx[%d]=sum;\n",i+k);

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=By[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sy[%d]=sum;\n",i+k);

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=Bz[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sz[%d]=sum;\n",i+k);

				}

			}
		}else

		for(c=0; c < maxc+1; c++)
		for(d=0; d < maxd+1; d++)
			if(c!=0  && d!=maxd){
				i =4*(MAXL+1)*( a*na + b*nb +     c*nc +  d   );
				j =4*(MAXL+1)*( a*na + b*nb + (c-1)*nc +  d+1 );
				k =4*(MAXL+1)*( a*na + b*nb + (c-1)*nc +  d   );
				n = a+b+c+d;

				printf(
				"    Bx[%d] = Bx[%d];" "By[%d] = By[%d];" "Bz[%d] = Bz[%d];\n"
				,i+n,j+n,  i+n,j+n,  i+n,j+n
				);

				                              printf("    sum=xc-xd;\n");
				for(n=0; n<=(a+b+c+d-1); n++) printf("    Bx[%d] = Bx[%d] - sum*Bx[%d];\n",i+n,j+n,k+n);

				                              printf("    sum=yc-yd;\n");
				for(n=0; n<=(a+b+c+d-1); n++) printf("    By[%d] = By[%d] - sum*By[%d];\n",i+n,j+n,k+n);

				                              printf("    sum=zc-zd;\n");
				for(n=0; n<=(a+b+c+d-1); n++) printf("    Bz[%d] = Bz[%d] - sum*Bz[%d];\n",i+n,j+n,k+n);

				n = a+b+c+d;
				for(k=0; k<=maxa+maxb+maxc+maxd-n; k++){

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=Bx[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sx[%d]=sum;\n",i+k);

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=By[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sy[%d]=sum;\n",i+k);

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=Bz[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sz[%d]=sum;\n",i+k);

				}

			}else{
				i =4*(MAXL+1)*(a*na + b*nb + c*nc + d);

				printf(
				"    eriB_1d(Bx+%d, %d, %d, %d, %d, xp,xa,xb,xq,xc,xd,eta1,eta2,zeta);\n"
				"    eriB_1d(By+%d, %d, %d, %d, %d, yp,ya,yb,yq,yc,yd,eta1,eta2,zeta);\n"
				"    eriB_1d(Bz+%d, %d, %d, %d, %d, zp,za,zb,zq,zc,zd,eta1,eta2,zeta);\n"
				,i,a,b,c,d,  i,a,b,c,d,  i,a,b,c,d);

				n = a+b+c+d;
				for(k=0; k<=maxa+maxb+maxc+maxd-n; k++){

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=Bx[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sx[%d]=sum;\n",i+k);

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=By[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sy[%d]=sum;\n",i+k);

					                  printf("    sum=0.0;\n");
					for(j=0;j<=n;j++) printf("    sum+=Bz[%d]*F[%d];\n",i+j,k+j);
					                  printf("    Sz[%d]=sum;\n",i+k);

				}

			}

	printf(
	"    return t;\n"
	"}\n"
	);
}


// gen : generate the "unrolled" version the computeQuartetEE code
//
// Feb 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void gen(
	char *str,                // string name
	struct GTOShell_t *P,     // pointer to P shell
	struct GTOShell_t *Q,     // pointer to Q shell
	struct GTOShell_t *I,     // pointer to I shell
	struct GTOShell_t *J){    // pointer to J shell

	int kX, kY, kZ;                // summation loop
	int p,q,i,j;                   // basis function index
	int nEE;                       // integral index

	// index preparation
	nEE = 0;
	for(p=0; p < P->nBasis; p++)
	for(q=0; q < Q->nBasis; q++)
	for(i=0; i < I->nBasis; i++)
	for(j=0; j < J->nBasis; j++){

		// precompute index for efficientcy
		kZ =    (J->maxL+1);
		kY = kZ*(I->maxL+1);
		kX = kY*(Q->maxL+1);
		iBx[nEE] = 4*(MAXL+1)*(P->l[p]*kX + Q->l[q]*kY + I->l[i]*kZ + J->l[j]);
		iBy[nEE] = 4*(MAXL+1)*(P->m[p]*kX + Q->m[q]*kY + I->m[i]*kZ + J->m[j]);
		iBz[nEE] = 4*(MAXL+1)*(P->n[p]*kX + Q->n[q]*kY + I->n[i]*kZ + J->n[j]);
#define MXYZZERO -4
#define MXYZERO -3
#define MXZZERO -2
#define MYZZERO -1
#define MXZERO 0
#define MYZERO 1
#define MZZERO 2
#define MXMAX  3
#define MYMAX  4
#define MZMAX  5		
		mX[nEE]  = P->l[p] + Q->l[q] + I->l[i] + J->l[j];
		mY[nEE]  = P->m[p] + Q->m[q] + I->m[i] + J->m[j];
		mZ[nEE]  = P->n[p] + Q->n[q] + I->n[i] + J->n[j];
		if(mX[nEE]==0 && mY[nEE]==0 && mZ[nEE]==0)      mT[nEE] = MXYZZERO;
		else if(mX[nEE]==0 && mY[nEE]==0)               mT[nEE] = MXYZERO;
		else if(mX[nEE]==0 && mZ[nEE]==0)               mT[nEE] = MXZZERO;
		else if(mY[nEE]==0 && mZ[nEE]==0)               mT[nEE] = MYZZERO;
		else if(mX[nEE]==0)                             mT[nEE] = MXZERO;
		else if(mY[nEE]==0)                             mT[nEE] = MYZERO;
		else if(mZ[nEE]==0)                             mT[nEE] = MZZERO;
		else if(mX[nEE] > mY[nEE] && mX[nEE] > mZ[nEE]) mT[nEE] = MXMAX;
		else if(mY[nEE] > mZ[nEE])                      mT[nEE] = MYMAX;
		else                                            mT[nEE] = MZMAX;

		// increase number of EE index
		nEE++;
	}

	printf("void compute%sQuartetEE(\n"
	       "const struct GTOShell_t *P,\n"
	       "const struct GTOShell_t *Q,\n"
	       "const struct GTOShell_t *I,\n"
	       "const struct GTOShell_t *J,\n"
	       "double *EEStore){\n"
	,str);
	printf("int pC,qC,iC,jC;\n");
	printf("double r,rpqi,EE,tSum;\n");
	printf("memset(EEStore,0,sizeof(double)*%d);\n"
	,P->nBasis * Q->nBasis * I->nBasis * J->nBasis);

	printf("for(pC=0; pC < P->nContract; pC++)                           \n"
	       "for(qC=0; qC < Q->nContract; qC++)                           \n"
	       "for(iC=0; iC < I->nContract; iC++)                           \n"
	       "for(jC=0; jC < J->nContract; jC++){                          \n"
	       "    r=genSet%sBxyzSxyzF(                                     \n"
	       "                      P->x,P->y,P->z,P->exps[pC],            \n"
	       "                      Q->x,Q->y,Q->z,Q->exps[qC],            \n"
	       "                      I->x,I->y,I->z,I->exps[iC],            \n"
	       "                      J->x,J->y,J->z,J->exps[jC],            \n"
	       "                      Bx,By,Bz,Sx,Sy,Sz,F);                  \n"
	       "    if(r < PRIMITIVE_CUTOFF) continue;                       \n"
	,str);

	nEE = 0;
	for(p=0; p < P->nBasis; p++)
	for(q=0; q < Q->nBasis; q++)
	for(i=0; i < I->nBasis; i++){
		printf("    rpqi = r * P->coef[pC][%d] * Q->coef[qC][%d] * I->coef[iC][%d];\n",p,q,i);
	for(j=0; j < J->nBasis; j++){

		switch(mT[nEE]){

		case MXYZZERO:
			printf("    EE = F[0];\n");
		break;

		case MXYZERO:
			printf("    EE = Sz[%d];\n",iBz[nEE]);
		break;

		case MXZZERO:
			printf("    EE = Sy[%d];\n",iBy[nEE]);
		break;

		case MYZZERO:
			printf("    EE = Sx[%d];\n",iBx[nEE]);
		break;

		case MXZERO:
			printf("    EE = 0.0;\n");
			for(kY=mY[nEE];kY>=0;kY--)
				printf("    EE += By[%d]*Sz[%d];\n",iBy[nEE]+kY,iBz[nEE]+kY);
		break;

		case MYZERO:
			printf("    EE = 0.0;\n");
			for(kZ=mZ[nEE];kZ>=0;kZ--)
				printf("    EE += Bz[%d]*Sx[%d];\n",iBz[nEE]+kZ,iBx[nEE]+kZ);
		break;

		case MZZERO:
			printf("    EE = 0.0;\n");
			for(kX=mX[nEE];kX>=0;kX--)
				printf("    EE += Bx[%d]*Sy[%d];\n",iBx[nEE]+kX,iBy[nEE]+kX);
		break;

		case MXMAX: 
			printf("    EE = 0.0;\n");
			for(kY=mY[nEE];kY>=0;kY--){
				printf("    tSum = 0.0;\n");
				for(kZ=mZ[nEE];kZ>=0;kZ--)
					printf("    tSum+=Bz[%d]*Sx[%d];\n",iBz[nEE]+kZ,iBx[nEE]+kY+kZ);
				printf("    EE += By[%d]*tSum;\n",iBy[nEE]+kY);
			}
		break;

		case MYMAX:
			printf("    EE = 0.0;\n");
			for(kZ=mZ[nEE];kZ>=0;kZ--){
				printf("    tSum = 0.0;\n");
				for(kX=mX[nEE];kX>=0;kX--)
					printf("    tSum+=Bx[%d]*Sy[%d];\n",iBx[nEE]+kX,iBy[nEE]+kX+kZ);
				printf("    EE += Bz[%d]*tSum;\n",iBz[nEE]+kZ);
			}	
		break;

		case MZMAX:
			printf("    EE = 0.0;\n");
			for(kX=mX[nEE];kX>=0;kX--){
				printf("    tSum = 0.0;\n");
				for(kY=mY[nEE];kY>=0;kY--)
					printf("    tSum+=By[%d]*Sz[%d];\n",iBy[nEE]+kY,iBz[nEE]+kX+kY);
				printf("    EE += Bx[%d]*tSum;\n",iBx[nEE]+kX);
			}
		break;
		}

		printf("    EEStore[%d] += EE * rpqi * J->coef[jC][%d];\n",
		nEE,j);

		nEE++;
	}
	}

	printf("}                                                            \n");
	printf("}                                                            \n");

}


// genStore: generate the "unrolled" version of storeQuartetEE code
//
// Feb 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
// Feb 24, 2013 - Teepanis Chachiyo
//     Change the name to UstoreQuartetEE
//
void genUStore(
	char *str,                   // name string
	struct GTOShell_t *P,        // pointer to P shell
	struct GTOShell_t *Q,        // pointer to Q shell
	struct GTOShell_t *I,        // pointer to I shell
	struct GTOShell_t *J){       // pointer to J shell

	int p,q,i,j;
	int pp,qq,ii;

	ii = J->nBasis;
	qq = I->nBasis * J->nBasis;
	pp = Q->nBasis * I->nBasis * J->nBasis;

	printf(
	"void Ustore%sQuartetEE(\n"
	"    const struct GTOShell_t *P,\n"
	"    const struct GTOShell_t *Q,\n"
	"    const struct GTOShell_t *I,\n"
	"    const struct GTOShell_t *J,\n"
	"    int nBasis,\n"
	"    double *GA, double *GB,\n"
	"    const double *PT, const double *PA, const double *PB,\n"
	"    const double *EEStore){\n"

	"    const double *PTptr, *PAptr, *PBptr;\n"
	"    register double Ga,Gb;\n"

	"    double *GAptr, *GBptr;\n"
	,str);


/*
//
// First generation
//
#define STORE_COULOMB(P,Q,I,J,p,q,i,j,nEE)                                     \
	printf("    PTptr = PT + " #I "->min * nBasis + " #J "->min;\n");          \
	for(i=0; i < I->nBasis; i++){                                              \
		for(j=0; j < J->nBasis; j++)                                           \
			printf("    P1[%d] = PTptr[%d];\n",i*J->nBasis+j,j);               \
		printf("    PTptr+=nBasis;\n");                                        \
	}                                                                          \
	printf("    memset(G1,0,%d*sizeof(double));\n",P->nBasis*Q->nBasis);       \
	printf("    memset(G2,0,%d*sizeof(double));\n",P->nBasis*Q->nBasis);       \
                                                                               \
	for(p=0; p < P->nBasis; p++)                                               \
	for(q=0; q < Q->nBasis; q++){                                              \
		printf("    Ga=0.0; Gb=0.0;\n");                                       \
		for(i=0; i < I->nBasis; i++)                                           \
		for(j=0; j < J->nBasis; j++){                                          \
			printf("    EE  = EEStore[%d]; EE += EE;\n",nEE);                  \
			printf("    Ga += P1[%d]*EE;\n",i*J->nBasis+j);                    \
			printf("    Gb += P1[%d]*EE;\n",i*J->nBasis+j);                    \
		}                                                                      \
		printf("    G1[%d] += Ga;\n",p*Q->nBasis+q);                           \
		printf("    G2[%d] += Gb;\n",p*Q->nBasis+q);                           \
	}                                                                          \
                                                                               \
	for(p=0; p < P->nBasis; p++){                                              \
		printf("    GAptr = GA + (%d+" #P "->min)*nBasis + " #Q "->min;\n",p); \
		printf("    GBptr = GB + (%d+" #P "->min)*nBasis + " #Q "->min;\n",p); \
		for(q=0; q < Q->nBasis; q++){                                          \
			printf("    GAptr[%d] += G1[%d];\n",q,p*Q->nBasis+q);              \
			printf("    GBptr[%d] += G2[%d];\n",q,p*Q->nBasis+q);              \
		}                                                                      \
	}

#define STORE_EXCHANGE(P,Q,I,J,p,q,i,j,nEE)                                    \
	printf("    PAptr = PA + " #Q "->min * nBasis + " #J "->min;\n");          \
	printf("    PBptr = PB + " #Q "->min * nBasis + " #J "->min;\n");          \
	for(q=0; q < Q->nBasis; q++){                                              \
		for(j=0; j < J->nBasis; j++){                                          \
			printf("    P1[%d] = PAptr[%d];\n",q*J->nBasis+j,j);               \
			printf("    P2[%d] = PBptr[%d];\n",q*J->nBasis+j,j);               \
		}                                                                      \
		printf("    PAptr+=nBasis;\n");                                        \
		printf("    PBptr+=nBasis;\n");                                        \
	}                                                                          \
	printf("    memset(G1,0,%d*sizeof(double));\n",P->nBasis*I->nBasis);       \
	printf("    memset(G2,0,%d*sizeof(double));\n",P->nBasis*I->nBasis);       \
                                                                               \
	for(p=0; p < P->nBasis; p++)                                               \
	for(i=0; i < I->nBasis; i++){                                              \
		printf("    Ga=0.0; Gb=0.0;\n");                                       \
		for(q=0; q < Q->nBasis; q++)                                           \
		for(j=0; j < J->nBasis; j++){                                          \
			printf("    EE  = EEStore[%d];\n",nEE);                            \
			printf("    Ga += P1[%d]*EE;\n",q*J->nBasis+j);                    \
			printf("    Gb += P2[%d]*EE;\n",q*J->nBasis+j);                    \
		}                                                                      \
		printf("    G1[%d] -= Ga;\n",p*I->nBasis+i);                           \
		printf("    G2[%d] -= Gb;\n",p*I->nBasis+i);                           \
	}                                                                          \
                                                                               \
	for(p=0; p < P->nBasis; p++){                                              \
		printf("    GAptr = GA + (%d+" #P "->min)*nBasis + " #I "->min;\n",p); \
		printf("    GBptr = GB + (%d+" #P "->min)*nBasis + " #I "->min;\n",p); \
		for(i=0; i < I->nBasis; i++){                                          \
			printf("    GAptr[%d] += G1[%d];\n",i,p*I->nBasis+i);              \
			printf("    GBptr[%d] += G2[%d];\n",i,p*I->nBasis+i);              \
		}                                                                      \
	}                                                                          \
	for(i=0; i < I->nBasis; i++){                                              \
		printf("    GAptr = GA + (%d+" #I "->min)*nBasis + " #P "->min;\n",i); \
		printf("    GBptr = GB + (%d+" #I "->min)*nBasis + " #P "->min;\n",i); \
		for(p=0; p < P->nBasis; p++){                                          \
			printf("    GAptr[%d] += G1[%d];\n",p,p*I->nBasis+i);              \
			printf("    GBptr[%d] += G2[%d];\n",p,p*I->nBasis+i);              \
		}                                                                      \
	}

*/

#define STORE_COULOMB(P,Q,I,J,p,q,i,j,nEE)                                     \
	printf("    PTptr = PT + " #I "->min * nBasis + " #J "->min;\n");          \
	for(i=0; i < I->nBasis; i++){                                              \
		for(j=0; j < J->nBasis; j++)                                           \
			printf("    P12[%d] = PTptr[%d];\n",i*J->nBasis+j,j);              \
		printf("    PTptr+=nBasis;\n");                                        \
	}                                                                          \
                                                                               \
	for(p=0; p < P->nBasis; p++){                                              \
		for(q=0; q < Q->nBasis; q++){                                          \
			printf("    Ga=0.0;\n");                                           \
			for(i=0; i < I->nBasis; i++)                                       \
			for(j=0; j < J->nBasis; j++){                                      \
				printf("    Ga += P12[%d]*EEStore[%d];\n",i*J->nBasis+j,nEE);  \
			}                                                                  \
			printf("    G12[%d]  = Ga+Ga;\n",p*Q->nBasis+q);                   \
		}                                                                      \
    }                                                                          \
                                                                               \
	printf("    if(" #P "->min >= " #Q "->min){\n");                           \
	printf("        GAptr = GA + " #P "->min*nBasis + " #Q "->min;\n");        \
	printf("        GBptr = GB + " #P "->min*nBasis + " #Q "->min;\n");        \
	for(p=0; p < P->nBasis; p++){                                              \
		for(q=0; q < Q->nBasis; q++){                                          \
			printf("        GAptr[%d] += G12[%d];",  q,p*Q->nBasis+q);         \
			printf(       " GBptr[%d] += G12[%d];\n",q,p*Q->nBasis+q);         \
		}                                                                      \
		printf("        GAptr+=nBasis; GBptr+=nBasis;\n");                     \
	}                                                                          \
	printf("    }\n");                                                         \
	printf("    if(" #Q "->min >= " #P "->min){\n");                           \
	printf("        GAptr = GA + " #Q "->min*nBasis + " #P "->min;\n");        \
	printf("        GBptr = GB + " #Q "->min*nBasis + " #P "->min;\n");        \
	for(q=0; q < Q->nBasis; q++){                                              \
		for(p=0; p < P->nBasis; p++){                                          \
			printf("        GAptr[%d] += G12[%d];",  p,p*Q->nBasis+q);         \
			printf(       " GBptr[%d] += G12[%d];\n",p,p*Q->nBasis+q);         \
		}                                                                      \
		printf("        GAptr+=nBasis; GBptr+=nBasis;\n");                     \
	}                                                                          \
	printf("    }\n");


#define STORE_EXCHANGE(P,Q,I,J,p,q,i,j,nEE)                                    \
	printf("    PAptr = PA + " #Q "->min * nBasis + " #J "->min;\n");          \
	printf("    PBptr = PB + " #Q "->min * nBasis + " #J "->min;\n");          \
	for(q=0; q < Q->nBasis; q++){                                              \
		for(j=0; j < J->nBasis; j++){                                          \
			printf("    P12[%d] = PAptr[%d];"  ,2*(q*J->nBasis+j)+0,j);        \
			printf(   " P12[%d] = PBptr[%d];\n",2*(q*J->nBasis+j)+1,j);        \
		}                                                                      \
		printf("    PAptr+=nBasis; PBptr+=nBasis;\n");                         \
	}                                                                          \
                                                                               \
	for(p=0; p < P->nBasis; p++)                                               \
	for(i=0; i < I->nBasis; i++){                                              \
		printf("    Ga=0.0; Gb=0.0;\n");                                       \
		for(q=0; q < Q->nBasis; q++)                                           \
		for(j=0; j < J->nBasis; j++){                                          \
			printf("    Ga += P12[%d]*EEStore[%d];",  2*(q*J->nBasis+j)+0,nEE);\
			printf(   " Gb += P12[%d]*EEStore[%d];\n",2*(q*J->nBasis+j)+1,nEE);\
		}                                                                      \
		printf("    G12[%d]  = Ga;",  2*(p*I->nBasis+i)+0);                    \
		printf(   " G12[%d]  = Gb;\n",2*(p*I->nBasis+i)+1);                    \
	}                                                                          \
                                                                               \
	printf("    if(" #P "->min >= " #I "->min){\n");                           \
	printf("        GAptr = GA + " #P "->min*nBasis + " #I "->min;\n");        \
	printf("        GBptr = GB + " #P "->min*nBasis + " #I "->min;\n");        \
	for(p=0; p < P->nBasis; p++){                                              \
		for(i=0; i < I->nBasis; i++){                                          \
			printf("        GAptr[%d] -= G12[%d];",  i,2*(p*I->nBasis+i)+0);   \
			printf(       " GBptr[%d] -= G12[%d];\n",i,2*(p*I->nBasis+i)+1);   \
		}                                                                      \
		printf("        GAptr+=nBasis; GBptr+=nBasis;\n");                     \
	}                                                                          \
	printf("    }\n");                                                         \
	printf("    if(" #I "->min >= " #P "->min){\n");                           \
	printf("        GAptr = GA + " #I "->min*nBasis + " #P "->min;\n");        \
	printf("        GBptr = GB + " #I "->min*nBasis + " #P "->min;\n");        \
	for(i=0; i < I->nBasis; i++){                                              \
		for(p=0; p < P->nBasis; p++){                                          \
			printf("        GAptr[%d] -= G12[%d];",  p,2*(p*I->nBasis+i)+0);   \
			printf(       " GBptr[%d] -= G12[%d];\n",p,2*(p*I->nBasis+i)+1);   \
		}                                                                      \
		printf("        GAptr+=nBasis; GBptr+=nBasis;\n");                     \
	}                                                                          \
	printf("    }\n");                                        

	STORE_COULOMB (P,Q,I,J,p,q,i,j,p*pp+q*qq+i*ii+j);
	STORE_COULOMB (I,J,P,Q,i,j,p,q,p*pp+q*qq+i*ii+j);

	STORE_EXCHANGE(P,Q,I,J,p,q,i,j,p*pp+q*qq+i*ii+j);
	STORE_EXCHANGE(Q,P,I,J,q,p,i,j,p*pp+q*qq+i*ii+j);
	STORE_EXCHANGE(P,Q,J,I,p,q,j,i,p*pp+q*qq+i*ii+j);
	STORE_EXCHANGE(Q,P,J,I,q,p,j,i,p*pp+q*qq+i*ii+j);

	printf(
	"}\n"
	);
#undef STORE_COULOMB
#undef STORE_EXCHANGE
}


// genRStore: generate the "unrolled" version of RstoreQuartetEE code
//
// Feb 24, 2013 - Teepanis Chachiyo
//     Initially copied from genStore
//
void genRStore(
	char *str,                   // name string
	struct GTOShell_t *P,        // pointer to P shell
	struct GTOShell_t *Q,        // pointer to Q shell
	struct GTOShell_t *I,        // pointer to I shell
	struct GTOShell_t *J){       // pointer to J shell

	int p,q,i,j;
	int pp,qq,ii;

	ii = J->nBasis;
	qq = I->nBasis * J->nBasis;
	pp = Q->nBasis * I->nBasis * J->nBasis;

	printf(
	"void Rstore%sQuartetEE(\n"
	"    const struct GTOShell_t *P,\n"
	"    const struct GTOShell_t *Q,\n"
	"    const struct GTOShell_t *I,\n"
	"    const struct GTOShell_t *J,\n"
	"    int nBasis,\n"
	"    double *GA,\n"
	"    const double *PT, const double *PA,\n"
	"    const double *EEStore){\n"

	"    const double *PTptr, *PAptr;\n"
	"    register double Ga;\n"

	"    double *GAptr;\n"
	,str);


#define STORE_COULOMB(P,Q,I,J,p,q,i,j,nEE)                                     \
	printf("    PTptr = PT + " #I "->min * nBasis + " #J "->min;\n");          \
	for(i=0; i < I->nBasis; i++){                                              \
		for(j=0; j < J->nBasis; j++)                                           \
			printf("    P12[%d] = PTptr[%d];\n",i*J->nBasis+j,j);              \
		printf("    PTptr+=nBasis;\n");                                        \
	}                                                                          \
                                                                               \
	for(p=0; p < P->nBasis; p++){                                              \
		for(q=0; q < Q->nBasis; q++){                                          \
			printf("    Ga=0.0;\n");                                           \
			for(i=0; i < I->nBasis; i++)                                       \
			for(j=0; j < J->nBasis; j++){                                      \
				printf("    Ga += P12[%d]*EEStore[%d];\n",i*J->nBasis+j,nEE);  \
			}                                                                  \
			printf("    G12[%d]  = Ga+Ga;\n",p*Q->nBasis+q);                   \
		}                                                                      \
    }                                                                          \
                                                                               \
	printf("    if(" #P "->min >= " #Q "->min){\n");                           \
	printf("        GAptr = GA + " #P "->min*nBasis + " #Q "->min;\n");        \
	for(p=0; p < P->nBasis; p++){                                              \
		for(q=0; q < Q->nBasis; q++){                                          \
			printf("        GAptr[%d] += G12[%d];",  q,p*Q->nBasis+q);         \
		}                                                                      \
		printf("        GAptr+=nBasis;\n");                                    \
	}                                                                          \
	printf("    }\n");                                                         \
	printf("    if(" #Q "->min >= " #P "->min){\n");                           \
	printf("        GAptr = GA + " #Q "->min*nBasis + " #P "->min;\n");        \
	for(q=0; q < Q->nBasis; q++){                                              \
		for(p=0; p < P->nBasis; p++){                                          \
			printf("        GAptr[%d] += G12[%d];",  p,p*Q->nBasis+q);         \
		}                                                                      \
		printf("        GAptr+=nBasis;\n");                                    \
	}                                                                          \
	printf("    }\n");


#define STORE_EXCHANGE(P,Q,I,J,p,q,i,j,nEE)                                    \
	printf("    PAptr = PA + " #Q "->min * nBasis + " #J "->min;\n");          \
	for(q=0; q < Q->nBasis; q++){                                              \
		for(j=0; j < J->nBasis; j++){                                          \
			printf("    P12[%d] = PAptr[%d];"  ,q*J->nBasis+j,j);              \
		}                                                                      \
		printf("    PAptr+=nBasis;\n");                                        \
	}                                                                          \
                                                                               \
	for(p=0; p < P->nBasis; p++)                                               \
	for(i=0; i < I->nBasis; i++){                                              \
		printf("    Ga=0.0;\n");                                               \
		for(q=0; q < Q->nBasis; q++)                                           \
		for(j=0; j < J->nBasis; j++){                                          \
			printf("    Ga += P12[%d]*EEStore[%d];",  q*J->nBasis+j,nEE);      \
		}                                                                      \
		printf("    G12[%d]  = Ga;",  p*I->nBasis+i);                          \
	}                                                                          \
                                                                               \
	printf("    if(" #P "->min >= " #I "->min){\n");                           \
	printf("        GAptr = GA + " #P "->min*nBasis + " #I "->min;\n");        \
	for(p=0; p < P->nBasis; p++){                                              \
		for(i=0; i < I->nBasis; i++){                                          \
			printf("        GAptr[%d] -= G12[%d];",  i,p*I->nBasis+i);         \
		}                                                                      \
		printf("        GAptr+=nBasis;\n");                                    \
	}                                                                          \
	printf("    }\n");                                                         \
	printf("    if(" #I "->min >= " #P "->min){\n");                           \
	printf("        GAptr = GA + " #I "->min*nBasis + " #P "->min;\n");        \
	for(i=0; i < I->nBasis; i++){                                              \
		for(p=0; p < P->nBasis; p++){                                          \
			printf("        GAptr[%d] -= G12[%d];",  p,p*I->nBasis+i);         \
		}                                                                      \
		printf("        GAptr+=nBasis;\n");                                    \
	}                                                                          \
	printf("    }\n");                                        

	STORE_COULOMB (P,Q,I,J,p,q,i,j,p*pp+q*qq+i*ii+j);
	STORE_COULOMB (I,J,P,Q,i,j,p,q,p*pp+q*qq+i*ii+j);

	STORE_EXCHANGE(P,Q,I,J,p,q,i,j,p*pp+q*qq+i*ii+j);
	STORE_EXCHANGE(Q,P,I,J,q,p,i,j,p*pp+q*qq+i*ii+j);
	STORE_EXCHANGE(P,Q,J,I,p,q,j,i,p*pp+q*qq+i*ii+j);
	STORE_EXCHANGE(Q,P,J,I,q,p,j,i,p*pp+q*qq+i*ii+j);

	printf(
	"}\n"
	);
#undef STORE_COULOMB
#undef STORE_EXCHANGE
}

int main(){
	struct GTOShell_t L;   // sp type shell
	struct GTOShell_t S;   // s  type shell
	struct GTOShell_t D;   // d  type shell

	// initialze sp type shell
	L.nBasis = 4;
	L.maxL   = 1;
	L.l[0] = 0; L.l[1] = 1; L.l[2] = 0; L.l[3] = 0;
	L.m[0] = 0; L.m[1] = 0; L.m[2] = 1; L.m[3] = 0;
	L.n[0] = 0; L.n[1] = 0; L.n[2] = 0; L.n[3] = 1;

	// initialize s type shell
	S.nBasis = 1;
	S.maxL   = 0;
	S.l[0] = 0; 
	S.m[0] = 0; 
	S.n[0] = 0;

	// initialize d type shell
	D.nBasis = 6;
	D.maxL   = 2;
	D.l[0] = 2; D.l[1] = 0; D.l[2] = 0; D.l[3] = 1; D.l[4] = 1; D.l[5] = 0;
	D.m[0] = 0; D.m[1] = 2; D.m[2] = 0; D.m[3] = 1; D.m[4] = 0; D.m[5] = 1;
	D.n[0] = 0; D.n[1] = 0; D.n[2] = 2; D.n[3] = 0; D.n[4] = 1; D.n[5] = 1;

#define GEN(A,B,C,D)                                 \
	genSetB (#A #B #C #D,&A,&B,&C,&D); printf("\n"); \
	gen     (#A #B #C #D,&A,&B,&C,&D); printf("\n"); \
	genRStore(#A #B #C #D,&A,&B,&C,&D); printf("\n");\
	genUStore(#A #B #C #D,&A,&B,&C,&D); printf("\n"); 

	// S L D 
	GEN(S,S,S,S);
	GEN(L,L,L,L);
	GEN(D,D,D,D);

	// SL pair
	GEN(S,L,L,L);
	GEN(L,S,L,L);
	GEN(L,L,S,L);
	GEN(L,L,L,S);

	GEN(L,S,S,S);
	GEN(S,L,S,S);
	GEN(S,S,L,S);
	GEN(S,S,S,L);

	GEN(S,S,L,L);
	GEN(S,L,S,L);
	GEN(S,L,L,S);
	GEN(L,S,S,L);
	GEN(L,S,L,S);
	GEN(L,L,S,S);

	// DS pair
	GEN(D,S,S,S);
	GEN(S,D,S,S);
	GEN(S,S,D,S);
	GEN(S,S,S,D);

	GEN(S,D,D,D);
	GEN(D,S,D,D);
	GEN(D,D,S,D);
	GEN(D,D,D,S);

	GEN(S,S,D,D);
	GEN(S,D,S,D);
	GEN(S,D,D,S);
	GEN(D,S,S,D);
	GEN(D,S,D,S);
	GEN(D,D,S,S);

	// DL pair
	GEN(D,L,L,L);
	GEN(L,D,L,L);
	GEN(L,L,D,L);
	GEN(L,L,L,D);

	GEN(L,D,D,D);
	GEN(D,L,D,D);
	GEN(D,D,L,D);
	GEN(D,D,D,L);

	GEN(L,L,D,D);
	GEN(L,D,L,D);
	GEN(L,D,D,L);
	GEN(D,L,L,D);
	GEN(D,L,D,L);
	GEN(D,D,L,L);

	// 2S and LD pair
	GEN(S,S,L,D);
	GEN(S,S,D,L);
	GEN(S,L,S,D);
	GEN(S,D,S,L);
	GEN(S,L,D,S);
	GEN(S,D,L,S);
	GEN(L,S,S,D);
	GEN(D,S,S,L);
	GEN(L,S,D,S);
	GEN(D,S,L,S);
	GEN(L,D,S,S);
	GEN(D,L,S,S);

	// 2L and SD pair
	GEN(L,L,S,D);
	GEN(L,L,D,S);
	GEN(L,S,L,D);
	GEN(L,D,L,S);
	GEN(L,S,D,L);
	GEN(L,D,S,L);
	GEN(S,L,L,D);
	GEN(D,L,L,S);
	GEN(S,L,D,L);
	GEN(D,L,S,L);
	GEN(S,D,L,L);
	GEN(D,S,L,L);

	// 2D and SL pair
	GEN(D,D,S,L);
	GEN(D,D,L,S);
	GEN(D,S,D,L);
	GEN(D,L,D,S);
	GEN(D,S,L,D);
	GEN(D,L,S,D);
	GEN(S,D,D,L);
	GEN(L,D,D,S);
	GEN(S,D,L,D);
	GEN(L,D,S,D);
	GEN(S,L,D,D);
	GEN(L,S,D,D);

}
