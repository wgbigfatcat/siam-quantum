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
#include <string.h>
#include "option.h"

// option_help : print out all available options
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
// Oct 20, 2010 - Theerapon Khumla
//  Add force option
//
// Oct 22, 2010 - Teepanis Chachiyo
//  Add geometry optimization and SCFConv option
//
// May 21, 2011 - Teepanis Chachiyo
//  Add DIIS option and set it to the default option
//
// May 21, 2011 - Teepanis Chachiyo
//  Add maxmem option
//
// June 4, 2011 - Nanta Sophonrat
//  Add potential option
//
// Oct 4, 2012 - Teepanis Chachiyo
//  Add and DIIS checkpoint options
//
// Oct 20, 2012 - Teepanis Chachiyo
//  Add MP2 options
//
// Jan 26, 2013 - Teepanis Chachiyo
//  Add DIIS3 options and set default to DIIS4
//
// March, 17, 2013 - Teepanis Chachiyo
//  Change default of SCFDrag to 0.25
//
void option_help(){
	printf(
	"[OPTIONS] : \n"
	"                                                     \n"
	"Ab Initio Method:                                    \n"
	"-Q=INT         Set total molecular charge (default=0)\n"
	"-M=INT         Set molecular spin multiplicity (M=2S+1)\n"
	"-RHF           Restricted Hartree-Fock (default for singlet state)\n"
	"-UHF           Unrestricted Hartree-Fock (default if M > 1)\n"
	"-FORCE         Calculate force acting on nuclei\n"
	"-OPT           Request geometry optimization\n"
	"-MP2           Request MP2 energy calculations\n"
	"-MECP=INT,INT  Request MECP between two spin multiplicities\n"
	"                                                                     \n"
	"SCF Cycle:                                                           \n"
	"-GUESS=DIAG    Use identity density matrix as initial guess (default)\n"
	"-GUESS=CORE    Use density from core hamiltonian as initial guess\n"
	"-GUESS=CHECK   Use density from checkpoint as initial guess\n"
	"-SCFDIIS       Use 4-point DIIS method for convergence (default)\n"
	"-SCFDIIS3      Use 3-point DIIS method for convergence\n"
	"-SCFDIIS2      Use 2-point DIIS method for convergence\n"
	"-SCFDAMP       Use simple weighting method for convergence\n"
	"-SCFDRAG=REAL  Set SCF drag coefficient in range 0.0-1.0 (default=0.25)\n"
	"-SCFCONV=REAL  Set SCF convergence threshold (default=1.0E-6)\n"
	"-SCFMAX=INT    Set maximum number of scf cycle (default=80)\n"
	"-SCFACC=3STEP  Use increasing integral accuracy in 3 steps (default)\n"
	"-SCFACC=1STEP  Use fixed integral accuracy\n"
	"-MAXMEM=INT    Set maximum memory per CPU in Mbyte (default=250)\n"
	"                                                        \n"
	"Checkpoint File:\n"
	"-LCHECK        Do not perform SCF but load info from checkpoint\n"
	"-SCHECK        Save checkpoint file at the end (default=no)\n"
	"-SCHECK=ALL    Save checkpoint file every scf cycle (default=no)\n"
	"-FCHECK=STR    Specify file name for checkpoint (default=checkpoint.txt)\n"
	"-LDMATRIX      Load density matrix at the beginning (default=no)\n"
	"-SDMATRIX      Save density matrix at the end (default=no)\n"
	"-FDMATRIX=STR  Specify file name for density matrix (default=dmatrix.txt)\n"
	"                                                         \n"
	"Volume Output:                                           \n"
	"-DENSITY       Print electron density  volume information\n"
	"-POTENTIAL     Print electric potential volume information\n"
	"-MOUP=INT      Print spin up mo. volume info (index starts at 1)\n"
	"-MODN=INT      Print spin dn mo. volume info (index starts at 1)\n"
	"-VOLCUT=REAL   Set accuracy for computing volume info (default=1.0E-4)\n"
	"-XSF           Volume info. will be in XSF  format to 'volume.xsf'\n"
	"-CUBE          Volume info. will be in CUBE format to 'volume.cube'\n"
	"-GAUSSIAN      Emulate Gaussian output to 'gaussian.log' (for GabeEdit)\n"
	"                                                         \n"
	"Parallel Run:                                            \n"
	"-NCPU=INT      Set the number of CPUs (default=1)\n"
	"-PREFIX=STR    Set prefix string for the job (default=SQ)\n"
	"                                                         \n"
	"Geometry Optimization:                                   \n"
	"-OPTMAX=INT    Maximum number of iterations (default=30) \n"
	"                                                         \n"                                
	"Minimum Energy Crossing Point (MECP):                    \n"
	"-MECPMAX=INT   Maximum number of iterations (default=30) \n"
	"-FCHECKA=STR   State A checkpoint file name (default=checkpointA.txt)\n"
	"-FCHECKB=STR   State B checkpoint file name (default=checkpointB.txt)\n"
	"-FDMATRIXA=STR State A density matrix file name (default=dmatrixA.txt)\n"
	"-FDMATRIXB=STR State B density matrix file name (default=dmatrixB.txt)\n"
	"-GAUSSINA=STR  State A Gaussian input file name (excluding .com)\n"
	"-GAUSSINB=STR  State B Gaussian input file name (excluding .com)\n"
	"-GAUSSEXE=STR  Gaussian program execution string\n"
	);
}


// parse_option : parse argv program argument into options
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//  Add printing molecular orbital option
//
// July 11, 2010 - Teepanis Chachiyo
//  Unrestricted calculation options
//
// Oct 20, 2010 - Theerapon Khumla
//  Add force option
//
// Oct 21, 2010 - Teepanis Chachiyo
//  Add geoemtry optimization and optin SCFConv
//
// May 21, 2011 - Teepanis Chachiyo
//  Add DIIS and make it a default option
//  Bugfix, compare strcmp the whole string if the option needs no argument
//
// May 21, 2011 - Teepanis Chachiyo
//  Add MAXMEM option
//
// June 4, 2011 - Nanta Sophonrat
//  Add potential option
//
// Sep 30, 2012 - Teepanis Chachiyo
//  Add DIIS2 option
//
// Oct 4, 2012 - Teepanis Chachiyo
//  Add checkpoint options
//
// Nov 26, 2012 - Teepanis Chachiyo
//  Add checkpoint at every scf and multiple integral accuracy cutoff
//
// Jan 26, 2013 - Teepanis Chachiyo
//  Add DIIS-4 convergence method
//
// Feb 24, 2013 - Teepanis Chachiyo
//  SCFMax is now 80 by default
//
// March 17, 2013 - Teepanis Chachiyo
//  Add MECP
//
// Jan 16, 2014 - Teepanis Chachiyo
//	Handle whichMO index (which the users think that the index start at 1)
//
void parse_option(struct option_t *opt, int argc, char *argv[]){

	int i;

	// set default
	opt->molCharge       = 0;
	opt->multiplicity    = 0;
	opt->RHF             = 0;
	opt->MP2             = 0;
	opt->UHF             = 0;
	opt->force           = 0;
	opt->opt             = 0;
	opt->outVolumeType   = VOLUME_NONE;
	opt->outWhichMO      = 0;
	opt->outFormat       = VOLUME_XSF;
	opt->outVolumeCut    = 1.0E-4;
	opt->outGAUSSIAN     = 0;
	opt->SCFGuess        = SCFGUESS_DIAG;
	opt->SCFConv         = 1.0E-6;
	opt->SCFCutoff       = 1.0E-15;
	opt->SCFDrag         = 0.25;
	opt->SCFMax          = 80;
	opt->convMethod      = CONVMETHOD_DIIS4;
	opt->SCFAccuracy     = SCFACCURACY_3STEP;
	opt->maxMem          = 250;
	opt->loadDMatrix     = 0;
	opt->saveDMatrix     = 0;
	strcpy(opt->DMatrixFile,"dmatrix.txt");
	opt->saveCheck       = 0;
	opt->saveCheckAll    = 0;
	opt->loadCheck       = 0;
	strcpy(opt->CheckFile,"checkpoint.txt");
	opt->nCPU            = 1;
	strcpy(opt->prefixStr,"SQ");
	opt->MECP            = 0;
	opt->mecpMax         = 30;
	opt->mecpMA          = 0;
	opt->mecpMB          = 0;
	strcpy(opt->DMatrixFileA,"dmatrixA.txt");
	strcpy(opt->DMatrixFileB,"dmatrixB.txt");
	strcpy(opt->CheckFileA,"checkpointA.txt");
	strcpy(opt->CheckFileB,"checkpointB.txt");
	strcpy(opt->gaussEXE,"\0");
	strcpy(opt->gaussINA,"\0");
	strcpy(opt->gaussINB,"\0");
	opt->optMax          = 30;

	// loop throu all options
	for(i=3; i < argc; i++){
		if(strncmp(argv[i],"-Q=",3)==0){
			opt->molCharge = atoi(argv[i]+3);
			continue;
		}
		if(strncmp(argv[i],"-M=",3)==0){
			opt->multiplicity = atoi(argv[i]+3);
			continue;
		}
		if(strcmp(argv[i],"-RHF")==0){
			opt->RHF = 1;
			continue;
		}
		if(strcmp(argv[i],"-MP2")==0){
			opt->MP2 = 1;
			continue;
		}
		if(strcmp(argv[i],"-UHF")==0){
			opt->UHF = 1;
			continue;
		}
		if(strcmp(argv[i],"-FORCE")==0){
			opt->force = 1;
			continue;
		}
		if(strcmp(argv[i],"-OPT")==0){
			opt->opt = 1;
			continue;
		}
		if(strncmp(argv[i],"-OPTMAX=",8)==0){
			opt->optMax = atoi(argv[i]+8);
			continue;
		}
		if(strcmp(argv[i],"-DENSITY")==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_DENSITY_TOTAL;
			continue;
		}
		if(strcmp(argv[i],"-POTENTIAL")==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_POTENTIAL;
			continue;
		}
		if(strncmp(argv[i],"-MOUP=",6)==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_MO_ALPHA;
			opt->outWhichMO = atoi(argv[i]+6);

			// users think that the index starts at 1
			opt->outWhichMO = opt->outWhichMO - 1;
			continue;
		}
		if(strncmp(argv[i],"-MODN=",6)==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_MO_BETA;
			opt->outWhichMO = atoi(argv[i]+6);

			// users think that the index starts at 1
			opt->outWhichMO = opt->outWhichMO - 1;
			continue;
		}
		if(strncmp(argv[i],"-VOLCUT=",8)==0){
			opt->outVolumeCut = atof(argv[i]+8);
			continue;
		}
		if(strcmp(argv[i],"-GAUSSIAN")==0){
			opt->outGAUSSIAN  = 1;
			continue;
		}
		if(strcmp(argv[i],"-XSF")==0){
			opt->outFormat  = VOLUME_XSF;
			continue;
		}
		if(strcmp(argv[i],"-CUBE")==0){
			opt->outFormat = VOLUME_CUBE;
			continue;
		}
		if(strcmp(argv[i],"-GUESS=DIAG")==0){
			opt->SCFGuess = SCFGUESS_DIAG;
			continue;
		}
		if(strcmp(argv[i],"-GUESS=CORE")==0){
			opt->SCFGuess = SCFGUESS_CORE;
			continue;
		}
		if(strcmp(argv[i],"-GUESS=CHECK")==0){
			opt->SCFGuess = SCFGUESS_CHECK;
			continue;
		}
		if(strcmp(argv[i],"-SCFACC=3STEP")==0){
			opt->SCFAccuracy = SCFACCURACY_3STEP;
			continue;
		}
		if(strcmp(argv[i],"-SCFACC=1STEP")==0){
			opt->SCFAccuracy = SCFACCURACY_1STEP;
			continue;
		}
		if(strcmp(argv[i],"-SCFDIIS")==0){
			opt->convMethod = CONVMETHOD_DIIS4;
			continue;
		}
		if(strcmp(argv[i],"-SCFDIIS3")==0){
			opt->convMethod = CONVMETHOD_DIIS3;
			continue;
		}
		if(strcmp(argv[i],"-SCFDIIS2")==0){
			opt->convMethod = CONVMETHOD_DIIS2;
			continue;
		}
		if(strcmp(argv[i],"-SCFDAMP")==0){
			opt->convMethod = CONVMETHOD_DAMPING;
			continue;
		}
		if(strncmp(argv[i],"-MAXMEM=",8)==0){
			opt->maxMem = atoi(argv[i]+8);
			continue;
		}
		if(strncmp(argv[i],"-SCFCONV=",9)==0){
			opt->SCFConv = atof(argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-SCFDRAG=",9)==0){
			opt->SCFDrag = atof(argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-SCFMAX=",8)==0){
			opt->SCFMax = atoi(argv[i]+8);
			continue;
		}
		if(strcmp(argv[i],"-LDMATRIX")==0){
			opt->loadDMatrix = 1;
			continue;
		}
		if(strcmp(argv[i],"-SDMATRIX")==0){
			opt->saveDMatrix = 1;
			continue;
		}
		if(strncmp(argv[i],"-FDMATRIX=",10)==0){
			strcpy(opt->DMatrixFile, argv[i]+10);
			continue;
		}
		if(strcmp(argv[i],"-SCHECK")==0){
			opt->saveCheck = 1;
			continue;
		}
		if(strcmp(argv[i],"-SCHECK=ALL")==0){
			opt->saveCheckAll = 1;
			continue;
		}
		if(strcmp(argv[i],"-LCHECK")==0){
			opt->loadCheck = 1;
			continue;
		}
		if(strncmp(argv[i],"-FCHECK=",8)==0){
			strcpy(opt->CheckFile, argv[i]+8);
			continue;
		}
		if(strncmp(argv[i],"-NCPU=",6)==0){
			opt->nCPU = atoi(argv[i]+6);
			continue;
		}
		if(strncmp(argv[i],"-PREFIX=",8)==0){
			strcpy(opt->prefixStr, argv[i]+8);
			continue;
		}
		if(strncmp(argv[i],"-MECP=",6)==0){
			if(sscanf(argv[i]+6,"%d,%d",&opt->mecpMA,&opt->mecpMB)!=2){
				printf("parse_option - error cannot recognize option %s\n",argv[i]);
				exit(-1);
			}
			opt->MECP = 1;
			continue;
		}
		if(strncmp(argv[i],"-MECPMAX=",9)==0){
			opt->mecpMax = atoi(argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-FDMATRIXA=",11)==0){
			strcpy(opt->DMatrixFileA, argv[i]+11);
			continue;
		}
		if(strncmp(argv[i],"-FDMATRIXB=",11)==0){
			strcpy(opt->DMatrixFileB, argv[i]+11);
			continue;
		}
		if(strncmp(argv[i],"-FCHECKA=",9)==0){
			strcpy(opt->CheckFileA, argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-FCHECKB=",9)==0){
			strcpy(opt->CheckFileB, argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-GAUSSEXE=",10)==0){
			strcpy(opt->gaussEXE, argv[i]+10);
			continue;
		}
		if(strncmp(argv[i],"-GAUSSINA=",10)==0){
			strcpy(opt->gaussINA, argv[i]+10);
			continue;
		}
		if(strncmp(argv[i],"-GAUSSINB=",10)==0){
			strcpy(opt->gaussINB, argv[i]+10);
			continue;
		}

		// cannot recognize parameter
		printf("parse_option - error cannot recognize option %s\n", argv[i]);
		exit(-1);
	}

	// validate 
	if(opt->outVolumeType==VOLUME_MO_ALPHA || opt->outVolumeType==VOLUME_MO_BETA)
	if(opt->outWhichMO <= 0){
		printf("parse_option - error molecular orbital index should be greater than zero\n");
		exit(-1);
	}

	// validate SCFConv
	if(opt->SCFConv <= 0.0){
		printf("parse_option - error invalid SCFConv range\n");
		exit(-1);
	}

	// validate volumeCut
	if(opt->outVolumeCut <= 0.0){
		printf("parse_option - error invalid volumeCut range\n");
		exit(-1);
	}

	// validate SCFDrag
	if(opt->SCFDrag <= 0.0 || opt->SCFDrag > 1){
		printf("parse_option - error invalid SCFDrag range\n");
		exit(-1);
	}

	// validate SCFMax
	if(opt->SCFMax < 0){
		printf("parse_option - error invalid SCFMax range\n");
		exit(-1);
	}

	// validate maxMem
	if(opt->maxMem < 0){
		printf("parse_option - error invalid MAXMEM range\n");
		exit(-1);
	}

	// validate RHF and UHF choice
	if(opt->RHF + opt->UHF > 1){
		printf("parse_option - error can not choose both RHF and UHF\n");
		exit(-1);
	}

	// validate checkpoint related options
	if(opt->loadCheck && opt->opt){
		printf("parse_option - error OPT cannot be used with LCHECK\n");
		exit(-1);
	}
	if(opt->loadCheck && (opt->saveCheck || opt->saveCheckAll)){
		printf("parse_option - error LCHECK cannot be used with saving checkpoints\n");
		exit(-1);
	}
	if(opt->saveCheck && opt->saveCheckAll){
		printf("parse_option - error SCHECK cannot be used with SCHECK=ALL");
		exit(-1);
	}

	// optimization does not support mp2 yet
	if(opt->opt && opt->MP2){
		printf("parse_option - error OPT does not support MP2 at the moment\n");
		exit(-1);
	}

	// validate nCPU
	if(opt->nCPU <= 0){
		printf("parse_option - error invalid number of cpus\n");
		exit(-1);
	}

	// validate MECP
	if(opt->MECP){
		if(opt->opt){
			printf("parse_option - error MECP cannot be used with OPT\n");
			exit(-1);
		}
		if(opt->MP2){
			printf("parse_option - error MECP cannot be used with MP2\n");
			exit(-1);
		}
	}

}

