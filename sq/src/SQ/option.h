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

#ifndef OPTION_H
#define OPTION_H
struct option_t{
int molCharge;        // total molecular charge
int multiplicity;     // molecule spin multiplicity = 2*s + 1

int outVolumeType;    // type of volume information to compute and save
#define VOLUME_NONE          0
#define VOLUME_DENSITY_TOTAL 1
#define VOLUME_DENSITY_SPIN  2
#define VOLUME_MO_ALPHA      3
#define VOLUME_MO_BETA       4
#define VOLUME_POTENTIAL     5

int outWhichMO;       // molecular orbital index to compute volume information
int outFormat;        // file format to save the volume information to
#define VOLUME_CUBE   1
#define VOLUME_XSF    2

double outVolumeCut;  // cutoff accurary for computing volume information
int outGAUSSIAN;      // flag to export output in Gaussian log format

double SCFConv;       // SCF convergence threshold
double SCFDrag;       // SCF drag coefficient
double SCFCutoff;     // SCF integral schwarz cut-off
int SCFMax;           // maximum number of scf cycle
#define CONVMETHOD_DAMPING   0
#define CONVMETHOD_DIIS2     1
#define CONVMETHOD_DIIS3     2
#define CONVMETHOD_DIIS4     3
int convMethod;       // scf convergence method
int SCFAccuracy;      // how to handle integral accuracy during scf
#define SCFACCURACY_3STEP 0
#define SCFACCURACY_1STEP 1

#define SCFACCURACY_3STEP_A 1.0E-3
#define SCFACCURACY_3STEP_B 1.0E-6
#define SCFACCURACY_3STEP_C 1.0E-9

int SCFGuess;         // guess method
#define SCFGUESS_DIAG  0
#define SCFGUESS_CORE  1
#define SCFGUESS_CHECK 2

int maxMem;            // maximum memory in Mbytes to use as storage 
int loadDMatrix;       // flag to load density matrix
int saveDMatrix;       // flag to save density matrix
int saveCheck;         // flag to save checkpoint file
int saveCheckAll;      // flag to save checkpoint file every scf cycle
int loadCheck;         // flag to read molecular orbital only without SCF
int opt;               // float to request geometry optimization
int optMax;            // maximum number of geometry optimization cycle
int RHF;               // flag to use RHF
int MP2;               // flag to use MP2
int UHF;               // flag to use UHF
int force;             // flag to compute force acting on nuclei
char DMatrixFile[256]; // file name to process density matrix
char CheckFile[256];   // file name to process checkpoint file
int nCPU;              // number of cpu to perform parallel calculations
char prefixStr[256];   // prefix string for job identification
int MECP;              // flag to request MECP
int mecpMax;           // maximum number of MECP cycle
int mecpMA;            // multiplicity of the first state
int mecpMB;            // multiplicity of the second state
char DMatrixFileA[256];// file name to process density matrix of the first state
char DMatrixFileB[256];// file name to process density matrix of the second state
char CheckFileA[256];  // file name to process checkpoint file
char CheckFileB[256];  // file name to process checkpoint file 
char gaussEXE[256];    // Gaussian program execution string
char gaussINA[256];   // Gaussian input file for the first state
char gaussINB[256];   // Gaussian input file for the second state
};
void parse_option(struct option_t *opt, int argc, char *argv[]);
void option_help();
#endif

