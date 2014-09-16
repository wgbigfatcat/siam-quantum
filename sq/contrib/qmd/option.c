#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "option.h"

void parseHelp(int argc, char *argv[]){
	printf(
	"QMD - Quantum Molecular Dynamics Simulation Program      \n"
	"                                                         \n"
	"by Chutchawan Jaisuk                                     \n"
	"                                                         \n"
	"usage: %s <system.xyz> [OPTIONS]                         \n"
	"                                                         \n"
	"OPTIONS:                                                 \n"
	"    -HELP             Print out help information         \n"
	"    -MAXSTEP=<int>    Set maximum number of steps        \n"
	"    -DT=<real>        Set time step in second            \n"
	"    -V0=<int>         Set initial velocity               \n"
	"                      0 - all velocities zero            \n"
	"                      1 - random velocities              \n"
	"    -R0=<int>         Set initial displacement           \n"
	"                      0 - all displacements zero         \n"
	"                      1 - random displacements           \n"
	"    -TEM=<real>       Set temperature in Kelvin          \n"
	"    -QM=<string>      Choose G03 or SQ                   \n"
	"                                                         \n"
	,argv[0]
	);
}

void parseOption(int argc,char *argv[],struct option_t *opt){
	int i;          // loop index

	// set default value
	opt->maxStep = 1000;
	opt->dt      = 1.0E-15;
	opt->R0      = 0;
	opt->V0      = 0;
	opt->tem     = 0;
	strncpy(opt->inpKW,"sq",1024);
	
	// check required option
	if(argc < 2){
		parseHelp(argc,argv);
		exit(-1);
	}
	
	// read file name
	strncpy(opt->inpName,argv[1],1024);
	
	
	// process each option
	for(i=2; i < argc; i++){
		
		if(strncmp(argv[i],"-HELP",5)==0){
			parseHelp(argc,argv);
			exit(-1);
			continue;
		}
		
		if(strncmp(argv[i],"-MAXSTEP=",9)==0){
			opt->maxStep = atoi(argv[i]+9);
			if(opt->maxStep <= 0){
				printf("parseOption: Invalid MaxStep range\n");
				exit(-1);
			}
			continue;
		}
		
		if(strncmp(argv[i],"-DT=",4)==0){
			opt->dt = atof(argv[i]+4);
			if(opt->dt <= 0){
				printf("parseOption: Invalid dt range\n");
				exit(-1);
			}
			continue;
		}
		
		if(strncmp(argv[i],"-R0=",4)==0){
			opt->R0 = atoi(argv[i]+4);
			if(opt->R0!=0&&opt->R0!=1){
				printf("parseOption: Invalid R0 range\n");
				exit(-1);
			}
			continue;
		}
		
		if(strncmp(argv[i],"-V0=",4)==0){
			opt->V0 = atoi(argv[i]+4);
			if(opt->V0!=0&&opt->V0!=1){
				printf("parseOption: Invalid V0 range\n");
				exit(-1);
			}
			continue;
		}
		
		if(strncmp(argv[i],"-TEM=",5)==0){
			opt->tem = atof(argv[i]+5);
			continue;
		}
		
		if(strncmp(argv[i],"-QM=",4)==0){
			strncpy(opt->inpKW,argv[i]+4,1024);
			if(strcmp(opt->inpKW,"g03")!=0&&strcmp(opt->inpKW,"G03")!=0
				&&strcmp(opt->inpKW,"sq")!=0&&strcmp(opt->inpKW,"SQ")!=0){
				printf("parseOption: Invalid key word\n");
				exit(-1); 
			}
			continue;
		}
		
		parseHelp(argc,argv);
		exit(-1);
	}
}
