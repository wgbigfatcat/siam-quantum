#ifndef OPTION_H
#define OPTION_H
struct option_t{
	int maxStep;        // maximum number of time step to simulate
	double dt;          // time step in seconds
	double tem;         // temperature of system in Kelvin
	int V0;             // initial velocity 
	int R0;             // initial diplacement
	char inpName[1024]; // file name for reading the system info in xyz format
	char inpKW[1024];   // key word of ab initio quantum chemistry program  
};
#endif
