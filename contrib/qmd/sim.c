#include <stdlib.h>
#include <stdio.h>

#include "qmd.h"
#include "option.h"
#include "sim.h"

// set initial displacement
void set_displacement(struct molecule_t *ptr,   // pointer to molecule
					  int R0){                  // key word
	int i;                                      // loop index 
	
	if(R0==0){};
	if(R0==1){
		for(i=0; i < ptr->nAtom; i++){
			ptr->x[i]=ptr->x[i]+(pow(-1,rand())*(0.01*(rand()%3))*1.0e-10);
			ptr->y[i]=ptr->y[i]+(pow(-1,rand())*(0.01*(rand()%3))*1.0e-10);
			ptr->z[i]=ptr->z[i]+(pow(-1,rand())*(0.01*(rand()%3))*1.0e-10);
		}
	}
}

// set initail velocity
void set_velocity(struct molecule_t *ptr,        // pointer to molecule
				  int V0){                       // key word
	int i;                                       // loop index 
	
	for ( i=0; i < ptr->nAtom; i++){
		if(V0==0){
			ptr->vx[i]=0;
			ptr->vy[i]=0;
			ptr->vz[i]=0;
		}
		if(V0==1){
			ptr->vx[i]=pow(-1,rand())*(300+(rand()%40));
			ptr->vy[i]=pow(-1,rand())*(300+(rand()%40));
			ptr->vz[i]=pow(-1,rand())*(300+(rand()%40));
			printf("\nvx= %lf vy= %lf vz= %lf",ptr->vx[i],ptr->vy[i],ptr->vz[i]);
		}
	}
}

// update velocity of all atom in system
void updateNewton_velocity(
	struct molecule_t *ptr,   // pointer to molecule
	struct option_t *opt){    // pointer to option
	int i;                    // loop index
	double dt;                // time step
	
	// set time step
	dt = opt->dt;
		
	for ( i=0; i < ptr->nAtom; i++){
		
		// update velocity
		ptr->vx[i] = ptr->Fx[i]*(dt/ptr->m[i])+ptr->vx[i];
		ptr->vy[i] = ptr->Fy[i]*(dt/ptr->m[i])+ptr->vy[i];
		ptr->vz[i] = ptr->Fz[i]*(dt/ptr->m[i])+ptr->vz[i];
	}
}

// update position of all atom in system
void updateNewton_position(
	struct molecule_t *ptr,   // pointer to molecule
	struct option_t *opt){    // pointer to option
	int i;                    // loop index
	double dt;                // time step
	
	// set time step
	dt = opt->dt;
	
	for ( i=0; i < ptr->nAtom; i++){

		// update position
		ptr->x[i] = ptr->vx[i]*dt+ptr->x[i];
		ptr->y[i] = ptr->vy[i]*dt+ptr->y[i];
		ptr->z[i] = ptr->vz[i]*dt+ptr->z[i];
	}
}

// control temperature by rescale velocities
void control_velocity(int totalAtom,           // total atom
					  struct system_t *sys,    // pointer to syatem
					  struct molecule_t *ptr,  // pointer to molecule
					  double tem){             // temperature
	double Ek=0;                               // total kinetic energy
	double lamda;                              // scaling factor
	int i,j;                                   // loop index
	
	// Do not rescale velocity
	if(tem<=0){};
	
	// Rescale velocity 
	if(tem>0){
		
		// calulate total kinetic energy
		for(i=0;i<sys->nMolecule;i++){
		get_kinetic(sys->mols+i,&Ek);
		}

		// calculate scaling factor
		lamda = sqrt((3/2)*BOLZMANN*tem*totalAtom/Ek);
				
		// rescale velocity
		for(i=0;i<sys->nMolecule;i++){
			for(j=0;j<ptr->nAtom;j++){
				ptr->vx[j] = ptr->vx[j]*lamda;
				ptr->vy[j] = ptr->vy[j]*lamda;
				ptr->vz[j] = ptr->vz[j]*lamda;
				printf("\nvx= %lf vy= %lf vz= %lf",ptr->vx[j],ptr->vy[j],ptr->vz[j]);
			}
			ptr = ptr+1;
		}
	}
}

// calulate total kinetic energy
void get_kinetic(struct molecule_t *ptr,double *Ek){
	int i;
	double E;
	
	for(i=0;i<ptr->nAtom;i++){
		E = 0.5*(ptr->m[i])*(pow(ptr->vx[i],2)+pow(ptr->vy[i],2)
		    +pow(ptr->vz[i],2));
	   *Ek = E+(*Ek);       
	}
}
