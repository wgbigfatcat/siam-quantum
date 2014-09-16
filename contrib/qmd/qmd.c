#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qmd.h"
#include "option.h"

int main(int argc, char *argv[]){
	struct option_t opt;                        // option data
	struct system_t *sys;                       // pointer to system
	struct molecule_t *mols;                    // pointer to molecule
	int totalAtom=0;                            // total atom
	int i,j;                                    // loop index 
	FILE *fd;                                   // file pointer
	FILE *fd1;                                  // file pointer
	FILE *fd2;                                  // file pointer
	FILE *fd3;                                  // file pointer

	parseOption(argc, argv, &opt);

	// openfile
	fd=fopen(opt.inpName,"r");
	if(fd==NULL){
		printf("main: cannot read file %s\n", opt.inpName);
		exit(-1);
	}
	
	// open .xyz file
	fd1=fopen("vmd.xyz","w");
	
	// open .txt file
	fd2=fopen("Coordinate.txt","w");
	fd3=fopen("O-O_distance.txt","w");
	
	// read system
	sys = read_system(fd,&totalAtom);
	
	for(i=0; i < sys->nMolecule; i++){
		print_molecule(sys->mols + i);            // print molecule in system
		Z2kg(sys->mols + i);                      // store mass in kg 
		set_displacement(sys->mols+i,opt.R0);     // set initial displacement
		set_velocity(sys->mols + i,opt.V0);       // set initial velocity
	}
	
	// control velocity
	control_velocity(totalAtom,sys,sys->mols,opt.tem);
	
	// print initial position to .xyz file
	fprintf(fd1,"%d\n\n",totalAtom);
	for(i=0; i < sys->nMolecule; i++)	outXYZ(sys->mols+i,fd1); 
	
	// print initial position to .txt file
	for(i=0 ; i < sys->nMolecule;i++)	outTXT(sys->mols+i,fd2);
	fprintf(fd2,"\n");
	
	// calculate O-O radial distance
	radial_distance(sys,sys->nMolecule,fd3);	
	
	printf("\ncontinue to perform simulation %d cycle",opt.maxStep);
	
	// perform simulation
	for(j=1;j<=opt.maxStep;j++){
		printf("\n\tstep = %d",j);
		
		// getforce
		getforce_total(sys,sys->nMolecule,opt.inpKW);
		
		// update velocity 
		for(i=0; i < sys->nMolecule; i++)	
		updateNewton_velocity(sys->mols+i,&opt);
		
		// control velocity
		control_velocity(totalAtom,sys,sys->mols,opt.tem);
		
		// update position
		for(i=0; i < sys->nMolecule; i++)	
		updateNewton_position(sys->mols+i,&opt);
		
		// print position to .xyz file format
		fprintf(fd1,"%d\n\n",totalAtom);
		for(i=0; i < sys->nMolecule; i++)	outXYZ(sys->mols+i,fd1);
		
		// print position to .txt file format  
		for(i=0 ; i < sys->nMolecule;i++)	outTXT(sys->mols+i,fd2);
		fprintf(fd2,"\n");	
		
		// calculate O-O radial distance
		radial_distance(sys,sys->nMolecule,fd3);	
	}	
		
	// close file
	fclose(fd);
	fclose(fd1);
	fclose(fd2);
	fclose(fd3);			

	printf("\nfinish");
	return 0;
}

// read_system: read file in xyz format and construct system information
//
// Oct 21, 2010 - Chutchawan Jaisuk
//   Initial implementation
//
struct system_t *read_system(FILE *fd, int *totalAtom){
	struct system_t *SysPtr;                    // pointer to system
	struct molecule_t *ptr;                     // pointer to molecule
	int i,j,k;                                  // index 
	char str[2];                                // store symbol of atom
	
	// allocate data slot for system
	SysPtr = calloc(1,sizeof(struct system_t));
	SysPtr->nMolecule = 0;
	
	// allocte memory to store molecule_t
	SysPtr->mols = calloc(MAXMOLECULE, sizeof(struct molecule_t));
	
	for(k=0; k < MAXMOLECULE; k++){
	
		ptr = SysPtr->mols + k;
	
		if(fscanf(fd,"%d",&i)!=1){
			printf("\nTerminate Reading System\n\n");
			break;
		}
		
		ptr->nAtom = i;
	
		// allocate memory for array of structure
		ptr->x = calloc(i,sizeof(double));
		ptr->y = calloc(i,sizeof(double));
		ptr->z = calloc(i,sizeof(double));
		ptr->vx = calloc(i,sizeof(double));
		ptr->vy = calloc(i,sizeof(double));
		ptr->vz = calloc(i,sizeof(double));
		ptr->Fx = calloc(i,sizeof(double));
		ptr->Fy = calloc(i,sizeof(double));
		ptr->Fz = calloc(i,sizeof(double));
		ptr->m = calloc(i,sizeof(double));
		ptr->Z = calloc(i,sizeof(int));
	
		// scan position 
		for(j=0;j<i;j++){
		fscanf(fd,"%s %lf %lf %lf",str,&(ptr->x[j]),&(ptr->y[j]),&(ptr->z[j]));
			
			//change Angstrom to SI Unit (meter)
			ptr->x[j]=ptr->x[j]*ANGSTROM2METER;
			ptr->y[j]=ptr->y[j]*ANGSTROM2METER;
			ptr->z[j]=ptr->z[j]*ANGSTROM2METER;
			
			// change symbol of atom to atomic number
			ptr->Z[j]=sym2Z(str,0);
			}
			
		SysPtr->nMolecule = SysPtr->nMolecule + 1;
		
		*totalAtom = *totalAtom+ptr->nAtom;
	}
	return SysPtr;
}

// print molecule info
void print_molecule(struct molecule_t *ptr){   
	int i;                                   // loop index
	char sym[3];                             // store symbol of atom
	
	printf("\t%d\n",ptr->nAtom);
	
	for(i=0;i< ptr->nAtom ;i++){
		Z2sym(ptr->Z[i],sym);
		printf("\t%s %lf %lf %lf \n",sym,ptr->x[i]*1.0e10,ptr->y[i]*1.0e10,
		       ptr->z[i]*1.0e10);
	}
	printf("\n");
}

// print atom info on .xyz file
void outXYZ(struct molecule_t *ptr,           // molecule pointer
			FILE *fd2){                       // file pointer
	int i;                                    // loop index
	char sym[3];                              // store symbol of atom
	
	// write atom info
	for(i=0; i< ptr->nAtom;i++){
		Z2sym(ptr->Z[i],sym);
		fprintf(fd2,"%s %lf %lf %lf \n",sym,ptr->x[i]*1.0e10,ptr->y[i]*1.0e10,
		        ptr->z[i]*1.0e10);
	}
}

// print atom info on .txt file
void outTXT(struct molecule_t *ptr,           // molecule pointer
			FILE *fd3){                       // file pointer
	int i;                                    // loop index
	char sym[3];                              // store symbol of atom
		
	//write atom info
	for(i=0 ; i < ptr->nAtom ; i++){
		fprintf(fd3,"%e %e %e ",ptr->x[i],ptr->y[i],ptr->z[i]);
	}
}

// analyse radial distance
void radial_distance(struct system_t *ptr,      // pointer to system
					 int nMolecule,             // number of molecule
					 FILE *fd3){                // file pointer
	int i,j;                                    // loop index

	for(i=0;i<nMolecule;i++){
		for(j=0;j<nMolecule;j++){
			if(i<j){
				OO_distance(ptr->mols+i,ptr->mols+j,fd3);
			}
		}
	}
	fprintf(fd3,"\n");
}

// analyse O-O distance 
void OO_distance(struct molecule_t *ptr1,       // pointer to molecule
			   struct molecule_t *ptr2,         // pointer to molecule
			   FILE *fd3){                      // file pointer
	int i,j;                // loop index
	double Ox,Oy,Oz,OO;     // O-O distance in x,y and z direction
	
	for(i=0;i<ptr1->nAtom;i++){
		for(j=0;j<ptr2->nAtom;j++){
			if(ptr1->Z[i]==8&&ptr2->Z[j]==8){
				Ox = ptr1->x[i]-ptr2->x[j];
				Oy = ptr1->y[i]-ptr2->y[j];
				Oz = ptr1->z[i]-ptr2->z[j];
				OO = sqrt(Ox*Ox+Oy*Oy+Oz*Oz);
				fprintf(fd3,"%lf\t",OO*1.0e10);
			}
		}
	}
}
