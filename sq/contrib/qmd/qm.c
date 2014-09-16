// run ab initio quantum chemistry program (G03 or SQ) 
#include <stdlib.h>
#include <stdio.h>

#include "qmd.h"
#include "qm.h"

// runG03 function for a molecule
void runG03_1(struct molecule_t *ptr){
	FILE *fd;                   // file pointer
	char sym[3];                // store symbol of atom
	int i;                      // loop index

	// open file
	fd=fopen("g03.xyz","w+");
	
	// write header
	fprintf(fd, "#P RHF/3-21G NOSYMM Force\n");
	fprintf(fd, "\ntitle\n\n");
	fprintf(fd, "0 1\n");
	
	// write atom info
	for(i=0; i < ptr->nAtom; i++){
		
		//change atomic number to symbol of atom 
		Z2sym(ptr->Z[i],sym);
		
		// write atom name
		fprintf(fd,"%s",sym);
		
		// write atom position (in angstrom)
		fprintf(fd,"%f %f %f\n",ptr->x[i]*1.0e10,ptr->y[i]*1.0e10,
		ptr->z[i]*1.0e10);
	}
	fprintf(fd, "\n");
	
	// close file
	fclose(fd);
	
	// execute gaussian program
	system("g03 g03.xyz");
}

// runG03 function for 2 molecule
void runG03_2(struct molecule_t *ptr1,struct molecule_t *ptr2){
	FILE *fd;                   // file pointer
	char sym[3];                // store symbol of atom
	int i;                      // loop index
	
	// open file
	fd=fopen("g03.xyz","w+");
	
	// write header
	fprintf(fd, "#P RHF/3-21G NOSYMM Force\n");
	fprintf(fd, "\ntitle\n\n");
	fprintf(fd, "0 1\n");
	
	// write atom info of 1st molecule
	for(i=0; i < ptr1->nAtom; i++){
		
		//change atomic number to symbol of atom 
		Z2sym(ptr1->Z[i],sym);
		
		// write atom name
		fprintf(fd, "%s ",sym);
		
		// write atom position (in angstrom)
		fprintf(fd, "%f %f %f\n", ptr1->x[i]*1.0e10,ptr1->y[i]*1.0e10,ptr1->z[i]*1.0e10);
	}		
	
	// write atom info of 2nd molecule
	for(i=0; i < ptr2->nAtom; i++){
		
		//change atomic number to symbol of atom 
		Z2sym(ptr2->Z[i],sym);
		
		// write atom name
		fprintf(fd, "%s ",sym);
		
		// write atom position (in angstrom)
		fprintf(fd, "%f %f %f\n", ptr2->x[i]*1.0e10,ptr2->y[i]*1.0e10,ptr2->z[i]*1.0e10);
	}
	fprintf(fd, "\n");
	
	// close file
	fclose(fd);
	
	// execute gaussian program
	system("g03 g03.xyz");
}

// run SQ function for a molecule
void runSQ_1(struct molecule_t *ptr){
	FILE *fd;                   // file pointer
	char sym[3];                // store symbol of atom
	int i;                      // loop index
	
	//open file
	fd=fopen("sq.xyz","w+");
	
	// write number of atom in molecule
	fprintf(fd,"%d\n\n",ptr->nAtom);
	
	// write atom info
	for(i=0; i < ptr->nAtom; i++){
		
		//change atomic number to symbol of atom 
		Z2sym(ptr->Z[i],sym);
		
		// write atom name
		fprintf(fd, "%s ",sym);
		
		// write atom position (in angstrom)
		fprintf(fd, "%f %f %f\n", ptr->x[i]*1.0e10,ptr->y[i]*1.0e10,ptr->z[i]*1.0e10);
	}
	// close file
	fclose(fd);
	
	//execute Siam Quantum programe
	system("sq sq.xyz 321g.txt -FORCE > sq.out");
}	

// run SQ function for 2 molecule
void runSQ_2(struct molecule_t *ptr1,struct molecule_t *ptr2){
	FILE *fd;                   // file pointer
	char sym[3];                // store symbol of atom
	int i;                      // loop index
	
	// open file
	fd = fopen("sq.xyz","w+");
	
	// write number of atom in molecule
	fprintf(fd,"%d\n\n",(ptr1->nAtom)+(ptr2->nAtom));
	
	// write atom info of 1st molecule
	for(i=0; i < ptr1->nAtom; i++){
		
		//change atomic number to symbol of atom 
		Z2sym(ptr1->Z[i],sym);
		
		// write atom name
		fprintf(fd, "%s ",sym);
		
		// write atom position (in angstrom)
		fprintf(fd, "%f %f %f\n", ptr1->x[i]*1.0e10,ptr1->y[i]*1.0e10,ptr1->z[i]*1.0e10);
	}
	
	// write atom info 2nd molecule
	for(i=0; i < ptr2->nAtom; i++){
		
		//change atomic number to symbol of atom 
		Z2sym(ptr2->Z[i],sym);
		
		// write atom name
		fprintf(fd, "%s ",sym);
		
		// write atom position (in angstrom)
		fprintf(fd, "%f %f %f\n", ptr2->x[i]*1.0e10,ptr2->y[i]*1.0e10,ptr2->z[i]*1.0e10);
	}
	// close file
	fclose(fd);
	
	//execute Siam Quantum programe
	system("sq sq.xyz 321g.txt -FORCE > sq.out");
}

// get force 1 molecule
void getforce(char *fname,              // name of output file from G03 or SQ
			  struct molecule_t *ptr,   // pointer to molecule
			  int nMolecule,            // number of molecule
			  char *option){            // QM keyword
	int i;                              // loop index
	FILE *fd;                           // file pointer
	double Fx;                          // store force in x direction
	double Fy;                          // store force in y direction
	double Fz;                          // store force in z direction
	
	// open .out file (G03 or SQ)
	fd = fopen(fname,"r");
	
	// search to force keyword
	if(strcmp(option,"g03")==0||strcmp(option,"G03")==0){
		if(findf(fd,6,"Number","Number","X","Y","Z","*")==EOF){
			printf("Error cannot find keyword\n");
			exit(-1);
		}
	}
	if(strcmp(option,"sq")==0||strcmp(option,"SQ")==0){
		if(findf(fd,6,"Number","Number","Fx","Fy","Fz","*")==EOF){
			printf("Error cannot find keyword\n");
			exit(-1);
		}
	}
	
	for(i=0; i< ptr->nAtom; i++){
		
		// scan force
		fscanf(fd,"%*s %*s %lf %lf %lf",&Fx,&Fy,&Fz);
		
		// calculate total force
		ptr->Fx[i]= (ptr->Fx[i])-(nMolecule-2)*Fx;
		ptr->Fy[i]= (ptr->Fy[i])-(nMolecule-2)*Fy;
		ptr->Fz[i]= (ptr->Fz[i])-(nMolecule-2)*Fz;
		
		// change unit from hartres/bohr to joule/m
		ptr->Fx[i] = ptr->Fx[i]*HARTREE2JOULE/(BOHR2ANGSTROM*ANGSTROM2METER);	
		ptr->Fy[i] = ptr->Fy[i]*HARTREE2JOULE/(BOHR2ANGSTROM*ANGSTROM2METER);	
		ptr->Fz[i] = ptr->Fz[i]*HARTREE2JOULE/(BOHR2ANGSTROM*ANGSTROM2METER); 		
	}
	
	// close file
	fclose(fd);
}

// get force 2 molecule
void getforce_pair(char *fname,           // name of out put file from G03 or SQ
			  	   struct molecule_t *ptr,   // pointer to molecule
			  	   char *option){            // QM keyword

	int i;                                   // loop index
	FILE *fd;                                // file pointer
	double Fx;                               // store force in x direction
	double Fy;                               // store force in y direction
	double Fz;                               // store force in z direction
	
	// open .out file (G03 or SQ)
	fd = fopen(fname,"r");
	
	// search to force keyword
	if(strcmp(option,"g03")==0||strcmp(option,"G03")==0){
		if(findf(fd,6,"Number","Number","X","Y","Z","*")==EOF){
			printf("Error cannot find keyword\n");
			exit(-1);
		}
	}
	if(strcmp(option,"sq")==0||strcmp(option,"SQ")==0){
		if(findf(fd,6,"Number","Number","Fx","Fy","Fz","*")==EOF){
			printf("Error cannot find keyword\n");
			exit(-1);
		}
	}
	
	for(i=0; i< ptr->nAtom; i++){
		
		// scan force
		fscanf(fd,"%*s %*s %lf %lf %lf",&Fx,&Fy,&Fz);
		
		ptr->Fx[i]=(ptr->Fx[i])+Fx;
		ptr->Fy[i]=(ptr->Fy[i])+Fy;
		ptr->Fz[i]=(ptr->Fz[i])+Fz;
		
	}
	
	// close file
	fclose(fd);
}

// get total force (in joule/m unit)
void getforce_total(struct system_t *sys,        // pointer to system
					int nMolecule,               // number of molecule
					char *option){               // QM keyword		
	int i,j;                                     // loop index
		
	// run and get force 2 molecule
	for(i=0;i<nMolecule;i++){
		for(j=0;j<nMolecule;j++){
			if(i!=j){
				if(strcmp(option,"g03")==0||strcmp(option,"G03")==0){
					runG03_2(sys->mols+i,sys->mols+j);
					getforce_pair("g03.out",sys->mols+i,option);
				}
				if(strcmp(option,"sq")==0||strcmp(option,"SQ")==0){
					runSQ_2(sys->mols+i,sys->mols+j);
					getforce_pair("sq.out",sys->mols+i,option);
				}
			}
		}
	}
	// run and get force only 1 molecule
	for(i=0;i<nMolecule;i++){
		if(strcmp(option,"g03")==0||strcmp(option,"G03")==0){
			runG03_1(sys->mols+i);
			getforce("g03.out",sys->mols+i,nMolecule,option);
		}
		if(strcmp(option,"sq")==0||strcmp(option,"SQ")==0){
			runSQ_1(sys->mols+i);
			getforce("sq.out",sys->mols+i,nMolecule,option);
		}
	}
}
