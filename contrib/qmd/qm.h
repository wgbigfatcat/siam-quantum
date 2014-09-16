#ifndef QM_H
#define QM_H

void runG03_1(struct molecule_t *ptr);
void runG03_2(struct molecule_t *ptr1,struct molecule_t *ptr2);
void runSQ_1(struct molecule_t *ptr);
void runSQ_2(struct molecule_t *ptr1,struct molecule_t *ptr2);
void getforce(char *fname,                    // name of output file from G03
			  struct molecule_t *ptr,         // molecule pointer
			  int nMolecule,                  // pointer to molecule
			  char *option);                  // QM key word
void getforce_pair(char *fname,               // name of output file from G03
			  	   struct molecule_t *ptr,    // pointer to molecule
				   char *option);             // QM keyword
void getforce_total(struct system_t *sys,     // pointer to system
					int nMolecule,            // number of molecule
					char *option);            // QM keword
#endif
