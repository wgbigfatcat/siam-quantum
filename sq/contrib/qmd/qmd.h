#ifndef QMD_H
#define QMD_H
#define MAXMOLECULE 200
#define PROTONMASS 1.67262158e-27
#define HARTREE2JOULE 4.3597482e-18
#define BOHR2ANGSTROM 0.529177249 
#define ANGSTROM2METER 1.0e-10
#define BOLZMANN 1.3806504e-23

// molecular data structure definition
struct molecule_t{
	int nAtom;     // number of atom in the molecule
	double *x;     // x coordinate
	double *y;     // y coordinate
	double *z;     // z coordinate
	double *vx;    // velocity in x direction
	double *vy;    // velocity in y direction
	double *vz;    // velocity in z direction
	double *Fx;    // force in x direction
	double *Fy;    // force in y direction
	double *Fz;    // force in z direction
	double *m;     // mass in kg
	int *Z;        // atomic number of the nuclei
};

// the entire system data structure
struct system_t{                
	   int nMolecule;           // number of molecule in the system
	   struct molecule_t *mols; // pointer to each molecule data
};

struct system_t *read_system(FILE *fd, int *totalAtom);
void print_molecule(struct molecule_t *ptr);
void outXYZ(struct molecule_t *ptr,FILE *fd2);
void outTXT(struct molecule_t *ptr,FILE *fd3);
void radial_distance(struct system_t *ptr,int nMolecule,FILE *fd3);
void OO_distance(struct molecule_t *ptr1,struct molecule_t *ptr2,FILE *fd3);
#endif
