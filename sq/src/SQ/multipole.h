#ifndef MULTIPOLE_H
#define MULTIPOLE_H

#include "basis.h"

#define MULTIPOLE_CHARGE_CUTOFF 1.0E-18
#define MULTIPOLE_RADII2_CUTOFF 1.0E-00

// multipole expansion of a product between two basis functions
struct Multipole_t{
	double q;                    // total charge
	double x,y,z;                // center
};

struct Multipole_t *genMultipole(
	int nBasis,                     // the number of basis functions
	const struct GTOBasis_t *gto);  // pointer to basis function information

double coulombMultipole(
	int nBasis,                       // number of basis function
	int p, int q,                     // index to first multipole
	int i, int j,                     // index to second multipole
	const struct Multipole_t *mpole); // pointer to multipole info

#endif // MULTIPOLE_H
