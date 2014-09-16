#ifndef SIM_H
#define SIM_H

void set_displacement(struct molecule_t *ptr,int R0);
void set_velocity(struct molecule_t *ptr,int V0);
void updateNewton_velocity(
	struct molecule_t *ptr,           // pointer to molecule 
	struct option_t *opt);            // pointer to option
void updateNewton_position(
	struct molecule_t *ptr,           // pointer to molecule
	struct option_t *opt);            // pointer to option
void control_velocity(
	int totalAtom,                    // total atom
	struct system_t *sys,             // pointer to system
    struct molecule_t *ptr,           // pointer to molecule
	double tem);                      // temperature
void get_kinetic(struct molecule_t *ptr, // pointer to molecule
	double *Ek);                         // pointer to kinetic energy
#endif
