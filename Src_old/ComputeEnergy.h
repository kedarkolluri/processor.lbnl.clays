/*
 * ComputeEnergy.h
 *
 *  Created on: Aug 31, 2009
 *      Author: kedar
 */

#ifndef COMPUTEENERGY_H_
#define COMPUTEENERGY_H_
#include <global.h>
#include <Utilities.h>

double compute_energy_difference( atomic_dat *atom_ref1, atomic_dat *atom_now1, int n_now, double lower_limit, bool replace,int *number_of_atoms);
double compute_energy( atomic_dat *atom_now, int n_now);

#endif /* COMPUTEENERGY_H_ */
