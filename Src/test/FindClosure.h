/*
 * FindClosure.h
 *
 *  Created on: Jul 30, 2009
 *      Author: kedar
 */

#ifndef FINDCLOSURE_H_
#define FINDCLOSURE_H_
#include <global.h>
#include <Utilities.h>

void FindClosure (atomic_dat *atom_ref, atomic_dat *atom_curr, int n, double Hcry_ref[3][3], double Hcry_curr[3][3], bool transform, char *filename, int number_of_atoms, double closure_vector[3]);


#endif /* FINDCLOSURE_H_ */
