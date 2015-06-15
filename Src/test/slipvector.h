/*
 * slipvector.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef SLIPVECTOR_H_
#define SLIPVECTOR_H_
#include <global.h>
#include <Utilities.h>
#include <BurgersvectorAnalysis.h>
//void init_geometry();
void init_geometry(char *STRING1, char *STRING2);
void compute_slipvector(atomic_dat *, atomic_dat *, double (*) [3],double (*)[3],bool);
void compute_inplaneslipvector(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform);
void compute_absolute_displacement(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform);

void compute_absolute_displacement(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform,int interface_type, int atom_type);
void compute_absolute_displacement_interface(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3]);




#endif /* SLIPVECTOR_H_ */
