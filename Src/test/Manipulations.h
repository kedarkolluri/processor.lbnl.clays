/*
 * Manipulations.h
 *
 *  Created on: Nov 18, 2009
 *      Author: kedar
 */

#ifndef MANIPULATIONS_H_
#define MANIPULATIONS_H_
#include <global.h>
#include <Utilities.h>
#include <ReadData.h>
#include <WriteData.h>
#include <ChangeAtoms.h>
#include <InitiateCNAandOthers.h>
#include <CNAandOthers.h>
int insert_atoms_recur(int coord_num_limit, atomic_dat *atom_now, int n_now, double H_now[3][3], bool consider, bool first_iter);
void shift(atomic_dat *atom_now,int n_now,double Hcry_now[3][3],double arr_shift[3]);
void save_atoms(atomic_dat *atom_now, int n_now, double H_now[3][3]);


#endif /* MANIPULATIONS_H_ */
