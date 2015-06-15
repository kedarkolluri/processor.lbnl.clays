/*
 * ChangeAtoms.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef CHANGEATOMS_H_
#define CHANGEATOMS_H_
#include <global.h>
#include <Utilities.h>
#include <CopyAtom.h>


void insert_atoms(double (*arr)[4], int len);
void insert_atoms_specific(double (*arr)[5], int len);

void delete_atoms(int *arr, int len);
void delete_duplicates();
void delete_atom_single(int tag);

int change_atom_types_randomly(int type_now, int type_needed, double fraction_or_natoms, double min_dist=0);

void delete_atoms_randomly(int type, double fraction_or_natoms, double dist=0);

int change_bridging_oxygen(int root_type, int type_of_neighbor, int change_type_to);


#endif /* CHANGEATOMS_H_ */
