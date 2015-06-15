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
void delete_atom_single(int tag);

#endif /* CHANGEATOMS_H_ */
