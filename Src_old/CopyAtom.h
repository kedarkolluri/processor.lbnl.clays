/*
 * CopyAtom.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef COPYATOM_H_
#define COPYATOM_H_
#include <global.h>
#include <Utilities.h>
void copy_atomstruct( atomic_dat *main_atom, atomic_dat *copy_atom, int atoms_number);
void copy_atom_indv( atomic_dat *main_atom, atomic_dat *copy_atom, int main_number,int copy_number);

void swap_atom( atomic_dat *main_atom, int main_number,int copy_number);
#endif /* COPYATOM_H_ */
