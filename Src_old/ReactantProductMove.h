/*
 * ReactantProductMove.h
 *
 *  Created on: Jul 20, 2009
 *      Author: kedar
 */

#ifndef REACTANTPRODUCTMOVE_H_
#define REACTANTPRODUCTMOVE_H_
#include <global.h>
#include <Utilities.h>
#include <ReadData.h>
#include <WriteData.h>
#include <CNAandOthers.h>
#include <CopyAtom.h>
#include "InitiateCNAandOthers.h"
extern bool ADDITIONAL_ATOMS;
extern bool RINGS_NO;
extern char ADDITIONAL_ATOMS_FILE[512] ;
int reactant_product_config_changes(char *filename_reactant, char *filename_product,int file_number, bool reverse, bool more, bool cfg_save, bool lammps_save, int neighs_to_stick);
//int move_partconfig(atomic_dat *atom_reactant,atomic_dat *atom_product, double H_reactant[3][3],double H_product[3][3],int n,bool reverse, bool more);
int prepare_for_chain( int start_filenumber, int end_filenumber, int interval_filenumber, int start_number, int *end_number, bool minimize_yes_no, char *script_title, bool cfg_save, bool lammps_save, int neighs_to_stick);


#endif /* REACTANTPRODUCTMOVE_H_ */
