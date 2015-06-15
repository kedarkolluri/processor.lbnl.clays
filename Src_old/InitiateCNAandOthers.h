/*
 * InitiateCNAandOthers.h
 *
 *  Created on: Jul 20, 2009
 *      Author: kedar
 */

#ifndef INITIATECNAANDOTHERS_H_
#define INITIATECNAANDOTHERS_H_
#include <global.h>
#include <Utilities.h>
#include <CNAandOthers.h>
#include <NbrList.h>
#include <FindRings.h>

extern double CNA_RATIO_1;
extern double CNA_RATIO_2;
extern double CNA_RATIO_12;

int prepare_nbrlist(double H1[3][3], int max_number_atoms);
void set_neighbordistances(char *filename);

void compute_CNA_and_others( atomic_dat *atom_now, int n_now, double H_now[3][3]);
void compute_CNA_and_others( atomic_dat *atom_now, int n_now, double H_now[3][3], bool extra);



#endif /* INITIATECNAANDOTHERS_H_ */
