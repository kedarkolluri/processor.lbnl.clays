/*
 * CNAandOthers.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef CNAANDOTHERS_H_
#define CNAANDOTHERS_H_
#include <global.h>
#include <Utilities.h>

void compute_acklandnotation( atomic_dat *atom_now, int n_now, double H_now[3][3]);

//int common_neighbors(int i,int j, int atom_type,int *return_value,double len_z);
void compute_CNA( atomic_dat *atom_now, int n_now, double H_now[3][3]);

void coord_number(double H1[3][3]);




#endif /* CNAANDOTHERS_H_ */
