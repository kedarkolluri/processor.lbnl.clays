/*
 * FindRings.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef FINDRINGS_H_
#define FINDRINGS_H_
#include <global.h>
#include <Utilities.h>

//bool is_distance_acceptable(atomic_dat *atom_now,int i,int ring_seq[MAX_LOOP_SIZE],int len, double H_now[3][3]);

//int complete_loop(int i, int j,int *counter, int ring_seq[MAX_LOOP_SIZE], atomic_dat *atom_now, double H_now[3][3]);
void compute_rings(atomic_dat *atom_now, int n_now, double H_now[3][3]);


#endif /* FINDRINGS_H_ */
