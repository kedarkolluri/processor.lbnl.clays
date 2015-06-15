/*
 * CreateChainOfStates.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef CREATECHAINOFSTATES_H_
#define CREATECHAINOFSTATES_H_
#include <global.h>
#include <Utilities.h>
#include <WriteData.h>
#include <ReadData.h>
#include <CopyAtom.h>
int create_chain_of_states(char *filename_reactant, char *filename_product, int number_of_states, int start_chain_number, bool save_first, bool save_last);
int create_recurring_chain(int start_filenumber, int end_filenumber, int interval_filenumber, int begin, int chain_length);
int create_recurring_chain_equidistant(int start_filenumber, int end_filenumber, int interval_filenumber, int begin, int chain_length);
double compute_distance_scalar(atomic_dat *initial_dat, atomic_dat *final_dat, int n, double HCry_init[3][3],double Hcry_final[3][3]);

#endif /* CREATECHAINOFSTATES_H_ */
