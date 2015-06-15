/*
 * Disrigistry.h
 *
 *  Created on: Jun 6, 2011
 *      Author: kedar
 */

#ifndef DISRIGISTRY_H_
#define DISRIGISTRY_H_
#include <global.h>
#include <Utilities.h>
#include <BurgersvectorAnalysis.h>
#include<finite_difference.h>

//#include "finite_difference.h"
#include "slipvector.h"
#include "NbrList.h"
#include "CNAandOthers.h"
#include "CopyAtom.h"
#include "WriteData.h"
#include "ReadData.h"
#include "InitiateCNAandOthers.h"
extern bool DISRIGISTRY_A;
void compute_disrigistry(atomic_dat *atom_now, int n_now,double H_now[3][3]);
void compute_disrigistry(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform);
void perform_disrigistry(char *firstfile,char *secondfile,int dis_type, int perform_mode,int s_atom_type,int s_atom_interface);


#endif /* DISRIGISTRY_H_ */
