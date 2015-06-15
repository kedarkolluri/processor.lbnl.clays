/*
 * NEB.h
 *
 *  Created on: Sep 14, 2009
 *      Author: kedar
 */

#ifndef NEB_H_
#define NEB_H_
#include <global.h>
#include <Utilities.h>
#include <ReadData.h>
#include <WriteData.h>
#include <CopyAtom.h>
#include <ComputeEnergy.h>
#include <OPToptions.h>

extern int NEB_SAVE_INFO;
extern bool NEB_CLIMB;
extern int NEB_CLIMB_AFTER;
extern int NEB_CLIMB_TOL;
extern int NEB_CLIMB_CONDITION;
extern int NEB_MAX_ITER;
extern double NEB_FORCE_TOL;
extern bool NEB_MODIFY_LINE_SEARCH;
extern bool NEB_SAVE_EACHSTEP;
extern double NEB_dt;
extern double NEB_springK;
extern int NEB_SLOT;
extern bool NEB_BACKTRACK_M;


int NEB(int *filename_array,int chain_length);


#endif /* NEB_H_ */
