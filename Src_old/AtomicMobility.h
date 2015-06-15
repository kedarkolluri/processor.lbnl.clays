/*
 * AtomicMobility.h
 *
 *  Created on: Jul 22, 2009
 *      Author: kedar
 */

#ifndef ATOMICMOBILITY_H_
#define ATOMICMOBILITY_H_
#include <global.h>
#include <Utilities.h>
#include <ReactantProductMove.h>
void determine_mobility(atomic_dat *atom_ref,atomic_dat *atom_curr,int n,double Href[3][3],double H_curr[3][3],char *filename_save, int tag_reference, int tag_target, bool append);


#endif /* ATOMICMOBILITY_H_ */
