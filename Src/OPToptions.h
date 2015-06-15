/*
 * OPToptions.h
 *
 *  Created on: Oct 6, 2009
 *      Author: kedar
 */

#ifndef OPTOPTIONS_H_
#define OPTOPTIONS_H_
#include <global.h>
#include <Utilities.h>
#include <ReadData.h>
#include <CopyAtom.h>
#include <ComputeEnergy.h>


extern double OPTOPTIONS_NEB_LIMIT;
enum OPTOPTIONS_NEB_STATIC_MODE_OPTIONS{
  ENERGY,
  DISPLACEMENT,
  NONE,

};
extern OPTOPTIONS_NEB_STATIC_MODE_OPTIONS OPTOPTIONS_NEB_STATIC_MODE;

int getdynamic(int *chain_array, int chain_length, int *dynamic_list, int *single_len, bool *bool_list);

#endif /* OPTOPTIONS_H_ */
