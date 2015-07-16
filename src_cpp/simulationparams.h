/*
 * global.h
 *
 *  Created on: Jul 18, 2009
 *      Author: kedar
 */

#ifndef SIMULATIONPARAMS_H_
#define SIMULATIONPARAMS_H_
#include <global.h>
#include <utilities.h>

atom_types_info set_atom_types_info(const char *filename);
void set_types_to_atom_from_element_name(simcell &MDcell);

#endif /* SIMULATIONPARAMS_H_ */
