/*
 * ReadData.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef READDATA_H_
#define READDATA_H_
#include <global.h>
#include <utilities.h>

simcell read_lammps(char *filename);
simcell read_xyz_VESTA(const char *filename);



#endif /* READDATA_H_ */
