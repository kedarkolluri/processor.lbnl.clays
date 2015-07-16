/*
 * WriteData.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef WRITEDATA_H_
#define WRITEDATA_H_
#include <global.h>
#include <utilities.h>
extern bool LAMMPS_CHARGE;
extern bool LAMMPS_MOLECULE;
//extern char GULP_SHELL_ATOM[80];
//extern char GULP_SHELL_TYPE[80];
//extern double GULP_CORE_CHARGE;
//extern double GULP_SHELL_CHARGE;

int save_lammps(int pp,simcell &MDcell);

int save_xyz_VESTA(int pp, simcell MDcell);




#endif /* WRITEDATA_H_ */
