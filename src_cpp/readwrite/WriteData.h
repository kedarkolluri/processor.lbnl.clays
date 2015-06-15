/*
 * WriteData.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef WRITEDATA_H_
#define WRITEDATA_H_
#include <global.h>
#include <Utilities.h>
extern bool LAMMPS_CHARGE;
extern bool LAMMPS_MOLECULE;
extern char GULP_SHELL_ATOM[80];
extern char GULP_SHELL_TYPE[80];
extern double GULP_CORE_CHARGE;
extern double GULP_SHELL_CHARGE;
int save_cfg_interface(int pp, double H[3][3], int interface,int atom_type);

int save_cfg(int pp, double H[3][3]);
int save_cfg_test(int pp, double H[3][3]);
int save_lammps(int pp,double H[3][3]);
int save_gulp(int pp, double H[3][3]);
int save_lammps_specific(char *prefix,int pp,atomic_dat *atom_, double H[3][3]);
int save_cfg_cut(int pp, double H[3][3], double x_low,double x_high, double y_low, double y_high, double z_low, double z_high, char *suffix);

int save_xyz_VESTA(int pp, double H[3][3]);

int write_A(const char *, int, atomic_dat *, double H_here[3][3]);



#endif /* WRITEDATA_H_ */
