/*
 * ReadData.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef READDATA_H_
#define READDATA_H_
#include <global.h>
#include <Utilities.h>

int read_cfg(char *filename, atomic_dat *atom_fill);
int read_A(void);
int read_A_exact(char *filename);
int read_xian_ming_blas(char *filename);
int read_xian_ming_blas2(char *filename);
int read_lammps(char *filename, atomic_dat *atom_fill, bool create_H, bool new_format);
int read_lammps_general(char *filename);



#endif /* READDATA_H_ */
