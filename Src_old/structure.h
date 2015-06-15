/*
 * structure.h
 *
 *  Created on: Jul 18, 2009
 *      Author: kedar
 */

#ifndef STRUCTURE_H_
#define STRUCTURE_H_

#include <fileport.h>
// Some minor tests but perhaps removed
#ifndef MAX_COORD
#define MAX_COORD    24
#endif

typedef struct atomic_dat {
	double rx,ry,rz,vx,vy,vz,fx,fy,fz,ke,pe,sx,sy,sz,ux,uy,uz;
	int type;
	double ma;
	int coord;
	int coord_id[MAX_COORD];
	int neigh_bonds[MAX_COORD];
	int interface;
	double disrigistry[3];
	int disrigistry_number_atoms;
	int disrigistry_number_atoms_negligible;
	int CNA;
	int ackN;
	int neigh_config;
	double drig;
	double delr[4];
	int BV;
};


#endif /* STRUCTURE_H_ */
