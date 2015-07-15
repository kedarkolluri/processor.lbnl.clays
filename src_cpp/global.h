/*
 * global.h
 *
 *  Created on: Jul 18, 2009
 *      Author: kedar
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_
//#include <mainprog.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <tuple>
//#include <typeinfo>
//#include <structure.h>
using namespace std;

#ifndef PI
#define PI    (3.14159265358979323846264338327950288419716939937510582097494)
#endif

#ifndef DISRIG_LIMIT
#define DISRIG_LIMIT 0.4
#endif

#ifndef STRING_LIMIT
#define STRING_LIMIT 512
#endif

#ifndef MAX_LOOP_SIZE
#define MAX_LOOP_SIZE   24
#endif
#ifndef ANGLE_DEVIATION
#define ANGLE_DEVIATION 0.955
#endif

extern int MAX_LOOP_VALUE;
extern FILE *fptr;

extern double ref_vectors[42][3];

extern double H0[3][3],H0_inv[3][3],Hcry[3][3],Hcry_inv[3][3],H[3][3],H_inv[3][3];
extern double crystal0[6];

extern double **BVS;
extern int VECTORS_ENUM;
extern int PLANES_ENUM;

enum{FCC,BCC,HCP,FCC_UNKNOWN,BCC_UNKNOWN,ICO,UNKNOWN, UNKNOWN_RING3,UNKNOWN_RING4,UNKNOWN_RING5,
	UNKNOWN_RING6,UNKNOWN_RING7,UNKNOWN_RING8,UNKNOWN_RING9,UNKNOWN_RING10, UNKNOWN_RING11,
	UNKNOWN_RING12,UNKNOWN_RING13,UNKNOWN_RING14,UNKNOWN_RING15};

extern int *nbr_ptr,*nbr_ptr1,*nbr_lst,*map1,*map_full,*head,*list,mx,my,mz,ncell;

extern double lx,ly,lz,xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz;

extern double rlst;
extern double rlstsq;
extern double **rcoordsq;

extern int format,target_format;
extern std::string inputfilename, outputfilename;

extern int interface_atoms;
extern int *tag_array_interface_atoms;
extern int total_number_of_files;

extern double *store_snapshots;
extern double *rms_displacement;
extern int filenumber_start,filenumber_end,filenumber_interval;


extern int n;
extern int n_types;
extern double **latt_cutoff;
extern double *element_weights;
extern double *element_charges;
//each array in the vector has 3 items (bond_type, atom_type_1, atom_type_2)
extern std::vector<int> element_bonds;

//each array in the vector has 3 items
//(bond_type, atom_type_1, atom_type_2, atom_type_3)
extern std::vector<int> element_angles;

#ifndef MAX_ELEMS
#define MAX_ELEMS    20
#endif
extern string *element_names;


#ifndef MAX_COORD
#define MAX_COORD    24
#endif

typedef struct atomic_dat_extra{
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

typedef struct atomic_dat {
	double rx,ry,rz,vx,vy,vz,fx,fy,fz,ke,pe,sx,sy,sz,ux,uy,uz,charge;
	int type;
	int molID;
	double ma;
	char elem[80];
	int coord;
	int coord_id[MAX_COORD];
	int neigh_bonds[MAX_COORD];
	atomic_dat_extra *localinfo;
	atomic_dat()
	{
		molID = 1;
	}
	//store a `dummy` to use for throw-away cases
	double dummy[3];
};

extern atomic_dat *atom;
extern std::vector<atomic_dat> atomvector;
//each set of 3 values (type, id1, id2) in the vector correspond to a bond
extern std::vector<int> bonds;
//each set of 4 values (type, id1, id2, id2) in the vector correspond to an angle
extern std::vector<int> angles;

extern bool TEST_BOOL;
extern int MAX_COORDNUM;

#endif /* GLOBAL_H_ */
