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
#include <map>
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

enum{FCC,BCC,HCP,FCC_UNKNOWN,BCC_UNKNOWN,ICO,UNKNOWN, UNKNOWN_RING3,UNKNOWN_RING4,UNKNOWN_RING5,
	UNKNOWN_RING6,UNKNOWN_RING7,UNKNOWN_RING8,UNKNOWN_RING9,UNKNOWN_RING10, UNKNOWN_RING11,
	UNKNOWN_RING12,UNKNOWN_RING13,UNKNOWN_RING14,UNKNOWN_RING15};

extern int *nbr_ptr,*nbr_ptr1,*nbr_lst,*map1,*map_full,*head,*list,mx,my,mz,ncell;

extern int format,target_format;

extern std::string inputfilename, outputfilename;

extern int total_number_of_files;


extern int filenumber_start,filenumber_end,filenumber_interval;

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
	double r[3], v[3], f[3], s[3], u[3], charge, ke, pe;
	int type;
	int molID;
	atomic_dat()
	{
		molID = 1;
	}
	double ma; //this is actually redundant but helps if there are errors
	char elem[80]; //this is redundant but helps if there are errors
	int coord;
	int coord_id[MAX_COORD];
	int neigh_bonds[MAX_COORD];
	std::vector<atomic_dat_extra> localinfo;

 //store a vector called `data` that contains extra data
 std::vector<double> data;
 //store a vector of strings that contains the names of the data
 std::vector<string> datanames;
};

typedef struct atom_types_info {
	int n_types;
	double rlst = 4.0; //defaults to 4.0
	double rlstsq = rlst*rlst;
	std::map<std::tuple<int, int>, double> rcoordsq;
	// we are using maps here so that element type numbers need not be consequtive
	std::map< int, double> latt_cutoff;
	std::map<int, string> element_names;
	std::map<int, double> element_weights;
	std::map<int, double> element_charges;
	//each array in the vector has 3 items (bond_type, atom_type_1, atom_type_2)
	std::vector< std::tuple<int, int, int> > element_bonds;
	//(angle_type, atom_type_1, atom_type_2, atom_type_3)
	std::vector< std::tuple<int, int, int, int> > element_angles;
};

typedef struct simcell{
	//simulation cell dimensions
	double Hcry[3][3], Hinv[3][3], crystal[6];
	double lx,ly,lz,xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz;
	//number of atoms in the simulation cell
	int n;
	atom_types_info atoms_info;
	//we are using map here so that atom ids need not be consequtive
	std::map<int, atomic_dat> atomdata;
	// bonds in the system - a vector of tuples
	std::vector<std::tuple<int, int, int> > bonds;
	//angles in the system - a vector of tuples
	std::vector<std::tuple<int, int, int, int> > angles;

	//arbitrary containers to take all other system parameters
	std::vector<double> systemprops;
	std::vector<string> systemprops_names;
};
extern bool TEST_BOOL;
extern int MAX_COORDNUM;

#endif /* GLOBAL_H_ */
