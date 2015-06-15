/*
 * global.cpp
 *
 *  Created on: Jul 18, 2009
 *      Author: kedar
 */

#include <global.h>


//void set_global_variables()
//{
	int MAX_LOOP_VALUE;
	double a_Cu = 3.615, a_Nb = 3.3008;

	double Cu_alpha_matrix[3][3];
	double Nb_alpha_matrix[3][3];
	double KS1_trans_mas[3][3];

	double ref_vectors[42][3];
	int basic_n[4];

	int count_disrig_dev=0;

	double H0[3][3],H0_inv[3][3],Hcry[3][3],Hcry_inv[3][3],H[3][3],H_inv[3][3];
	double crystal0[6];

	double H0_geo[3][3], H0_geo_inv[3][3];

	double **BVS;
	int VECTORS_ENUM;
	int PLANES_ENUM;


	int *nbr_ptr=NULL,*nbr_ptr1=NULL,*nbr_lst=NULL,*map1=NULL,*map_full=NULL,*head=NULL,*list=NULL,mx,my,mz,ncell;

	double lx,ly,lz,xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz;

	double H_Cu[3][3],H_Nb[3][3],H_total[3][3];
	double n_Cu, n_Nb, lx_Cu, ly_Cu, lz_Cu, lx_Nb, ly_Nb, lz_Nb;

	double rlst = 4.0; //hard coded cut off value of rlist -- need to change according to need
	double rlstsq = rlst*rlst;
	//double z_start=0.548,z_end = 0.624; // for KSmin

	//double z_start = 0.422, z_end = 0.484; //KS1

	double z_start = 0.42, z_end = 0.49; //KS2

	double Cu_massf = 0.406168318;

	double **rcoordsq;

	FILE *fptr;
	int format,target_format;
	//char *inputfilename, *outputfilename;
	std::string inputfilename, outputfilename;

	int interface_atoms;
	int *tag_array_interface_atoms;
	int total_number_of_files=0;

	double *store_snapshots;
	double *rms_displacement;
	int filenumber_start,filenumber_end,filenumber_interval;

	std::vector<int*> bonds;
	std::vector<int*> angles;


	int interface_Cu_n;
	int n;
	int n_types =0;;
	double **latt_cutoff;
	double *element_charges;
	double *element_weights;
	string *element_names;
	std::vector<int> element_bonds;
	std::vector<int> element_angles;


	atomic_dat *atom;
	bool TEST_BOOL;
	int MAX_COORDNUM;

//}
