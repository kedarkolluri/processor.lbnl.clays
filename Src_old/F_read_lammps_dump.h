//	=================================================================	//
// Name        : F_read_lammps_dump.h
// Author      : Richard E. Baumer
// Copyright   : Copyright 2009 Richard E. Baumer. All rights Reserved
// Description : System command functions
// Date		   : September 17, 2010
// Date		   : Uknown; Modifications by Kedar
//	=================================================================	//

#ifndef F_READ_LAMMPS_DUMP_H
#define F_READ_LAMMPS_DUMP_H

#include <string>
#include "gzstream.h"
#include <fstream>
#include <gsl/gsl_matrix.h>

//	AAAT Include Files
#include "global.h"
#include "F_system_funcs.h"
//#include "F_gsl_mat_manips.h"

using namespace std;

//	Function Declaration
void F_read_lammps_dump();

//	Function Definition
void F_read_lammps_dump(char *filename)
{
	clock_t start_time;	//	Gets the starting time of the program

	//	Declare local variables
	//	=============================================================	//
    igzstream datafile;						//	Datastream in the gz format
	//ifstream datafile;						//	Data input stream
	int i = 0, j = 0;						//	Loop variables
	int ct = 0;
	int ATOM_INDEX = 0, ATOM_TYPE = 0;		//	Holds Atom index and atomic type (both integer values)
	double temp_val = 0.0;					//	Temporary double values
	double shift_axis_val[3] = {0.0, 0.0, 0.0};	//	Amount to shift SIM_CELL vectors so origin is at {0 0 0}
	string temp_line;						//	String container to hold line
	string* ptr_temp_str;					//	Temp pointer to string
	int ID_COL_NUM = 0, TYPE_COL_NUM = 0;	//	ID and TYPE column indices
	int X_COL_NUM = 0, Y_COL_NUM = 0, Z_COL_NUM = 0;	//	Index corresponding to the x y and z data
	//	=============================================================	//

	//	If FLAG_DUMP == 1, deallocate previously allocated memory
	//	Always clear number of atoms
	//	=============================================================	//
	if (FLAG_DUMP == 1) {
		gsl_matrix_free(R0);
		gsl_matrix_free(SIM_CELL);
		gsl_matrix_free(cfg_props);
		delete [] PTR_TYPE;
		delete [] ptr_prop_str; }
	n = 0;
//	R_CR_LINKED = 0.0;		//	Bin size of the current linked list
//	R_CR_VERLET = 0.0;		//	Cutoff radius of the current verlet neighbor list
	//	=============================================================	//

/*
	//	Unzip the input file
	//	=============================================================	//
	size_t found = INPUT_FILENAME.find(".gz");
	if (found!=string::npos) {
		unzip (INPUT_FILENAME);
		INPUT_FILENAME = INPUT_FILENAME.substr(0, int(found));}
	//	=============================================================	//
*/

	//Open the datafile
// modified by Kedar: input is not string but  *char
//	datafile.open(INPUT_FILENAME.c_str());
	datafile.open(filename);

	//Read the time step
	getline(datafile, temp_line);
	datafile >> TIMESTEP;
	cout << "made it here " << TIMESTEP << endl;

	//Read the number of atoms
	getline(datafile, temp_line);
	getline(datafile, temp_line);
	datafile >> NUM_ATOMS;
	cout << "made it here " << NUM_ATOMS << endl;

	//	Allocate memory for the datafile based on number of atoms
	//	=============================================================	//
	R0 = gsl_matrix_alloc(NUM_ATOMS, 3);
	SIM_CELL = gsl_matrix_alloc(3,3);
	gsl_matrix_set_identity(SIM_CELL);
	PTR_TYPE = new int[NUM_ATOMS];
	//	=============================================================	//

	//Read simulation cell vectors; determine system density
	//	=============================================================	//
	getline(datafile, temp_line);
	getline(datafile, temp_line);
	for (i=0; i<3; i++)	{
		datafile >> shift_axis_val[i];
		datafile >> temp_val;	
		gsl_matrix_set(SIM_CELL, i, i, temp_val-shift_axis_val[i]);}
	DENSITY = F_density(SIM_CELL, NUM_ATOMS);
	//	=============================================================	//
	getline(datafile, temp_line);

	//Read the atomic type, coordinates, velocities, and auxiliary properties
	//	=============================================================	//

	//	Read the line with property names; split into string
	getline(datafile, temp_line);
	ptr_temp_str = F_split_string(temp_line, NUM_PROPS);

	if (ptr_temp_str[0] == "ITEM:" && ptr_temp_str[1] == "ATOMS")
	{
		//	Since the data appears to be formatted correctly, allocate the property string
		NUM_PROPS=NUM_PROPS-2;
		ptr_prop_str = new string[NUM_PROPS];

		//	Store the properties into the ptr_prop_str array
		for (i = 2; i<NUM_PROPS+2; i++)
		{
			//	Get the property string array value
			ptr_prop_str[i-2] = ptr_temp_str[i];

			//	Check for the id, type, and {x,y,z} column numbers
			if (ptr_prop_str[i-2] == "id")
				ID_COL_NUM = i-2;
			if (ptr_prop_str[i-2] == "type")
				TYPE_COL_NUM = i-2;
			if (ptr_prop_str[i-2] == "x" || ptr_prop_str[i-2] == "xu" || ptr_prop_str[i-2] == "sx")
				X_COL_NUM = i-2;
			else if (ptr_prop_str[i-2] == "y" || ptr_prop_str[i-2] == "yu" || ptr_prop_str[i-2] == "sy")
				Y_COL_NUM = i-2;
			else if (ptr_prop_str[i-2] == "z" || ptr_prop_str[i-2] == "zu" || ptr_prop_str[i-2] == "sz")
				Z_COL_NUM = i-2;
		}
	}
	else {
		F_print_errors(4);}

	cout << NUM_PROPS << endl;
	
	//	Allocate the configuration matrix: The number of auxiliary properties is NUM_PROPS - 5;
	if (NUM_PROPS-5 > 0) {
		cfg_props = gsl_matrix_alloc(NUM_ATOMS, NUM_PROPS-5);}

	for (i=0; i<NUM_ATOMS; i++)
	{
		ct=0;
		for (j=0; j<NUM_PROPS; j++)
		{
			//	Get the atomic index
			if (j == ID_COL_NUM){
				datafile >> ATOM_INDEX;}

			//	Get the atomic type
			else if (j == TYPE_COL_NUM) {
				datafile >> ATOM_TYPE;
				PTR_TYPE[ATOM_INDEX-1] = ATOM_TYPE;
				if (ATOM_TYPE == 1)
					NUM_A_ATOMS++;
				else if (ATOM_TYPE == 2)
					NUM_B_ATOMS++;
			}

			//	Get the x-coordinate
			else if (j == X_COL_NUM) {
				datafile >> temp_val;
				gsl_matrix_set(R0, ATOM_INDEX-1, 0, temp_val);}

			//	Get the y-coordinate
			else if (j == Y_COL_NUM) {
				datafile >> temp_val;
				gsl_matrix_set(R0, ATOM_INDEX-1, 1, temp_val);}

			//	Get the z-coordinate
			else if (j == Z_COL_NUM) {
				datafile >> temp_val;
				gsl_matrix_set(R0, ATOM_INDEX-1, 2, temp_val);}

			//	Get the atomic properties
			else if (j != ID_COL_NUM && j != TYPE_COL_NUM && j != X_COL_NUM && j != Y_COL_NUM && j != Z_COL_NUM){
				datafile >> temp_val;
				gsl_matrix_set(cfg_props, ATOM_INDEX-1, ct, temp_val);
				ct++;}
		}
	}
	//	=============================================================	//

	//	Create the property string array
	//	=============================================================	//

	//	Clear pointer to the ptr_temp_str; reallocate with NUM_PROPS length; and
	//	give it the values of ptr_prop_str
	delete [] ptr_temp_str;
	ptr_temp_str = new string[NUM_PROPS];
	for (j = 0; j<NUM_PROPS; j++){
		ptr_temp_str[j]=ptr_prop_str[j];}

	//	Clear pointer to the ptr_prop_str; reallocate with NUM_PROPS-5, which is the true number of auxiliary properties (Total - id - type - x -y -z)
	delete [] ptr_prop_str;
	ptr_prop_str = new string[NUM_PROPS-5];

	//	 Assign values to the ptr_prop_str
	ct=0;
	for (j = 0; j<NUM_PROPS; j++)
	{
		if (j != ID_COL_NUM && j != TYPE_COL_NUM && j != X_COL_NUM && j != Y_COL_NUM && j != Z_COL_NUM)
		{
			ptr_prop_str[ct]=ptr_temp_str[j];
			ct++;
		}
	}
	delete [] ptr_temp_str;
	NUM_PROPS = ct;
	//	=============================================================	//

	//	Print a summary of the read operation
	//	=============================================================	//
	cout << "---------------------------------------------------------------" << endl;
	cout << "read_dump summary:" << endl;
	cout << NUM_ATOMS << " atoms read from input file: " << INPUT_FILENAME << endl;
	cout << NUM_A_ATOMS << " A_type atoms, " << NUM_B_ATOMS << " B_type atoms" << endl;
	cout << "atomic type is stored in PTR_TYPE array" << endl;
	cout << "atomic coordinates are stored in R0 array" << endl;
	cout << "From data file, the auxiliary properties are stored in cfg_prop matrix as follows:" << endl;
	for (i = 0; i<NUM_PROPS; i++){
		cout << i << " " << ptr_prop_str[i] << endl;}
	cout << "---------------------------------------------------------------" << endl;
	//	=============================================================	//

	//	Zip the dump file
//	zip(INPUT_FILENAME);
//	cout << "Zip time: " << (clock() - start_time)/double(CLOCKS_PER_SEC) << " secs" << endl;

	//Close the datafile	
	datafile.close();

	//	Set FLAG_DUMP to true
	FLAG_DUMP = 1;
	FLAG_S0 = 0;
//	FLAG_VERLET_LIST = 0;
//	FLAG_LINKED_LIST = 0;
//	FLAG_PBC = 0;
}
#endif
