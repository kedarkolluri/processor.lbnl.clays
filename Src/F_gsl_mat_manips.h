//======================================================================================
// Name        : F_mat_manip.h
// Author      : Richard E. Baumer
// Version     : 5
// NOTES:	   : Final version prepared July 23, 2009; 
//				 July 29, 2009 - Added F_rdf_r_cr for computing the cutoff radius
//				 July 29, 2009 - Added F_density for computing system density
//				 July 30, 2009 - Added the histogram calculator
//				 August 31, 2009 - Added trajectory saver
// Copyright   : Copyright 2009 Richard E. Baumer. All rights Reserved
// Description : Matrix inversion, multiplication, reduced coordinates, and dot product.
//				 Vector magnitude, 1x3 vector cross product
// Date		   : August 31, 2009
//======================================================================================

#ifndef F_GSL_MAT_MANIPS_H
#define F_GSL_MAT_MANIPS_H

//Included libraries and files | Namespaced Declaration
//======================================================================================
#include <iostream>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_linalg.h>       // Header file for Linear Algebra
#include <gsl/gsl_matrix.h>       // Header file for Matrix definitions
#include <gsl/gsl_permutation.h>  // Header file for Permutation definitions
#include <gsl/gsl_blas.h>         // Header file for BLAS
#include <gsl/gsl_vector.h>       // Header for Vector definitions
#include <cmath>
using namespace std;			  // Namespaced declaration
//======================================================================================

//FUNCTION DECLARATION - Print GSL Matrix to Terminal
//======================================================================================
void F_gsl_mat_print(gsl_matrix* I,			/*GSL Matrix*/
					 int NUM_ROWS,			/*Number of Rows*/
					 int NUM_COLS,			/*Number of Columns*/
					 int Out_Prec			/*Number of digits to display for the dec*/
					 );

//FUNCTION DECLARATION - Write GSL Matrix to File
//======================================================================================
void F_gsl_mat_fsave(string file_name,		/*String filename*/
					 gsl_matrix* I,			/*GSL Matrix*/
					 int NUM_ROWS,			/*Num of Rows*/
					 int NUM_COLS,			/*Num of Cols*/
					 int Out_Prec			/*Number of digits to display for the dec*/
					 );

//FUNCTION DECLARATION - Inversion
//======================================================================================
gsl_matrix* gsl_mat_invrs(gsl_matrix *M		/* Input matrix; inverse form returned */
						  ); 

//FUNCTION DECLARATION - Matrix Multiplication
//======================================================================================
gsl_matrix* gsl_mat_mult(gsl_matrix *A,		/*Input Matrix*/
						 gsl_matrix *B,		/*Input Matrix*/
						 int A_ROWS,		/*Number of rows in A*/
						 int A_COLS,		/*Number of cols in A*/
						 int B_ROWS,		/*Number of rows in B*/
						 int B_COLS			/*Number of cols in B*/
						 ); 

//FUNCTION DECLARATION - Reduced Coordinate Matrix
//======================================================================================
gsl_matrix* reduced_coords_matrix(gsl_matrix *lattice,		/*Atomic Coords*/
								  gsl_matrix *system_cell,	/*System Cell Matrix*/
								  int NUM_ATOMS
								  ); 

void reduced_coords_vector(gsl_vector* ri,		/*Atomic Coords*/
						   gsl_matrix* sim_cell	/*System Cell Matrix*/
						   );

//FUNCTION DECLARATION - 3x1 vector-vector cross product, C = AxB
//======================================================================================
gsl_vector* gsl_vector3_cross_prod(gsl_vector *A,		/*Vector 1*/
								   gsl_vector *B		/*Vector 2*/
								   ); 

//FUNCTION DECLARATION - cutoff distance for RDF calculation
//======================================================================================
double F_rdf_r_cr(gsl_matrix* Sim_Cell		/*Simulation cell axes vector matrix*/
				  );

//FUNCTION - system volume
//======================================================================================
double F_volume(gsl_matrix* Sim_Cell		/*Simulation cell axes vector matrix*/
				);

//FUNCTION DECLARATION - system density
//======================================================================================
double F_density(gsl_matrix* Sim_Cell,		/*Simulation cell axes vector matrix*/
				 int NUM_ATOMS				/*Number of atoms in the system*/
				 );

//FUNCTION DECLARATION - remove zero (< or > 10^-16) values from ycol in xy (2 col) matrix
//======================================================================================
gsl_matrix* F_remove_y_0s(gsl_matrix *XY_Matrix,	/*Input two col matrix*/
						  int NUM_ROWS,				/*number of rows */
						  int &NUM_PRINT_ROWS		/*returned: number in print_mat*/
						  );

//FUNCTION DECLARATION - read a non-format, delimited, gls-matrix (NUM_ROWS x NUM_COLS)
//======================================================================================
gsl_matrix* F_read_gsl_file(string filename,	/*input filename*/
							int NUM_COLS,		/*Number of input columns*/
							int &NUM_ATOMS
							);

//FUNCTION DECLARATION - Calculate the histogram for a given data set
//======================================================================================
gsl_matrix* F_histogram(gsl_matrix* m,		/*Data matrix*/
						int COL,			/*Column from which to get the histogram*/
						int NUM_ROWS,		/*Number of rows in the matrix*/
						double bin,			/*bin size for histogram*/
						int& NUM_BINS	/*Number of print rows*/
						);

//FUNCTION DECLARATION - Save Trajectory File (0 through 5 order derivatives)
//======================================================================================
void F_save_trajectory(string filename,		/*String filename*/
					   gsl_matrix* r0,		/*GSL Matrix*/
					   gsl_matrix* r1,		/*GSL Matrix*/
					   gsl_matrix* r2,		/*GSL Matrix*/
					   gsl_matrix* r3,		/*GSL Matrix*/
					   gsl_matrix* r4,		/*GSL Matrix*/
					   gsl_matrix* r5,		/*GSL Matrix*/
					   gsl_matrix* sim_cell,	/*GSL Matrix*/
					   int NUM_ATOMS,			/*Num of Rows*/
					   int Out_Prec			/*Number of digits to display for the dec*/
					   );

//FUNCTION DECLARATION - Distance from atom I to J, magnitude
//======================================================================================
double F_dist(gsl_matrix* s0,
			  gsl_matrix* Sim_Cell,
			  int ATOM_I,
			  int ATOM_J
			  );

//FUNCTION DECLARATION - Wrap a coordinate to between (-0.5 and 0.5]
//======================================================================================
double F_reduced_round_pbc(double x,		/*	Any number between [0 1]	*/
						   int TOGGLE		/* TOGGLE = 0, use round(); TOGGLE = 1, (-0.5, 0.5] */
						   );

//FUNCTION - Inversion
//======================================================================================
// Function to compute the Inverse of a square 3x3 matrix using LU Decomposition.
// The matrices need to be initialized in the function that they are called from.
// M is the input matrix and M_inverse is the matrix that will hold the inverse.
//======================================================================================
gsl_matrix* gsl_mat_invrs(gsl_matrix *M) 
{
	gsl_matrix *M_inverse = gsl_matrix_alloc(3,3);
	
    int s; // Sign of the permutation for the LU decomposition
	
    //Defining the permutation vector
    gsl_permutation *p = gsl_permutation_alloc(3);
	
    gsl_linalg_LU_decomp(M, p, &s); // Performing the LU decomposition
	
    gsl_linalg_LU_invert(M, p, M_inverse); // Performing the matrix inversion
	
	return M_inverse;
}

//FUNCTION - Matrix Multiplication
//======================================================================================
// Function to multiply 2 matrices:
// gsl_blas_dgemm is called to execute the multiplication
// C = alpha*op(A)*op(B) + beta * C
// CblasNoTrans: The matrices to be multiplied should not be transposed
// alpha = 1.0
// beta = 0.0
// A and B are the input matrices and C is the product matrix
//======================================================================================
gsl_matrix* gsl_mat_mult(gsl_matrix *A,		/*Input Matrix*/
						 gsl_matrix *B,		/*Input Matrix*/
						 int A_ROWS,		/*Number of rows in A*/
						 int A_COLS,		/*Number of cols in A*/
						 int B_ROWS,		/*Number of rows in B*/
						 int B_COLS			/*Number of cols in B*/
						 ) 
{
	if (A_COLS != B_ROWS)
	{
		cout << "ERROR: Number of input columns (Matrix 1) does";
		cout << " not match number of input rows (Matrix 2)" << endl;
	}
			
    gsl_matrix* C = gsl_matrix_alloc(A_ROWS, B_COLS);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, C);
	return C;
}

//FUNCTION - Reduced Coordinate Matrix
//======================================================================================
//Function to compute the reduced coordinates of a lattice system.
// Reduced lattice =  Lattice matrix * Inverse of the System Cell
// The matrices should be initialized in the function they are called in.
// Arguments: lattice ---> Coordinates of atoms in the lattice
//            system_cell ---> Matrix holding the size of the system
//           reduced_lattice ---> Output matrix that hold the reduced 
//                                 coordinates of the atoms
//
//======================================================================================
gsl_matrix* reduced_coords_matrix(gsl_matrix *lattice,		/*Atomic Coords*/
								  gsl_matrix *system_cell,	/*System Cell Matrix*/
								  int NUM_ATOMS
								  ) 
{
	// Initializing the System cell inverse matrix
    gsl_matrix *system_cell_inverse = gsl_matrix_alloc(3, 3);
	gsl_matrix_memcpy(system_cell_inverse, system_cell);

    // Computing the inverse of the system cell
    system_cell_inverse = gsl_mat_invrs(system_cell_inverse);

    // Computing [reduced lattice] = [lattice] * [system cell inverse]
	gsl_matrix* reduced_lattice = gsl_mat_mult(lattice, system_cell_inverse, NUM_ATOMS, 3, 3, 3);
	
    // Deallocating the block [system_cell_inverse]
    gsl_matrix_free(system_cell_inverse);
	
	return reduced_lattice;
}


//FUNCTION - Reduced Coordinate Vector
//======================================================================================
void reduced_coords_vector(gsl_vector* ri,		/*Atomic Coords*/
						   gsl_matrix* sim_cell	/*System Cell Matrix*/
						   ) 
{
	// Initializing the System cell inverse matrix
    gsl_matrix *system_cell_inverse = gsl_matrix_alloc(3, 3);
	gsl_matrix_memcpy(system_cell_inverse, sim_cell);
	
    // Computing the inverse of the system cell
    system_cell_inverse = gsl_mat_invrs(system_cell_inverse);
	
    // Computing [reduced lattice] = [lattice] * [system cell inverse]
	gsl_vector* si = gsl_vector_alloc(3);
	gsl_blas_dgemv(CblasTrans, 1.0, system_cell_inverse, ri, 0, si);
	gsl_vector_memcpy(ri, si);
	
	
    // Deallocating the block [system_cell_inverse]
    gsl_matrix_free(system_cell_inverse);
	gsl_vector_free(si);
}


//FUNCTION - Write GSL Matrix to Terminal
//I: GSL Matrix containing the values to be displayed
//NUM_ROWS: number of rows of Matrix I (Or the number to be printed)
//NUM_COLS: number of columns in Matrix I (Or the number to be printed)
//Out_Prec: number of decimnal digits to be displayed
//======================================================================================
void F_gsl_mat_print(gsl_matrix* I,			/*GSL Matrix*/
					 int NUM_ROWS,			/*Number of Rows*/
					 int NUM_COLS,			/*Number of Columns*/
					 int Out_Prec			/*Number of digits to display for the dec*/
					 )
{
	
	for (int i = 0; i<NUM_ROWS; i++)
	{
		for (int j = 0; j<NUM_COLS; j++)
		{			
			cout << scientific << setprecision(Out_Prec) << gsl_matrix_get(I, i, j) << "  " ;
		}
		cout << endl;
	}
}
//======================================================================================

//FUNCTION - Save GSL Matrix to File
//filename: name of the output file (include file extension)
//I: GSL Matrix containing the values to be displayed
//NUM_ROWS: number of rows of Matrix I (Or the number to be printed)
//NUM_COLS: number of columns in Matrix I (Or the number to be printed)
//Out_Prec: number of decimnal digits to be displayed
//======================================================================================
void F_gsl_mat_fsave(string filename,		/*String filename*/
					 gsl_matrix* I,			/*GSL Matrix*/
					 int NUM_ROWS,			/*Num of Rows*/
					 int NUM_COLS,			/*Num of Cols*/
					 int Out_Prec			/*Number of digits to display for the dec*/
					 )
{
	//Open the structure file
	ofstream datafile;
	datafile.open(filename.c_str());
	
	//Write the real coordinate values
	for (int i=0; i< NUM_ROWS; i++)
	{
		for (int j=0; j<NUM_COLS; j++)
		{
			datafile << scientific << setprecision(Out_Prec) << "\t" << gsl_matrix_get(I, i, j);	
		}
		datafile << endl;
	}
	
	//Close the Datafile
	datafile.close();
}

//FUNCTION - 3x1 vector-vector cross product, C = AxB
//======================================================================================
gsl_vector* gsl_vector3_cross_prod(gsl_vector *A,		/*Vector 1*/
								   gsl_vector *B		/*Vector 2*/
								   )
{
	//Allocate the output vector
	gsl_vector* C = gsl_vector_alloc(3);
	
	//Perform the vector cross product: A[ax, ay, az] x B[bx, by, bz]
	//AxB = C = [(Ay*Bz-By*Az), -(Ax*Bz-Bx*Az), (Ax*By-Bx*Ay)]
	gsl_vector_set(C, 0, (gsl_vector_get(A, 1)*gsl_vector_get(B, 2) - gsl_vector_get(B, 1)*gsl_vector_get(A, 2) ));
	gsl_vector_set(C, 1, -(gsl_vector_get(A, 0)*gsl_vector_get(B, 2) - gsl_vector_get(B, 0)*gsl_vector_get(A, 2) ));
	gsl_vector_set(C, 2, (gsl_vector_get(A, 0)*gsl_vector_get(B, 1) - gsl_vector_get(B, 0)*gsl_vector_get(A, 1) ));
	
	//Return C
	return C;
	
}

//FUNCTION - cutoff distance for RDF calculation
//======================================================================================
double F_rdf_r_cr(gsl_matrix* Sim_Cell		/*Simulation cell axes vector matrix*/
				  )
{
	//Declare global variables
	double r_cr = 0;							/*Cutoff distance*/
	gsl_vector *a = gsl_vector_alloc(3);		/*axis 1*/
	gsl_vector *b = gsl_vector_alloc(3);		/*axis 2*/
	gsl_vector *c = gsl_vector_alloc(3);		/*axis 3*/
	double r1 = 0;
	double r2 = 0;
	double r3 = 0;
	gsl_vector *N1 = gsl_vector_alloc(3);		/*cross-product result 1*/
	gsl_vector *N2 = gsl_vector_alloc(3);		/*cross-product result 2*/
	gsl_vector *N3 = gsl_vector_alloc(3);		/*cross-product result 3*/
	
	//Read in values to the a, b, c vectors
	for (int i=0; i<3; i++)
	{
		gsl_vector_set(a, i, gsl_matrix_get(Sim_Cell, 0, i));
		gsl_vector_set(b, i, gsl_matrix_get(Sim_Cell, 1, i));
		gsl_vector_set(c, i, gsl_matrix_get(Sim_Cell, 2, i));
	}
	
	//Compute the radius values
	N1 = gsl_vector3_cross_prod(a,b);	/*Get cross product*/
	gsl_blas_ddot(c, N1, &r1);			/*Get dot product*/
	r1 = r1/(2.0*gsl_blas_dnrm2(N1));	/*Get projection, dividing by 2*/

	N2 = gsl_vector3_cross_prod(b,c);	/*Get cross product*/
	gsl_blas_ddot(a, N2, &r2);			/*Get dot product*/
	r2 = r2/(2.0*gsl_blas_dnrm2(N2));	/*Get projection, dividing by 2*/
	
	N3 = gsl_vector3_cross_prod(c,a);	/*Get cross product*/
	gsl_blas_ddot(b, N3, &r3);			/*Get dot product*/
	r3 = r3/(2.0*gsl_blas_dnrm2(N3));	/*Get projection, dividing by 2*/
	
	//Determine the minimum radius:
	r_cr = r1;
	
	if (r2 < r_cr)
	{
		r_cr = r2;
	}
	
	if (r3 < r_cr)
	{	
		r_cr = r3;
	}
	
	//Return the mimumum radius that was computed
	return r_cr;
}

//FUNCTION - system volume
//======================================================================================
double F_volume(gsl_matrix* Sim_Cell		/*Simulation cell axes vector matrix*/
				 )
{
	//Variable Declaration
	double vol = 0.0;
	gsl_vector *a = gsl_vector_alloc(3);		/*axis 1*/
	gsl_vector *b = gsl_vector_alloc(3);		/*axis 2*/
	gsl_vector *c = gsl_vector_alloc(3);		/*axis 3*/
	gsl_vector *N1 = gsl_vector_alloc(3);		/*axis 3*/
	
	//Read in values to the a, b, c vectors
	for (int i=0; i<3; i++)
	{
		gsl_vector_set(a, i, gsl_matrix_get(Sim_Cell, 0, i));
		gsl_vector_set(b, i, gsl_matrix_get(Sim_Cell, 1, i));
		gsl_vector_set(c, i, gsl_matrix_get(Sim_Cell, 2, i));
	}
	
	//Compute the density
	N1 = gsl_vector3_cross_prod(a,b);
	gsl_blas_ddot(N1, c, &vol);
	
	gsl_vector_free(N1);
	gsl_vector_free(a);
	gsl_vector_free(b);
	gsl_vector_free(c);
	
	//Exit and return calculated system density
	return vol;
}


//FUNCTION - system density
//======================================================================================
double F_density(gsl_matrix* Sim_Cell,		/*Simulation cell axes vector matrix*/
				 int NUM_ATOMS				/*Number of atoms in the system*/
				 )
{
	//Variable Declaration
	double density = 0.0;
	double vol = 0.0;
	
	//Get system volume
	vol = F_volume(Sim_Cell);
	
	density = double(NUM_ATOMS) / vol;

	//Exit and return calculated system density
	return density;
}

//FUNCTION - remove zero (< or > 10^-16) values from ycol in xy (2 col) matrix
//======================================================================================
gsl_matrix* F_remove_y_0s(gsl_matrix *XY_Matrix,	/*Input two col matrix*/
						  int NUM_ROWS,				/*number of rows */
						  int &NUM_PRINT_ROWS		/*returned: number in print_mat*/
						  )
{
	//Declare Variables
	int i = 0;
	int j = 0;
	
	for (i=0; i<NUM_ROWS; i++)
	{
		if (gsl_matrix_get(XY_Matrix,i,1) > 1e-14 || gsl_matrix_get(XY_Matrix, i, 1) < -1e-14)
		{
			j++;
		}
	}
	
	gsl_matrix *XY_Matrix_Plot = gsl_matrix_alloc(j, 2);
	
	j=0;
	for (i=0; i<NUM_ROWS; i++)
	{
		if (gsl_matrix_get(XY_Matrix,i,1) > 1e-14 || gsl_matrix_get(XY_Matrix, i, 1) < -1e-14)
		{
			gsl_matrix_set(XY_Matrix_Plot, j, 0, gsl_matrix_get(XY_Matrix, i, 0));
			gsl_matrix_set(XY_Matrix_Plot, j, 1, gsl_matrix_get(XY_Matrix, i, 1));
			j++;
		}
	}
	NUM_PRINT_ROWS = j;
	return	XY_Matrix_Plot;
}

//FUNCTION - read a non-format, delimited, gls-matrix (NUM_ROWS x NUM_COLS)
//======================================================================================
gsl_matrix* F_read_gsl_file(string filename,	/*input filename*/
							int NUM_COLS,		/*Number of input columns*/
							int &NUM_ATOMS
							)
{
	//Declare global variables
	ifstream datafile;
	int i = 0;
	int j = 0;
	gsl_matrix* data;
	int NUM_ROWS = 0;
	double temp_val = 0;
	
	//Determine the number of rows in the file
	datafile.open(filename.c_str());
	while (!datafile.fail())
	{
		//Read the first column value
		datafile >> temp_val;
		if (datafile.fail())
		{
			break;
		}
		
		for (j =1; j<NUM_COLS; j++)
		{
			
			datafile >> temp_val;
		}	
		NUM_ROWS++;
	}
	datafile.close();

	//Error check - make sure that file has data in it, else quit
	if (NUM_ROWS == 0)
	{
		cout << "ERROR: EMPTY FILE or FILE OPENING FAILURE" << endl;
		data = gsl_matrix_alloc(1, NUM_COLS);
		return data;
	}
	
	//Allocate data matrix
	data = gsl_matrix_alloc(NUM_ROWS, NUM_COLS);

	//Reopen the file and read in data
	datafile.open(filename.c_str());
	while (!datafile.fail())
	{
		//Read the first column value
		datafile >> temp_val;
		if (datafile.fail())
		{
			break;
		}
		gsl_matrix_set(data, i, 0, temp_val);
		
		for (j =1; j<NUM_COLS; j++)
		{
			
			datafile >> temp_val;
			gsl_matrix_set(data, i, j, temp_val);
		}	
		i++;
	}
	datafile.close();
	
	//Exit and return values
	NUM_ATOMS = NUM_ROWS;
	return data;
	
}

// Description : 
// Date		   : July 30, 2009


//FUNCTION - Calculate the histogram for a given data set
//======================================================================================
gsl_matrix* F_histogram(gsl_matrix* m,		/*Data matrix*/
						int COL,			/*Column from which to get the histogram*/
						int NUM_ROWS,		/*Number of rows in the matrix*/
						double bin,			/*bin size for histogram*/
						int& NUM_BINS		/*Number of print rows*/
						)
{
	//Declare global variables
	double min_val = 0.0;
	double max_val = 0.0;
	int bin_int = 0;
	int i = 0;
	
	//Determine the number of bins
	gsl_matrix_minmax(m, &min_val, &max_val);
	NUM_BINS = floor((max_val - min_val)/bin)+1;
	
	//Shift the data so that it starts at zero
	gsl_matrix* data_shifted = gsl_matrix_alloc(NUM_ROWS, 2);
	for (i=0; i<NUM_ROWS; i++)
	{
		gsl_matrix_set(data_shifted, i, 0, gsl_matrix_get(m, i, COL-1)-min_val);
	}
	
	//Allocate the histogram matrix w.r.t. m and bin data
	gsl_matrix* hist = gsl_matrix_alloc(NUM_BINS, 2);
	for (i = 0; i<NUM_BINS; i++)
	{
		gsl_matrix_set(hist, i, 0, i*bin+min_val);
	}
	
	for (i=0; i<NUM_ROWS; i++)
	{
		bin_int = floor( gsl_matrix_get(data_shifted, i, 0)/(bin) );
		gsl_matrix_set( hist, bin_int, 1, (gsl_matrix_get(hist, bin_int, 1)+2) );
	}
	
	gsl_matrix* hist_return = F_remove_y_0s(hist, NUM_BINS, NUM_BINS);
	
	//Exit and return histogram
	return hist_return;
	
}

//FUNCTION - Save Trajectory File (0 through 5 order derivatives)
//filename: name of the output file (include file extension)
//I: GSL Matrix containing the values to be displayed
//NUM_ROWS: number of rows of Matrix I (Or the number to be printed)
//NUM_COLS: number of columns in Matrix I (Or the number to be printed)
//Out_Prec: number of decimnal digits to be displayed
//======================================================================================
void F_save_trajectory(string filename,		/*String filename*/
					   gsl_matrix* r0,		/*GSL Matrix*/
					   gsl_matrix* r1,		/*GSL Matrix*/
					   gsl_matrix* r2,		/*GSL Matrix*/
					   gsl_matrix* r3,		/*GSL Matrix*/
					   gsl_matrix* r4,		/*GSL Matrix*/
					   gsl_matrix* r5,		/*GSL Matrix*/
					   gsl_matrix* sim_cell,	/*GSL Matrix*/
					   int NUM_ATOMS,			/*Num of Rows*/
					   int Out_Prec			/*Number of digits to display for the dec*/
					   )
{
	int i, j;
	int NUM_BLANK_LINES;
	
	//Open the structure file
	ofstream datafile;
	datafile.open(filename.c_str());
	
	//Write the first line
	for (i=0; i<10; i++)
	{
		datafile << " " << 4;
	}
	datafile << " " << NUM_ATOMS << endl;
	
	//Write blank lines
	NUM_BLANK_LINES = 8;
	for (i=0; i < NUM_BLANK_LINES; i++)
	{
		datafile << endl;
	}
	
	//Write the simulation cell matrix values
	for (i=0; i< 3; i++)
	{
		datafile << " " << gsl_matrix_get(sim_cell, i, 0);	
		datafile << " " << gsl_matrix_get(sim_cell, i, 1);	
		datafile << " " << gsl_matrix_get(sim_cell, i, 2);
		datafile << endl;
	}
	
	//Write blank lines
	NUM_BLANK_LINES = 6;
	for (i=0; i < NUM_BLANK_LINES; i++)
	{
		datafile << endl;
	}
	
	//Write the real coordinate values
	for (i=0; i< NUM_ATOMS; i++)
	{
		//Write r0
		for (j=0; j<3; j++)
		{
			datafile << scientific << setprecision(Out_Prec) << "\t" << gsl_matrix_get(r0, i, j);	
		}
		
		//Write r1
		for (j=0; j<3; j++)
		{
			datafile << scientific << setprecision(Out_Prec) << "\t" << gsl_matrix_get(r1, i, j);	
		}
		
		//Write r2
		for (j=0; j<3; j++)
		{
			datafile << scientific << setprecision(Out_Prec) << "\t" << gsl_matrix_get(r2, i, j);	
		}
		
		//Write r3
		for (j=0; j<3; j++)
		{
			datafile << scientific << setprecision(Out_Prec) << "\t" << gsl_matrix_get(r3, i, j);	
		}
		
		//Write r4
		for (j=0; j<3; j++)
		{
			datafile << scientific << setprecision(Out_Prec) << "\t" << gsl_matrix_get(r4, i, j);	
		}
		
		//Write r5
		for (j=0; j<3; j++)
		{
			datafile << scientific << setprecision(Out_Prec) << "\t" << gsl_matrix_get(r5, i, j);	
		}
		
		datafile << endl;
	}
	
	//Close the Datafile
	datafile.close();
}

//FUNCTION DECLARATION - Distance from atom I to J, magnitude
//======================================================================================
double F_dist(gsl_matrix* s0,
			  gsl_matrix* Sim_Cell,
			  int ATOM_I,
			  int ATOM_J
			  )
{
	double S_ij_x = 0.0;
	double S_ij_y = 0.0;
	double S_ij_z = 0.0;
	double R_ij_mag=0;
	gsl_vector *S_ij = gsl_vector_alloc(3);
	gsl_vector *R_ij = gsl_vector_alloc(3);	
	
	//Get Nearest Periodic Image in Reduced Lattice Space & store to vec
	S_ij_x =  gsl_matrix_get(s0, ATOM_J, 0) - gsl_matrix_get(s0, ATOM_I, 0);
	gsl_vector_set(S_ij, 0, ( S_ij_x - round(S_ij_x)) );
	
	S_ij_y =  gsl_matrix_get(s0, ATOM_J, 1) - gsl_matrix_get(s0, ATOM_I, 1);
	gsl_vector_set(S_ij, 1, ( S_ij_y - round(S_ij_y)) );
	
	S_ij_z =  gsl_matrix_get(s0, ATOM_J, 2) - gsl_matrix_get(s0, ATOM_I, 2);
	gsl_vector_set(S_ij, 2, ( S_ij_z - round(S_ij_z)) );
	
	//Get the Real space version of S_ij and find magnitude
	gsl_blas_dgemv(CblasNoTrans, 1, Sim_Cell, S_ij, 0, R_ij);	//S_ij -> R_ij
	R_ij_mag = gsl_blas_dnrm2(R_ij);							//Norm of R_ij	
	
	//Free Memory
	gsl_vector_free(S_ij);
	gsl_vector_free(R_ij);
	
	//Return distance and exit
	return R_ij_mag;
	
	
}

//FUNCTION - Wrap a coordinate to between (-0.5 and 0.5]
//======================================================================================
double F_reduced_round_pbc(double x,		/*	Any number between [0 1]	*/
						   int TOGGLE		/* TOGGLE = 0, use round(); TOGGLE = 1, (-0.5, 0.5] */
							)
{
	double epsilon = 1E-10;
	double intpart = 0.0;
	double fracpart = 0.0;
	
	if (TOGGLE == 0)	/*(-0.5, 0.5) */
	{
		return(round(x));
	}
	
	if (TOGGLE == 1)	/*(-0.5, 0.5] */
	{
		fracpart = modf(x, &intpart);
		if (x >= 0.0)
		{
			if (x - (intpart + 0.5) < epsilon)
			{
				return(intpart);
			}
			else
			{
				return(floor(x+0.5));
			}
		}
		if (x < 0.0)
		{
			if (x - (intpart + 0.5) < epsilon)
			{
				return(floor(x+0.5-epsilon));
			}
			else
			{
				return(floor(x+0.5));
			}
		}
	}
	
	if (TOGGLE == 2)	/*[-0.5, 0.5] */
	{
		fracpart = modf(x, &intpart);
		if (abs(x - (intpart + 0.5)) < epsilon)
		{
			return(intpart);
		}
		else
		{
			return(floor(x+0.5));
		}
	}
		
}


#endif
