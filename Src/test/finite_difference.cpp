/*
 *  finite_difference.cpp
 *
 *
 *  Created by Kedarnath Kolluri on 6/7/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "finite_difference.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <fstream>
#include <math.h>

double M2inv(double **A,double **B)
{
	double determinant;
	B[0][0] = A[1][1];
	B[1][1] = A[0][0];
	B[0][1] = -1.0*A[0][1];
	B[1][0] = -1.0*A[1][0];
	determinant = A[0][0]*A[1][1]-A[0][1]*A[1][0];
	if(determinant==0)
	{
		return 0;
	}
	B[0][0] = B[0][0]/determinant;
	B[1][0] = B[1][0]/determinant;
	B[0][1] = B[0][1]/determinant;
	B[1][1] = B[1][1]/determinant;
	return determinant;
}

double M3inv (double **A, double **B)
{
    double determinant;
    B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1];
    B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2];
    B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0];
    B[1][0] = A[1][2]*A[2][0]-A[1][0]*A[2][2];
    B[2][1] = A[2][0]*A[0][1]-A[2][1]*A[0][0];
    B[0][2] = A[0][1]*A[1][2]-A[0][2]*A[1][1];
    B[2][0] = A[1][0]*A[2][1]-A[2][0]*A[1][1];
    B[0][1] = A[2][1]*A[0][2]-A[0][1]*A[2][2];
    B[1][2] = A[0][2]*A[1][0]-A[1][2]*A[0][0];
    determinant = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    if (determinant == 0.)
    {

		return 0;
    }
    B[0][0] /= determinant;
    B[1][1] /= determinant;
    B[2][2] /= determinant;
    B[1][0] /= determinant;
    B[2][1] /= determinant;
    B[0][2] /= determinant;
    B[2][0] /= determinant;
    B[0][1] /= determinant;
    B[1][2] /= determinant;
    return (determinant);
} /* end M3inv() */

int getderivative3(double **y,double **x, int rows, double **dydx) //assumes that the number of columns are 3
{
	if(rows<3) return -1;
	double x1[rows*3];
	double y1[rows*3];

//	std::cout<<"entered  dfdf \t"<<rows<<"\n";
/*
	for(int i=0;i<rows;i++)
	{

		std::cout<<x[i][0]<<"\t"<<x[i][1]<<"\t"<<x[i][2]<<" x \t";
		std::cout<<y[i][0]<<"\t"<<y[i][1]<<"\t"<<y[i][2]<<" y \n";
	}
*/
	int a=0;
	for(int i=0;i<rows;i++)
	{
		x1[i*3] = x[i][0];x1[i*3+1]=x[i][1];x1[i*3+2]=x[i][2];
		y1[i*3] = y[i][0];y1[i*3+1]=y[i][1];y1[i*3+2]=y[i][2];

	}

//	std::cout<<"entered t\n";
	double c[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	gsl_matrix_view X;
	gsl_matrix_view C;
	try{
		X = gsl_matrix_view_array(x1,rows,3);
		C = gsl_matrix_view_array(c,3,3);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&X.matrix,&X.matrix,0.0,&C.matrix);
	}
	catch (int e)
	{
		return 2;
	}
//	std::cout<<"entered t\n";
	double **c1,**c1inv;
	c1 = (double **) malloc(3*sizeof(double));
	c1inv = (double **) malloc(3*sizeof(double));
	for(int i=0;i<3;i++)
	{
		c1[i] = (double *) malloc(3*sizeof(double)); c1inv[i]= (double *) malloc(3*sizeof(double));
	}
	int k=0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
	//		std::cout << c[k]<<"\t";
			c1[i][j] = c[k];
			k++;
		}
	}

	double det = M3inv(c1,c1inv);
//	std::cout<<"entered t\t"<<det<<"\n";
	if(det<=0)
	{
		free(c1);
		free(c1inv);
		return 1;
	}
	k=0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			c[k] = c1inv[i][j];
			k++;
		}
	}
//	std::cout<<"entered t\n";
	double c2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double c3[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	try{
		gsl_matrix_view C1 = gsl_matrix_view_array(c,3,3);
		gsl_matrix_view C2 = gsl_matrix_view_array(c2,3,3);
		gsl_matrix_view C3 = gsl_matrix_view_array(c3,3,3);
		gsl_matrix_view Y = gsl_matrix_view_array(y1,rows,3);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&X.matrix,&Y.matrix,0.0,&C2.matrix);

		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&C1.matrix,&C2.matrix,0.0,&C3.matrix);
//		std::cout<<"entered t\n";
	}
	catch (int e)
	{
//		std::cout<<"entered inside catch t\n";
		free(c1);
		free(c1inv);
		return 2;
	}

	k=0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			dydx[i][j] = c3[k];
			k++;
		}
	}
//	std::cout <<"finished\n";
		free(c1);
		free(c1inv);
//			std::cout <<"finished a;ifgha;gh'adlfdfkljdffffffffffffffffffffffffffffffffffffffffffffffffffffff\n";

	return 0;
}



int getderivative2(double **y,double **x, int rows, double **dydx) //assumes that the number of columns are 3
{
	double x1[rows*2];
	double y1[rows*2];

//	std::cout<<"entered  dfdf \t"<<rows<<"\n";
/*
	for(int i=0;i<rows;i++)
	{

		std::cout<<x[i][0]<<"\t"<<x[i][1]<<"\t"<<x[i][2]<<" x \t";
		std::cout<<y[i][0]<<"\t"<<y[i][1]<<"\t"<<y[i][2]<<" y \n";
	}
*/
	int a=0;
	for(int i=0;i<rows;i++)
	{
		x1[i*2] = x[i][0];x1[i*2+1]=x[i][1];
		y1[i*2] = y[i][0];y1[i*2+1]=y[i][1];

	}

//	std::cout<<"entered t\n";
	double c[] = {0.0, 0.0, 0.0, 0.0};
	gsl_matrix_view X;
	gsl_matrix_view C;
	try{
		X = gsl_matrix_view_array(x1,rows,2);
		C = gsl_matrix_view_array(c,2,2);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&X.matrix,&X.matrix,0.0,&C.matrix);
	}
	catch (int e)
	{
		return 2;
	}
//	std::cout<<"entered t\n";
	double **c1,**c1inv;
	c1 = (double **) malloc(2*sizeof(double));
	c1inv = (double **) malloc(2*sizeof(double));
	for(int i=0;i<3;i++)
	{
		c1[i] = (double *) malloc(2*sizeof(double)); c1inv[i]= (double *) malloc(2*sizeof(double));
	}
	int k=0;
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
	//		std::cout << c[k]<<"\t";
			c1[i][j] = c[k];
			k++;
		}
	}

	double det = M2inv(c1,c1inv);
//	std::cout<<"entered t\t"<<det<<"\n";
	if(det<=0)
	{
		free(c1);
		free(c1inv);
		return 1;
	}
	k=0;
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			c[k] = c1inv[i][j];
			k++;
		}
	}
//	std::cout<<"entered t\n";
	double c2[] = {0.0, 0.0, 0.0, 0.0};
	double c3[] = {0.0, 0.0, 0.0, 0.0};
	try{
		gsl_matrix_view C1 = gsl_matrix_view_array(c,2,2);
		gsl_matrix_view C2 = gsl_matrix_view_array(c2,2,2);
		gsl_matrix_view C3 = gsl_matrix_view_array(c3,2,2);
		gsl_matrix_view Y = gsl_matrix_view_array(y1,rows,2);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&X.matrix,&Y.matrix,0.0,&C2.matrix);

		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&C1.matrix,&C2.matrix,0.0,&C3.matrix);
//		std::cout<<"entered t\n";
	}
	catch (int e)
	{
//		std::cout<<"entered inside catch t\n";s
		free(c1);
		free(c1inv);
		return 2;
	}

	k=0;
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			dydx[i][j] = c3[k];
			k++;
		}
	}
//	std::cout <<"finished a;ifgha;gh'adlfdfkljdffffffffffffffffffffffffffffffffffffffffffffffffffffff\n";
		free(c1);
		free(c1inv);
	return 0;
}

int getderivative(double **y,double **x, int rows, double **dydx, int dimension) //assumes that the number of columns are 3
{
	double x1[rows*dimension];
	double y1[rows*dimension];

//	std::cout<<"entered  dfdf \t"<<rows<<"\n";
/*
	for(int i=0;i<rows;i++)
	{

		std::cout<<x[i][0]<<"\t"<<x[i][1]<<"\t"<<x[i][2]<<" x \t";
		std::cout<<y[i][0]<<"\t"<<y[i][1]<<"\t"<<y[i][2]<<" y \n";
	}
*/
	int a=0;
	for(int i=0;i<rows;i++)
	{
		x1[i*dimension] = x[i][0];x1[i*dimension+1]=x[i][1];x1[i*dimension+2]=x[i][2];
		y1[i*dimension] = y[i][0];y1[i*dimension+1]=y[i][1];y1[i*dimension+2]=y[i][2];

	}

//	std::cout<<"entered t\n";
	double c[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	gsl_matrix_view X;
	gsl_matrix_view C;
	try{
		X = gsl_matrix_view_array(x1,rows,dimension);
		C = gsl_matrix_view_array(c,dimension,dimension);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&X.matrix,&X.matrix,0.0,&C.matrix);
	}
	catch (int e)
	{
		return 2;
	}
//	std::cout<<"entered t\n";
	double **c1,**c1inv;
	c1 = (double **) malloc(dimension*sizeof(double));
	c1inv = (double **) malloc(dimension*sizeof(double));
	for(int i=0;i<dimension;i++)
	{
		c1[i] = (double *) malloc(dimension*sizeof(double)); c1inv[i]= (double *) malloc(dimension*sizeof(double));
	}
	int k=0;
//	std::cout << "a'*a is*******\n";
	for(int i=0;i<dimension;i++)
	{
		for(int j=0;j<dimension;j++)
		{
//			std::cout << c[k]<<"\t";
			c1[i][j] = c[k];
			k++;
		}
	}
//	std::cout <<"\n*******\n";
	double det;
	if(dimension==3)
	{
		det= M3inv(c1,c1inv);

	}else if(dimension==2)
	{
		det =M2inv(c1,c1inv);
	}
//	std::cout<<"entered t\t"<<det<<"\n";
	if(det<=0)
	{
		free(c1);
		free(c1inv);
		return 1;
	}
	k=0;
//	std::cout << "c inverse is\n";
	for(int i=0;i<dimension;i++)
	{
		for(int j=0;j<dimension;j++)
		{
			c[k] = c1inv[i][j];
//			std::cout << c[k]<<"\t";

			k++;
		}
	}
//std::cout <<"\n*******\n";
//	std::cout<<"entered t\n";
	double c2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double c3[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	try{
		gsl_matrix_view C1 = gsl_matrix_view_array(c,dimension,dimension);
		gsl_matrix_view C2 = gsl_matrix_view_array(c2,dimension,dimension);
		gsl_matrix_view C3 = gsl_matrix_view_array(c3,dimension,dimension);
		gsl_matrix_view Y = gsl_matrix_view_array(y1,rows,dimension);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&X.matrix,&Y.matrix,0.0,&C2.matrix);

		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&C1.matrix,&C2.matrix,0.0,&C3.matrix);
//		std::cout<<"entered t\n";
	}
	catch (int e)
	{
//		std::cout<<"entered inside catch t\n";
		free(c1);
		free(c1inv);
		return 2;
	}

	k=0;
	for(int i=0;i<dimension;i++)
	{
		for(int j=0;j<dimension;j++)
		{
			dydx[i][j] = c3[k];
			k++;
		}
	}
//	std::cout <<"finished\n";
		free(c1);
		free(c1inv);
//			std::cout <<"finished a;ifgha;gh'adlfdfkljdffffffffffffffffffffffffffffffffffffffffffffffffffffff\n";

	return 0;
}



/*
int main (int argc, char *argv[])
{
	std::cout<<"hi\n";
	double **x;
	double **y;
	//(struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
	x = (double **) malloc(4*sizeof(double));
	y = (double **) malloc(4*sizeof(double));
	for(int i=0;i<4;i++)
	{
		x[i] = (double *) malloc(3*sizeof(double));
		y[i]= (double *) malloc(3*sizeof(double));
	}
	double **dydx;
	dydx = (double **) malloc(3*sizeof(double));
	for(int i=0;i<3;i++)
	{
		dydx[i] = (double *) malloc(3*sizeof(double));
	}

	x[0][0] = 1;x[0][1] = 1; x[0][2]= 1; x[1][0]= 1;x[1][1]= 2;x[1][2]= 3; x[2][0] =3;x[2][1]= 4;x[2][2]= 5; x[3][0] =2;x[3][1]= 7;x[3][2]= 8;
	y[0][0] = 9;y[0][1] = 10; y[0][2]= 11; y[1][0]= 22;y[1][1]= 22;y[1][2]= 22; y[2][0] =40;y[2][1]= 42;y[2][2]= 44; y[3][0] =63;y[3][1]=64;y[3][2]= 65;
	int ret_val = getderivative(y,x,4,dydx);

	if(ret_val==0)
	{
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			std::cout <<dydx[i][j]<<"\t";
		}
		std::cout <<"\n";
	}
	}else
	{
		std::cout<<ret_val<<"\n";

	}

}
*/
