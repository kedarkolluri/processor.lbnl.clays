/*
 * Utilities.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */
#include "global.h"

#ifndef UTILITIES_H_
#define UTILITIES_H_


extern bool UTILITIES_ZIP;

int min(int a, int b);


double min(double a, double b);

int max(int a, int b);


double max(double a, double b);

bool check_repeat(int j,int *arr,int len);

void zip (char *filename);


void execute_system_command( char *command);

void unzip (char *filename);

//sort ascending order
void sort (int count,int *data);


void V3norm(double *B);

double V3dot (double *a, double *b);

double V3mulM3 (double a[3], double B[3][3], double c[3]);

void V3cross (double a[3], double b[3], double c[3]);

double M3mulV3(double A[3][3], double b[3], double c[3]);

double M3inv (double A[3][3], double B[3][3]);

void M3mul (double A[3][3], double B[3][3], double C[3][3]);
void M2mul (double A[2][2], double B[2][2], double C[2][2]);
void M2nms (double A[2][2], double C[2][2], double *trace, double *vmises);

void H_crystal(double H[3][3],double crystal[6]);
void crystal_H(double crystal[6],double H[3][3]);

bool check_file_exists(char *filename);
string* split_string(string tmp_str, int &num_sub_strings);
string* get_next_splits( ifstream &inputfile_h, int &numval);

#endif /* UTILITIES_H_ */
