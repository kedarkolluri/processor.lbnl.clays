
/*
 *  linkedlist_custom.h
 *
 *
 *  Created by Kedarnath Kolluri on 5/18/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <time.h>

#define TWOPI           6.283185307179586476925287
#define PI              3.141592653589793238462643
#define PID2            1.570796326794896619231322

typedef struct {
  double x,y,z;
} XYZ;
double Determinant(double **,int);
int find_sphere(double x[4],double y[4],double z[4], double *x0,double *y0,double *z0,double *r);
