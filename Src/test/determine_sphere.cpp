/*
 *  linkedlist_custom.cpp
 *
 *
 *  Created by Kedarnath Kolluri on 5/18/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "determine_sphere.h"

/*
   Recursive definition of determinate using expansion by minors.
*/
double Determinant(double **a,int n)
{
  int i,j,j1,j2;
  double det = 0;
  double **m = NULL;

  if (n < 1) { /* Error */

  } else if (n == 1) { /* Shouldn't get used */
    det = a[0][0];
  } else if (n == 2) {
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  } else {
    det = 0;
    for (j1=0;j1<n;j1++) {
      m = (double**) malloc((n-1)*sizeof(double *));
      for (i=0;i<n-1;i++)
	m[i] = (double*) malloc((n-1)*sizeof(double));
      for (i=1;i<n;i++) {
	j2 = 0;
	for (j=0;j<n;j++) {
	  if (j == j1)
	    continue;
	  m[i-1][j2] = a[i][j];
	  j2++;
	}
      }
      det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
      for (i=0;i<n-1;i++)
	free(m[i]);
      free(m);
    }
  }
  return(det);
}



int find_sphere(double x[4],double y[4],double z[4], double *x0,double *y0,double *z0,double *r)
{

  int i;
  double r1,theta,phi;
  XYZ c,p[4];
  double **a=NULL,m11,m12,m13,m14,m15;

 /* Create 4 random points on the surface */
  for (i=0;i<4;i++) {

    p[i].x = x[i] ;
    p[i].y = y[i] ;
    p[i].z = z[i] ;
  }

  /* Malloc the array for the minor arrays */
  a = (double**) malloc(4*sizeof(double *));
  for (i=0;i<4;i++)
    a[i] = (double*) malloc(4*sizeof(double));

  /* Find determinant M11 */
  for (i=0;i<4;i++) {
    a[i][0] = p[i].x;
    a[i][1] = p[i].y;
    a[i][2] = p[i].z;
    a[i][3] = 1;
  }
  m11 = Determinant(a,4);

  /* Find determinant M12 */
  for (i=0;i<4;i++) {
    a[i][0] = p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z;
    a[i][1] = p[i].y;
    a[i][2] = p[i].z;
    a[i][3] = 1;
  }
  m12 = Determinant(a,4);

  /* Find determinant M13 */
  for (i=0;i<4;i++) {
    a[i][0] = p[i].x;
    a[i][1] = p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z;
    a[i][2] = p[i].z;
    a[i][3] = 1;
  }
  m13 = Determinant(a,4);

  /* Find determinant M14 */
  for (i=0;i<4;i++) {
    a[i][0] = p[i].x;
    a[i][1] = p[i].y;
    a[i][2] = p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z;
    a[i][3] = 1;
  }
  m14 = Determinant(a,4);

  /* Find determinant M15 */
  for (i=0;i<4;i++) {
    a[i][0] = p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z;
    a[i][1] = p[i].x;
    a[i][2] = p[i].y;
    a[i][3] = p[i].z;
  }
  m15 = Determinant(a,4);

  if (fabs(m11) < 1e-5) {
   // fprintf(stderr,"The points don't define a sphere!\n");
    //fprintf(stderr,"Determinants: %g %g %g %g %g\n",m11,m12,m13,m14,m15);

    //    exit(-1);
    return(-1);
  }

  *x0 = 0.5 * m12 / m11;
  *y0 = 0.5 * m13 / m11;
  *z0 = 0.5 * m14 / m11;
  double alpha = *x0;double beta = *y0;double gamma=*z0;
  (*r) = sqrt((alpha)*(alpha) + (beta)*(beta) + (gamma)*(gamma) - m15/m11);
//  fprintf(stderr,"Sphere: (%g,%g,%g), %g\n",*x0,*y0,*z0,*r);

  for (i=0;i<4;i++)
    free(a[i]);
  free(a);
  return(0);
}









