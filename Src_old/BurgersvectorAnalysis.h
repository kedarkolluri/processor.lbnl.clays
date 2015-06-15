/*
 * BurgersvectorAnalysis.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef BURGERSVECTORANALYSIS_H_
#define BURGERSVECTORANALYSIS_H_
#include <global.h>
#include <Utilities.h>

double find_BV_plane(double *bv, int *plane, int *bvl);

void determine_BV_plane(double rx,double ry,double rz, int *plane, double *angle);

void determine_BV_plane(int index, int num_slip,double rx,double ry,double rz);

#endif /* BURGERSVECTORANALYSIS_H_ */
