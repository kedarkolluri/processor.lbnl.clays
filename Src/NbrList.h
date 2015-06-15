/*
 * NbrList.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef NBRLIST_H_
#define NBRLIST_H_
#include <global.h>
#include <Utilities.h>


/* --------------------------------------------------------------------- */


int cellno(int ix, int iy, int iz);


void mapcells();


void make_nbr_lst(double H1[3][3]);



#endif /* NBRLIST_H_ */
