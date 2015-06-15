/*
 * MakeSystem.h
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#ifndef MAKESYSTEM_H_
#define MAKESYSTEM_H_
#include <global.h>
#include <Utilities.h>
#include <ReadData.h>
#include <WriteData.h>
void create_system(char *Cu_filename, char *Nb_filename, bool keep_xy, bool keepcms_together, bool keep_max);

void create_alpha_system(int what_alpha, bool remove_pbc);

#endif /* MAKESYSTEM_H_ */
