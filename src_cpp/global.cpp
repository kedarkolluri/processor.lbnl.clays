/*
 * global.cpp
 *
 *  Created on: Jul 18, 2009
 *      Author: kedar
 */

#include <global.h>


int *nbr_ptr=NULL,*nbr_ptr1=NULL,*nbr_lst=NULL,*map1=NULL,*map_full=NULL,*head=NULL,*list=NULL,mx,my,mz,ncell;

FILE *fptr;
int format,target_format;
//char *inputfilename, *outputfilename;
std::string inputfilename, outputfilename;

int interface_atoms;
int *tag_array_interface_atoms;
int total_number_of_files=0;

int filenumber_start,filenumber_end,filenumber_interval;


bool TEST_BOOL;
int MAX_COORDNUM;
