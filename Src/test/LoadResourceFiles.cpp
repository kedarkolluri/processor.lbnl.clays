/*
 * LoadResourceFiles.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#include <LoadResourceFiles.h>

void load_reference_vectors(char *REF_STRING)
{
	//reference vectors are in the local folder by the name "ref.vectors" with their distances which are discarded here
	// changed for the strucutres given by Demkowicz - X MADE TO Y AND Y MADE TO X
	FILE *fptr;
	double dummy;
	int counter=0;
	basic_n[0] = 12;
	basic_n[1] = 8;
	basic_n[2] = 12;
	basic_n[3] = 10;
	//fptr = fopen("ref.vectors","r");
	fptr = fopen(REF_STRING,"r");
	for(int i=0;i<basic_n[0];i++)
	{
		fscanf(fptr,"%lf %lf %lf %lf",&ref_vectors[counter][0],&ref_vectors[counter][1],&ref_vectors[counter][2],&dummy);
		cout << counter<<"\t"<<dummy<<"\n";
		counter++;
	}
	for(int i=0;i<basic_n[1];i++)
	{
		fscanf(fptr,"%lf %lf %lf %lf",&ref_vectors[counter][0],&ref_vectors[counter][1],&ref_vectors[counter][2],&dummy);
		cout << counter<<"\t"<<dummy<<"\n";
		counter++;
	}

	for(int i=0;i<basic_n[2];i++)
	{
		fscanf(fptr,"%lf %lf %lf %lf",&ref_vectors[counter][0],&ref_vectors[counter][1],&ref_vectors[counter][2],&dummy);
		cout << counter<<"\t"<<dummy<<"\n";
		counter++;
	}

	for(int i=0;i<basic_n[3];i++)
	{
		fscanf(fptr,"%lf %lf %lf %lf",&ref_vectors[counter][0],&ref_vectors[counter][1],&ref_vectors[counter][2],&dummy);
		cout << counter<<"\t"<<dummy<<"\n";
		counter++;
	}

	fclose(fptr);

}
