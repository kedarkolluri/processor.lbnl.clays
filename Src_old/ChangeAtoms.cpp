/*
 * ChangeAtoms.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#include <ChangeAtoms.h>

void delete_atoms(int *arr, int len)
{

	int n_new = n-len;
	atomic_dat *atom_new;
	atom_new = (struct atomic_dat *) malloc((n_new+3)*sizeof(struct atomic_dat));
	int j = 0;
	for(int i=0;i<n;i++)
	{
		if(!check_repeat(i,arr,len))
		{
			copy_atom_indv(atom,atom_new, i,j);
			j++;
		}
	}
	n = n_new;
	free(atom);
	atom = atom_new;



}

void delete_atom_single(int tag)
{
	int n_new = n-1;
	atomic_dat *atom_new;
	atom_new = (struct atomic_dat *) malloc((n_new+3)*sizeof(struct atomic_dat));
	int j = 0;
	for(int i=0;i<n;i++)
	{
		if(i!=tag)
		{
			copy_atom_indv(atom,atom_new, i,j);
			j++;
		}
	}
	n = n_new;
	free(atom);
	atom = atom_new;
}

void insert_atoms(double (*arr)[4], int len)
{

	int n_new = n+len;
	atomic_dat *atom_new;
	atom_new = (struct atomic_dat *) malloc((n_new+3)*sizeof(struct atomic_dat));
	int j = 0;
	for(int i=0;i<n;i++)
	{

		copy_atom_indv(atom,atom_new, i,j);
		j++;
	}

	for(int i =n;i<n_new;i++)
	{
		atom_new[i].sx = arr[(i-n)][0];
		atom_new[i].sy = arr[(i-n)][1];
		atom_new[i].sz = arr[(i-n)][2];
		atom_new[i].type = arr[(i-n)][3];
	}
	n = n_new;
	free(atom);
	atom = atom_new;



}

void insert_atoms_specific(double (*arr)[5], int len)
{

	int n_new = n+len;
	atomic_dat *atom_new;
	atom_new = (struct atomic_dat *) malloc((n_new+3)*sizeof(struct atomic_dat));
int j = 0;
	for(int i=0;i<n;i++)
	{

		for(int k=0;k<len;k++)
		{
			if(arr[k][4]==i)
			{
				atom_new[j].sx = arr[k][0];
				atom_new[j].sy = arr[k][1];
				atom_new[j].sz = arr[k][2];
				atom_new[j].type = arr[k][3];
				j++;
			}
		}
		copy_atom_indv(atom,atom_new, i,j);
		j++;
	}
	n = n_new;
	free(atom);
	atom = atom_new;



}
