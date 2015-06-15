/*
 * FindClosure.cpp
 *
 *  Created on: Jul 30, 2009
 *      Author: kedar - some simple tests on sept 30 2014
 */

#include <FindClosure.h>

void FindClosure (atomic_dat *atom_ref, atomic_dat *atom_curr, int n, double Hcry_ref[3][3], double Hcry_curr[3][3], bool transform, char *filename, int number_of_atoms, double closure_vector[3])
{

	double s_curr[3],s_ref[3],r_curr[3],r_ref[3],s_prev_curr[3],s_prev_ref[3];

	FILE *fptr;
	int atom_list[number_of_atoms];
	for(int i=0;i<number_of_atoms;i++) atom_list[i]=-1;
	fptr = fopen(filename,"r");
	if(fptr==NULL)
	{
		cout << "atom file list not found, exiting\n";
		return ;
	}
	for(int i=0;i<number_of_atoms;i++)
	{
		fscanf(fptr,"%d",&atom_list[i]);

	}
	bool start=true;
	for(int i=0;i<n;i++)
	{
		if(check_repeat(i,atom_list,number_of_atoms))
		{
			if(start)
			{
				start=false;

				for(int j=0;j<3;j++)
				{
					s_curr[j]=s_ref[j]=r_curr[j]=r_ref[j]=0.0;
				}
			}else
			{
				s_prev_curr[0] += atom_curr[i].sx;s_prev_curr[1]+=atom_curr[i].sy;s_prev_curr[2]+=atom_curr[i].sz;
				s_prev_ref[0] += atom_ref[i].sx;s_prev_ref[1]+=atom_ref[i].sy;s_prev_ref[2]+=atom_ref[i].sz;
				for(int j=0;j<3;j++)
				{
					s_prev_curr[j]=s_prev_curr[j]-(int)(s_prev_curr[j]*2);
					s_prev_ref[j]=s_prev_ref[j]-(int)(s_prev_ref[j]*2);
					s_curr[j] +=s_prev_curr[j];
					s_ref[j] +=s_prev_ref[j];
				}


			}
			s_prev_curr[0] = -1.0*atom_curr[i].sx;s_prev_curr[1]=-1.0*atom_curr[i].sy;s_prev_curr[2]=-1.0*atom_curr[i].sz;
			s_prev_ref[0] = -1.0*atom_ref[i].sx;s_prev_ref[1]=-1.0*atom_ref[i].sy;s_prev_ref[2]=-1.0*atom_ref[i].sz;


		}

		V3mulM3(s_curr,Hcry_curr,r_curr);
		V3mulM3(s_ref,Hcry_ref,r_ref);
		double r_mod_curr[3],r_mod_ref[3],r_mod_del[3];
		V3mulM3(r_curr,KS1_trans_mas,r_mod_curr);
		V3mulM3(r_ref,KS1_trans_mas,r_mod_ref);
		for(int j=0;j<3;j++)
		{
			r_mod_del[j]=r_mod_curr[j]-r_mod_ref[j];
		}
		V3mulM3(r_mod_del,H0_geo,closure_vector);

	}

}
