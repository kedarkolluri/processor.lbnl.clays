/*
 * OPToptions.cpp
 *
 *  Created on: Oct 6, 2009
 *      Author: kedar
 */

#include "OPToptions.h"

OPTOPTIONS_NEB_STATIC_MODE_OPTIONS OPTOPTIONS_NEB_STATIC_MODE = NONE;
double OPTOPTIONS_NEB_LIMIT = 0.0009;

int getdynamic(int *chain_array, int chain_length, int *dynamic_list, int *single_len, bool *bool_list)
{
	if(dynamic_list == NULL)
	{
		*single_len = 0;
	}else
	{
		free(dynamic_list); dynamic_list = NULL; *single_len = 0;

	}
	if(OPTOPTIONS_NEB_STATIC_MODE == ENERGY)
	{
		int count = 0;
		atomic_dat *initial_state;
		for(int ii=0;ii<chain_length;ii++)
		{
			char filename[80]="dat.",str[80];
			sprintf(str,"%d",chain_array[ii]);
			strcat(filename,str);

			int is_read = read_lammps(filename,atom,true,true);

			if(ii==0)
			{
				initial_state = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

				copy_atomstruct(atom, initial_state,n);
				for(int i =0;i<n;i++)
				{
					bool_list[i] = false;
				}
			}else
			{
				for(int i=0;i<n;i++)
				{
					if(fabs(atom[i].pe-initial_state[i].pe)>OPTOPTIONS_NEB_LIMIT)
					{
						if(!bool_list[i]) count++;
						bool_list[i] = true;

					}
				}
			}

		}
		dynamic_list = (int *)malloc(count*sizeof(int));
		*single_len = count;
		cout << count<<" total dynamic atoms are\n";
		count = 0;
		for(int i = 0;i<n;i++)
		{
			if(bool_list[i]){
				dynamic_list[count] = i;
				count++;
			}
		}
	}else if(OPTOPTIONS_NEB_STATIC_MODE == NONE)
	{
		cout << "all atoms being considered for minimization\n";

	}

}
