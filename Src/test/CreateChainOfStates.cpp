/*
 * CreateChainOfStates.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */
#include <CreateChainOfStates.h>
/*
void determine_distances(int iter, int *filename_array, bool all)
{

	FILE *fptr;
	atomic_dat *atom_prev;
	double initial = energies[0];
	if(all)
	{

		double f_mag_1[total_length];

		for(int i=0;i<total_length;i++)
		{
			char filename[80]="dat.",str[80];
			sprintf(str,"%d",filename_array[i]);
			strcat(filename,str);
			int is_read = read_lammps(filename,atom,true,true);
			if(i == 0)
			{
				atom_prev = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
				copy_atomstruct(atom, atom_prev,n);
			}

			reac_coord[i] =0.0;
			f_mag_1[i] = 0.0;
			for(int ii=0;ii<n;ii++)
			{
				double r[3],s[3];
				double delx = (atom[ii].sx-atom_prev[ii].sx);
				double dely = (atom[ii].sy-atom_prev[ii].sy);
				double delz = (atom[ii].sz-atom_prev[ii].sz);
				delx = delx-(int)(delx*2);
				dely = dely-(int)(dely*2);
				delz = delz-(int)(delz*2);
				s[0] = delx;s[1] =dely;s[2] = delz;
				V3mulM3(s,Hcry,r);
				reac_coord[i] +=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
				double ting = sqrt(atom[ii].fx*atom[ii].fx+atom[ii].fy*atom[ii].fy+atom[ii].fz*atom[ii].fz);
				if(ting > f_mag_1[i]) f_mag_1[i]= ting;

			}

			reac_coord[i] = sqrt(reac_coord[i]);

			if(i!=0)
			{
				reac_coord[i]+=reac_coord[i-1];
				copy_atomstruct(atom, atom_prev,n);
			}

		}
		free(atom_prev);
		char fn[80]="energy.",str[80];
		sprintf(str,"%d",iter);
					strcat(fn,str);
		fptr = fopen(fn,"w");
		fprintf(fptr, "\n");
		double R1_max = reac_coord[total_length-1]/R_max;
		for(int i =0;i< total_length;i++)
		{
			fprintf(fptr, "%d %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf\n",i,energies[i]-initial,neb_forces[i],true_forces[i],spring_forces[i],t_mag_all[i],t_l_all[i],t_r_all[i],iter, reac_coord[i]/reac_coord[total_length-1],reac_coord[i], reac_coord[total_length-1],sqrt(max_neb_forces[i]),f_mag_1[i]);
		}
		fclose(fptr);
	}
	char fn1[80]="neb_energy.info";
	fptr = fopen(fn1,"w");
	for(int i =0;i< total_length;i++)
	{
		fprintf(fptr, "%d %lf %lf %lf %lf %d %lf\n",i,energies[i]-initial,neb_forces[i],true_forces[i],spring_forces[i],iter,sqrt(max_neb_forces[i]));
	}
	fclose(fptr);


}
*/
int create_chain_of_states(char *filename_reactant, char *filename_product, int number_of_states, int start_chain_number, bool save_first, bool save_last)
{
//assumes constant H matrix and isoatomic configurations
	atomic_dat *atom_start,*atom_end;
        atom = NULL;
	double H_start[3][3]; double H_end[3][3];
	cout << filename_reactant<<"\t"<<filename_product<<"\n";
	int is_read = read_lammps(filename_reactant, atom,true,true);
	if(is_read!=0) return -1;
	atom_start = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

	int n1 = n;

	copy_atomstruct(atom, atom_start,n);
	//	free(atom);
	is_read = read_lammps(filename_product, atom,true,true);
	if(is_read!=0) return -2;
	atom_end = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
	cout << "both files read\n";
	if(n1!=n) { cout << "not isoatomic start and end files. Can not build chain of states\n"; return -1;}

	copy_atomstruct(atom, atom_end,n);
	//	free(atom);

	int save_number=0;
	for(int k=0;k<number_of_states;k++)
	{
	  if(atom!=NULL) free(atom);
		atom = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
		copy_atomstruct(atom_start, atom,n);
		for(int i=0;i<n;i++)
		{

			double delx = (atom_end[i].sx-atom_start[i].sx);
			double dely = (atom_end[i].sy-atom_start[i].sy);
			double delz = (atom_end[i].sz-atom_start[i].sz);
			delx = delx-(int)(delx*2);
			dely = dely-(int)(dely*2);
			delz = delz-(int)(delz*2);
			atom[i].sx +=(k*1.0/((number_of_states-1)*1.0))*delx;
			atom[i].sy +=(k*1.0/((number_of_states-1)*1.0))*dely;
			atom[i].sz +=(k*1.0/((number_of_states-1)*1.0))*delz;
			atom[i].vx = atom[i].vy=atom[i].vz = 0.0;


		}
	//	compute_CNA_and_others(atom,n, Hcry);
	//	compute_disrigistry(atom,atom_start,Hcry,Hcry,false);
	//	compute_slipvector(atom,atom_start,Hcry,Hcry,false);
		if(((save_first)&&(k==0))||((k!=0)&&(k!=(number_of_states-1)))||((save_last)&&(k==number_of_states-1)))
		{

		//	save_cfg(save_number+start_chain_number,Hcry);
			save_lammps(save_number+start_chain_number,Hcry);
			save_number++;
		}
		free(atom); atom=NULL;
	}

	free(atom_start); free(atom_end);
	return 0;
}



int create_recurring_chain(int start_filenumber, int end_filenumber, int interval_filenumber, int begin, int chain_length)
{
cout <<"hey, entered here\n";
int change_prev=start_filenumber;
int change_next=change_prev;
int chain_number_begin = begin;
int chain_number_now = chain_number_begin;
cout << start_filenumber<<"\t"<<end_filenumber<<"\t"<<interval_filenumber<<"\n";
for(int ii=start_filenumber;ii<=end_filenumber;ii+=interval_filenumber)
{
	char filename[80]="dat.",str[80];
	char filename2[80]="dat.",str2[80];
	sprintf(str,"%d",ii);
	strcat(filename,str);
	if((ii-start_filenumber)==0)
	{
		sprintf(str2,"%d",ii);
		strcat(filename2,str2);

	}else
	{
		sprintf(str2,"%d",change_prev);
		strcat(filename2,str2);
	}
	sprintf(str2,"%d",change_next);
	cout << str2<<"new chain***************\n";
	cout << filename<<"\t"<<filename2<<"\t"<<chain_number_now<<"\n";
	int done = 10;
	if((ii-start_filenumber)!=0)
	{
		if(chain_number_begin==chain_number_now)
		{
			done = create_chain_of_states(filename2, filename,chain_length,chain_number_now,true,true);

		}else
		{
			done = create_chain_of_states(filename2, filename,chain_length+1,chain_number_now,false,true);

		}

	}
	if(done ==0)
	{
		chain_number_now+=chain_length;
		change_prev = ii;
		//change_next=ii;
	}
}
	return 0;
}

int create_recurring_chain_equidistant(int start_filenumber, int end_filenumber, int interval_filenumber, int begin, int total_chain_length)
{
cout <<"hey, entered here\n";
int change_prev=start_filenumber;
int change_next=change_prev;
int chain_number_begin = begin;
int chain_number_now = chain_number_begin;
int chain_length=0;
cout << start_filenumber<<"\t"<<end_filenumber<<"\t"<<interval_filenumber<<"\n";
int actual_number = 0;
for(int ii=start_filenumber;ii<=end_filenumber;ii+=interval_filenumber)
{
	char filename[80]="dat.",str[80];
	char filename2[80]="dat.",str2[80];
	sprintf(str,"%d",ii);
	strcat(filename,str);

	if ( check_file_exists(filename) ) actual_number++;
}
cout << actual_number<<" is the number of files available\n";

double del_reac_coord[actual_number];

actual_number = 0;
atomic_dat *atom_prev;
for(int i=start_filenumber;i<=end_filenumber;i+=interval_filenumber)
{
		char filename[80]="dat.",str[80];
		sprintf(str,"%d",i);
		strcat(filename,str);
		if ( check_file_exists(filename) )
		{
			int is_read = read_lammps(filename,atom,true,true);
			if((i-start_filenumber) == 0)
			{
				atom_prev = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
				copy_atomstruct(atom, atom_prev,n);
			}
			del_reac_coord[actual_number] =0.0;
			for(int ii=0;ii<n;ii++)
			{
				double r[3],s[3];
				double delx = (atom[ii].sx-atom_prev[ii].sx);
				double dely = (atom[ii].sy-atom_prev[ii].sy);
				double delz = (atom[ii].sz-atom_prev[ii].sz);
				delx = delx-(int)(delx*2);
				dely = dely-(int)(dely*2);
				delz = delz-(int)(delz*2);
				s[0] = delx;s[1] =dely;s[2] = delz;
				V3mulM3(s,Hcry,r);
//				cout << "here"<<"\t"<<ii<<"\n";
				del_reac_coord[actual_number] +=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];

			}
			del_reac_coord[actual_number] = sqrt(del_reac_coord[actual_number]);
			if((i-start_filenumber)!=0)
			{
				del_reac_coord[actual_number] += del_reac_coord[actual_number-1];
				copy_atomstruct(atom, atom_prev,n);
			}
			actual_number++;
		}


}

for(int i =0;i<actual_number;i++)
{
	del_reac_coord[i] = del_reac_coord[i]/del_reac_coord[actual_number-1];
	//cout << i<<"\t"<<del_reac_coord[i]<<"\n";
}

int counter = 1;
	for(int ii=start_filenumber;ii<=end_filenumber;ii+=interval_filenumber)
	{
		char filename[80]="dat.",str[80];
		char filename2[80]="dat.",str2[80];
		sprintf(str,"%d",ii);
		strcat(filename,str);
		if((ii-start_filenumber)==0)
		{
			sprintf(str2,"%d",ii);
			strcat(filename2,str2);

		}else
		{
			sprintf(str2,"%d",change_prev);
			strcat(filename2,str2);
		}
		sprintf(str2,"%d",change_next);
		cout << str2<<"new chain***************\n";
		cout << filename<<"\t"<<filename2<<"\t"<<chain_number_now<<"\n";

		int done = 10;
		if ( (check_file_exists(filename)) && (check_file_exists(filename2))&&((ii-start_filenumber)!=0))
		{

			chain_length = (int)(((double)total_chain_length)*(del_reac_coord[counter]-del_reac_coord[counter-1])+0.5);
			if(counter==actual_number-1)
			{

				cout << "LAST LEG\t"<< chain_length+chain_number_now<<"\t"<<(chain_number_now-begin)<<"\n";
				if((chain_length+chain_number_now)!=(total_chain_length+begin))
				{
					chain_length = total_chain_length+begin-chain_number_now;
					cout << "came here\n";
				}
			}
			cout << filename<<"\t"<<filename2<<"\t"<<chain_length<<"\t"<<counter<<"\t"<<actual_number<<" CHAIN INFORMATION\n";

			if(chain_number_begin==chain_number_now)
				{
					done = create_chain_of_states(filename2, filename,chain_length,chain_number_now,true,true);

				}else
				{
					cout << "\t here\n";
					done = create_chain_of_states(filename2, filename,chain_length+1,chain_number_now,false,true);

				}

				if(done ==0)
				{
					chain_number_now+=chain_length;
					change_prev = ii;
					//change_next=ii;
					counter++;
				}else
				{
					cout << "some error in create_chain_of_states_recurring_equidistant\n";exit(1);
				}
		}
	}
	return 0;
}


double compute_distance_scalar(atomic_dat *atom_start, atomic_dat *atom_end, int n_now, double HCry_init[3][3],double Hcry_final[3][3])
{
//for now, assuming that Hcry_init = Hcry_final
	double r=0.0;

	for(int i=0;i<n_now;i++)
	{
		double delx = (atom_end[i].sx-atom_start[i].sx);
		double dely = (atom_end[i].sy-atom_start[i].sy);
		double delz = (atom_end[i].sz-atom_start[i].sz);
		delx = delx-(int)(delx*2);
		dely = dely-(int)(dely*2);
		delz = delz-(int)(delz*2);
		r = r+ delx*delx+dely*dely+delz*delz;

	}

	return sqrt(r);
}
