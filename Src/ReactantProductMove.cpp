/*
 * ReactantProductMove.cpp
 *
 *  Created on: Jul 20, 2009
 *      Author: kedar
 */

#include <ReactantProductMove.h>
bool ADDITIONAL_ATOMS=false;
bool RINGS_NO=false;
char ADDITIONAL_ATOMS_FILE[512] ="";

int extra_atoms = 0;
int *extra_atom_list;
bool read_atoms(char *filename_e_atoms)
{
	bool return_val = false;
	FILE *fptr;
	fptr = fopen(filename_e_atoms,"r");
	if(fptr!=NULL)
	{
		int dummy1=-1;
		fscanf(fptr,"%d",&dummy1);
		if((dummy1>0))
		{
			extra_atoms = dummy1;
			extra_atom_list = (int *) malloc((extra_atoms)*sizeof(int));

			for(int i=0;i<extra_atoms;i++)
			{
				extra_atom_list[i] = -1;
				fscanf(fptr,"%d",&dummy1);
				if(dummy1>0) extra_atom_list[i] = dummy1;
			}
			return_val = true;
		}else
		{
			cout << "error in the input file, ignoring the file\n";
		}
	}

	return return_val;
}

void save_atom_indices_in(int num)
{
	FILE *fptr1;

	char filename[80]="atom_indices_part.",str[80];
	sprintf(str,"%d",num);
	strcat(filename,str);
	fptr1 = fopen(filename,"w");

	for(int i =0;i<n;i++)
	{
		if(atom[i].ackN == -1)
		{
			fprintf(fptr1,"%d\n",1);
		}else
		{
			fprintf(fptr1,"%d\n",0);
		}
	}
	fclose(fptr1);
}

int mark_recurse(atomic_dat *atom_reactant,atomic_dat *atom_product, int atom_number,int nbr_pass)
{
	//cout << "entered here\t"<< atom_number<<"\t"<<nbr_pass<<"\n";
	for(int j=0;j<MAX_COORD;j++)
	{

		int nbr_pass1 = nbr_pass;
		int k = atom_product[atom_number].coord_id[j];
		if((k>-1))
		{
			if(!((atom_product[k].CNA>6)&&(atom_product[k].interface==1)&&(atom_product[k].type==1)))
			{
				if((atom_product[k].interface==1)&&(atom_product[k].type==1))
		//		fprintf(fptr,"%d %d %d %d %lf24.12 %lf24.12 %lf24.12 \n",k,atom_product[k].type,atom_reactant[k].CNA,atom_product[k].CNA,atom_reactant[k].sx-atom_product[k].sx,atom_reactant[k].sy-atom_product[k].sy,atom_reactant[k].sz-atom_product[k].sz);
                               //if(atom_product[k].type!=2)
                                // {
                                  atom_reactant[k].sx = atom_product[k].sx;
                                  atom_reactant[k].sy = atom_product[k].sy;
                                  atom_reactant[k].sz = atom_product[k].sz;
                                  atom_reactant[k].ackN = -1;
                                 //}
			}

			nbr_pass1 = nbr_pass1-1;
			if(nbr_pass1>0)
			{
				mark_recurse(atom_reactant,atom_product, k,nbr_pass1);

			}
		}

	}
}

int move_partconfig_recursive(atomic_dat *atom_reactant,atomic_dat *atom_product, double H_reactant[3][3],double H_product[3][3],int n,bool reverse, int neighbors_toconsider)
{
	// for now atoms in 4- and 5-member rings and atoms that are neighbors to these atoms
	int here=0;
	//general if want to view - Begin
	for(int i=0;i<n;i++)
	{
		atom_reactant[i].ackN = 0;
	}
	//general if want to view - End

	for(int i=0;i<n;i++)
	{
//		if((atom_product[i].CNA>6)&&(atom_product[i].interface==1)&&(atom_product[i].type==1))
	  if((atom_product[i].CNA>6))
		{

			atom_reactant[i].sx = atom_product[i].sx;
			atom_reactant[i].sy = atom_product[i].sy;
			atom_reactant[i].sz = atom_product[i].sz;
					atom_reactant[i].ackN = -1;
			here++; cout << here<<" SDFSDFSFSDFSDFS\n";
			if(i==0) cout <<" I was detected as CNA 6\n";
			int nbrs_pass = neighbors_toconsider;
			mark_recurse(atom_reactant,atom_product, i,nbrs_pass);
		}
	}



}



int move_partconfig(atomic_dat *atom_reactant,atomic_dat *atom_product, double H_reactant[3][3],double H_product[3][3],int n,bool reverse, bool more)
{
	// for now atoms in 4- and 5-member rings and atoms that are neighbors to these atoms
	FILE *fptr;
	fptr = fopen("change_pos","w");
	int here=0;
	for(int i=0;i<n;i++)
	{
		atom_reactant[i].ackN = 0;
	}
	for(int i=0;i<n;i++)
	{
		if((atom_product[i].CNA>6)&&(atom_product[i].interface==1)&&(atom_product[i].type==1))
		{
			if(atom_reactant[i].CNA<=6)
			fprintf(fptr,"%d %d %d %d %24.12lf %24.12lf %24.12lf \n",i,atom_product[i].type,atom_reactant[i].CNA,atom_product[i].CNA,atom_reactant[i].sx-atom_product[i].sx,atom_reactant[i].sy-atom_product[i].sy,atom_reactant[i].sz-atom_product[i].sz);

			atom_reactant[i].sx = atom_product[i].sx;
			atom_reactant[i].sy = atom_product[i].sy;
			atom_reactant[i].sz = atom_product[i].sz;
					atom_reactant[i].ackN = -1;
			here++; cout << here<<" SDFSDFSFSDFSDFS\n";
			if(i==0) cout <<" I was detected as CNA 6\n";
			for(int j=0;j<MAX_COORD;j++)
			{
				int k = atom_product[i].coord_id[j];
				if((k>-1))
				{
					if(!((atom_product[k].CNA>6)&&(atom_product[k].interface==1)&&(atom_product[k].type==1)))
					{
						if((atom_product[k].interface==1)&&(atom_product[k].type==1))
				//		fprintf(fptr,"%d %d %d %d %lf24.12 %lf24.12 %lf24.12 \n",k,atom_product[k].type,atom_reactant[k].CNA,atom_product[k].CNA,atom_reactant[k].sx-atom_product[k].sx,atom_reactant[k].sy-atom_product[k].sy,atom_reactant[k].sz-atom_product[k].sz);

						atom_reactant[k].sx = atom_product[k].sx;
						atom_reactant[k].sy = atom_product[k].sy;
						atom_reactant[k].sz = atom_product[k].sz;
						atom_reactant[k].ackN = -1;
						//atom_reactant[k].CNA=11;
						if(k==0) cout <<i<<"\t K was detected as CNA 6\n";
					}
					if(more){
					// - begin here for more nbrs
					for(int j1=0;j1<MAX_COORD;j1++)
					{
						int k1 = atom_product[k].coord_id[j1];
						if((k1>-1))
						{
							if(!((atom_product[k1].CNA>6)&&(atom_product[k1].interface==1)&&(atom_product[k1].type==1)))
							{
								fprintf(fptr,"%d %d %lf %lf %lf \n",k1,atom[k1].type,atom_reactant[k1].sx-atom_product[k1].sx,atom_reactant[k1].sy-atom_product[k1].sy,atom_reactant[k1].sz-atom_product[k1].sz);
								atom_reactant[k1].sx = atom_product[k1].sx;
								atom_reactant[k1].sy = atom_product[k1].sy;
								atom_reactant[k1].sz = atom_product[k1].sz;
								atom_reactant[k1].ackN = -1;
							//	atom_reactant[k1].CNA=11;
								if(k1==0) cout <<i<<"\t"<<k1<<"\t K1 was detected as CNA 6\n";
							}

							for(int j2=0;j2<MAX_COORD;j2++)
							{
								int k2 = atom_product[k1].coord_id[j2];
								if((k2>-1))
								{
									if(!((atom_product[k2].CNA>6)&&(atom_product[k2].interface==1)&&(atom_product[k2].type==1)))
									{
										fprintf(fptr,"%d %d %lf24.12 %lf24.12 %lf24.12 \n",k2,atom[k2].type,atom_reactant[k2].sx-atom_product[k2].sx,atom_reactant[k2].sy-atom_product[k2].sy,atom_reactant[k2].sz-atom_product[k2].sz);
										atom_reactant[k2].sx = atom_product[k2].sx;
										atom_reactant[k2].sy = atom_product[k2].sy;
										atom_reactant[k2].sz = atom_product[k2].sz;
										atom_reactant[k2].ackN = -1;
									//	atom_reactant[k2].CNA=11;
										if(k2==0) cout <<i<<"\t"<<k2<<"\t K2 was detected as CNA 6\n";
									}


									for(int j3=0;j3<MAX_COORD;j3++)
									{
										int k3 = atom_product[k2].coord_id[j3];
										if((k3>-1))
										{
											if(!((atom_product[k3].CNA>6)&&(atom_product[k3].interface==1)&&(atom_product[k3].type==1)))
											{
												//											fprintf(fptr,"%d %d %lf %lf %lf \n",k2,atom[k2].type,atom_reactant[k2].sx-atom_product[k2].sx,atom_reactant[k2].sy-atom_product[k2].sy,atom_reactant[k2].sz-atom_product[k2].sz);
												atom_reactant[k3].sx = atom_product[k3].sx;
												atom_reactant[k3].sy = atom_product[k3].sy;
												atom_reactant[k3].sz = atom_product[k3].sz;
												atom_reactant[k3].ackN = -1;
											//	atom_reactant[k3].CNA=11;
												if(k3==0) cout <<i<<"\t"<<k3<<"\t K3 was detected as CNA 6\n";
											}
										}
									}
								}
							}
						}
					}
					}
				}
			}
		}
	}

	fclose(fptr);

}


int reactant_product_config_changes(char *filename_reactant, char *filename_product,int file_number, bool reverse, bool more, bool cfg_save, bool lammps_save, int neighs_to_stick)
{
	atomic_dat *atom_ref;

	double H_ref[3][3];

	//int is_read = read_lammps(filename_reactant, atom,true,true);
	int is_read = read_lammps_general(filename_reactant);
	if(is_read!=0) return -1;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			H_ref[i][j] = Hcry[i][j];
			cout << H_ref[i][j]<<"\t";
		}
		cout <<" Reference H*****\n";
	}
	compute_CNA_and_others(atom,n, Hcry);
	atom_ref = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

	copy_atomstruct(atom, atom_ref,n);

	is_read = read_lammps_general(filename_product);
	if(is_read!=0)
	{
			cout << "HERE INSIDE OUT\n";
			free(atom_ref);
			return -2;
	}
//	else{
//	  // this is setting the thing by changing a little bit; for now.
//		cout << "there is some atom change going on, is this correct?\n\n\n";
//	  for(int i=0;i<n;i++)
//	    {
//	      double sij[3];
//	      double rij[3];
//	      sij[0] =atom[i].sx;
//	      sij[1] =atom[i].sy;
//	      sij[2] =atom[i].sz;
//
//	      V3mulM3(sij,Hcry,rij);
//	      if(atom[i].type==1)
//	        {
//	          rij[0] = rij[0]-3;
//	          rij[1]= rij[1]-2;
//	        }else if(atom[i].type==2)
//	          {
//	            rij[0] = rij[0]+3;
//	            rij[1]= rij[1]+2;
//	          }
//	      sij[0]=sij[1]=sij[2] = 0;
//	      V3mulM3(rij,Hcry_inv,sij);
//	      atom[i].sx = sij[0];atom[i].sy=sij[1];atom[i].sz=sij[2];
//	      if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
//	      if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
//	      if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;
//
//	      if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
//	      if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
//	      if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;
//	    }
//	}
	compute_CNA_and_others(atom,n, Hcry);

//	FILE *fptr;
//	fptr = fopen("check_positions","w");

	for(int i=0;i<n;i++)
	{

	//	cout << "WARNING WARNING WARNING\n\n\n ***************\n product CNA set to -1 for some debugging\n";
	//	cout << "change code in ReactantProductMove.cpp if this is not desirable\n";
	//	cout << "***********************************\n";

		if(RINGS_NO)
		{
			if(atom[i].CNA>6)
			{
			atom[i].CNA = 6;
			atom_ref[i].CNA = 6;

			}

		}else
		{
			if(atom_ref[i].CNA>6)
				{
				atom[i].CNA=atom_ref[i].CNA;

				}

		}
	}

//	ADDITIONAL_ATOMS = true;
	cout <<"\t\t\t\t Entering check\n";
	if(read_atoms(ADDITIONAL_ATOMS_FILE))
	{
	  cout << "Able to read file \n\n\n\n\n";
		for(int i=0;i<extra_atoms;i++)
		{
			atom[extra_atom_list[i]].CNA = 10;
			atom_ref[extra_atom_list[i]].CNA = 10;
	//		cout << extra_atom_list[i]<<"\t read this atoms number\n\n\n\n";
		}
	}else
	{
		ADDITIONAL_ATOMS = false;
	}



	//move_partconfig(atom_ref,atom,H_ref,Hcry,n,reverse,more);
	move_partconfig_recursive(atom_ref,atom,H_ref,Hcry,n,reverse,neighs_to_stick);
	copy_atomstruct(atom_ref, atom,n);
	free(atom_ref);
	//atom = atom_ref;

//	compute_CNA_and_others(atom,n, Hcry);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			Hcry[i][j] = H_ref[i][j];
			cout << Hcry[i][j]<<"\t";
		}
		cout <<" Reference H*****\n";
	}
//			compute_CNA_and_others(atom,n, Hcry);
//		compute_disrigistry(atom,atom_start,Hcry,Hcry,false);
//		compute_slipvector(atom,atom_start,Hcry,Hcry,false);
	if(cfg_save) save_cfg(file_number,H_ref);
	if(lammps_save) save_lammps(file_number,H_ref);
	//save_atom_indices_in(file_number);
	cout << "here outside out getting out\n";
	//free(atom_ref);
	return 0;
}



int prepare_for_chain( int start_filenumber, int end_filenumber, int interval_filenumber, int start_number, int *end_number, bool minimize_yes_no, char *script_title, bool cfg_save, bool lammps_save, int neighs_to_stick)
{

int change_prev=start_filenumber;
int change_next=start_number;
for(int ii=start_filenumber;ii<=end_filenumber;ii+=interval_filenumber)
{
	char filename[80]="dat.",str[80];
	char filename2[80]="dat.",str2[80];
	sprintf(str,"%d",ii);
	strcat(filename,str);
	if((ii-filenumber_start)==0)
	{
		sprintf(str2,"%d",ii);
		strcat(filename2,str2);

	}else
	{
		sprintf(str2,"%d",change_prev);
		strcat(filename2,str2);
	}
	sprintf(str2,"%d",change_next);
	cout << str2<<"new modified***************\n";
	cout << filename<<"\t"<<filename2<<" here here"<<"\t"<<ii<<"\n";
	int done = -5;
	if((ii-filenumber_start)!=0)
		{
			done = reactant_product_config_changes(filename2,filename,change_next,false,true, cfg_save, lammps_save, neighs_to_stick);
		}
	if(done==0)
	{
		if(minimize_yes_no)
		{
			char system_command[80]="./",str3[80];
			strcat(system_command, script_title);
			strcat(system_command, " ");
			strcat(system_command,str2);
			strcat(system_command," ");
			sprintf(str3,"%d",change_next+1);
			strcat(system_command,str3);
			strcat(system_command," 1");
			cout << system_command<<"new modified***************\n";
			execute_system_command(system_command);
		}

	change_prev = change_next;
	change_next=change_next+1;
	}
}
*end_number = change_prev;
}


