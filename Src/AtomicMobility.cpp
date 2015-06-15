/*
 * AtomicMobility.cpp
 *
 *  Created on: Jul 22, 2009 dfd
 *      Author: kedar
 */

#include "AtomicMobility.h"

void find_relative_distance(atomic_dat *atom_curr,atomic_dat *atom_ref,int n,int atom_id,double r12[3], double H_here[3][3])
{
/*
	static int **check_r;
	static bool start=false;
	if(!start)
	{
		check_r= (int **) malloc((n)*sizeof(int));
		for(int i=0;i<n;i++)
		{
			check_r[i]= (int *) malloc((n)*sizeof(int));
		}
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				check_r[i][j]=-1;
			}
		}
	}
*/
	int jbeg = nbr_ptr[atom_id];
	int jend = nbr_ptr1[atom_id];
	double sij[3]; for(int h=0;h<3;h++)sij[h]=0.0;
	double sij_net[3];for(int h=0;h<3;h++)sij_net[h]=0.0;
	for(int h=0;h<3;h++)r12[h]=0.0;
	int count=0;

	for(int jnab = jbeg; jnab<=jend;jnab++)
	{
		//int j = atom_curr[atom_id].coord_id[jnab];
	 int j = nbr_lst[jnab];
		if(j>-1)
		{
			if((atom_curr[j].interface==1)&&(atom_curr[j].type==1))
			{
//				(check_r[atom_id][j]==-1)&&
				//check_r[atom_id][j]=1;
//				check_r[j][atom_id]=1;
				sij[0] = atom_curr[atom_id].sx - atom_curr[j].sx;
				sij[1] = atom_curr[atom_id].sy - atom_curr[j].sy;
				sij[2] = atom_curr[atom_id].sz - atom_curr[j].sz;
				sij[0] = sij[0]-(int)(sij[0]*2);
				sij[1] = sij[1]-(int)(sij[1]*2);
				sij[2] = sij[2]-(int)(sij[2]*2);
				double rij[3];for(int h=0;h<3;h++)rij[h]=0.0;
				V3mulM3(sij,H_here,rij);
		double rijsq =rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
		if (rijsq <= rcoordsq[atom_curr[atom_id].type-1][atom_curr[j].type-1])
		{
			double sdel[3];for(int h=0;h<3;h++) sdel[h]=0.0;
			double rdel[3];for(int h=0;h<3;h++) rdel[h]=0.0;
			sdel[0]=sij[0];
			sdel[0]=sij[1];
			sdel[2]=sij[2];

			sij[0] = atom_ref[atom_id].sx - atom_ref[j].sx;
			sij[1] = atom_ref[atom_id].sy - atom_ref[j].sy;
			sij[2] = atom_ref[atom_id].sz - atom_ref[j].sz;
			sij[0] = sij[0]-(int)(sij[0]*2);
			sij[1] = sij[1]-(int)(sij[1]*2);
			sij[2] = sij[2]-(int)(sij[2]*2);


			sdel[0]-=sij[0];
			sdel[0]-=sij[1];
			sdel[2]-=sij[2];
			V3mulM3(sdel,H_here,rdel);
			if(atom_id==20333)
			{
				cout << count <<" for atom is \t"<<r12[0]<<"\t"<<r12[1]<<"\t"<<r12[2]<<"\n";
			}
			if(((sqrt(rdel[0]*rdel[0]+rdel[1]*rdel[1]+rdel[2]*rdel[2]))) > 0.4)
			{
				count++;
				for(int h=0;h<3;h++) r12[h]+=rdel[h];
			}
		}
			}
		}
	}
	if(count<0)
	{
		for(int h=0;h<3;h++) r12[h]/=count;
	}

}

void neighbors_check(atomic_dat *data_toconsider, atomic_dat *refdata, int n_here, int *markasflagged)
{
	for(int i=0;i<n_here;i++)
	{
		data_toconsider[i].ackN = 0;
		if(data_toconsider[i].CNA<7)
		data_toconsider[i].CNA = 0;
		for(int j=0;j<MAX_COORD;j++)
		{
			int k = data_toconsider[i].coord_id[j];

			if ((!(check_repeat(k,refdata[i].coord_id,MAX_COORD)))&&(data_toconsider[k].interface==1)&&(data_toconsider[k].type==1))
			{

				data_toconsider[i].ackN++;
			}
			int k1 = refdata[i].coord_id[j];
			if ((!(check_repeat(k1,data_toconsider[i].coord_id,MAX_COORD)))&&(data_toconsider[k].interface==1)&&(data_toconsider[k].type==1))
			{
				data_toconsider[i].ackN++;
			}


		}


		if(data_toconsider[i].ackN>1)
		{
			markasflagged[i] = data_toconsider[i].ackN;
		}else
		{
			data_toconsider[i].ackN=0;
			markasflagged[i] = 0;
		}
	}
}
void recursive_neighbors(int levels, atomic_dat *data_to_consider,int n, int *markasflagged)
{

	cout << "Levels here are \n \n"<< levels<<"\n";
	int levels_here = levels-1;
	static int start_count=0;
	if(levels_here>=0)
	{
		for(int i=0;i<n;i++)
		{
			if(markasflagged[i]==1)
			{
				for(int j=0;j<MAX_COORD;j++)
				{
					int k = data_to_consider[i].coord_id[j];
					if((k>-1)&&(markasflagged[k]<1))
					{

						markasflagged[k] = 2;
						start_count++;
					}
				}
				markasflagged[i]=3;
			}

		}

		for(int i=0;i<n;i++)
		{
			if(markasflagged[i]==2) markasflagged[i]=1;
		}
		cout << start_count << "\t is the count at the end of level\t"<<levels<<"\n";
		recursive_neighbors(levels_here, data_to_consider,n,markasflagged);
	}

}

void find_atoms_to_consider_formobility(atomic_dat *atom_ref,atomic_dat *atom_curr,int n, int *consider_formobility)
{
		// for now atoms in 4- and 5-member rings and atoms that are neighbors to these atoms
		FILE *fptr;
		int here=0;
		atomic_dat *atom_product;
		for(int count=0;count<2;count++)
		{
			if(count==0) atom_product = atom_ref;
			if(count==1) atom_product = atom_curr;

		for(int i=0;i<n;i++)
		{
			if((atom_product[i].CNA>6)&&(atom_product[i].interface==1)&&(atom_product[i].type==1))
			{

				consider_formobility[i] = 1;
				atom_curr[i].ackN = -15;
				for(int j=0;j<MAX_COORD;j++)
				{
					int k = atom_product[i].coord_id[j];
					if((k>-1))
					{
						if(!((atom_product[k].CNA>6)&&(atom_product[k].interface==1)&&(atom_product[k].type==1)))
						{

							consider_formobility[k] = 1;
							atom_curr[k].ackN = -15;
						}
						if(true){
						// - begin here for more nbrs
						for(int j1=0;j1<MAX_COORD;j1++)
						{
							int k1 = atom_product[k].coord_id[j1];
							if((k1>-1))
							{

								consider_formobility[k1] = 1;
								atom_curr[k1].ackN = -15;
								for(int j2=0;j2<MAX_COORD;j2++)
								{
									int k2 = atom_product[k1].coord_id[j2];
									if((k2>-1))
									{
										/*
										consider_formobility[k2] = 1;
										atom_curr[k2].ackN = -15;

										for(int j3=0;j3<MAX_COORD;j3++)
										{
											int k3 = atom_product[k2].coord_id[j3];
											if((k3>-1))
											{
												consider_formobility[k3] = 1;
												atom_curr[k3].ackN = -15;
											}
										}*/
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

void determine_mobility1(atomic_dat *atom_ref,atomic_dat *atom_curr,int n,double Href[3][3],double H_curr[3][3],char *filename_save, int tag_reference,int tag_target, bool append, bool relative)
{
	int *consider_formobility;
	consider_formobility = (int *) malloc((n)*sizeof(int));
	for(int i=0;i<n;i++)
	{
		consider_formobility[i]=0;
		atom_curr[i].delr[0]=atom_curr[i].delr[1]=atom_curr[i].delr[2]=atom_curr[i].delr[3]=0.0;
	}

	static double prev_sign=1;
	double sign_to_consider=0;
	double mag_sign_to_consider=0;
//	find_atoms_to_consider_formobility(atom_ref,atom_curr,n, consider_formobility);
/*
	atomic_dat *atom_product;
			for(int count=0;count<2;count++)
			{
				if(count==0) atom_product = atom_ref;
				if(count==1) atom_product = atom_curr;
				int start_count=0;
				for(int i=0;i<n;i++)
				{
					if((atom_product[i].CNA>6)&&(atom_product[i].interface==1)&&(atom_product[i].type==1))
					{

						consider_formobility[i] = 1;
						start_count++;
					}
				}

				cout << start_count<<"\t"<<"start_count is\n";
				recursive_neighbors(tag, atom_product,n, consider_formobility);
			}
*/
	neighbors_check(atom_curr,atom_ref,n,consider_formobility);



	//assuming no change in box
	double r1[3],s1[3],r2[3],s2[3],s12[3],r12[3];
	for(int i=0;i<3;i++){r1[i]=0.0;s1[i]=0.0;r2[i]=0.0;s2[i]=0.0;s12[i]=0.0;r12[i]=0.0;}
	double r13[3];
	//double r1rel[3],s1rel[3],r2rel[3],s2rel[3],s12rel[3],r12rel[3];
	//for(int i=0;i<3;i++){r1[i]=0.0;s1[i]=0.0;r2[i]=0.0;s2[i]=0.0;s12[i]=0.0;r12[i]=0.0;}

	double square_d=0;
	double defect_square_d;
	double defect_temp_d[3];
	for(int kk=0;kk<3;kk++) defect_temp_d[kk]=0;
	double count_square_d=0;
	int abcd=0;
	for(int i=0;i<n;i++)
	{
		if(consider_formobility[i]>0)
		{
			abcd++;
			//atom_curr[i].ackN = -15;
		}
	}
	cout << abcd<<"\t end count is \n";
	/*
	char fname[80]="",str[20];
		strcat(fname,filename_save);
		strcat(fname,"_indv_");
		sprintf(str,"%d",tag_reference);
		strcat(fname,str);
		strcat(fname,"_");
		sprintf(str,"%d",tag_target);
		strcat(fname,str);
		cout << fname<<"\n";
		FILE *fptr_indv;
		fptr_indv=fopen(fname,"w");
		cout << "opened file to write\n";
*/
	FILE *fptr;
	char append_val[6]="";
	if(append)
		{
			strcat(append_val,"a");
		}else
		{
			strcat(append_val,"w");
		}

	fptr=fopen(filename_save,append_val);
	if(relative)
	{
		fprintf(fptr," \n\nRelative Displacements \n");
	}else
	{
		fprintf(fptr," \n\nAbsolute Displacements  \n");
	}

	for(int i=0;i<n;i++)
	{
		atom_curr[i].delr[3] = (double) consider_formobility[i]*1.0;
		atom_curr[i].BV =-1;

		if((consider_formobility[i]>0)&&(atom_curr[i].interface==1)&&(atom_curr[i].type==1))
		{
			if(i==20752)
			{
				cout << "HOW DID THIS COME HERE\n\n\n\n\n\n";
			}
			s1[0]=atom_ref[i].sx;s1[1]=atom_ref[i].sy;s1[2]=atom_ref[i].sz;
			s2[0]=atom_curr[i].sx;s2[1]=atom_curr[i].sy;s2[2]=atom_curr[i].sz;
			for(int kk=0;kk<3;kk++)
			{
				s12[kk]=s1[kk]-s2[kk];
				s12[kk] = s12[kk]-(int)(s12[kk]*2);
			}

			double H_t_now[3][3];
			M3mul(H_curr,KS1_trans_mas,H_t_now);

			/*for(int j=0;j<3;j++)
			{
				cout << H_t_now[j][0]<<"\t"<<H_t_now[j][1]<<"\t"<<H_t_now[j][2]<<"\n";
			}*/

			V3mulM3(s1,Href,r1);
			V3mulM3(s2,H_curr,r2);
			V3mulM3(s12,H_t_now,r12);


				if(relative)
				{
					find_relative_distance(atom_curr,atom_ref,n,i,r12,H_t_now);
				}else
				{
					find_relative_distance(atom_curr,atom_ref,n,i,r13,H_t_now);

				}


			if(i==20333)
			{
				cout << "before"<<"\t";
				cout << r12[0]<<"\t"<<r12[1]<<"\t"<<r12[2]<<"\n";

			}
			double temp_square_d=0;

			//temp_square_d += (atom_ref[i].ux-atom_curr[i].ux)*(atom_ref[i].ux-atom_curr[i].ux);
			//temp_square_d += (atom_ref[i].uy-atom_curr[i].uy)*(atom_ref[i].uy-atom_curr[i].uy);
			//temp_square_d += (atom_ref[i].uz-atom_curr[i].uz)*(atom_ref[i].uz-atom_curr[i].uz);

			//project the displacements onto the [01-1] plane
			double project_p[3][3],project_p_inv[3][3];
			project_p[0][0]=0.0;project_p[0][1]=1/sqrt(2);project_p[0][2]=-1/sqrt(2);
			project_p[1][0]=-2/sqrt(6);project_p[1][1]=1/sqrt(6);project_p[1][2]=1/sqrt(6);

			V3cross(project_p[0],project_p[1],project_p[2]);
		//	cout <<project_p[2][0]<<"\t"<<project_p[2][1]<<"\t"<<project_p[2][2]<<"\n";
			M3inv(project_p,project_p_inv);

			//temp_square_d += (r1[0]-r2[0])*(r1[0]-r2[0]);
			//temp_square_d += (r1[1]-r2[1])*(r1[1]-r2[1]);
			//temp_square_d += (r1[2]-r2[2])*(r1[2]-r2[2]);

			double v[3];
			double v13[3];

			V3mulM3(r12,H0_geo,v);
			V3mulM3(v,project_p_inv,r12);

			V3mulM3(r13,H0_geo,v13);
			V3mulM3(v13,project_p_inv,r13);

			for(int kk=0;kk<1;kk++) temp_square_d +=r12[kk]*r12[kk];

			atom_curr[i].drig = r12[0];
			if(i==20333)
			{
				cout << sqrt(temp_square_d)<<"\n";
				cout << "after"<<"\t";
				cout << r12[0]<<"\t"<<r12[1]<<"\t"<<r12[2]<<"\n";

			}
			/*
			if(i==39119)
			{
				cout << atom_ref[i].ux<<"\t"<<atom_ref[i].uy<<"\t"<<atom_ref[i].uz<<"\n";

				cout << atom_ref[i].sx<<"\t"<<atom_ref[i].sy<<"\t"<<atom_ref[i].sz<<"\n";

				cout << atom_curr[i].ux<<"\t"<<atom_curr[i].uy<<"\t"<<atom_curr[i].uz<<"\n";
				cout << atom_curr[i].sx<<"\t"<<atom_curr[i].sy<<"\t"<<atom_curr[i].sz<<"\n";
			}*/
			if(((temp_square_d > 0.5)&&(atom_curr[i].ackN>3))||(temp_square_d > 0.64))
			{
				for(int kk=0;kk<3;kk++)atom_curr[i].delr[kk]=r12[kk];
				atom_curr[i].delr[3] = 1;
				square_d+=temp_square_d;//*atom_ref[i].ma;
				count_square_d +=1;//atom_ref[i].ma;

				cout << i<<"\t"<<count_square_d<< " HEREREREHERER\t"<< temp_square_d<<"\t";
				for(int kk=0;kk<1;kk++) {defect_temp_d[kk]+=fabs(r12[kk]);}//cout <<i<<"\t"<< kk<<"\t"<<defect_temp_d[kk]<<"\t"<<r12[kk]<<"\n";}
				if(r12[0]>0)
				{
					sign_to_consider++;
				}else
				{
					sign_to_consider--;
				}
				mag_sign_to_consider+=r12[0];
				double dis = sqrt(r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2]);
				cout << defect_temp_d[0]<<"\t"<<atom_curr[i].ackN<<"\n";
				fprintf(fptr,"%d %d %d %d %lf %lf %lf %lf %d\n",tag_reference, tag_target,(int)count_square_d,i,sqrt(temp_square_d),dis,r12[0],defect_temp_d[0],atom_curr[i].ackN);


			}

		}
	}
	square_d=0.0;
	count_square_d=0;
	if(mag_sign_to_consider*sign_to_consider<0)
	{
		if(fabs(mag_sign_to_consider)>0.9) sign_to_consider = -1*sign_to_consider;
	}
	if(sign_to_consider==0)
	{
		sign_to_consider = mag_sign_to_consider/(fabs(mag_sign_to_consider));
	}

	for(int kk=0;kk<3;kk++) {defect_temp_d[kk]=0.0;}
	for(int i=0;i<n;i++)
	{
		if((atom_curr[i].delr[3]>0.1)&&(atom_curr[i].delr[0]*sign_to_consider)>0)
		{
			square_d+=atom_curr[i].delr[0]*atom_curr[i].delr[0];
			for(int kk=0;kk<3;kk++) {defect_temp_d[kk]+=atom_curr[i].delr[kk];cout <<i<<"\t"<< kk<<"\t"<<defect_temp_d[kk]<<"\t"<<r12[kk]<<"\n";}
			count_square_d++;
			atom_curr[i].BV = count_square_d;
			atom_curr[i].CNA = 12;
		}
	}
	defect_square_d=0.0;
	for(int kk=0;kk<3;kk++)
	{
		defect_square_d+=defect_temp_d[kk]*defect_temp_d[kk];
	}
	//defect_square_d = sqrt(defect_temp_d);
	int sign = 1.0;
	if(sign_to_consider<0) sign= -1.0;

	fprintf(fptr,"\n%d %d %lf %lf %lf %lf %d %d\n",tag_reference, tag_target,square_d,count_square_d,sqrt(defect_square_d),sqrt(square_d/count_square_d),abcd,sign);
/*
	if(relative)
	{
		fprintf(fptr," Relative Displacements End \n\n");
	}else
	{
		fprintf(fptr," Absolute Displacements End \n\n");
	}*/
	fclose(fptr);
}

void determine_mobility(atomic_dat *atom_ref,atomic_dat *atom_curr,int n,double Href[3][3],double H_curr[3][3],char *filename_save, int tag_reference,int tag_target, bool append)
{
//	determine_mobility1(atom_ref,atom_curr, n,Href, H_curr,filename_save,  tag_reference, tag_target,  append,true);
	determine_mobility1(atom_ref,atom_curr, n,Href, H_curr,filename_save,  tag_reference, tag_target,  append,false);



}
