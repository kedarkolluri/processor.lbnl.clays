/*
 * Disrigistry.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: kedar
 */

#include "Disrigistry.h"

int ref_vector_n = 0;
double **ref_vector_list;
int closest_atom_mode=1;
int bias_list_1=-1;
int bias_list_2=-1;

void set_angle_pattern (double *angle)
{
	if(*angle>-3)
	{
		if((*angle<=15)||(*angle>165)) *angle = 1;
		if((*angle>15)||(*angle<=45)) *angle =2;
		if((*angle>45)||(*angle<=75)) *angle =3;
		if((*angle>75)||(*angle<=105)) *angle =4;
		if((*angle>105)||(*angle<=135)) *angle =5;
		if((*angle>135)||(*angle<=165)) *angle =6;
	}

}

void compute_tensor_mean_square(double **reference, double **target,int rows,int columns, double **dydx)
{
	//cout << rows<<"\n";
	int max_arr_size=20;
	int dimension_t = columns;
	bool pad = false;
//	if(columns==2)
//		{
//			columns =3;
//			pad = true;
//		}
	int rejected = 0;
	int processed = 0;
	double **y,**x;
	double **dydx_here;

	y = (double **) malloc(rows*sizeof(double));
	x = (double **) malloc(rows*sizeof(double));
	dydx_here = (double **) malloc(columns*sizeof(double));
	for(int i=0;i<rows;i++)
	{
		y[i] = (double *) malloc(columns*sizeof(double));
		x[i] = (double *) malloc(columns*sizeof(double));
		if(i<columns) dydx_here[i] = (double *) malloc(columns*sizeof(double));
	}
	//cout << "set\n";
	for (int i =0;i<rows;i++)
	{
		for (int j=0;j<columns;j++)
		{
			if(j<dimension_t)
			{
			x[i][j] = target[i][j];
			y[i][j] = reference[i][j];

			}else
			{
				x[i][j] = y[i][j] =1.0;
			}
			if(i<columns)dydx_here[i][j]=0.0;
		}
	}
	int r =0;
	//cout << "here r is 0\n";
	if(columns==3)
	{
		//cout << "r value is in 3\n";
		r = getderivative3(y,x, rows, dydx_here);
		//cout << r << "\n";
	 }else if(columns==2)
	 {
		//cout << "r value is \n";
		r = getderivative2(y,x, rows, dydx_here);
		//cout << r << "\n";
	 }


		for(int i=0;i<dimension_t;i++)
		{
			for(int j=0;j<dimension_t;j++)
			{
				dydx[i][j] = dydx_here[i][j];
			}
		}

//	cout << " at the end before memory free\n";
	for(int k=0;k<rows;k++)
	{
			free(x[k]);
			free(y[k]);
			if(k<columns) free(dydx_here[k]);
	}
	free(x);free(y);free(dydx_here);

}

void compute_diff_ten(double **diff_ten,int rows,int columns )
{
	cout << rows<<"\n";
	int max_arr_size=20;
	int dimension_t = columns;
	int rejected = 0;
	int processed = 0;
	FILE *fptr2;
	fptr2 = fopen("test.dis","w");
	fclose(fptr2);
	FILE *fptr3;
	fptr3 = fopen("aux.dis","w");
	fclose(fptr3);
	for(int i=0;i<n;i++)
	{
		if(atom[i].interface>0)
		{
			//cout <<"here again\t"<<i<<"\n";
			int collate_counter=0;
			double **y,**x,**dydx;
			y = (double **) malloc(max_arr_size*sizeof(double));
			x = (double **) malloc(max_arr_size*sizeof(double));
			dydx = (double **) malloc(max_arr_size*sizeof(double));
			for(int k=0;k<dimension_t;k++)
			{
				dydx[k] = (double *) malloc(dimension_t*sizeof(double));
				for(int dd=0;dd<dimension_t;dd++)
				{
					dydx[k][dd] = 0;
				}
			}
			for(int k=0;k<max_arr_size;k++)
			{
				y[k] = (double *) malloc(dimension_t*sizeof(double));
				x[k] = (double *) malloc(dimension_t*sizeof(double));
				for(int dd=0;dd<dimension_t;dd++)
				{
					x[k][dd] = 0.0;
					y[k][dd] = 0.0;
				}

			}
			int diff_ten_counter=0;
			int sign = 1;
			//cout <<"here again\t"<<i<<"\n";
			while((diff_ten_counter<rows)&&(collate_counter<max_arr_size))
			{
				bool collate=false;
				if((int)diff_ten[diff_ten_counter][0]==i)
				{

					sign = 1.0;
					collate = true;
				}else if((int) diff_ten[diff_ten_counter][1]==i)
				{
					//sign=-1.0;
					collate = true;
				}else
				{

					double r1[3],s1[3],s2[3],s3[3],r2[3],r3[3],s12[3],s13[3];
					s1[0] = atom[i].sx;s1[1] = atom[i].sy;s1[2]=atom[i].sz;
					int c0 = diff_ten[diff_ten_counter][0];
					int c1 = diff_ten[diff_ten_counter][1];

					s2[0] = atom[c0].sx;s2[1] = atom[c0].sy;s2[2]=atom[c0].sz;
					s3[0] = atom[c1].sx;s3[1] = atom[c1].sy;s3[2]=atom[c1].sz;

					V3mulM3(s1,Hcry,r1);

					for(int h1=0;h1<3;h1++)
					{
						s12[h1] = s1[h1]-s2[h1];
						s12[h1] = s12[h1]-(int)(s12[h1]*2);

						s13[h1] = s1[h1]-s3[h1];
						s13[h1] = s13[h1]-(int)(s13[h1]*2);
					}
					V3mulM3(s12,Hcry,r2);
					V3mulM3(s13,Hcry,r3);


					double distance=0.0, distance1 = 0.0,distance2=0.0;
					for(int d=0;d<dimension_t;d++)
					{
						distance+=(r1[d]-diff_ten[diff_ten_counter][2+d])*(r1[d]-diff_ten[diff_ten_counter][2+d]);
						distance1+= r2[d]*r2[d];
						distance2+= r3[d]*r3[d];

					}
					if((sqrt(distance2)<3.1)&&(sqrt(distance1)<3.1))
					{
						collate = true;
//						cout <<i<<"\t"<<sqrt(distance)<<"\t"<<sqrt(distance1)<<"\t"<<sqrt(distance2)<<"\n";
//						if((c0==i)||(c1==i)) cout <<i<<"\t"<<c0<<"\t"<<c1<<" HH\n";
					}else
					{
						diff_ten_counter++;
					}
					//diff_ten_counter++;
				}

				if((collate)&&(collate_counter<max_arr_size))
				{
					for(int d=0;d<dimension_t;d++)
					{
						x[collate_counter][d] = diff_ten[diff_ten_counter][2+d];
						y[collate_counter][d] = diff_ten[diff_ten_counter][5+d];//*sign;

					}

					diff_ten_counter++;
					collate_counter++;
				//	cout <<"here again\t"<<i<<"\t"<<collate_counter<<"\n";
				}
			}

			int r=0;
			if(dimension_t==3)
			{
				r = getderivative3(y,x, max_arr_size, dydx);
			 }else if(dimension_t==2)
			 {
				r = getderivative2(y,x, max_arr_size, dydx);
			 }

			if((i==27054)||(i==25331))
			{
				for(int h=0;h<atom[i].coord;h++)
				{

					for(int g=0;g<dimension_t;g++)
					{
						cout <<y[h][g]<<"\t";
					}
					cout <<";\n";
				}
				cout <<"\n";

				for(int h=0;h<atom[i].coord;h++)
				{

					for(int g=0;g<dimension_t;g++)
					{
						cout <<x[h][g]<<"\t";
					}
					cout <<";\n";
				}

				cout <<"\n";
				for(int h=0;h<dimension_t;h++)
				{
					for(int g=0;g<dimension_t;g++)
					{
						cout << dydx[h][g]<<"\t";
					}
					cout <<"\n";
				}
				cout <<"\n\n";
			}

			atom[i].delr[0] = atom[i].delr[1]=atom[i].delr[2]=atom[i].delr[3]=0.0;
		//	atom[i].delr[0] = r*1.0;

			if(r==0)
			{

				for(int h=0;h<dimension_t-1;h++)
				{
					for(int g=0;g<dimension_t-1;g++)
					{
				//		cout << dydx[h][g]<<"\t";
						atom[i].delr[0] += 0.5*dydx[h][g]*dydx[h][g];
					}
				//	cout <<"\n";
				}


				atom[i].delr[0] = sqrt(atom[i].delr[0]);
				if((i==27054)||(i==25331)) cout << "norm is \t"<<atom[i].delr[0]<<"\n";
			//	if(atom[i].delr[0]>0.5) atom[i].delr[0] = 0;
			//	if(atom[i].delr[0]<0.1) atom[i].delr[0] = 0;
			//	if(atom[i].type==2) atom[i].delr[0] = -1;
				processed++;
			}else
			{
				//cout <<i<<"\t"<<atom[i].CNA<<"\t"<<atom[i].type<<"\n";
				atom[i].delr[0]=-1;
				rejected++;
			}
			/*
			if(atom[i].interface!=0)
			{
				fptr2 = fopen("test.dis","a");
				fprintf(fptr2,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",i,atom[i].sx,atom[i].sy,atom[i].sz,atom[i].delr[0],dydx[0][0],dydx[0][1],dydx[0][2],dydx[1][0],dydx[1][1],dydx[1][2],dydx[2][0],dydx[2][1],dydx[2][2],atom[i].CNA);
				fclose(fptr2);
			}
			if((atom[i].interface==1)&&(atom[i].type==1))
			{
				fptr3 = fopen("aux.dis","a");
				fprintf(fptr3,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",i,dydx[0][0],dydx[0][1],dydx[0][2],dydx[1][0],dydx[1][1],dydx[1][2],dydx[2][0],dydx[2][1],dydx[2][2],atom[i].CNA);
				fclose(fptr3);
			}
*/
			for(int k=0;k<dimension_t;k++)
			{
					free(dydx[k]);

			}
			for(int k=0;k<max_arr_size;k++)
			{
					free(x[k]);
					free(y[k]);
			}
			free(x);free(y);free(dydx);
		}
	}

	cout << "total atoms reject are\t"<< rejected<<"\t and processed are\t"<<processed<<"\n";
}

int determine_correspondence(int ref_type,double *rij)
{
	int start = 0;
	for(int i=0;i<ref_type;i++)
	{
		start = start+basic_n[i];
	}

	int ref_id=-1;
	double degree_deviation_prev = 180000;
	double deviation_sign = 1;
	for(int i=start;i<(start+basic_n[ref_type]-1);i++)
	{
		double ref_r[3];
		ref_r[0] = ref_vectors[i][0];ref_r[1]=ref_vectors[i][1];ref_r[2]=ref_vectors[i][2];
		double deviation = V3dot(ref_r,rij);
		double deviation_degrees = acos(fabs(deviation))*180/PI;
		if(deviation_degrees<degree_deviation_prev)
		{
			degree_deviation_prev = deviation_degrees;
			ref_id = i;
			if (deviation<0) deviation_sign = -1;
		}
	}
	rij[0] = rij[0]*deviation_sign; rij[1] = rij[1]*deviation_sign;rij[2] = rij[2]*deviation_sign;

	if(degree_deviation_prev>30)
	{
		count_disrig_dev++;
		cout <<"\t and least angle is\t"<<degree_deviation_prev<<"\t"<<count_disrig_dev<<"\n";
		//ref_id = -1;
	}

	return(ref_id);
}

void determine_correspondence_3(atomic_dat *atom_now, int n_now, double H_now[3][3],int atom_i,int atom_j,int *corr_vect_ref,int *corr_vect,double **corr_vect_ref_d,double **corr_vect_d, double distt[3])
{
	double angles_list[6][6];
	int angles_list_single[6];
	for(int i=0;i<6;i++)
	{
		angles_list_single[i]=-1;
		for(int j=0;j<6;j++)
		{
			angles_list[i][j]=180000000;
		}
	}
	double sij_0[3],rij_0[3];
	sij_0[0]=atom_now[atom_i].sx-atom_now[atom_j].sx;
	sij_0[1]=atom_now[atom_i].sy-atom_now[atom_j].sy;
	sij_0[2]=0;//atom_now[atom_i].sz-atom_now[atom_j].sz;
	sij_0[0] = sij_0[0]-(int)(sij_0[0]*2);
	sij_0[1] = sij_0[1]-(int)(sij_0[1]*2);
	sij_0[2] = sij_0[2]-(int)(sij_0[2]*2);
	V3mulM3(sij_0,H_now,rij_0);
	distt[0]= rij_0[0];
	distt[1] = rij_0[1];
	distt[2] = rij_0[2];
//	cout << atom_i<<"\t"<<atom_j<<"\t"<<rij_0[0]<<"\t"<<rij_0[1]<<"\t"<<rij_0[2]<<"\n";
	for (int i=0;i<6;i++)
	{
		int lowest_angle=180000;
		for(int j=0;j<6;j++)
		{

			int j_ref = corr_vect_ref[i];
			int j_main = corr_vect[j];
		//	cout << "I am starting fine\t"<<j_ref<<"\t"<<j_main<<"\n";
			double sij[3],rij_ref[3],rij_main[3];
			sij[0] = atom_now[atom_i].sx- atom_now[j_ref].sx;
			sij[1] = atom_now[atom_i].sy- atom_now[j_ref].sy;
			sij[2] = 0;//atom_now[atom_i].sz- atom_now[j_ref].sz;
			sij[0] = sij[0]-(int)(sij[0]*2);
			sij[1] = sij[1]-(int)(sij[1]*2);
			sij[2] = sij[2]-(int)(sij[2]*2);
			V3mulM3(sij,H_now,rij_ref);
			corr_vect_ref_d[i][0] = rij_ref[0];
			corr_vect_ref_d[i][1] = rij_ref[1];
			corr_vect_ref_d[i][2] = rij_ref[2];
			sij[0] = atom_now[atom_j].sx- atom_now[j_main].sx;
			sij[1] = atom_now[atom_j].sy- atom_now[j_main].sy;
			sij[2] = 0;//atom_now[atom_i].sz- atom_now[j_ref].sz;
			sij[0] = sij[0]-(int)(sij[0]*2);
			sij[1] = sij[1]-(int)(sij[1]*2);
			sij[2] = sij[2]-(int)(sij[2]*2);
			V3mulM3(sij,H_now,rij_main);
			//angles_list[i][j] = V3dot(ref_ref,rij_main);
			//cout << "before deviation determination ";
			double deviation = V3dot(rij_ref,rij_main);
			angles_list[i][j] = acos((deviation))*180/PI;

			//if(atom_i==7059)
		//	cout << j_ref<<"\t"<<j_main<<"\t"<<angles_list[i][j]<<"\n";

			if(angles_list[i][j]<lowest_angle)
			{
				angles_list_single[i]=j_main;
				lowest_angle = angles_list[i][j];
				corr_vect_d[i][0] = rij_main[0];//+rij_0[0];
				corr_vect_d[i][1] = rij_main[1];//+rij_0[1];
				corr_vect_d[i][2] = rij_main[2];//+rij_0[2];
			}

		}
//		cout << "finally it is\t"<<angles_list_single[i]<<"\n";
	}

	for(int i=0;i<6;i++)
	{
		corr_vect[i]=angles_list_single[i];
//		cout <<corr_vect_ref_d[i][0]<<"\t"<<corr_vect_ref_d[i][1]<<"\t"<<corr_vect_ref_d[i][2]<<"\t\t";
//		cout <<corr_vect_d[i][0]<<"\t"<<corr_vect_d[i][1]<<"\t"<<corr_vect_d[i][2]<<"\n";
	}

	//compare the angles to determine the smallest sets of angles
}

void determine_correspondence_2(atomic_dat *atom_now, int n_now, double H_now[3][3],int atom_i,int *corr_vect_ref,int *corr_vect)
{
	double angles_list[6][6];
	int angles_list_single[6];
	for(int i=0;i<6;i++)
	{
		angles_list_single[i]=-1;
		for(int j=0;j<6;j++)
		{
			angles_list[i][j]=180000000;
		}
	}

	for (int i=0;i<6;i++)
	{
		int lowest_angle=180000;
		for(int j=0;j<6;j++)
		{
			int j_ref = corr_vect_ref[i];
			int j_main = corr_vect[j];
			double sij[3],rij_ref[3],rij_main[3];
			sij[0] = atom_now[atom_i].sx- atom_now[j_ref].sx;
			sij[1] = atom_now[atom_i].sy- atom_now[j_ref].sy;
			sij[2] = 0;//atom_now[atom_i].sz- atom_now[j_ref].sz;
			sij[0] = sij[0]-(int)(sij[0]*2);
			sij[1] = sij[1]-(int)(sij[1]*2);
			sij[2] = sij[2]-(int)(sij[2]*2);
			V3mulM3(sij,H_now,rij_ref);

			sij[0] = atom_now[atom_i].sx- atom_now[j_main].sx;
			sij[1] = atom_now[atom_i].sy- atom_now[j_main].sy;
			sij[2] = 0;//atom_now[atom_i].sz- atom_now[j_ref].sz;
			sij[0] = sij[0]-(int)(sij[0]*2);
			sij[1] = sij[1]-(int)(sij[1]*2);
			sij[2] = sij[2]-(int)(sij[2]*2);
			V3mulM3(sij,H_now,rij_main);
			//angles_list[i][j] = V3dot(ref_ref,rij_main);
			double deviation = V3dot(rij_ref,rij_main);
			angles_list[i][j] = acos(fabs(deviation))*180/PI;
			//cout << j_ref<<"\t"<<j_main<<"\t"<<angles_list[i][j]<<"\n";

			if(angles_list[i][j]<lowest_angle)
			{
				angles_list_single[i]=j_main;
				lowest_angle = angles_list[i][j];
			}

		}
		//cout << "finally it is\t"<<angles_list_single[i]<<"\n";
	}

	for(int i=0;i<6;i++)
	{
		corr_vect[i]=angles_list_single[i];
	}

	//compare the angles to determine the smallest sets of angles
}

void update_closest(int ii, int jj, atomic_dat *atom_now, double H_now[3][3],int *id, double *dist )
{
	double sij[3],rij[3];
	double rijsq=0;
	sij[0] = atom_now[ii].sx- atom_now[jj].sx;
	sij[1] = atom_now[ii].sy- atom_now[jj].sy;
	sij[2]=0;
//	sij[2] = atom_now[ii].sz- atom_now[jj].sz;
	sij[0] = sij[0]-(int)(sij[0]*2);
	sij[1] = sij[1]-(int)(sij[1]*2);
	sij[2] = sij[2]-(int)(sij[2]*2);

	if((closest_atom_mode!=1)&&(bias_list_1>-1)&&(bias_list_2>-1))
	{
		double s_xt = (atom_now[bias_list_1].sx- atom_now[bias_list_2].sx);
		double s_yt = (atom_now[bias_list_1].sy- atom_now[bias_list_2].sy);
		s_xt = s_xt-(int)(s_xt*2);
		s_yt = s_yt-(int)(s_yt*2);
		sij[0]=sij[0]+s_xt;
		sij[1]=sij[1]+s_yt;
		sij[0] = sij[0]-(int)(sij[0]*2);
		sij[1] = sij[1]-(int)(sij[1]*2);

	}

	V3mulM3(sij,H_now,rij);
	rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
	if(rijsq < *dist){ *id =jj; *dist = rijsq;}
}

void compute_disregistry_3(atomic_dat *atom_now, int n_now, double H_now[3][3], int suffix,int s_atom_type,int s_atom_interface)
{
	// find the closest bcc atom to a given fcc atom
	// find 6 in-plane nearest bcc atoms to this bcc atom
	// and 6 in-plane fcc atoms to the chose fcc atom
	// (a) 6 in-plane copper neighbors to the considered copper atom
	// find disregistry between the 6 cu-nb bonds and 6 cu-cu bonds are finding correspondence
	// save disregistry

	// for now, it is assumed that the atoms coord_id thing will contain the desired neighbors
	// this has to be accomplished by changing the cut-off radius so that we will get it


	char filename1[80]="disregistry_data_3_",str[80];
	char filename2[80]="mean_disregistry_data_3_";
	char filename3[80]="mean_disregistry_data_3.diff_";
	sprintf(str,"%d",suffix);
	strcat(filename1,str);
	strcat(filename2,str);
	strcat(filename3,str);
		FILE *fptr;
		fptr = fopen(filename1,"w");
		FILE *fptr2;
		fptr2=fopen(filename2,"w");

		double dist[3];

		double **corr_vect_ref_d;
		double **corr_vect_d;
		double **dydx;
		dydx=(double **) malloc(2*sizeof(double));
		for(int i=0;i<2;i++)
		{
			dydx[i] = (double *) malloc(2*sizeof(double));
		}

		corr_vect_ref_d = (double **) malloc(6*sizeof(double));
		corr_vect_d = (double **) malloc(6*sizeof(double));
		for(int i=0;i<6;i++)
		{
			corr_vect_ref_d[i]  = (double *) malloc(3*sizeof(double));
			corr_vect_d[i]  = (double *) malloc(3*sizeof(double));
			for(int j=0;j<3;j++)
			{
				corr_vect_ref_d[i][j] = -1000000000.0;
				corr_vect_d[i][j] = -1000000000.0;
			}
		}

		int slgi_counter=0;
		for(int i=0;i<n_now;i++)
		{
			atom_now[i].disrigistry[0] = -1e6;
			atom_now[i].disrigistry[1] = -1e6;
			atom_now[i].disrigistry[2] = -1e6;
			atom_now[i].delr[0] = -1e6;
			atom_now[i].delr[1]  = -1e6;
			atom_now[i].delr[2]  = -1e6;
			atom_now[i].delr[3]  = -1e6;
			atom_now[i].drig = -1;

			int corr_vect_ref[6];
			int corr_vect[6];
			bool discard = false;
			for(int ii=0;ii<6;ii++)
			{
				corr_vect_ref[ii]=-1;
				corr_vect[ii]=-1;
			}
			int c_atom=i;
			int cc_atom=i;
			if((atom_now[i].interface==1)&&(atom_now[i].type==1))
			{
				// determine the closest type=2 atom and at the same time fill corr_vect_ref
				double nbnbr_init_dist = 10000000;
				int nbnbr_init_id = -1;
				int count_inplane_same=0;
				for(int jnab=0;jnab<MAX_COORD;jnab++)
				{
					int j = atom_now[i].coord_id[jnab];
					if(j>-1)
					{
						// have to change here if we want to do between cu and cu_alpha
						if((atom_now[j].type==s_atom_type)&&(atom_now[j].interface==s_atom_interface))
						{
							//cout << "before\t"<<i<<"\t"<<j<<"\n";
							update_closest(i,j,atom_now, H_now,&nbnbr_init_id,&nbnbr_init_dist);
							//cout <<"after\t"<<i<<"\t"<<j<<"\n";
						}else if ((atom_now[j].interface==1)&&(atom_now[j].type==1))
						{
							if(count_inplane_same<6)
							{
							corr_vect_ref[count_inplane_same]=j;
							count_inplane_same++;
							}else
							{
							//	cout << "more than 6 neighbors in cu, exiting"<<"\n";
								discard = true;
								//exit(1);
							}
						}
					}
				}
				if(count_inplane_same<6)discard = true;
				// fills the nb neighbors to the closest nb atom
				count_inplane_same =0;
				for(int jnab =0; jnab<MAX_COORD;jnab++)
				{
					int j = atom_now[nbnbr_init_id].coord_id[jnab];
					if(j>-1)
					{
						// have to change here if we want to do between cu and cu_alpha
						if((atom_now[j].type==s_atom_type)&&(atom_now[j].interface==s_atom_interface))
						{
							if(count_inplane_same<6)
							{
							corr_vect[count_inplane_same]=j;
							count_inplane_same++;
							}else
							{
							//	cout << "more than 6 neighbors in Nb, exiting "<<nbnbr_init_id<<"\n";
								discard = true; //exit(1);
							}
						}
					}
				}
				if(count_inplane_same<6)discard = true;
				//cout <<"hi\n";
				//if(discard) cout << "decided to discard\n";
				if(!discard)
				{
					//cout << "before determining the correspondence\n";
					determine_correspondence_3(atom_now, n_now, H_now,i,nbnbr_init_id,corr_vect_ref,corr_vect,corr_vect_ref_d,corr_vect_d,dist);
					//cout << "determined the correspondence\n";
					// currently no checks are performed to see if the correspondence
					// is good. Have to add checks but that is a training


					compute_tensor_mean_square(corr_vect_ref_d,corr_vect_d, 6, 2,dydx);
				//	cout << "computed the tensor\n";
					if(i==13229)
					{
						for(int j_j=0;j_j<6;j_j++)
						{
							for(int j_k=0;j_k<2;j_k++)
							{
								cout << corr_vect_ref_d[j_j][j_k]<<" ";
							}
							cout <<"; ";
						}
						cout <<"\n";
						for(int j_j=0;j_j<6;j_j++)
						{
							for(int j_k=0;j_k<2;j_k++)
							{
								cout << corr_vect_d[j_j][j_k]<<" ";
							}
							cout <<"; ";
						}


						cout << dydx[0][0]<<" "<<dydx[0][1]<<" "<<dydx[1][0]<<" "<<dydx[1][1];
					}
					double sij_mean[3];
					sij_mean[0]=sij_mean[1]=sij_mean[2]=0.0;
					int av_nm=0;
					for(int i_here=0;i_here<6;i_here++)
					{
						if((corr_vect[i_here]>-1)&&(corr_vect_ref[i_here]>-1))
						{

							double sij_now[3],rij_now[3],rijsq_now;
							double sij_ref[3],rij_ref[3],rijsq_ref;
							double sxmid_now, symid_now,szmid_now,smid_now[3],rmid_now[3];
							// Added the difference between the centers ---CHECK
							//cout << "****************** Added the difference between the centers ---CHECK\n";
							//-(atom_now[nbnbr_init_id].sx-atom_now[i].sx)
							sij_now[0] = atom_now[corr_vect[i_here]].sx-atom_now[corr_vect_ref[i_here]].sx;
							sij_now[1] = atom_now[corr_vect[i_here]].sy-atom_now[corr_vect_ref[i_here]].sy;
							sij_now[2] = 0.0;//atom_now[corr_vect[i]].sz-atom_now[corr_vect_ref[i]].sz;
							for(int kk=0;kk<3;kk++) sij_now[kk] = sij_now[kk]-(int)(sij_now[kk]*2);

							sij_ref[0] = atom_now[i].sx - atom_now[corr_vect_ref[i_here]].sx;
							sij_ref[1] = atom_now[i].sy - atom_now[corr_vect_ref[i_here]].sy;
							sij_ref[2] = atom_now[i].sz - atom_now[corr_vect_ref[i_here]].sz;

							sij_ref[0] = sij_ref[0]-(int)(sij_ref[0]*2);
							sij_ref[1] = sij_ref[1]-(int)(sij_ref[1]*2);
							sij_ref[2] = sij_ref[2]-(int)(sij_ref[2]*2);

							smid_now[0]=atom_now[i].sx+sij_ref[0]/2.0;
							smid_now[1]=atom_now[i].sy+sij_ref[1]/2.0;
							smid_now[2]=atom_now[i].sz+sij_ref[2]/2.0;
							for(int ii2=0;ii2<3;ii2++)
							{
								if(smid_now[ii2]>1) smid_now[ii2]=smid_now[ii2]-1;
							if(smid_now[ii2]<0) smid_now[ii2]=smid_now[ii2]+1;
							}
							V3mulM3(smid_now,H_now,rmid_now);

							V3mulM3(sij_now,H_now,rij_now);
							double dx = (rij_now[0]);
							double dy = (rij_now[1]);
							double dz = (rij_now[2]);

							bool discard1 = false;
							if((rij_ref[0]*rij_ref[0]+rij_ref[1]*rij_ref[1]+rij_ref[2]*rij_ref[2])>0.5*(H_now[0][0]*H_now[0][0]+H_now[1][1]*H_now[1][1]+H_now[2][2]*H_now[2][2]))
							{
								dx=0.0;dy=0.0;dz=0.0;
								discard1 = true;
							}
							if(fabs(dx) > 0.5*(H_now[0][0])) {dz = dy = dx = 0.0; discard = true;}
							if(fabs(dy) > 0.5*(H_now[1][1])) {dz = dy = dx = 0.0; discard = true;}
							if(fabs(dz) > 0.5*(H_now[2][2])) {dz = dy = dx = 0.0; discard = true;}
							if(!discard1)
							{
								av_nm++;
								sij_mean[0]+=sij_now[0];
								sij_mean[1]+=sij_now[1];
								sij_mean[2]+=sij_now[2];
								fprintf(fptr,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,corr_vect[i_here],atom[i].type,atom[corr_vect[i_here]].type,dx,dy,dz,nbnbr_init_dist,atom_now[i].sx,atom_now[i].sy,atom_now[i].sz,smid_now[0],smid_now[1],smid_now[2]);
								//fprintf(fptr,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,corr_vect[i_here],atom[i].type,atom[corr_vect[i_here]].type,corr_vect_ref_d[i_here][0],corr_vect_ref_d[i_here][1],corr_vect_ref_d[i_here][2],corr_vect_d[i_here][0],corr_vect_d[i_here][1],corr_vect_d[i_here][2],nbnbr_init_dist,atom_now[i].sx,atom_now[i].sy,atom_now[i].sz);

							}
						}
					}
				//	cout << "out of almost evreything for a single atom\n";
					double l_strain[2][2];
					double dydx_2[2][2];
					dydx_2[0][0] = dydx[0][0];dydx_2[0][1] = dydx[0][1];dydx_2[1][0] = dydx[1][0];dydx_2[1][1] = dydx[1][1];
					double trace=0;
					double vmises=0;
					M2nms(dydx_2,l_strain,&trace,&vmises);
				//	cout << "just before saving\n";
					if(i==2984)
					{
						cout << "herer \t\t"<<l_strain[0][0]<<"\t"<<l_strain[0][1]<<"\t"<<l_strain[1][0]<<"\t"<<l_strain[1][1]<<"\n";
					}
//					av_nm=1;
					fprintf(fptr2,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,atom_now[i].sx,atom_now[i].sy,atom_now[i].sz,sij_mean[0]/(av_nm*1.0),sij_mean[1]/(av_nm*1.0),sij_mean[2]/(av_nm*1.0),dist[0],dist[1],dist[2],dydx[0][0],dydx[0][1],dydx[1][0],dydx[1][1],trace,vmises);

					atom_now[i].drig = 1;
					atom_now[i].disrigistry[0] = sij_mean[0]/(av_nm*1.0);
					atom_now[i].disrigistry[1] = sij_mean[1]/(av_nm*1.0);
					atom_now[i].disrigistry[2] = sij_mean[2]/(av_nm*1.0);
					atom_now[i].delr[0] = dydx[0][0];
					atom_now[i].delr[1]  = dydx[0][1];
					atom_now[i].delr[2]  = dydx[1][0];
					atom_now[i].delr[3]  = dydx[1][1];



				}

			}

		}
		fclose(fptr);
		fclose(fptr2);
		double **dR;
		double **ddisreg;
//		double **ddisreg10;
//		double **ddisreg11;
//		double **ddisreg12;
//		double **ddisreg22;
		double **dRddisreg;
		dRddisreg=(double **) malloc(2*sizeof(double));
		dR = (double **) malloc(6*sizeof(double));
		ddisreg = (double **) malloc(6*sizeof(double));
//		ddisreg10 = (double **) malloc(6*sizeof(double));
//		ddisreg11 = (double **) malloc(6*sizeof(double));
//		ddisreg12 = (double **) malloc(6*sizeof(double));
//		ddisreg22 = (double **) malloc(6*sizeof(double));
		for(int i=0;i<2;i++)
		{
			dRddisreg[i] = (double *) malloc(2*sizeof(double));
		}


		for(int i=0;i<6;i++)
		{
			dR[i]  = (double *) malloc(3*sizeof(double));
			ddisreg[i]  = (double *) malloc(3*sizeof(double));
			for(int j=0;j<3;j++)
			{
				dR[i][j] = -1000000000.0;
				ddisreg[i][j] = -1000000000.0;
			}
		}

		fptr2=fopen(filename3,"w");
		//cout <<"in the second, start\n";
		for(int i=0;i<n_now;i++)
		{
			double norm=0;
			double direction[3];
			double direction_s[3];
			direction[0] = direction[1]=direction[2]=0.0;
			direction_s[0] = direction_s[1]=direction_s[2]=0.0;
			if((atom_now[i].interface==1)&&(atom_now[i].type==1)&&(atom_now[i].drig==1))
			{
				int coord_c=0;
				for(int jnab=0;jnab<MAX_COORD;jnab++)
				{
					//cout << "jnab value is "<< jnab<<"\t";
					int j = atom_now[i].coord_id[jnab];
					//cout << "and j is "<<j<<"\n";
					if(j>-1)
					{
						//cout <<"since j is not -1, I am writing again\t" <<j<<"\n";
						if((atom_now[j].interface==1)&&(atom_now[j].type==1)&&(atom_now[j].drig==1)&&(coord_c<6))
						{
							//cout << i<<"\t"<<j<<" just entered here\t"<<coord_c<<"\n";
							double sij_now[3],rij_now[3],delsij[3],delrij[3];
							sij_now[0] = atom_now[i].sx-atom_now[j].sx;
							sij_now[1] = atom_now[i].sy-atom_now[j].sy;
							sij_now[2] = 0.0;//atom_now[corr_vect[i]].sz-atom_now[corr_vect_ref[i]].sz;
							for(int kk=0;kk<3;kk++) sij_now[kk] = sij_now[kk]-(int)(sij_now[kk]*2);
							V3mulM3(sij_now,H_now,rij_now);
							delsij[0] = atom_now[i].disrigistry[0]-atom_now[j].disrigistry[0];
							delsij[1] = atom_now[i].disrigistry[1]-atom_now[j].disrigistry[1];
							delsij[2] = atom_now[i].disrigistry[2]-atom_now[j].disrigistry[2];
							V3mulM3(delsij,H_now,delrij);
							direction_s[0]+=delsij[0];direction_s[1]+=delsij[1];
							if(coord_c<6)
							{
							dR[coord_c][0]=rij_now[0];
							dR[coord_c][1]=rij_now[1];
							ddisreg[coord_c][0]=delrij[0];
							ddisreg[coord_c][1] = delrij[1];
							}
							coord_c++;
							double sqr_v = sqrt(delrij[0]*delrij[0]+delrij[1]*delrij[1]);
							norm = norm+sqr_v;
							//direction[0]+=delrij[0]/sqrv;direction[1]+=delrij[1]/sqrv;
							direction[0]+=delrij[0];direction[1]+=delrij[1];
							//if(i==8171)
				//			cout <<j<<"\n"<<rij_now[0]<<"\t"<<rij_now[1]<<";\nd "<<delrij[0]<<"\t"<<delrij[1]<<";\n";
						}
					}
				}
			//	cout << "here?\t"<<coord_c<<"\n";
				if(coord_c!=0)
				{
					//cout <<"here2\n";
					int coord_max_val=coord_c;
					//cout <<"here3\n";
					if(coord_c>6){ coord_max_val = 6;}

					//cout << "final tensor compute\n";
					if(coord_max_val>2)
					{
						//cout << "before\n";
						compute_tensor_mean_square(ddisreg,dR, coord_max_val, 2,dRddisreg);
						//cout <<"and done\n";
					}
					//cout << "abcd\n";
					double indicator = -1;
					double angle = atan(direction[1]/direction[0])*180/PI;
					if(angle<0) angle = angle+180;

					if(((sqrt(direction[1]*direction[1]+direction[0]*direction[0])/(coord_c*1.0))<0.7))
					{
						angle = -30.0;
					}
/*
					if(angle>-3)
						{
							if((angle<=15)||(angle>165)) angle = 1;
							if((angle>15)&&(angle<=45)) angle =2;
							if((angle>45)&&(angle<=75)) angle =3;
							if((angle>75)&&(angle<=105)) angle =4;
							if((angle>105)&&(angle<=135)) angle =5;
							if((angle>135)&&(angle<=165)) angle =6;
						}
*/

					fprintf(fptr2,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,atom_now[i].sx,atom_now[i].sy,atom_now[i].sz,norm/coord_c,sqrt(direction[1]*direction[1]+direction[0]*direction[0])/coord_c,angle,direction_s[0]/(coord_c*1.0),direction_s[1]/(coord_c*1.0),dRddisreg[0][0],dRddisreg[0][1],dRddisreg[1][0],dRddisreg[1][1]);
					//cout <<"finished writing\n";
				}
			}
		}
		fclose(fptr2);

}



void compute_disregistry_2(atomic_dat *atom_now, int n_now, double H_now[3][3])
{
	// find the closest bcc atom to a given fcc atom
	// find 6 in-plane nearest bcc atoms to this bcc atom
	// find correspondence between these 6 atoms and
	// (a) 6 in-plane copper neighbors to the considered copper atom
	// find disregistry between the 6 cu-nb bonds and 6 cu-cu bonds are finding correspondence
	// save disregistry

	// for now, it is assumed that the atoms coord_id thing will contain the desired neighbors
	// this has to be accomplished by changing the cut-off radius so that we will get it

		FILE *fptr;
		fptr = fopen("disregistry_data_2","w");
		FILE *fptr2;
		fptr2=fopen("mean_disregistry_data","w");

		for(int i=0;i<n_now;i++)
		{
			int corr_vect_ref[6];
			int corr_vect[6];
			for(int ii=0;ii<6;ii++)
			{
				corr_vect_ref[ii]=-1;
				corr_vect[ii]=-1;
			}
			int c_atom=i;
			int cc_atom=i;
			if((atom_now[i].interface==1)&&(atom_now[i].type==1))
			{
				// determine the closest type=2 atom and at the same time fill corr_vect_ref
				double nbnbr_init_dist = 10000000;
				int nbnbr_init_id = -1;
				int count_inplane_same=0;
				for(int jnab=0;jnab<MAX_COORD;jnab++)
				{
					int j = atom_now[i].coord_id[jnab];
					if(j>-1)
					{
						if(atom_now[j].type==2)
						{
							update_closest(i,j,atom_now, H_now,&nbnbr_init_id,&nbnbr_init_dist);
						}else if (atom_now[j].interface==1)
						{
							if(count_inplane_same<6)
							{
							corr_vect_ref[count_inplane_same]=j;
							count_inplane_same++;
							}else
							{
								cout << "more than 6 neighbors in cu, exiting";
								exit(1);
							}
						}
					}
				}

				// fills the nb neighbors to the closest nb atom
				count_inplane_same =0;
				for(int jnab =0; jnab<MAX_COORD;jnab++)
				{
					int j = atom_now[nbnbr_init_id].coord_id[jnab];
					if(j>-1)
					{
						if((atom_now[j].type==2)&&(atom_now[j].interface==1))
						{
							if(count_inplane_same<6)
							{
							corr_vect[count_inplane_same]=j;
							count_inplane_same++;
							}else
							{
								cout << "more than 6 neighbors in Nb, exiting";
								exit(1);
							}
						}
					}
				}
				if(i==7059)
				{
					for(int ii3=0;ii3<6;ii3++)
					{
						cout << corr_vect_ref[ii3]<<"\t"<<corr_vect[ii3]<<"\n";
					}
				}
				determine_correspondence_2(atom_now, n_now, H_now,i,corr_vect_ref,corr_vect);
				double sij_mean[3];
				sij_mean[0]=sij_mean[1]=sij_mean[2]=0.0;
				int av_nm=0;
				for(int i_here=0;i_here<6;i_here++)
				{
					if((corr_vect[i_here]>-1)&&(corr_vect_ref[i_here]>-1))
					{

						double sij_now[3],rij_now[3],rijsq_now;
						double sij_ref[3],rij_ref[3],rijsq_ref;
						double sxmid_now, symid_now,szmid_now,smid_now[3],rmid_now[3];
						sij_now[0] = atom_now[corr_vect[i_here]].sx-atom_now[corr_vect_ref[i_here]].sx;
						sij_now[1] = atom_now[corr_vect[i_here]].sy-atom_now[corr_vect_ref[i_here]].sy;
						sij_now[2] = 0.0;//atom_now[corr_vect[i]].sz-atom_now[corr_vect_ref[i]].sz;
						for(int kk=0;kk<3;kk++) sij_now[kk] = sij_now[kk]-(int)(sij_now[kk]*2);

						sij_ref[0] = atom_now[i].sx - atom_now[corr_vect_ref[i_here]].sx;
						sij_ref[1] = atom_now[i].sy - atom_now[corr_vect_ref[i_here]].sy;
						sij_ref[2] = atom_now[i].sz - atom_now[corr_vect_ref[i_here]].sz;

						sij_ref[0] = sij_ref[0]-(int)(sij_ref[0]*2);
						sij_ref[1] = sij_ref[1]-(int)(sij_ref[1]*2);
						sij_ref[2] = sij_ref[2]-(int)(sij_ref[2]*2);

						smid_now[0]=atom_now[i].sx+sij_ref[0]/2.0;
						smid_now[1]=atom_now[i].sy+sij_ref[1]/2.0;
						smid_now[2]=atom_now[i].sz+sij_ref[2]/2.0;
						for(int ii2=0;ii2<3;ii2++)
						{
							if(smid_now[ii2]>1) smid_now[ii2]=smid_now[ii2]-1;
							if(smid_now[ii2]<0) smid_now[ii2]=smid_now[ii2]+1;
						}

						V3mulM3(sij_now,H_now,rij_now);
						V3mulM3(smid_now,H_now,rmid_now);
						double dx = (rij_now[0]);
						double dy = (rij_now[1]);
						double dz = (rij_now[2]);

						if(i==7113)
						{
							cout << "requested data\t"<<i<<"\t"<<corr_vect[i_here]<<"\t"<<"\t"<<atom_now[corr_vect[i_here]].type<<"\t"<<rij_ref[0]<<"\t"<<rij_ref[1]<<"\t"<<rij_ref[2];
							cout <<"\t"<<rij_now[0]<<"\t"<<rij_now[1]<<"\t"<<rij_now[2]<<"\n";
						}

						if(i==4262)
						{
							cout << rij_now[0]<<"\t"<<rij_now[1]<<"\t"<<rij_ref[0]<<"\t"<<rij_ref[1]<<"\t"<<dx<<"\t"<<dy<<"\n";
						}

						//dx=dx*sign_dir; dy=dy*sign_dir;dz=dz*sign_dir;
						double tot_disrig = sqrt(dx*dx+dy*dy+dz*dz);
						bool discard = false;
						if((rij_ref[0]*rij_ref[0]+rij_ref[1]*rij_ref[1]+rij_ref[2]*rij_ref[2])>0.5*(H_now[0][0]*H_now[0][0]+H_now[1][1]*H_now[1][1]+H_now[2][2]*H_now[2][2]))
						{
							dx=0.0;dy=0.0;dz=0.0;
							discard = true;
						}
						if(fabs(dx) > 0.5*(H_now[0][0])) {dz = dy = dx = 0.0; discard = true;}
						if(fabs(dy) > 0.5*(H_now[1][1])) {dz = dy = dx = 0.0; discard = true;}
						if(fabs(dz) > 0.5*(H_now[2][2])) {dz = dy = dx = 0.0; discard = true;}
						if(!discard)
						{
							av_nm++;
							sij_mean[0]+=sij_now[0];
							sij_mean[1]+=sij_now[1];
							sij_mean[2]+=sij_now[2];
							fprintf(fptr,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf \n", i,corr_vect[i_here],atom[i].type,atom[corr_vect[i_here]].type,rmid_now[0],rmid_now[1],rmid_now[2],dx,dy,dz);
						}
					}
				}
				fprintf(fptr2,"%d %lf %lf %lf %lf %lf %lf\n",i,atom_now[i].sx,atom_now[i].sy,atom_now[i].sz,sij_mean[0]/(av_nm*1.0),sij_mean[1]/(av_nm*1.0),sij_mean[2]/(av_nm*1.0));
				if(i==7059)
				{
					for(int ii3=0;ii3<6;ii3++)
					{
						cout << corr_vect_ref[ii3]<<"\t"<<corr_vect[ii3]<<"\n";
					}
				}
			}
		}
		fclose(fptr);
		fclose(fptr2);
}

void compute_disrigistry(atomic_dat *atom_now, int n_now,double H_now[3][3])
{

	double sxi,syi,szi,sij[3],rij[3],rijsq;
	double sij_mid[3],rij_mid[3];
	int ref_type;
	FILE *fptr;
	fptr = fopen("interface_disrigistry","w");
	int counter=0;
	int rejected_counter=0;

	int jbeg,jend,jnab;

	for(int i=0;i<n;i++)
	{
		jbeg = nbr_ptr[i];
		jend = nbr_ptr1[i];
		//
		double correspond_info[(jend-jbeg+1)][7];
		//	cout <<"entered herer\t"<<i<<"\n";
		sxi = atom_now[i].sx;
		syi = atom_now[i].sy;
		szi = atom_now[i].sz;
	    for (jnab = jbeg; jnab <= jend; jnab++)
		{
			int j = nbr_lst[jnab]; //((i==19115)||(i==9371)||(j==19115)||(j==9371))
			if(((atom_now[i].interface==1)&&(atom_now[i].type==1))||((atom_now[j].interface==1)&&(atom_now[j].type==1))){
				sij[0] = sxi - atom_now[j].sx;
				sij[1] = syi - atom_now[j].sy;
				sij[2] = szi - atom_now[j].sz;
				sij[0] = sij[0]-(int)(sij[0]*2);
				sij[1] = sij[1]-(int)(sij[1]*2);
				sij[2] = sij[2]-(int)(sij[2]*2);

				V3mulM3(sij,H_now,rij);
				rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
				//
				if (rijsq <= rcoordsq[atom_now[i].type-1][atom_now[j].type-1])
				{
					sij_mid[0] = (-1.0*sxi + atom_now[j].sx);sij_mid[1] = (-1.0*syi + atom_now[j].sy);sij_mid[2] = (-1.0*szi + atom_now[j].sz);

					sij_mid[0] = sxi+0.5*(sij_mid[0]-(int)(sij_mid[0]*2));sij_mid[1] = syi+0.5*(sij_mid[1]-(int)(sij_mid[1]*2));sij_mid[2] = szi+0.5*(sij_mid[2]-(int)(sij_mid[2]*2));
					V3mulM3(sij_mid,H_now,rij_mid);
					correspond_info[jnab-jbeg][0] = rij_mid[0]; correspond_info[jnab-jbeg][1] = rij_mid[1];correspond_info[jnab-jbeg][2] = rij_mid[2];

					if((atom_now[i].interface==1)&&(atom_now[i].type==1)) ref_type = 2;if((atom_now[i].interface==1)&&(atom_now[i].type==2)) ref_type = 3;
					if((atom_now[i].interface==2)&&(atom_now[i].type==1)) ref_type = 0;if((atom_now[i].interface==2)&&(atom_now[i].type==2)) ref_type = 1;
					if((atom_now[i].interface<1)&&(atom_now[i].type==1)) ref_type = 0;if((atom_now[i].interface<1)&&(atom_now[i].type==2)) ref_type = 1;
					//cout << ref_type<<"\t"<<i<<"\t"<<atom[i].coord<<"\t"<<j<<"\t"<<atom[j].coord<<"\n";
					correspond_info[jnab-jbeg][3] = 1.0*determine_correspondence(ref_type,rij);
					double slip_len=0.0;
					if(correspond_info[jnab-jbeg][3]>-1)
					{
						int k = (int) correspond_info[jnab-jbeg][3];
						double tpos[3],tneg[3];
						tpos[0] = ref_vectors[k][0]+rij[0];
						tpos[1] = ref_vectors[k][1]+rij[1];
						tpos[2] = ref_vectors[k][2]+rij[2];

						tneg[0] = ref_vectors[k][0]-rij[0];
						tneg[1] = ref_vectors[k][1]-rij[1];
						tneg[2] = ref_vectors[k][2]-rij[2];

						if((tpos[0]*tpos[0]+tpos[1]*tpos[1]+tpos[2]*tpos[2])>(tneg[0]*tneg[0]+tneg[1]*tneg[1]+tneg[2]*tneg[2]))
						{
							correspond_info[jnab-jbeg][4] = ref_vectors[k][0]-rij[0];
							correspond_info[jnab-jbeg][5] = ref_vectors[k][1]-rij[1];
							correspond_info[jnab-jbeg][6] = ref_vectors[k][2]-rij[2];
						}else
						{
							correspond_info[jnab-jbeg][4] = ref_vectors[k][0]+rij[0];
							correspond_info[jnab-jbeg][5] = ref_vectors[k][1]+rij[1];
							correspond_info[jnab-jbeg][6] = ref_vectors[k][2]+rij[2];
						}
						slip_len = correspond_info[jnab-jbeg][4]*correspond_info[jnab-jbeg][4]+correspond_info[jnab-jbeg][5]*correspond_info[jnab-jbeg][5];//+ correspond_info[jnab-jbeg][6]*correspond_info[jnab-jbeg][6];

					}else
					{
						correspond_info[jnab-jbeg][4] = 0;correspond_info[jnab-jbeg][5] = 0;correspond_info[jnab-jbeg][6] = 0;
						cout<<"***************NO VECTORS FOUND*******\n"<<i<<"\t"<<j<<"\n";
					}

					if(slip_len>0.3)
						fprintf(fptr,"%d %d %d %d %d %lf %lf %lf %lf %lf %lf\n",counter,i,j,atom[i].CNA,atom[j].CNA,correspond_info[jnab-jbeg][0],correspond_info[jnab-jbeg][1],correspond_info[jnab-jbeg][2],correspond_info[jnab-jbeg][4],correspond_info[jnab-jbeg][5],correspond_info[jnab-jbeg][6]);
					counter++;

				}else
				{
					correspond_info[jnab-jbeg][0] = -1000;correspond_info[jnab-jbeg][1] = -1000;correspond_info[jnab-jbeg][2] = -1000;
					correspond_info[jnab-jbeg][3]= -1;
					correspond_info[jnab-jbeg][4] = 0;correspond_info[jnab-jbeg][5] = 0;correspond_info[jnab-jbeg][6] = 0;
				}
			}
		}

	}


}


void compute_disrigistry(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform)
{
	// Will compute the disrigistry vector for all the atoms in the current configuration with respect to the reference configuration
	// Will compute for those atoms that are neighbors in the current configuration

	// assumes that the neighbor lists that are created are for the current configuration.

	// assumes also that the neighbors have been found and saved

	// The H_ref here is in crystal format - assumed - so use the transformation matrix to convert this H_ref to the one to be used

	double **diff_ten;
	int max_size=0;
	for(int i=0;i<n;i++)
	{
		max_size = max_size+atom[i].coord;
	}
	diff_ten = (double **) malloc(max_size*sizeof(double));
	for(int i=0;i<max_size;i++)
	{
		diff_ten[i] = (double *) malloc(8*sizeof(double));
		for(int j=0;j<8;j++)
		{
			diff_ten[i][j] = -1.0;
		}
	}

	double H_t_now[3][3];



	if(transform)
	{

		M3mul(H_now,KS1_trans_mas,H_t_now);
	}else
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				H_t_now[i][j] = H_now[i][j];
			}
			cout <<" H_ref\n";
		}

	}

	double sxi_now,syi_now,szi_now,sij_now[3],rij_now[3],rijsq_now;
	double sxi_ref,syi_ref,szi_ref,sij_ref[3],rij_ref[3],rijsq_ref;
	double sxmid_now, symid_now,szmid_now,smid_now[3],rmid_now[3];

	cout <<"in disrigistry\n";
	for(int i=0;i<n;i++)
	{

		for(int j=0;j<3;j++)
		{
			atom_now[i].disrigistry[j]=0.0;
		}
		atom[i].disrigistry_number_atoms = 0;
		atom[i].disrigistry_number_atoms_negligible = 0;
	}
	FILE *fptr;

	fptr = fopen("disregistry_data","w");
	int counter=1;
	int diff_ten_counter=0;
	for(int i=0;i<n;i++)
	{

		if((atom_now[i].interface==1)&&(atom_now[i].type==1))
		{
			atom_now[i].disrigistry[0]=0.0;
			atom_now[i].disrigistry[1]=0.0;
			atom_now[i].disrigistry[2]=0.0;
			atom_now[i].drig = 0.0;
			double av_nm=0;


			double sign_dir = 1.0;
			int jbeg = nbr_ptr[i];
			int jend = nbr_ptr1[i];
			sxi_now = atom_now[i].sx;syi_now = atom_now[i].sy;szi_now = atom_now[i].sz;
			sxi_ref = atom_ref[i].sx;syi_ref = atom_ref[i].sy; szi_ref = atom_ref[i].sz;

			for(int jnab=0;jnab<MAX_COORD;jnab++)
			{
				int j = atom_now[i].coord_id[jnab];
				if(j>-1)
				{
					if(((atom_now[i].interface==1)&&(atom_now[i].type==1))&&((atom_now[j].interface==1)&&(atom_now[j].type==2)))
	//				if(((atom_now[i].interface==1)&&(atom_now[i].type==1)&&(atom_now[j].type==1)))
					{
					if(atom_now[i].type==2) sign_dir = -1.0;

						sij_now[0] = sxi_now - atom_now[j].sx;sij_now[1] = syi_now - atom_now[j].sy;sij_now[2] = szi_now - atom_now[j].sz;

						sij_now[0] = sij_now[0]-(int)(sij_now[0]*2);sij_now[1] = sij_now[1]-(int)(sij_now[1]*2);sij_now[2] = sij_now[2]-(int)(sij_now[2]*2);

						smid_now[0]=sxi_now+sij_now[0]/2.0;smid_now[1]=syi_now+sij_now[1]/2.0;smid_now[2]=szi_now+sij_now[2]/2.0;
						for(int ii2=0;ii2<3;ii2++)
						{
							if(smid_now[ii2]>1) smid_now[ii2]=smid_now[ii2]-1;
							if(smid_now[ii2]<0) smid_now[ii2]=smid_now[ii2]+1;
						}

						V3mulM3(sij_now,H_t_now,rij_now);
						V3mulM3(smid_now,H_t_now,rmid_now);

						rijsq_now = rij_now[0]*rij_now[0]+rij_now[1]*rij_now[1]+rij_now[2]*rij_now[2];
						if (rijsq_now <= rcoordsq[atom_now[i].type-1][atom_now[j].type-1])
						{
							sij_ref[0] = sxi_ref - atom_ref[j].sx; sij_ref[1] = syi_ref - atom_ref[j].sy;sij_ref[2] = szi_ref - atom_ref[j].sz;
							sij_ref[0] = sij_ref[0]-(int)(sij_ref[0]*2);sij_ref[1] = sij_ref[1]-(int)(sij_ref[1]*2);sij_ref[2] = sij_ref[2]-(int)(sij_ref[2]*2);


							V3mulM3(sij_ref,H_ref,rij_ref);
							double sdr[3];
							for(int pp=0;pp<3;pp++)
							{
								sdr[pp] = sij_now[pp]-sij_ref[pp];
							}

							double dx = (rij_now[0]-rij_ref[0]);
							double dy = (rij_now[1]-rij_ref[1]);
							double dz = (rij_now[2]-rij_ref[2]);

							dx=dx*sign_dir; dy=dy*sign_dir;dz=dz*sign_dir;
							double tot_disrig = sqrt(dx*dx+dy*dy+dz*dz);
							bool discard = false;
							if((rij_ref[0]*rij_ref[0]+rij_ref[1]*rij_ref[1]+rij_ref[2]*rij_ref[2])>0.5*(H_now[0][0]*H_now[0][0]+H_now[1][1]*H_now[1][1]+H_now[2][2]*H_now[2][2]))
							{
								dx=0.0;dy=0.0;dz=0.0;
								discard = true;
							}
							if(fabs(dx) > 0.5*(H_now[0][0])) {dz = dy = dx = 0.0; discard = true;}
							if(fabs(dy) > 0.5*(H_now[1][1])) {dz = dy = dx = 0.0; discard = true;}
							if(fabs(dz) > 0.5*(H_now[2][2])) {dz = dy = dx = 0.0; discard = true;}
							if(!discard)
							{
								diff_ten[diff_ten_counter][0] = i;
								diff_ten[diff_ten_counter][1] = j;
								diff_ten[diff_ten_counter][2] = rmid_now[0];
								diff_ten[diff_ten_counter][3] = rmid_now[1];
								diff_ten[diff_ten_counter][4] = rmid_now[2];
								diff_ten[diff_ten_counter][5] = dx;
								diff_ten[diff_ten_counter][6] = dy;
								diff_ten[diff_ten_counter][7] = dz;
								diff_ten_counter++;
								fprintf(fptr,"%d %d %d %d %d %lf %lf %lf %lf %lf %lf\n",counter,i,j,atom[i].type,atom[j].type,smid_now[0],smid_now[1],smid_now[2],sdr[0],sdr[1],sdr[2]);

//								av_nm++;
//								atom_now[i].disrigistry[0]=
							}
						}
					}

				}
			}
		}

	}
	fclose(fptr);
//	cout << "diff_Tensor is \t"<<diff_ten_counter<<"\n";
//	compute_diff_ten(diff_ten,diff_ten_counter,3);

	for(int i=0;i<max_size;i++)
	{
		free(diff_ten[i]);
	}
		free(diff_ten);

}

int load_ref_vectors(char *ref_vector_file)
{
	//load reference vectors
	// for now from a file with only one set of vectors  i.e.
	// for now irrespective of where the dislocation is, this will use same ref vectors
	// need to extend it so that we can use different vectors for different regions
	// for example, one set for interface, one set for Cu-Cu and one set for Nb-Nb

	// the file will have vectors in the same frame as the target file matrix
	// here, we are only collecting them
	// for flexibility, first column will have the type number
	// as the first line, the reference file will contain the total number of rows

	FILE *fptr;
	double dummy;
	//fptr = fopen("ref.vectors","r");
	fptr = fopen(ref_vector_file,"r");
	if(fptr!=NULL)
	{
		fscanf(fptr,"%d",&ref_vector_n);
		if((ref_vector_n>0)&&(ref_vector_n<1000))
		{
			ref_vector_list = (double **) malloc(ref_vector_n*sizeof(double));
			for(int i=0;i<ref_vector_n;i++)
			{
				ref_vector_list[i] = (double*) malloc(4*sizeof(double));
			}
			for(int i=0;i<ref_vector_n;i++)
			{
				double dummy2, dummy3,dummy4;
				int dummy1;

				fscanf(fptr,"%d %lf %lf %lf",&dummy1,&dummy2,&dummy3,&dummy4);
				ref_vector_list[i][0] = (double)dummy1;
				ref_vector_list[i][1] = dummy2;
				ref_vector_list[i][2] = dummy3;
				ref_vector_list[i][3] = dummy4;
			}
		}else
		{
			cout << "too many vectors in reference file";
			fclose(fptr);
			exit(1);
		}
		fclose(fptr);
	}

	return 0;
}

void find_3_atoms(atomic_dat *atom_now, int n_now, double H_now[3][3],int *ref_list_p,int ref_a_list[6])
{
	// first find a fcc atom that is in the bulk and next to Cu alpha plane
	cout << "***** Entered finding the atom thing\n";
	bool found = false;
	int counter=0;
	do{
		if((atom_now[counter].sx>0.1)&&(atom_now[counter].sx<0.9)&&(atom_now[counter].sy>0.1)&&(atom_now[counter].sy<0.9))
		{

			if((atom_now[counter].type==1)&&(atom[counter].interface==0)&&(atom[counter].coord==12)&&(atom[counter].CNA==0))
			{
				int counter_inside=0;
				for(int jnab=0;jnab<MAX_COORD;jnab++)
				{
					int j = atom_now[counter].coord_id[jnab];
					if(j>-1)
					{
						if((atom[j].interface==2)&&(atom[j].CNA==0))counter_inside++;
					}
				}
				//cout << counter<<"\t"<<atom_now[counter].coord<<"\n";
				if(counter_inside==3)
				{
					*ref_list_p=counter;found=true;
				}else
				{
					counter++;
				}
			}else{
				//cout << counter<<"\t"<<atom_now[counter].type<<"\t"<<atom_now[counter].interface<<"\t"<<atom_now[counter].coord<<"\n";
				counter++;
			}
			}else
			{
				//cout << counter<<"\t"<<atom_now[counter].type<<"\t"<<atom_now[counter].interface<<"\t"<<atom_now[counter].coord<<"\n";
				counter++;
			}
	}while(!found);
	cout << *ref_list_p<<" is the atom I have found \n";
	//now find atoms in the next plane
	// this is done by discarding the previous plane which is tagged as interface=2
	// and then atoms in the same plane are rejected because they ar next to interface==2
	int counter_inside=0;
	int atom_list=0;
	for(int jnab=0;jnab<MAX_COORD;jnab++)
	{
		int j = atom_now[*ref_list_p].coord_id[jnab];
		if(j>-1)
		{
			bool done=false;
			do{
				if(atom_now[j].interface==2)
				{
					ref_a_list[atom_list]=j;
					atom_list++;
					done = true;
				}else
				{

					for(int jnab1=0;jnab1<MAX_COORD;jnab1++)
					{
						int k = atom_now[j].coord_id[jnab1];
						if(k>-1)
						{
							if(atom_now[k].interface==2){done=true;jnab1=MAX_COORD+1;}
						}
					}
					if(!done)
					{
						ref_a_list[atom_list]=j;
						atom_list++;
						done=true;
						cout <<j<<" is the neighbor\n";
					}
				}

			}while(!done);
		}
	}

}

void perform_disrigistry(char *firstfile,char *secondfile,int dis_type, int perform_mode,int s_atom_type,int s_atom_interface)
{
	atomic_dat *atom_ref;
	double H_ref[3][3];
	int is_read;

	cout << firstfile<<"\t"<<secondfile<<"\n";
	if(dis_type==1)
	{
		is_read = read_lammps(secondfile, atom,true,true);
		if(is_read==0)
		{
				compute_CNA_and_others(atom,n, Hcry);
		}else
		{
			cout << "secondary input file not found\n";
			exit(1);
		}
		atom_ref = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
		copy_atomstruct(atom, atom_ref,n);
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				H_ref[i][j] = Hcry[i][j];
				cout << H_ref[i][j]<<"\t";
			}
		}
		//compute_disrigistry(atom_main_now, atom, H_main_now, Hcry,fa);
	}else
	{
		load_ref_vectors(secondfile);
	}


	is_read = read_lammps(firstfile, atom,true,true);
	if(is_read==0)
	{
		compute_CNA_and_others(atom,n, Hcry);
	}else
	{
		cout << "main input file not found\n";
		exit(1);
	}

	if(dis_type==1)
	{
		compute_disrigistry(atom, atom_ref, Hcry, H_ref,false);
		if(perform_mode==1)
		{
			compute_disregistry_3(atom,n,Hcry,0,s_atom_type,s_atom_interface);
		}else if(perform_mode==2)
		{
			closest_atom_mode =2;
			compute_disregistry_3(atom,n,Hcry,0,s_atom_type,s_atom_interface);
			cout << "entering herer\n";
			int ref_a_list[6];
			int ref_list_p=-1;
			find_3_atoms(atom,n,Hcry,&ref_list_p,ref_a_list);
			cout << ref_list_p<<"\t"<<ref_a_list[0]<<"\t"<<ref_a_list[1]<<"\t"<<ref_a_list[2]<<"\n";
			//cout << ref_list_p<<"\t"<<ref_a_list[3]<<"\t"<<ref_a_list[4]<<"\t"<<ref_a_list[5]<<"\n";

			bias_list_1=ref_list_p;
			bias_list_2=ref_a_list[0];
			compute_disregistry_3(atom,n,Hcry,1,s_atom_type,s_atom_interface);
			bias_list_1=ref_list_p;
			bias_list_2=ref_a_list[1];
			compute_disregistry_3(atom,n,Hcry,2,s_atom_type,s_atom_interface);
			bias_list_1=ref_list_p;
			bias_list_2=ref_a_list[2];
			compute_disregistry_3(atom,n,Hcry,3,s_atom_type,s_atom_interface);
			bias_list_1=ref_list_p;
			bias_list_2=ref_a_list[3];
			compute_disregistry_3(atom,n,Hcry,4,s_atom_type,s_atom_interface);
			bias_list_1=ref_list_p;
			bias_list_2=ref_a_list[4];
			compute_disregistry_3(atom,n,Hcry,5,s_atom_type,s_atom_interface);
			bias_list_1=ref_list_p;
			bias_list_2=ref_a_list[5];
			compute_disregistry_3(atom,n,Hcry,6,s_atom_type,s_atom_interface);
		}
	}
}
