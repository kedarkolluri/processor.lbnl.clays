/*
 * slipvector.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */


/************************************************************************/
#include <slipvector.h>

void init_geometry(char *STRING1, char *STRING2)
{


	FILE *fptr;
	double x1,x2,x3;
	char str2[80];
	//fptr = fopen("/Users/kedar/Documents/include/geo-orient.dat","r");
	fptr = fopen(STRING1,"r");
	cout << "reading file info\n";

	fgets(str2,80,fptr);  fgets(str2,80,fptr);
	sscanf(str2,"%lf %lf %lf",&x1,&x2,&x3);

	cout << "H X values\t"<<x1<<"\t"<<x2<<"\t"<<x3<<"\n";

	H0_geo[0][0] = x1, H0_geo[0][1] =  x2, H0_geo[0][2] = x3;

	fgets(str2,80,fptr);
	sscanf(str2,"%lf %lf %lf",&x1,&x2,&x3);

	cout << "H X values\t"<<x1<<"\t"<<x2<<"\t"<<x3<<"\n";
	H0_geo[1][0] = x1, H0_geo[1][1] = x2, H0_geo[1][2] = x3;

	fgets(str2,80,fptr);
	sscanf(str2,"%lf %lf %lf",&x1,&x2,&x3);
	fclose(fptr);
	cout << "H X values\t"<<x1<<"\t"<<x2<<"\t"<<x3<<"\n";

	H0_geo[2][0] =  x1, H0_geo[2][1] =  x2, H0_geo[2][2] = x3;


	double H0_norm = sqrt(H0_geo[0][0]*H0_geo[0][0]+H0_geo[0][1]*H0_geo[0][1]+H0_geo[0][2]*H0_geo[0][2]);
	double H1_norm = sqrt(H0_geo[1][0]*H0_geo[1][0]+H0_geo[1][1]*H0_geo[1][1]+H0_geo[1][2]*H0_geo[1][2]);
	double H2_norm = sqrt(H0_geo[2][0]*H0_geo[2][0]+H0_geo[2][1]*H0_geo[2][1]+H0_geo[2][2]*H0_geo[2][2]);


	H0_geo[0][0] = H0_geo[0][0]/H0_norm;H0_geo[0][1] = H0_geo[0][1]/H0_norm;H0_geo[0][2] = H0_geo[0][2]/H0_norm;
	H0_geo[1][0] = H0_geo[1][0]/H1_norm;H0_geo[1][1] = H0_geo[1][1]/H1_norm;H0_geo[1][2] = H0_geo[1][2]/H1_norm;
	H0_geo[2][0] = H0_geo[2][0]/H2_norm;H0_geo[2][1] = H0_geo[2][1]/H2_norm;H0_geo[2][2] = H0_geo[2][2]/H2_norm;

	M3inv(H0_geo,H0_geo_inv);

	//fptr = fopen("/Users/kedar/Documents/include/fcc_bvs.dat","r");
	fptr = fopen(STRING2,"r");
	cout << "opened\n";
	fgets(str2,80,fptr);
	cout << str2<<"\n";
	fgets(str2,80,fptr);
	sscanf(str2,"%lf %lf",&x1,&x2);

	cout << x1<<"\t"<<x2<<"\n";
	int tot_Val = x1+x2;
	VECTORS_ENUM = x1;
	PLANES_ENUM = x2;
	cout <<   VECTORS_ENUM <<"\t"<< PLANES_ENUM <<"\n";
	BVS  = (double **)malloc((tot_Val)*sizeof(double));
	fgets(str2,80,fptr);
	cout << str2<<"\n";
	for (int id=0; id< tot_Val; id++)    BVS[id] = (double *)malloc(3*sizeof(double));
	for(int i=0;i<VECTORS_ENUM;i++)
    {
		fgets(str2,80,fptr);
		sscanf(str2,"%lf %lf %lf",&x1,&x2,&x3);
		BVS[i][0] = x1;   BVS[i][1] = x2;      BVS[i][2] = x3;
    }
	fgets(str2,80,fptr);
	for(int i=VECTORS_ENUM;i<(PLANES_ENUM+VECTORS_ENUM);i++)
    {
		fgets(str2,80,fptr);
		sscanf(str2,"%lf %lf %lf",&x1,&x2,&x3);
		BVS[i][0] = x1;   BVS[i][1] = x2;      BVS[i][2] = x3;
    }

	fclose(fptr);
	for(int i=0;i<tot_Val;i++)
    {
		cout<< "AAAAAAAAAAAAAAAAA\t"<<BVS[i][0] <<"\t"<<  BVS[i][1] <<"\t"<< BVS[i][2] <<"\n";
    }

}

void compute_slipvector(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform)
{

	// Will compute for those atoms that are neighbors in the current configuration

	// assumes that the neighbor lists that are created are for the current configuration.

	// assumes also that the neighbors have been found and saved

	// The H_ref here is in crystal format - assumed - so use the transformation matrix to convert this H_ref to the one to be used

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
		}

	}


	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout << H_now[i][j]<<"\t";
		}
	}


	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout << H_t_now[i][j]<<"\t";
		}
	}

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout << H_ref[i][j]<<"\t";
		}
	}


	double sxi_now,syi_now,szi_now,sij_now[3],rij_now[3],rijsq_now;
	double sxi_ref,syi_ref,szi_ref,sij_ref[3],rij_ref[3],rijsq_ref;
	double s_now[3],r_now[3];

	cout <<"in slip vector\n";
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

	fptr = fopen("slip_data","w");
	int counter=1;
	for(int i=0;i<n;i++)
	{

		sxi_now = atom_now[i].sx;syi_now = atom_now[i].sy;szi_now = atom_now[i].sz;
		s_now[0] = atom_now[i].sx;s_now[1] = atom_now[i].sy;s_now[2] = atom_now[i].sz;
		sxi_ref = atom_ref[i].sx;syi_ref = atom_ref[i].sy; szi_ref = atom_ref[i].sz;
		atom[i].disrigistry[2] = -2;
		atom[i].BV = -2;

		double r_sum[3];
		for(int c=0;c<3;c++) r_sum[c]=0.0;
		int n_slip=0;

		for (int jnab = 0; jnab < MAX_COORD; jnab++)
		{
			int j = atom_now[i].coord_id[jnab];
			if(j>-1)
			{

				sij_now[0] = sxi_now - atom_now[j].sx;sij_now[1] = syi_now - atom_now[j].sy;sij_now[2] = szi_now - atom_now[j].sz;

				sij_now[0] = sij_now[0]-(int)(sij_now[0]*2);sij_now[1] = sij_now[1]-(int)(sij_now[1]*2);sij_now[2] = sij_now[2]-(int)(sij_now[2]*2);

				V3mulM3(sij_now,H_t_now,rij_now);

				rijsq_now = rij_now[0]*rij_now[0]+rij_now[1]*rij_now[1]+rij_now[2]*rij_now[2];

				sij_ref[0] = sxi_ref - atom_ref[j].sx; sij_ref[1] = syi_ref - atom_ref[j].sy;sij_ref[2] = szi_ref - atom_ref[j].sz;
				sij_ref[0] = sij_ref[0]-(int)(sij_ref[0]*2);sij_ref[1] = sij_ref[1]-(int)(sij_ref[1]*2);sij_ref[2] = sij_ref[2]-(int)(sij_ref[2]*2);


				V3mulM3(sij_ref,H_ref,rij_ref);

				double dx = (rij_now[0]-rij_ref[0]);
				double dy = (rij_now[1]-rij_ref[1]);
				double dz = (rij_now[2]-rij_ref[2]);
				double tot_disrig = sqrt(dx*dx+dy*dy+dz*dz);

				if(tot_disrig>DISRIG_LIMIT)
				{
					r_sum[0] = r_sum[0]+dx;
					r_sum[1] = r_sum[1]+dy;
					r_sum[2] = r_sum[2]+dz;
					n_slip++;
				}




			}
		}

		if((n_slip>1))
		{
			int first_number,second_number,first_type,second_type;

			double slip = sqrt(r_sum[0]*r_sum[0]+r_sum[1]*r_sum[1]+r_sum[2]*r_sum[2])/(1.0*n_slip);

			atom_now[i].drig = slip;
			//	if((atom_now[i].interface==1)||(atom_now[i].interface==2))
			//	{
			V3mulM3(s_now,H_t_now,r_now);
			double r_sum_mod[3];
			V3mulM3(r_sum,KS1_trans_mas,r_sum_mod);
			determine_BV_plane(i, n_slip,-1.0*r_sum_mod[0]/n_slip,r_sum_mod[1]*-1.0/n_slip,r_sum_mod[2]*-1.0/n_slip);
			double *v;
			v   = (double *)malloc(3*sizeof(double));
			v[0] = -1.0*r_sum_mod[0]/n_slip;r_sum_mod[0]=0.0;
			v[1] = -1.0*r_sum_mod[1]/n_slip;r_sum_mod[1]=0.0;
			v[2] = -1.0*r_sum_mod[2]/n_slip;r_sum_mod[2]=0.0;
			V3norm(v);
			V3mulM3(v,H0_geo,r_sum_mod);
			free(v);
			if((atom_now[i].interface==1)&&(atom_now[i].type==1))
			{
//				fprintf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",counter,i,atom_now[i].type,r_now[0],r_now[1],r_now[2],r_sum[0]/n_slip,r_sum[1]/n_slip,r_sum[2]/n_slip,r_sum_mod[0],r_sum_mod[1],r_sum_mod[2],slip);
				fprintf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",counter,i,atom_now[i].type,r_now[0],r_now[1],r_now[2],r_sum[0],r_sum[1],r_sum[2],r_sum_mod[0],r_sum_mod[1],r_sum_mod[2],slip);

				counter++;

			}



			//				if(slip>0.8)

			//	}
		}else
		{
			atom_now[i].drig = 0.0;
		}


	}
	fclose(fptr);

}


void compute_inplaneslipvector(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform)
{

	// Will compute for those atoms that are neighbors in the current configuration

	// assumes that the neighbor lists that are created are for the current configuration.

	// assumes also that the neighbors have been found and saved

	// The H_ref here is in crystal format - assumed - so use the transformation matrix to convert this H_ref to the one to be used

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
		}

	}


	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout << H_now[i][j]<<"\t";
		}
	}


	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout << H_t_now[i][j]<<"\t";
		}
	}

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout << H_ref[i][j]<<"\t";
		}
	}


	double sxi_now,syi_now,szi_now,sij_now[3],rij_now[3],rijsq_now;
	double sxi_ref,syi_ref,szi_ref,sij_ref[3],rij_ref[3],rijsq_ref;
	double s_now[3],r_now[3];

	cout <<"in slip vector\n";
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

	fptr = fopen("slip_inplane_data","w");
	int counter=1;
	for(int i=0;i<n;i++)
	{

		sxi_now = atom_now[i].sx;syi_now = atom_now[i].sy;szi_now = atom_now[i].sz;
		s_now[0] = atom_now[i].sx;s_now[1] = atom_now[i].sy;s_now[2] = atom_now[i].sz;
		sxi_ref = atom_ref[i].sx;syi_ref = atom_ref[i].sy; szi_ref = atom_ref[i].sz;
		atom[i].disrigistry[2] = -2;
		atom[i].BV = -2;

		double r_sum[3];
		for(int c=0;c<3;c++) r_sum[c]=0.0;
		int n_slip=0;

		for (int jnab = 0; jnab < MAX_COORD; jnab++)
		{
			int j = atom_now[i].coord_id[jnab];
			if(j>-1)
			{
				if((atom_now[i].type==atom_now[j].type)&&(atom_now[i].interface==atom_now[j].interface))
				{

					sij_now[0] = sxi_now - atom_now[j].sx;sij_now[1] = syi_now - atom_now[j].sy;sij_now[2] = szi_now - atom_now[j].sz;

					sij_now[0] = sij_now[0]-(int)(sij_now[0]*2);sij_now[1] = sij_now[1]-(int)(sij_now[1]*2);sij_now[2] = sij_now[2]-(int)(sij_now[2]*2);

					V3mulM3(sij_now,H_t_now,rij_now);

					rijsq_now = rij_now[0]*rij_now[0]+rij_now[1]*rij_now[1]+rij_now[2]*rij_now[2];

					sij_ref[0] = sxi_ref - atom_ref[j].sx; sij_ref[1] = syi_ref - atom_ref[j].sy;sij_ref[2] = szi_ref - atom_ref[j].sz;
					sij_ref[0] = sij_ref[0]-(int)(sij_ref[0]*2);sij_ref[1] = sij_ref[1]-(int)(sij_ref[1]*2);sij_ref[2] = sij_ref[2]-(int)(sij_ref[2]*2);


					V3mulM3(sij_ref,H_ref,rij_ref);

					double dx = (rij_now[0]-rij_ref[0]);
					double dy = (rij_now[1]-rij_ref[1]);
					double dz = (rij_now[2]-rij_ref[2]);
					double tot_disrig = sqrt(dx*dx+dy*dy+dz*dz);

					if(tot_disrig>DISRIG_LIMIT)
					{
						r_sum[0] = r_sum[0]+dx;
						r_sum[1] = r_sum[1]+dy;
						r_sum[2] = r_sum[2]+dz;
						n_slip++;
					}



				}
			}
		}

		if((n_slip>1))
		{
			int first_number,second_number,first_type,second_type;

			double slip = sqrt(r_sum[0]*r_sum[0]+r_sum[1]*r_sum[1]+r_sum[2]*r_sum[2])/(1.0*n_slip);

			atom_now[i].drig = slip;
			//	if((atom_now[i].interface==1)||(atom_now[i].interface==2))
			//	{
			V3mulM3(s_now,H_t_now,r_now);
			double r_sum_mod[3];
			V3mulM3(r_sum,KS1_trans_mas,r_sum_mod);
			determine_BV_plane(i, n_slip,-1.0*r_sum_mod[0]/n_slip,r_sum_mod[1]*-1.0/n_slip,r_sum_mod[2]*-1.0/n_slip);
			double *v;
			v   = (double *)malloc(3*sizeof(double));
			v[0] = -1.0*r_sum_mod[0]/n_slip;r_sum_mod[0]=0.0;
			v[1] = -1.0*r_sum_mod[1]/n_slip;r_sum_mod[1]=0.0;
			v[2] = -1.0*r_sum_mod[2]/n_slip;r_sum_mod[2]=0.0;
			V3norm(v);
			V3mulM3(v,H0_geo,r_sum_mod);
			free(v);
			if((atom_now[i].interface==1)&&(atom_now[i].type==1))
			{
//				fprintf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",counter,i,atom_now[i].type,r_now[0],r_now[1],r_now[2],r_sum[0]/n_slip,r_sum[1]/n_slip,r_sum[2]/n_slip,r_sum_mod[0],r_sum_mod[1],r_sum_mod[2],slip);
				fprintf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",counter,i,atom_now[i].type,r_now[0],r_now[1],r_now[2],r_sum[0],r_sum[1],r_sum[2],r_sum_mod[0],r_sum_mod[1],r_sum_mod[2],slip);

				counter++;

			}



			//				if(slip>0.8)

			//	}
		}else
		{
			atom_now[i].drig = 0.0;
		}


	}
	fclose(fptr);

}

void compute_absolute_displacement(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform)
{
	// Will compute for those atoms that are neighbors in the current configuration

	// assumes that the neighbor lists that are created are for the current configuration.

	// assumes also that the neighbors have been found and saved

	// The H_ref here is in crystal format - assumed - so use the transformation matrix to convert this H_ref to the one to be used

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
		}

	}


	double rdel[3];
	double s_now[3],r_now[3];

	cout <<"in displacement vector vector\n";

	FILE *fptr;

	fptr = fopen("displacement_data","w");
	for(int i=0;i<n;i++)
	{
		double adjust[3];
		adjust[0]= -0.0;
		adjust[1] = 0.0;
		adjust[2]=0.0;
		s_now[0] = atom_now[i].sx;s_now[1] = atom_now[i].sy;s_now[2] = atom_now[i].sz;
		double s_temp[3];

		for(int j=0;j<3;j++)
		{
			s_temp[j]=s_now[j]+adjust[j];
	//		s_temp[j] = s_temp[j]-(int)(s_temp[j]*2);
			if(s_temp[j]>=1) s_temp[j]-=1;
			if(s_temp[j]<0) s_temp[j]+=1;

		}
		V3mulM3(s_temp,H_now,r_now);



		s_now[0] -= atom_ref[i].sx;s_now[1] -= atom_ref[i].sy; s_now[2] -= atom_ref[i].sz;
		s_now[0] = s_now[0]-(int)(s_now[0]*2);s_now[1] = s_now[1]-(int)(s_now[1]*2);s_now[2] = s_now[2]-(int)(s_now[2]*2);

		V3mulM3(s_now,H_now,rdel);

		double rdel_mag = sqrt(rdel[0]*rdel[0]+rdel[1]*rdel[1]+rdel[2]*rdel[2]);

		if((atom_now[i].interface==1)&&(atom_now[i].type==1))
//		if(rdel_mag>0.2)
		{
			fprintf(fptr,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",i,atom_now[i].type,r_now[0],r_now[1],r_now[2],rdel[0],rdel[1],rdel[2],rdel_mag,atom_now[i].pe);


		}


	}
	fclose(fptr);
}
void compute_absolute_displacement(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3],bool transform, int interface_type,int atom_type)
{

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
		}

	}


	double rdel[3];
	double s_now[3],r_now[3];

	cout <<"in displacement vector vector specific\n";
	   // cout << atom_type << "\t"<<interface_type<<"\n";

	FILE *fptr;

	fptr = fopen("displacement_data","w");
	for(int i=0;i<n;i++)
	{
		double adjust[3];
		adjust[0]= -0.0;
		adjust[1] = 0.0;
		adjust[2]=0.0;
		s_now[0] = atom_now[i].sx;s_now[1] = atom_now[i].sy;s_now[2] = atom_now[i].sz;
		double s_temp[3];

		for(int j=0;j<3;j++)
		{
			s_temp[j]=s_now[j]+adjust[j];
	//		s_temp[j] = s_temp[j]-(int)(s_temp[j]*2);
			if(s_temp[j]>=1) s_temp[j]-=1;
			if(s_temp[j]<0) s_temp[j]+=1;

		}
		V3mulM3(s_temp,H_now,r_now);



		s_now[0] -= atom_ref[i].sx;s_now[1] -= atom_ref[i].sy; s_now[2] -= atom_ref[i].sz;
		s_now[0] = s_now[0]-(int)(s_now[0]*2);s_now[1] = s_now[1]-(int)(s_now[1]*2);s_now[2] = s_now[2]-(int)(s_now[2]*2);

		V3mulM3(s_now,H_now,rdel);

		double rdel_mag = sqrt(rdel[0]*rdel[0]+rdel[1]*rdel[1]+rdel[2]*rdel[2]);

		if((atom_now[i].interface==interface_type)&&(atom_now[i].type==atom_type))
		{
		  //cout << "values are "<<atom_now[i].type<<"\t"<<atom_now[i].interface<<"\n";
		  fprintf(fptr,"%d %d %lf %lf %lf %lf %lf %lf %lf \n",i,atom_now[i].type,r_now[0],r_now[1],r_now[2],rdel[0],rdel[1],rdel[2],rdel_mag);


		}


	}
	fclose(fptr);
}

void compute_absolute_displacement_interface(atomic_dat *atom_now, atomic_dat *atom_ref, double H_now[3][3],double H_ref[3][3])
{
	// Will compute for those atoms that are neighbors in the current configuration

	// assumes that the neighbor lists that are created are for the current configuration.

	// assumes also that the neighbors have been found and saved

	// The H_ref here is in crystal format - assumed - so use the transformation matrix to convert this H_ref to the one to be used

	double H_t_now[3][3];

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			H_t_now[i][j] = H_now[i][j];
		}
	}




	double rdel[3];
	double s_now[3],r_now[3];

	cout <<"in displacement vector vector\n";

	FILE *fptr;

	fptr = fopen("displacement_data","w");
	for(int i=0;i<n;i++)
	{
		double adjust[3];
		adjust[0]= -0.0;
		adjust[1] = 0.0;
		adjust[2]=0.0;
		s_now[0] = atom_now[i].sx;s_now[1] = atom_now[i].sy;s_now[2] = atom_now[i].sz;
		double s_temp[3];

		for(int j=0;j<3;j++)
		{
			s_temp[j]=s_now[j]+adjust[j];
	//		s_temp[j] = s_temp[j]-(int)(s_temp[j]*2);
			if(s_temp[j]>=1) s_temp[j]-=1;
			if(s_temp[j]<0) s_temp[j]+=1;

		}
		V3mulM3(s_temp,H_now,r_now);



		s_now[0] -= atom_ref[i].sx;s_now[1] -= atom_ref[i].sy; s_now[2] -= atom_ref[i].sz;
		s_now[0] = s_now[0]-(int)(s_now[0]*2);s_now[1] = s_now[1]-(int)(s_now[1]*2);s_now[2] = s_now[2]-(int)(s_now[2]*2);

		V3mulM3(s_now,H_now,rdel);

		double rdel_mag = sqrt(rdel[0]*rdel[0]+rdel[1]*rdel[1]+rdel[2]*rdel[2]);

		if((atom_now[i].interface!=0))
//		if(rdel_mag>0.2)
		{
			fprintf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",i,atom_now[i].type,atom_now[i].interface,r_now[0],r_now[1],r_now[2],rdel[0],rdel[1],rdel[2],rdel_mag,atom_now[i].pe);


		}


	}
	fclose(fptr);
}
