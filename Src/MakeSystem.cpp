/*
 * MakeSystem.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */
#include <MakeSystem.h>

void create_system(char *Cu_filename, char *Nb_filename, bool keep_xy, bool keepcms_together, bool keep_max)
{
	atomic_dat *atomCu, *atomNb, *system;

	cout << Cu_filename<<"\n";
	read_lammps(Cu_filename,atomCu,true,true);

	H_total[0][0] = lx;
	H_total[1][1] = ly;
	H_total[2][2] = lz;
cout << lx<<"\t"<<ly<<"\t"<<lz<<" Copper lattice\n";
	H_total[0][1]=H_total[0][2]=H_total[1][0]=H_total[1][2]=H_total[2][0]=H_total[2][1]=0.0;



	atomCu = (struct atomic_dat *) malloc((n+5)*sizeof(struct atomic_dat));

	for (int i = 0; i < n; i++)
	{
		atomCu[i].sx = atom[i].sx;
		atomCu[i].sy = atom[i].sy;
		atomCu[i].sz = atom[i].sz;
		atomCu[i].pe = atom[i].pe;
		atomCu[i].type = 1;
		atomCu[i].ma = 63.546;

		atomCu[i].vx = atom[i].vx;
		atomCu[i].vy = atom[i].vy;
		atomCu[i].vz = atom[i].vz;
		atomCu[i].ux = atom[i].ux;
		atomCu[i].uy = atom[i].uy;
		atomCu[i].uz = atom[i].uz;

		double s[3],r[3];
		s[0] = 0.0; s[1]=0.0;s[2]=0.0;r[0]=0.0;r[1]=0.0;r[2]=0.0;
		s[0] = atomCu[i].sx;s[1]=atomCu[i].sy;s[2]=atomCu[i].sz;
		V3mulM3(s,H_total,r);

		atomCu[i].rx = r[0];atomCu[i].ry=r[1];atomCu[i].rz=r[2];
	}

	n_Cu = n;
	lx_Cu = lx, ly_Cu = ly; lz_Cu = lz;

	read_lammps(Nb_filename,atomNb,true,true);

	n_Nb = n;
	lx_Nb = lx, ly_Nb = ly; lz_Nb = lz;
	cout << lx<<"\t"<<ly<<"\t"<<lz<<" Niobium lattice\n";


	H_total[0][0] = lx_Nb;
	H_total[1][1] = ly_Nb;
	H_total[2][2] = lz_Nb;

	H_total[0][1]=H_total[0][2]=H_total[1][0]=H_total[1][2]=H_total[2][0]=H_total[2][1]=0.0;

	for (int i = 0; i < n_Nb; i++)
	{

		double s[3],r[3];
		s[0] = 0.0; s[1]=0.0;s[2]=0.0;r[0]=0.0;r[1]=0.0;r[2]=0.0;
		s[0] = atom[i].sx;s[1]=atom[i].sy;s[2]=atom[i].sz;

		V3mulM3(s,H_total,r);

		atom[i].rx = r[0];atom[i].ry=r[1];atom[i].rz=r[2];

	}




	if(keep_xy)
	{
		if(!keep_max)
		{
			H_total[0][0] = (lx_Cu+lx_Nb)/2.0;
			H_total[1][1] = (ly_Cu+ly_Nb)/2.0;
		}else
		{
			H_total[0][0] = max(lx_Cu,lx_Nb);
			H_total[1][1] = max(ly_Cu,ly_Nb);
		}
	}else
	{
		H_total[0][0] = max(lx_Cu,lx_Nb)*1.50;
		H_total[1][1] = max(ly_Cu,ly_Nb)*1.50;
	}

	double center_x = H_total[0][0]/2;
	double center_y = H_total[1][1]/2;

	H_total[2][2] = (lz_Cu+lz_Nb)*1.6;
	H_total[0][1]=H_total[0][2]=H_total[1][0]=H_total[1][2]=H_total[2][0]=H_total[2][1]=0.0;

	if(keep_xy)
	{
	double H_temp[3][3];
	H_temp[0][0] = H_total[0][0];H_temp[1][1]=H_total[1][1];H_temp[2][2] = lz_Cu;
	H_temp[0][1] = H_temp[0][2] = H_temp[1][0] = H_temp[1][2] = H_temp[2][0] = H_temp[2][1] = 0.0;

	for (int i = 0; i < n_Cu; i++)
		{
			double s[3],r[3];
			s[0] = 0.0; s[1]=0.0;s[2]=0.0;r[0]=0.0;r[1]=0.0;r[2]=0.0;
			s[0] = atomCu[i].sx;s[1]=atomCu[i].sy;s[2]=atomCu[i].sz;
			V3mulM3(s,H_temp,r);

			atomCu[i].rx = r[0];atomCu[i].ry=r[1];atomCu[i].rz=r[2];
		}
	H_temp[2][2] = lz_Nb;
	for (int i = 0; i < n_Nb; i++)
		{

			double s[3],r[3];
			s[0] = 0.0; s[1]=0.0;s[2]=0.0;r[0]=0.0;r[1]=0.0;r[2]=0.0;
			s[0] = atom[i].sx;s[1]=atom[i].sy;s[2]=atom[i].sz;

			V3mulM3(s,H_temp,r);
			atom[i].rx = r[0];atom[i].ry=r[1];atom[i].rz=r[2];

		}
	}


	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			Hcry[i][j]=H_total[i][j];
			//cout <<Hcry[i][j]<<"\t"<<H_total[i][j]<<"\n";


		}
	}

	double ratio_Cu = lz_Cu/H_total[2][2];
	double ratio_Nb = lz_Nb/H_total[2][2];

	double initial_spacing = 5/H_total[2][2];
	double reset_space = 0;//2.5885/H_total[2][2];
	double Cu_freesurface = (1-ratio_Cu-ratio_Nb-initial_spacing-reset_space)/2.0;
	double Hc_inv[3][3];

	initial_spacing = initial_spacing -reset_space;
	M3inv(Hcry,Hc_inv);

	//	cout << n<<"\t"<<atomCu[1].sx<<"\t"<<atom[1].sx<<"\t"<<atomNb[1].sx<<"\n";
	cout << ratio_Cu<<"\t"<<ratio_Nb<<"\t"<<initial_spacing<<" ting tong\n";

	system = (struct atomic_dat *) malloc(((n_Cu+n_Nb)+5)*sizeof(struct atomic_dat));
	double max_sz=-10.0;
	double cm_x = 0.0; double cm_y = 0.0;
	for (int i = 0; i < n_Cu; i++)
	{
		double s[3],r[3];

		s[0] = 0.0; s[1]=0.0;s[2]=0.0;r[0]=0.0;r[1]=0.0;r[2]=0.0;
		r[0] = atomCu[i].rx;r[1]=atomCu[i].ry;r[2]=atomCu[i].rz;
		V3mulM3(r,Hc_inv,s);

		system[i].sx = s[0];
		system[i].sy = s[1];
		system[i].sz = s[2]+Cu_freesurface;
		if(system[i].sz>max_sz)
		{
			max_sz = system[i].sz;
		}
		system[i].rx = atomCu[i].rx;
		system[i].ry = atomCu[i].ry;
		system[i].rz = atomCu[i].rz+Cu_freesurface*H_total[2][2];

		cm_x += system[i].sx;
		cm_y += system[i].sy;
		system[i].pe = atomCu[i].pe;
		system[i].ma = atomCu[i].ma;
		system[i].type = 1;
		system[i].vx = atomCu[i].vx;
		system[i].vy = atomCu[i].vy;
		system[i].vz = atomCu[i].vz;
		system[i].ux = atomCu[i].ux;
		system[i].uy = atomCu[i].uy;
		system[i].uz = atomCu[i].uz;
	}

	if(keepcms_together)
	{
		cm_x = cm_x/n_Cu;
		cm_y = cm_y/n_Cu;
		for (int i = 0; i < n_Cu; i++)
		{
			system[i].sx = system[i].sx-(cm_x-0.5);
			if(system[i].sx>=1) system[i].sx = system[i].sx-1;
			if(system[i].sx<0) system[i].sx = 1+system[i].sx;

			system[i].sy = system[i].sy-(cm_y-0.5);
			if(system[i].sy>=1) system[i].sy = system[i].sy-1;
			if(system[i].sy<0) system[i].sy = 1+system[i].sy;
		}
	}
	cm_x = 0;
	cm_y = 0;
	cout <<"out of Cu\n";
	cout << "max z value here is \t"<<max_sz<<"\t"<<Cu_freesurface<<"\n";
	for (int j = 0; j < n_Nb; j++)
	{
		int i = j+n_Cu;

		double s[3],r[3];

		s[0] = 0.0; s[1]=0.0;s[2]=0.0;r[0]=0.0;r[1]=0.0;r[2]=0.0;
		r[0] = atom[j].rx;r[1]=atom[j].ry;r[2]=atom[j].rz;


		V3mulM3(r,Hc_inv,s);


		system[i].sx = s[0];
		system[i].sy = s[1];
		cm_x = cm_x+system[i].sx;
		cm_y = cm_y+system[i].sy;
		system[i].sz = s[2]+max_sz+initial_spacing;
		system[i].rx = atom[j].rx;
		system[i].ry = atom[j].ry;
		system[i].rz = atom[j].rz+(max_sz+initial_spacing)*H_total[2][2];

		system[i].pe = atom[j].pe;
		system[i].ma = atom[j].ma;
		system[i].type = 2;
		system[i].vx = atom[j].vx;
		system[i].vy = atom[j].vy;
		system[i].vz = atom[j].vz;
		system[i].ux = atom[j].ux;
		system[i].uy = atom[j].uy;
		system[i].uz = atom[j].uz;
	}
	if(keepcms_together)
		{
			cm_x = cm_x/n_Nb;
			cm_y = cm_y/n_Nb;
			for (int i = n_Cu; i < n_Cu+n_Nb; i++)
			{
				system[i].sx = system[i].sx-(cm_x-0.5);
				if(system[i].sx>=1) system[i].sx = system[i].sx-1;
				if(system[i].sx<0) system[i].sx = 1+system[i].sx;

				system[i].sy = system[i].sy-(cm_y-0.5);
				if(system[i].sy>=1) system[i].sy = system[i].sy-1;
				if(system[i].sy<0) system[i].sy = 1+system[i].sy;
			}
		}
	n = n_Cu+n_Nb;
	/*
	if(keepcms_together)
	{
		cm_x = 0;cm_y = 0;
		for(int i = 0;i<n;i++)
		{
			cm_x = cm_x+system[i].sx;
			cm_y = cm_y+system[i].sy;
		}
		cm_x = cm_x/n;cm_y = cm_y/n;
		for (int i = 0; i < n; i++)
		{
			system[i].sx = system[i].sx-(cm_x-0.5);
			if(system[i].sx>=1) system[i].sx-1;
			if(system[i].sx<0) 1+system[i].sx;

			system[i].sy = system[i].sy-(cm_y-0.5);
			if(system[i].sy>=1) system[i].sy-1;
			if(system[i].sy<0) 1+system[i].sy;
		}
	}
	*/

	free(atom);
	free(atomCu);
	atom = system;
	n = n_Cu+n_Nb;
	cout << n_Cu+n_Nb<<"\t"<<atom[1].sx<<"\n";
	save_cfg(34,H_total);
	save_lammps(34,H_total);

}

void create_alpha_system(int what_alpha, bool remove_pbc)
{
	// what_alpha values 0 - Cu interface to Cu_alpha interface
	// what_alpha values 1 - Nb interface to Nb_alpha interface
	// what_alpha values 2 - Cu and Nb interfaces to Cu_alpha and Nb_alpha interfaces
	// what_alpha values 3 - Cu and Nb, whole system, to Cu_alpha and Nb_alpha
	// what_alpha values 6 - do nothing
	// what_alpha values 4 - Cu to Cu_alpha
	// what_alpha values 5 - Nb to Nb_alpha

	double hx=1.0,hy=1.0,hz=1.0;
	if(remove_pbc)
	{
		hx = 2.0;hy=2.0;hz=2.0;
	}
	double Hcry2[3][3];
	Hcry2[0][0] = Hcry[0][0]*hx;Hcry2[0][1] =0.0*hx;Hcry2[0][2] =0.0*hx;
	Hcry2[1][1] = Hcry[1][1]*hy;Hcry2[1][0] =0.0*hy;Hcry2[1][2] =0.0*hy;
	Hcry2[2][2] = Hcry[2][2]*hz;Hcry2[2][0] =0.0*hz;Hcry2[2][1] =0.0*hz;

	double Hi[3][3], Hi2[3][3];
	M3inv(Hcry,Hi);
	M3inv(Hcry2,Hi2);
	cout << Hcry2[0][0]<<"\t"<<Hcry2[1][1]<<"\t"<<Hcry2[2][2]<<" the new Hcry is \n";

	double max_x=0.0,max_y=0.0,max_z=0.0;

	for(int i=0;i<n;i++)
	{
		double r1[3],s1[3],r[3],s[3];
		s[0] = atom[i].sx;
		if(s[0]>=1) s[0]=s[0]-1;
		if(s[0]<0) s[0]=1+s[0];

		s[1] = atom[i].sy;
		if(s[1]>=1) s[1]=s[1]-1;
		if(s[1]<0) s[1]=1+s[1];

		s[2] = atom[i].sz;
		if(s[2]>=1) s[2]=s[2]-1;
		if(s[2]<0) s[2]=1+s[2];

		V3mulM3(s,Hcry,r);

		if(what_alpha==0)
		{
			if((atom[i].type==1)&&(atom[i].interface==1))
			{
				M3mulV3(Cu_alpha_matrix,r,r1);
				cout << "entered into what_alpha=\t"<<what_alpha<<"\t"<<i<<"\n";
				if(r1[0]> max_x) max_x=r1[0];
				if(r1[1]>max_y) max_y = r1[1];
				if(r1[2]>max_z) max_z = r1[2];
				if((i==5758)||(i==5759))
				{
					cout <<"\t"<<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r1[0]<<"\t"<<r1[1]<<"\t"<<r1[2]<<"\n";
				}

			}else
			{
				r1[0] = r[0];r1[1] = r[1];r1[2]=r[2];

			}
		}else if(what_alpha==1)
		{
			if((atom[i].type==2)&&(atom[i].interface==1))
			{
				M3mulV3(Nb_alpha_matrix,r,r1);
				cout << "entered into what_alpha=\t"<<what_alpha<<"\t"<<i<<"\n";
			}else
			{
				r1[0] = r[0];r1[1] = r[1];r1[2]=r[2];
			}
		}else if(what_alpha==2)
		{
			if((atom[i].type==1)&&(atom[i].interface==1))
			{
				M3mulV3(Cu_alpha_matrix,r,r1);
				cout << "entered into what_alpha=\t"<<what_alpha<<"\t"<<i<<"\n";

			}else if((atom[i].type==2)&&(atom[i].interface==1))
			{
				M3mulV3(Nb_alpha_matrix,r,r1);
				cout << "entered into what_alpha=\t"<<what_alpha<<"\t"<<i<<"\n";
			}else
			{
				r1[0] = r[0];r1[1] = r[1];r1[2]=r[2];
			}
		}else if(what_alpha==3)
		{
			if((atom[i].type==1))
			{
				M3mulV3(Cu_alpha_matrix,r,r1);
				//cout << "entered into what_alpha=\t"<<what_alpha<<"\t"<<i<<"\n";
			} else if((atom[i].type==2))
			{
				M3mulV3(Nb_alpha_matrix,r,r1);
				//cout << "entered into what_alpha=\t"<<what_alpha<<"\t"<<i<<"\n";
			}
		}else if(what_alpha==4)
		{
				if((atom[i].type==1))
				{
					M3mulV3(Cu_alpha_matrix,r,r1);
					cout << "entered into what_alpha=\t"<<what_alpha<<"\t"<<i<<"\n";
					if(r1[0]> max_x) max_x=r1[0];
					if(r1[1]>max_y) max_y = r1[1];
					if(r1[2]>max_z) max_z = r1[2];
					if((i==5758)||(i==5759))
					{
						cout <<"\t"<<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r1[0]<<"\t"<<r1[1]<<"\t"<<r1[2]<<"\n";
					}
				}else
				{
					r1[0] = r[0];r1[1] = r[1];r1[2]=r[2];
				}
		}else if(what_alpha==5)
		{

				if((atom[i].type==2))
				{
					M3mulV3(Nb_alpha_matrix,r,r1);
					//	cout << "entered into what_alpha=\t"<<what_alpha<<"\t"<<i<<"\n";
					if(r1[0]> max_x) max_x=r1[0];
					if(r1[1]>max_y) max_y = r1[1];
					if(r1[2]>max_z) max_z = r1[2];
					if((i==5758)||(i==5759))
					{
						cout <<"\t"<<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r1[0]<<"\t"<<r1[1]<<"\t"<<r1[2]<<"\n";
					}
				}else
				{
					r1[0] = r[0];r1[1] = r[1];r1[2]=r[2];
				}

		}else
		{
				r1[0] = r[0];r1[1] = r[1];r1[2]=r[2];
		}



		V3mulM3(r1,Hi2,s1); //with transformation

		if(s1[0]>=1) s1[0]=s1[0]-1;
		if(s1[0]<0) s1[0]=1+s1[0];

		if(s1[1]>=1) s1[1]=s1[1]-1;
		if(s1[1]<0) s1[1]=1+s1[1];

		if(s1[2]>=1) s1[2]=s1[2]-1;
		if(s1[2]<0) s1[2]=1+s1[2];


		atom[i].sx = s1[0];
		atom[i].sy = s1[1];
		atom[i].sz = s1[2];
		if((i==5758)||(i==5759))
		{
			cout <<"\t"<<s[0]<<"\t"<<s[1]<<"\t"<<s[2]<<"\t"<<s1[0]<<"\t"<<s1[1]<<"\t"<<s1[2]<<"\n";
		}
	}

	cout << max_x<<"\t"<<max_y<<"\t"<<max_z<<" MAX VALUES ARE \n";

	//	prepare_nbrlist(Hcry2,100);

	cout <<"neighbor list prepared\n";
	//	coord_number(Hcry2);
	//	compute_CNA_and_others(atom,n, Hcry2);
	save_lammps(1,Hcry2);
	save_cfg(1,Hcry2);

	for(int i = 0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			Hcry[i][j] = Hcry2[i][j];
		}
	}

	// THIS IS FOR KS2 CONSTRUCTION

	/*
	 int jbeg,jend,jnab;
	 int coord[n];
	 double sxi,syi,szi,sij[3],rij[3],rijsq;

	 for(int i=0;i<n;i++)
	 {
		 if((i==14660)||(i==11650)||(i==3650)||(i==5674))
		 {
			 sxi = atom[i].sx;
			 syi	= atom[i].sy;
			 szi = atom[i].sz;
			 for (jnab = 0; jnab < MAX_COORD; jnab++)
			 {

				 int j = atom[i].coord_id[jnab];
				 if(j>-1){
					 sij[0] = sxi - atom[j].sx;
					 sij[1] = syi - atom[j].sy;
					 sij[2] = szi - atom[j].sz;
					 sij[0] = sij[0]-(int)(sij[0]*2);
					 sij[1] = sij[1]-(int)(sij[1]*2);
					 sij[2] = sij[2]-(int)(sij[2]*2);

					 V3mulM3(sij,Hcry2,rij);
					 rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
					 //
					 if (rijsq <= rcoordsq[atom[i].type-1][atom[j].type-1])
					 {

						 cout << i<<"\t"<<j<<"\t"<<atom[i].type<<"\t"<<atom[j].type<<"\t"<<rij[0]<<"\t"<<rij[1]<<"\t"<<rij[2]<<"\n";
					 }
				 }
			 }
		 }
	 }
	 */



}
