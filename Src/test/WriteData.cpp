/*
 * WriteData.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */
#include <WriteData.h>
bool LAMMPS_CHARGE = false;
bool LAMMPS_MOLECULE = false;
int save_lammps(int pp,double H[3][3])
{
	cout << "entered save lammps\n";
	if(((H[0][1])>0.0000)||((H[0][2])>0.0000)||((H[1][2])>0.0000)||((H[0][1])<0.0000)||((H[0][2])<0.0000)||((H[1][2])<0.0000))
	{
		cout << "the lammps save value of pp\t"<<pp<<" not done because H matrix is not in the desired format\n";
		return(0);
	}


	double r1[3],s1[3];
	double limits[3][2];

	limits[0][0] = 1e31;
	limits[1][0] = 1e31;
	limits[2][0] = 1e31;
	limits[0][1] = -1e31;
	limits[1][1] = -1e31;
	limits[2][1] = -1e31;
/*
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			for(int k=0;k<2;k++)
			{
				s1[0] = k,s1[1]=j,s1[2]=i;
				V3mulM3(s1,H,r1);
				if(r1[0]<limits[0][0]) limits[0][0] = r1[0];
				if(r1[1]<limits[1][0]) limits[1][0] = r1[1];
				if(r1[2]<limits[2][0]) limits[2][0] = r1[2];

				if(r1[0]>limits[0][1]) limits[0][1] = r1[0];
				if(r1[1]>limits[1][1]) limits[1][1] = r1[1];
				if(r1[2]>limits[2][2]) limits[2][1] = r1[2];

				cout << "limits for\t"<< s1[0]<<"\t"<<s1[1]<<"\t"<<s1[2]<<"\t";
				cout << "LIMITS ARE\t";
				printf("%e %e %e %e %e %e\n",limits[0][0],limits[0][1],limits[1][0],limits[1][1],limits[2][0],limits[2][1]);

			}
		}
	}

*/


	int i ;
	FILE *fptr;
	char filename[30] = "dat_lammps.", str[20];
	sprintf(str,"%d",pp);
	strcat(filename,str);
	fptr=fopen(filename,"w");
	fprintf(fptr, "%s\n\n","input file for LAMMPS");
	fprintf(fptr, "%d %s\n\n",n, " atoms");
	fprintf(fptr, "%d %s\n\n", n_types, "atom types");
	fprintf(fptr, "%12.9f %12.9lf %s\n",0.0,H[0][0], " xlo xhi");
	fprintf(fptr, "%12.9f %12.9lf %s\n",0.0,H[1][1], " ylo yhi");
	fprintf(fptr, "%12.9f %12.9lf %s\n\n",0.0,H[2][2], " zlo zhi");
	fprintf(fptr, "%12.9f %12.9lf %12.9lf %s\n\n",H[1][0],H[2][0],H[2][1], " xy xz yz");

	fprintf(fptr, "%s\n\n", "Masses");
	for (int ii=1;ii<=n_types;ii++)
	{
		fprintf(fptr, "%d %lf\n", ii,element_weights[ii-1]);
	}
	fprintf(fptr, "%s\n\n", "Atoms");
	for (i = 0; i < n; i++)
	{

		double r[3],s[3];

		if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
		if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
		if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;

		if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
		if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
		if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;


		s[0] = atom[i].sx;s[1] = atom[i].sy;s[2] = atom[i].sz;

		V3mulM3(s,H,r);

		fprintf(fptr,"%d ",i+1);
		if(LAMMPS_MOLECULE)
		{
			fprintf(fptr,"%d ", atom[i].molID);
		}

		fprintf(fptr,"%d ",atom[i].type);

		if(LAMMPS_CHARGE)
		{
			fprintf(fptr,"%lf ",atom[i].charge);
		}

		fprintf(fptr,"%12.9lf %12.9lf %12.9lf\n",r[0],r[1],r[2]);
	}

	fprintf(fptr, "%s\n\n", "Velocities");

	for (i = 0; i < n; i++)
	{

		fprintf(fptr,"%d %12.9lf %12.9lf %12.9lf\n",i+1,atom[i].vx,atom[i].vy,atom[i].vz);
	}


	fclose(fptr);
	zip(filename);

	return(0);
}

int save_lammps_specific(char *prefix,int pp,atomic_dat *atom_, double H[3][3])
{

	if(((H[0][1])>0.0000)||((H[0][2])>0.0000)||((H[1][2])>0.0000)||((H[0][1])<0.0000)||((H[0][2])<0.0000)||((H[1][2])<0.0000))
	{
		cout << "the lammps save value of pp\t"<<pp<<" not done because H matrix is not in the desired format\n";
		return(0);
	}
	int num_main = 0;
	for (int i = 0; i < n; i++)
		{
			if(atom_[i].type==2)
			{
				num_main++;
			}
		}

	double r1[3],s1[3];

	int i ;
	FILE *fptr;
	char filename[80]="",str[20];
	strcat(filename,prefix);
	strcat(filename,".");

	sprintf(str,"%d",pp);
	strcat(filename,str);
	fptr=fopen(filename,"w");
	fprintf(fptr, "%s\n\n","input file for LAMMPS");
	fprintf(fptr, "%d %s\n\n",num_main, " atoms");
	fprintf(fptr, "%s\n\n", "2 atom types");
	fprintf(fptr, "%12.9f %12.9lf %s\n",0.0,H[0][0], " xlo xhi");
	fprintf(fptr, "%12.9f %12.9lf %s\n",0.0,H[1][1], " ylo yhi");
	fprintf(fptr, "%12.9f %12.9lf %s\n\n",0.0,H[2][2], " zlo zhi");
	fprintf(fptr, "%12.9f %12.9lf %12.9lf %s\n\n",H[1][0],H[2][0],H[2][1], " xy xz yz");

	fprintf(fptr, "%s\n\n", "Masses");
	fprintf(fptr, "%s\n", "1 63.546");
	fprintf(fptr, "%s\n\n", "2 92.90638");
	fprintf(fptr, "%s\n\n", "Atoms");
	int num = 1;
	for (i = 0; i < n; i++)
	{


		if(atom_[i].type==2)
		{
		double r[3],s[3];

		if(atom_[i].sx>1.0) atom_[i].sx = atom_[i].sx-1;
		if(atom_[i].sy>1.0) atom_[i].sy = atom_[i].sy-1;
		if(atom_[i].sz>1.0) atom_[i].sz = atom_[i].sz-1;

		if(atom_[i].sx<0.0) atom[i].sx = 1.0+atom_[i].sx;
		if(atom_[i].sy<0.0) atom_[i].sy = 1.0+atom_[i].sy;
		if(atom_[i].sz<0.0) atom_[i].sz = 1.0+atom_[i].sz;


		s[0] = atom_[i].sx;s[1] = atom_[i].sy;s[2] = atom_[i].sz;

		V3mulM3(s,H,r);

	//	if((atom[i].type==1)&&(atom[i].interface==1)) atom[i].type = 3;
		fprintf(fptr,"%d %d %12.9lf %12.9lf %12.9lf\n",num,atom_[i].type,r[0],r[1],r[2]);
		num++;
		}
	}
/*
	fprintf(fptr, "%s\n\n", "Velocities");

	for (i = 0; i < n; i++)
	{

		fprintf(fptr,"%d %12.9lf %12.9lf %12.9lf\n",i+1,atom_[i].vx,atom_[i].vy,atom_[i].vz);
	}
*/

	fclose(fptr);
	zip(filename);

	return(0);
}

int save_cfg_test(int pp, double H[3][3])
{


	FILE *fptr;
	char filename[80]="dat_atomeye.",str[20];
	if(pp==-1) pp=0;
	sprintf(str,"%06d",pp);
	strcat(filename,str);
	fptr=fopen(filename,"w");

	fprintf(fptr,"%s %d\n\n","Number of particles =",n);

	fprintf(fptr,"%s %lf %s\n","H0(1,1) =",H[0][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(1,2) =",H[0][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(1,3) =",H[0][2],"A");

	fprintf(fptr,"%s %lf %s\n","H0(2,1) =",H[1][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(2,2) =",H[1][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(2,3) =",H[1][2],"A");

	fprintf(fptr,"%s %lf %s\n","H0(3,1) =",H[2][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(3,2) =",H[2][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(3,3) =",H[2][2],"A");

	fprintf(fptr,"%s\n","eta(1,1) = 0.0\neta(1,2) = 0.0\neta(1,3) = 0.0\neta(2,2) = 0.0\neta(2,3) = 0.0\neta(3,3) = 0.0\n");

	fprintf(fptr,"%s\n",".NO_VELOCITY.");

	fprintf(fptr,"%s %d\n","entry_count =",15);
	for (int i = 0; i < n; i++)
    {

		char symb[80]="";
		if(atom[i].type==1)
		{
			if(atom[i].interface!=1)
			{
			strcat(symb,"A");
			cout <<"here\t"<<i<<"\n";
			}else
			{
				strcat(symb,"Cu");
			}
		}
		if(atom[i].type==2)
		{
		if(atom[i].interface!=1)
		{
		strcat(symb,"B");
		}else
		{
			strcat(symb,"Nb");
		}
		}
		if(atom[i].type>2) strcat(symb,"Cu");
		double r[3];
		double s[3];
		double ratio=1;
		s[0] = atom[i].sx; s[1]=atom[i].sy;s[2]=atom[i].sz;
		//	V3mulM3(s,Hcry,r);
		int neigh_config=0;
		if(atom[i].interface==0)
		{
			neigh_config =0;
		}else
		{
			if(atom[i].neigh_config==0)
			{
				neigh_config = atom[i].interface;
			}else
			{
				neigh_config = atom[i].neigh_config+2;
			}
		}
		double delrr = sqrt(atom[i].delr[0]*atom[i].delr[0]+atom[i].delr[1]*atom[i].delr[1]+atom[i].delr[2]*atom[i].delr[2]);

		fprintf(fptr,"%lf\n%s\n%12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %d %d %12.9lf %12.9lf %12.9lf %12.9lf %d %d\n",
				atom[i].ma,symb,atom[i].sx,atom[i].sy,atom[i].sz,atom[i].delr[0], delrr,atom[i].pe,i*1.0,atom[i].coord,atom[i].interface,atom[i].drig,atom[i].BV*1.0,atom[i].disrigistry[2],neigh_config*1.0,atom[i].ackN,atom[i].CNA);

    }

	fclose(fptr);
	zip(filename);

	return(0);

}



int save_cfg(int pp, double H[3][3])
{

	cout << "entered save cfg\n";

	FILE *fptr;
	char filename[80]="dat_atomeye.",str[20];
	if(pp==-1) pp=0;
	sprintf(str,"%06d",pp);
	strcat(filename,str);
	fptr=fopen(filename,"w");

	fprintf(fptr,"%s %d\n\n","Number of particles =",n);

	fprintf(fptr,"%s %lf %s\n","H0(1,1) =",H[0][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(1,2) =",H[0][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(1,3) =",H[0][2],"A");

	fprintf(fptr,"%s %lf %s\n","H0(2,1) =",H[1][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(2,2) =",H[1][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(2,3) =",H[1][2],"A");

	fprintf(fptr,"%s %lf %s\n","H0(3,1) =",H[2][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(3,2) =",H[2][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(3,3) =",H[2][2],"A");

	fprintf(fptr,"%s\n","eta(1,1) = 0.0\neta(1,2) = 0.0\neta(1,3) = 0.0\neta(2,2) = 0.0\neta(2,3) = 0.0\neta(3,3) = 0.0\n");

	fprintf(fptr,"%s\n",".NO_VELOCITY.");
	//cout << "before the thing about htis\n";
	fprintf(fptr,"%s %d\n","entry_count =",15);
	for (int i = 0; i < n; i++)
    {
		char symb[80]="";
		if(atom==NULL) cout << "atom is null\n";

		//strcat(symb,atom[i].elem);
	// << i << "\t"<<atom[i].elem << "\t"<< symb<<"fgfgf\n";
		double r[3];
		double s[3];
		double ratio=1;
		//cout << "\hi\n";
		s[0] = atom[i].sx; s[1]=atom[i].sy;s[2]=atom[i].sz;
		//	V3mulM3(s,Hcry,r);
		int neigh_config=0;
		//cout << "hi\n";
		if(atom[i].interface==0)
		{
			neigh_config =0;
		}else
		{
			if(atom[i].neigh_config==0)
			{
				neigh_config = atom[i].interface;
			}else
			{
				neigh_config = atom[i].neigh_config+2;
			}
		}
		//cout << "hi2\n";
		//cout << "before delrr\n";
		double delrr = sqrt(atom[i].delr[0]*atom[i].delr[0]+atom[i].delr[1]*atom[i].delr[1]+atom[i].delr[2]*atom[i].delr[2]);
		//cout << "after delrr\n";
		fprintf(fptr,"%lf\n%s\n%12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %d %d %12.9lf %12.9lf %12.9lf %12.9lf %d %d\n",
				atom[i].ma,atom[i].elem,atom[i].sx,atom[i].sy,atom[i].sz,atom[i].type*1.0, delrr,atom[i].pe,i*1.0,atom[i].coord,atom[i].interface,atom[i].drig,atom[i].BV*1.0,atom[i].disrigistry[2],neigh_config*1.0,atom[i].ackN,atom[i].CNA);

    }

	fclose(fptr);
	zip(filename);

	return(0);

}

int save_cfg_interface(int pp, double H[3][3], int interface,int atom_type)
{

	int n_here=0;
	for(int i=0;i<n;i++)
	{
		if(((atom[i].interface==interface)&&(atom[i].type==atom_type)&&(atom_type!=0))||((atom[i].interface==interface)&&(atom_type==0))) n_here++;
	}
	FILE *fptr;
	char filename[80]="dat_atomeye_interface",str[20],str1[20],str2[20];
	char filename2[80]="add_color.";
	if(pp==-1) pp=0;
	strcat(filename,"_");
	sprintf(str1,"%03d",interface);
	strcat(filename,str1);
	strcat(filename,"_");
	sprintf(str1,"%03d",atom_type);
	strcat(filename,str1);
	strcat(filename,".");
	sprintf(str,"%06d",pp);
	strcat(filename,str);
	strcat(filename2,filename);
//	FILE *fptr_addcolor;
//	fptr_addcolor = fopen(filename2,"w");
	fptr=fopen(filename,"w");

	fprintf(fptr,"%s %d\n\n","Number of particles =",n_here);

	fprintf(fptr,"%s %lf %s\n","H0(1,1) =",H[0][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(1,2) =",H[0][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(1,3) =",H[0][2],"A");

	fprintf(fptr,"%s %lf %s\n","H0(2,1) =",H[1][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(2,2) =",H[1][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(2,3) =",H[1][2],"A");

	fprintf(fptr,"%s %lf %s\n","H0(3,1) =",H[2][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(3,2) =",H[2][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(3,3) =",H[2][2],"A");

	fprintf(fptr,"%s\n","eta(1,1) = 0.0\neta(1,2) = 0.0\neta(1,3) = 0.0\neta(2,2) = 0.0\neta(2,3) = 0.0\neta(3,3) = 0.0\n");

	fprintf(fptr,"%s\n",".NO_VELOCITY.");

	fprintf(fptr,"%s %d\n","entry_count =",15);
	n_here=0;
	for (int i = 0; i < n; i++)
    {
		if(((atom[i].interface==interface)&&(atom[i].type==atom_type)&&(atom_type!=0))||((atom[i].interface==interface)&&(atom_type==0)))
		{
			char symb[80]="";
			if(atom[i].type==1) strcat(symb,"Mg");
			if(atom[i].type==2) strcat(symb,"O");
			double r[3];
			double s[3];
			double ratio=1;
			s[0] = atom[i].sx; s[1]=atom[i].sy;s[2]=atom[i].sz;
			//	V3mulM3(s,Hcry,r);
			//  double drig = sqrt(atom[i].disrigistry[0]*atom[i].disrigistry[0]+atom[i].disrigistry[1]*atom[i].disrigistry[1]+atom[i].disrigistry[2]*atom[i].disrigistry[2]);
			/*
			int neigh_config=0;
			if(atom[i].interface==0)
			{
				neigh_config =0;
			}else
			{
				if(atom[i].neigh_config==0)
				{
					neigh_config = atom[i].interface;
				}else
				{
					neigh_config = atom[i].neigh_config+2;
				}
			}
			*/
			double delrr = sqrt(atom[i].delr[0]*atom[i].delr[0]+atom[i].delr[1]*atom[i].delr[1]+atom[i].delr[2]*atom[i].delr[2]);

			fprintf(fptr,"%lf\n%s\n%12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %d %d %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %d\n",
					atom[i].ma,symb,atom[i].sx,atom[i].sy,atom[i].sz,atom[i].delr[3], delrr,atom[i].pe,i*1.0,atom[i].coord,atom[i].interface,atom[i].drig,atom[i].BV*1.0,atom[i].disrigistry[2],atom[i].delr[0],atom[i].ackN*1.0,atom[i].CNA);
	//		fprintf(fptr_addcolor,"%d %d\n",n_here,atom[i].CNA);

			n_here++;
		}

    }

	fclose(fptr);
//	fclose(fptr_addcolor);
	zip(filename);

	return(0);

}


int save_cfg_cut(int pp, double H[3][3], double x_low,double x_high, double y_low, double y_high, double z_low, double z_high, char *suffix)
{

	int n_here=0;
	for(int i=0;i<n;i++)
	{
		if((atom[i].sx>= x_low)&&(atom[i].sx<= x_high)&&(atom[i].sy>= y_low)&&(atom[i].sy<= y_high)&&(atom[i].sz>= z_low)&&(atom[i].sz<= z_high)) n_here++;
	}
	FILE *fptr;
	char filename[80]="dat_atomeye_partconfig",str[20],str1[20],str2[20];
//	char filename2[80]="add_color.";
	if(pp==-1) pp=0;
	strcat(filename,"_");
	strcat(filename,suffix);
	fptr=fopen(filename,"w");

	fprintf(fptr,"%s %d\n\n","Number of particles =",n_here);

	fprintf(fptr,"%s %lf %s\n","H0(1,1) =",H[0][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(1,2) =",H[0][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(1,3) =",H[0][2],"A");

	fprintf(fptr,"%s %lf %s\n","H0(2,1) =",H[1][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(2,2) =",H[1][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(2,3) =",H[1][2],"A");

	fprintf(fptr,"%s %lf %s\n","H0(3,1) =",H[2][0],"A");
	fprintf(fptr,"%s %lf %s\n","H0(3,2) =",H[2][1],"A");
	fprintf(fptr,"%s %lf %s\n\n","H0(3,3) =",H[2][2],"A");

	fprintf(fptr,"%s\n","eta(1,1) = 0.0\neta(1,2) = 0.0\neta(1,3) = 0.0\neta(2,2) = 0.0\neta(2,3) = 0.0\neta(3,3) = 0.0\n");

	fprintf(fptr,"%s\n",".NO_VELOCITY.");

	fprintf(fptr,"%s %d\n","entry_count =",15);
	n_here=0;
	for (int i = 0; i < n; i++)
    {
		if((atom[i].sx>= x_low)&&(atom[i].sx<= x_high)&&(atom[i].sy>= y_low)&&(atom[i].sy<= y_high)&&(atom[i].sz>= z_low)&&(atom[i].sz<= z_high))
		{
			char symb[80]="";
			if(atom[i].type==1) strcat(symb,"Cu");
			if(atom[i].type==2) strcat(symb,"Nb");
			double r[3];
			double s[3];
			double ratio=1;
			s[0] = atom[i].sx; s[1]=atom[i].sy;s[2]=atom[i].sz;
			//	V3mulM3(s,Hcry,r);
			//  double drig = sqrt(atom[i].disrigistry[0]*atom[i].disrigistry[0]+atom[i].disrigistry[1]*atom[i].disrigistry[1]+atom[i].disrigistry[2]*atom[i].disrigistry[2]);
			/*
			int neigh_config=0;
			if(atom[i].interface==0)
			{
				neigh_config =0;
			}else
			{
				if(atom[i].neigh_config==0)
				{
					neigh_config = atom[i].interface;
				}else
				{
					neigh_config = atom[i].neigh_config+2;
				}
			}
			*/
			double delrr = sqrt(atom[i].delr[0]*atom[i].delr[0]+atom[i].delr[1]*atom[i].delr[1]+atom[i].delr[2]*atom[i].delr[2]);

			fprintf(fptr,"%lf\n%s\n%12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %d %d %12.9lf %12.9lf %12.9lf %12.9lf %12.9lf %d\n",
					atom[i].ma,symb,atom[i].sx,atom[i].sy,atom[i].sz,atom[i].delr[3], delrr,atom[i].pe,i*1.0,atom[i].coord,atom[i].interface,atom[i].drig,atom[i].BV*1.0,atom[i].disrigistry[2],atom[i].delr[0],atom[i].ackN*1.0,atom[i].CNA);
	//		fprintf(fptr_addcolor,"%d %d\n",n_here,atom[i].CNA);

			n_here++;
		}

    }

	fclose(fptr);
//	fclose(fptr_addcolor);
	zip(filename);

	return(0);

}

int write_A_interface(char *prefix, int pp, atomic_dat *atom_here, double H_here[3][3])
{
	int n_interface =0;
	for(int i=0;i<n;i++)
	{
		if(atom_here[i].interface==1)
		{
			n_interface++;
		}
	}
	double dummy1;
	int dummy2=2;
	char *s;
	//cout << "\n\n\n\n\n\n"<<prefix<<"\n\n\n\n\n";
	FILE *fptr;
	char filename[80]="",str[20];
	if(pp==-1) pp=0;
	strcat(filename,prefix);
	strcat(filename,".");
	sprintf(str,"%d",pp);
	strcat(filename,str);
	cout << "\n\n\n\n\n\n"<<filename<<"\n\n\n\n\n";
	fptr = fopen(filename,"w");
	for(int k=0;k<10;k++)
	{
	fprintf(fptr,"     %d",dummy2);
	}
	fprintf(fptr,"   %d\n",n_interface);
	fprintf(fptr,"\n\n\n\n\n\n\n\n");
//	cout << n<<" n balue \n";

	fprintf(fptr,"%24.16lf %24.16lf %24.16lf\n",H_here[0][0],H_here[0][1],H_here[0][2]);
	fprintf(fptr,"%24.16lf %24.16lf %24.16lf\n",H_here[1][0],H_here[1][1],H_here[1][2]);
	fprintf(fptr,"%24.16lf %24.16lf %24.16lf\n",H_here[2][0],H_here[2][1],H_here[2][2]);
	fprintf(fptr,"\n\n\n\n\n\n");
	// FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2
	//for (int i=0;i<3;i++) H0[2][i] = H0[2][i];
	// FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2



	//cout << atom[5].rx<<"\t"<<atom[5].ry<<"\t"<<atom[5].rz<<"\n";
	//cout << "after collecting atoms\n";

	for(int i=0;i<n;i++)
	{
		if(atom_here[i].interface==1)
		{
		double r[3],s[3];
		s[0]=atom_here[i].sx;s[1]=atom_here[i].sy;s[2]=atom_here[i].sz;
		//atom[i].ux=atom[i].rx;atom[i].uy=atom[i].ry;atom[i].uz=atom[i].rz;

		V3mulM3(s,H_here,r);

		fprintf(fptr,"   %d %24.16lf %24.16lf %24.16lf\n",atom_here[i].type,r[0],r[1],r[2]);
		}

	}

	fclose(fptr);
	zip(filename);

	return(0);
}

int write_A(const char *prefix, int pp, atomic_dat *atom_here, double H_here[3][3])
{
	double dummy1;
	int dummy2=2;
	char *s;
	//cout << "\n\n\n\n\n\n"<<prefix<<"\n\n\n\n\n";
	FILE *fptr;
	char filename[80]="",str[20];
	if(pp==-1) pp=0;
	strcat(filename,prefix);
	char interfacefilename[80]="";
	strcat(interfacefilename,filename);
	strcat(interfacefilename,"_interface");
	strcat(filename,".");
	sprintf(str,"%d",pp);
	strcat(filename,str);
	cout << "\n\n\n\n\n\n"<<filename<<"\n\n\n\n\n";
	fptr = fopen(filename,"w");
	for(int k=0;k<10;k++)
	{
	fprintf(fptr,"     %d",dummy2);
	}
	fprintf(fptr,"   %d\n",n);
	fprintf(fptr,"\n\n\n\n\n\n\n\n");
//	cout << n<<" n balue \n";

	fprintf(fptr,"%24.16lf %24.16lf %24.16lf\n",H_here[0][0],H_here[0][1],H_here[0][2]);
	fprintf(fptr,"%24.16lf %24.16lf %24.16lf\n",H_here[1][0],H_here[1][1],H_here[1][2]);
	fprintf(fptr,"%24.16lf %24.16lf %24.16lf\n",H_here[2][0],H_here[2][1],H_here[2][2]);
	fprintf(fptr,"\n\n\n\n\n\n");
	// FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2
	//for (int i=0;i<3;i++) H0[2][i] = H0[2][i];
	// FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2



	//cout << atom[5].rx<<"\t"<<atom[5].ry<<"\t"<<atom[5].rz<<"\n";
	//cout << "after collecting atoms\n";

	for(int i=0;i<n;i++)
	{
		double r[3],s[3];
		s[0]=atom_here[i].sx;s[1]=atom_here[i].sy;s[2]=atom_here[i].sz;
		//atom[i].ux=atom[i].rx;atom[i].uy=atom[i].ry;atom[i].uz=atom[i].rz;

		V3mulM3(s,H_here,r);

		fprintf(fptr,"   %d %24.16lf %24.16lf %24.16lf\n",atom_here[i].type,r[0],r[1],r[2]);

	}

	fclose(fptr);
	zip(filename);
	write_A_interface(interfacefilename,pp, atom_here, H_here);
	return(0);
}


int save_gulp(int pp, double H[3][3])
{

	cout << "entered save gulp\n";

	FILE *fptr;
	char filename[80]="dat_gulp.",str[20];
	if(pp==-1) pp=0;
	sprintf(str,"%d",pp);
	strcat(filename,str);
	fptr=fopen(filename,"w");

	fprintf(fptr,"output xyz gulp_output\n\ntitle\n blah blah\nend\n\ncell\n");

 	fprintf(fptr,"%lf %lf %lf 90.0 90.0 90.0\n",H[0][0],H[1][1], H[2][2]);
	fprintf(fptr,"%s\n","Cartesian");

	for (int i = 0; i < n; i++)
    {
		char symb[80]="";
		if(atom==NULL) cout << "atom is null\n";

		//strcat(symb,atom[i].elem);
	// << i << "\t"<<atom[i].elem << "\t"<< symb<<"fgfgf\n";
		double r[3];
		double s[3];
		double ratio=1;
		for (i = 0; i < n; i++)
		{

			double r[3],s[3];

			if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
			if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
			if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;

			if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
			if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
			if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;


			s[0] = atom[i].sx;s[1] = atom[i].sy;s[2] = atom[i].sz;

			V3mulM3(s,H,r);
			double charge_core = 0.0;
			double charge_shell =0.0;

			// hardcoding charge values because some time I am not getting them back from lammps
			if(strcmp(atom[i].elem, "Mg")==0){
				atom[i].charge = 2.0;
			}else
			{
				atom[i].charge = -2.0;
			}

			charge_core = atom[i].charge;
			if(strcmp(atom[i].elem, GULP_SHELL_ATOM)==0)
			{
				fprintf(fptr, "%s %s %12.9lf %12.9lf %12.9lf %12.9lf\n", atom[i].elem, GULP_SHELL_TYPE, r[0], r[1],r[2], GULP_SHELL_CHARGE);
				charge_core = GULP_CORE_CHARGE;
			}

			fprintf(fptr,"%s %12.9lf %12.9lf %12.9lf %12.9lf\n",atom[i].elem, r[0],r[1],r[2], charge_core);
		}

    }

	fclose(fptr);
	zip(filename);

	return(0);

}

int save_xyz_VESTA(int pp, double H[3][3])
{

	cout << "entered save xyz VESTA\n";

	FILE *fptr;
	char filename[80]="dat_VESTA.",str[20];
	if(pp==-1) pp=0;
	sprintf(str,"%d",pp);
	strcat(filename,str);
	strcat(filename,".xyz");
	
	fptr=fopen(filename,"w");

	fprintf(fptr, "%d\n convert to xyz for viewing pleasure\n", n);


	for (int i = 0; i < n; i++)
	{
		char symb[80]="";
		if(atom==NULL) cout << "atom is null\n";

		//strcat(symb,atom[i].elem);
		// << i << "\t"<<atom[i].elem << "\t"<< symb<<"fgfgf\n";
		double r[3];
		double s[3];
		double ratio=1;
		for (i = 0; i < n; i++)
		{

			double r[3],s[3];

			if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
			if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
			if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;

			if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
			if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
			if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;


			s[0] = atom[i].sx;s[1] = atom[i].sy;s[2] = atom[i].sz;

			V3mulM3(s,H,r);
			fprintf(fptr,"%s %12.9lf %12.9lf %12.9lf \n",atom[i].elem, r[0],r[1],r[2]);
		}

	}

	fclose(fptr);
	zip(filename);

	return(0);

}
