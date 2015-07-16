/*
 * WriteData.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */
#include <writedata.h>
bool LAMMPS_CHARGE = false;
bool LAMMPS_MOLECULE = false;


int save_lammps(int pp,simcell &MDcell)
{
	if (MDcell.n <= 0)
	{
		cout << "not saved because there are no atoms in the system";
		return (-1);
	}
	if(((MDcell.Hcry[0][1])>0.0000)||((MDcell.Hcry[0][2])>0.0000)||((MDcell.Hcry[1][2])>0.0000)||((MDcell.Hcry[0][1])<0.0000)||((MDcell.Hcry[0][2])<0.0000)||((MDcell.Hcry[1][2])<0.0000))
	{
		cout << "the lammps save value of pp\t"<<pp<<" not done because H matrix is not in the desired format\n";
		return(-1);
	}

	double r1[3],s1[3];
	double limits[3][2];

	limits[0][0] = 1e31;
	limits[1][0] = 1e31;
	limits[2][0] = 1e31;
	limits[0][1] = -1e31;
	limits[1][1] = -1e31;
	limits[2][1] = -1e31;


	//int i ;
	FILE *fptr;
	char filename[30] = "dat_lammps.", str[20];
	sprintf(str,"%d",pp);
	strcat(filename,str);
	fptr=fopen(filename,"w");
	fprintf(fptr, "%s\n\n","input file for LAMMPS");
	fprintf(fptr, "%d %s\n",MDcell.n, " atoms");
	bool bonds_exist = false;
	bool angles_exist = false;
	if(MDcell.bonds.size()>0)
	{
		bonds_exist = true;
		fprintf(fptr, "%d %s\n", MDcell.bonds.size(), " bonds");
	}
	if(MDcell.angles.size()>0)
	{
		angles_exist = true;
		fprintf(fptr, "%d %s\n", MDcell.angles.size(), " angles");
	}
	fprintf(fptr, "\n");

	fprintf(fptr, "%d %s\n", MDcell.atoms_info.n_types, "atom types");
	if(bonds_exist) fprintf(fptr, "%d %s\n\n", MDcell.atoms_info.element_bonds.size(), "bond types");
	if(angles_exist) fprintf(fptr, "%d %s\n\n", MDcell.atoms_info.element_angles.size(), "angle types");
	fprintf(fptr, "\n");

	fprintf(fptr, "%12.9f %12.9lf %s\n",0.0,MDcell.Hcry[0][0], " xlo xhi");
	fprintf(fptr, "%12.9f %12.9lf %s\n",0.0,MDcell.Hcry[1][1], " ylo yhi");
	fprintf(fptr, "%12.9f %12.9lf %s\n\n",0.0,MDcell.Hcry[2][2], " zlo zhi");
	fprintf(fptr, "%12.9f %12.9lf %12.9lf %s\n\n",MDcell.Hcry[1][0],MDcell.Hcry[2][0],MDcell.Hcry[2][1], " xy xz yz");

	fprintf(fptr, "%s\n\n", "Masses");
	for (int ii=1;ii<=MDcell.atoms_info.n_types;ii++)
	{
		fprintf(fptr, "%d %lf\n", ii,MDcell.atoms_info.element_weights[ii]);
	}
	fprintf(fptr, "%s\n\n", "Atoms");
	for(auto iter = MDcell.atomdata.begin(); iter !=MDcell.atomdata.end(); ++iter)
	{
		atomic_dat localatom = iter->second;
		int tag = iter->first;
		for (int i =0 ; i < 3; i++)
		{
			if(localatom.s[i]>=1.0) localatom.s[i] = localatom.s[i]-1;
			if(localatom.s[i]<0.0) localatom.s[i] = localatom.s[i]+1;
		}
		V3mulM3(localatom.s, MDcell.Hcry, localatom.r);
		fprintf(fptr,"%d ",tag);
		if(LAMMPS_MOLECULE)
		{
			fprintf(fptr,"%d ", localatom.molID);
		}

		fprintf(fptr,"%d ",localatom.type);

		if(LAMMPS_CHARGE)
		{
			fprintf(fptr,"%lf ",localatom.charge);
		}

		fprintf(fptr,"%12.9lf %12.9lf %12.9lf\n",localatom.r[0],localatom.r[1],localatom.r[2]);
	}
	if(bonds_exist)
	{
		fprintf(fptr, "%s\n\n", "Bonds");
		int counter = 1;
		for(auto bonditer = MDcell.bonds.begin(); bonditer != MDcell.bonds.end(); ++bonditer)
		{
			fprintf(fptr, "%d %d %d %d\n", counter, std::get<0>(*bonditer), std::get<1>(*bonditer), std::get<2>(*bonditer));
			counter++;
		}
	}
	if(angles_exist)
	{
		fprintf(fptr, "%s\n\n", "Angles");
		int counter = 1;
		for(auto angleiter = MDcell.angles.begin(); angleiter != MDcell.angles.end(); ++angleiter)
		{
			fprintf(fptr, "%d %d %d %d %d\n", counter, std::get<0>(*angleiter), std::get<1>(*angleiter), std::get<2>(*angleiter),std::get<3>(*angleiter));
			counter++;
		}
	}

	fprintf(fptr, "%s\n\n", "Velocities");
	for(auto iter = MDcell.atomdata.begin(); iter !=MDcell.atomdata.end(); ++iter)
	{
		atomic_dat localatom = iter->second;
		int tag = iter->first;
		fprintf(fptr,"%d %12.9lf %12.9lf %12.9lf\n",tag,localatom.v[0],localatom.v[1],localatom.v[2]);
	}
	fclose(fptr);
	zip(filename);

	return(0);
}
/*
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
*/
