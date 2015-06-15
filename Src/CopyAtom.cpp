/*
 * CopyAtom.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#include <CopyAtom.h>

void copy_atomstruct( atomic_dat *main_atom, atomic_dat *copy_atom, int atoms_number)
{


//	cout << "copying...\t";
	for(int i=0;i<atoms_number;i++)
	{

		copy_atom_indv(main_atom,copy_atom,i,i);
	}


//	cout << "done\n";


}


void copy_atom_indv( atomic_dat *main_atom, atomic_dat *copy_atom, int main_number,int copy_number)
{


	int i = main_number;
	int j = copy_number;

	copy_atom[j].rx = main_atom[i].rx;copy_atom[j].ry = main_atom[i].ry;copy_atom[j].rz =main_atom[i].rz;

	copy_atom[j].vx = main_atom[i].vx;copy_atom[j].vy = main_atom[i].vy;copy_atom[j].vz = main_atom[i].vz;

	copy_atom[j].fx = main_atom[i].fx;copy_atom[j].fy = main_atom[i].fy;copy_atom[j].fz = main_atom[i].fz;

	copy_atom[j].sx = main_atom[i].sx;copy_atom[j].sy = main_atom[i].sy;copy_atom[j].sz = main_atom[i].sz;

	copy_atom[j].ux = main_atom[i].ux;copy_atom[j].uy = main_atom[i].uy;copy_atom[j].uz = main_atom[i].uz;

	copy_atom[j].ke = main_atom[i].ke;copy_atom[j].pe = main_atom[i].pe;copy_atom[j].type = main_atom[i].type;

	copy_atom[j].ma = main_atom[i].ma;copy_atom[j].coord = main_atom[i].coord;copy_atom[j].interface = main_atom[i].interface;
	copy_atom[j].disrigistry_number_atoms = main_atom[i].disrigistry_number_atoms;copy_atom[j].disrigistry_number_atoms_negligible = main_atom[i].disrigistry_number_atoms_negligible;

	copy_atom[j].charge = main_atom[i].charge;
	strncpy(copy_atom[j].elem,main_atom[i].elem,sizeof(copy_atom[j].elem));

	copy_atom[j].CNA = main_atom[i].CNA;copy_atom[j].ackN = main_atom[i].ackN;copy_atom[j].neigh_config = main_atom[i].neigh_config;

	copy_atom[j].BV = main_atom[i].BV;
	copy_atom[j].molID = main_atom[i].molID;

	for(int k=0;k<3;k++)
	{
		copy_atom[j].disrigistry[k] = main_atom[i].disrigistry[k];
		copy_atom[j].delr[k] = main_atom[i].delr[k];
	}
	copy_atom[j].delr[3] = main_atom[i].delr[3];

	for(int k=0;k<MAX_COORD;k++)
	{
		copy_atom[j].coord_id[k] =main_atom[i].coord_id[k];
	}



}

void swap_atom( atomic_dat *main_atom, int main_number,int copy_number)
{
	int i = main_number;
	int j = copy_number;
	double tmp_1[13]; for(int i=0;i<13;i++)tmp_1[i]=0;
	int tmp_2;

	tmp_1[0]=main_atom[j].rx;tmp_1[1]=main_atom[j].ry;tmp_1[2]=main_atom[j].rz;
	tmp_1[3]=main_atom[j].sx;tmp_1[4]=main_atom[j].sy;tmp_1[5]=main_atom[j].sz;
	tmp_1[6]=main_atom[j].ux;tmp_1[7]=main_atom[j].uy;tmp_1[8]=main_atom[j].uz;
	tmp_1[9]=main_atom[j].ke;tmp_1[10]=main_atom[j].pe;tmp_2=main_atom[j].type;
	tmp_1[11]=main_atom[j].ma;
	tmp_1[12]=main_atom[j].interface;

	main_atom[j].rx = main_atom[i].rx;main_atom[j].ry = main_atom[i].ry;main_atom[j].rz =main_atom[i].rz;
	main_atom[j].sx = main_atom[i].sx;main_atom[j].sy = main_atom[i].sy;main_atom[j].sz = main_atom[i].sz;
	main_atom[j].ux = main_atom[i].ux;main_atom[j].uy = main_atom[i].uy;main_atom[j].uz = main_atom[i].uz;
	main_atom[j].ke = main_atom[i].ke;main_atom[j].pe = main_atom[i].pe;main_atom[j].type = main_atom[i].type;
	main_atom[j].ma = main_atom[i].ma;
	main_atom[j].interface = main_atom[i].interface;

	main_atom[i].rx=tmp_1[0];main_atom[i].ry=tmp_1[1];main_atom[i].rz=tmp_1[2];
	main_atom[i].sx=tmp_1[3];main_atom[i].sy=tmp_1[4];main_atom[i].sz=tmp_1[5];
	main_atom[i].ux=tmp_1[6];main_atom[i].uy=tmp_1[7];main_atom[i].uz=tmp_1[8];
	main_atom[i].ke=tmp_1[9];main_atom[i].pe=tmp_1[10];main_atom[i].type=tmp_2;
	main_atom[i].ma=tmp_1[11];
	main_atom[i].interface = tmp_1[12];

}
