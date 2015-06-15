/*
 * Manipulations.cpp
 *
 *  Created on: Nov 18, 2009
 *      Author: kedar
 */
#include <Manipulations.h>

#ifndef LEN_H
#define LEN_H    4
#endif
double coords_to_insert[MAX_COORD][LEN_H];
int *inital_coords;

bool fill_coords_to_insert(atomic_dat *atom_h,int index, double H_now[3][3])//,double *toinsert_coords[4])
{
	int *arr;
	arr = (int *) malloc((MAX_COORD)*sizeof(int));
	for(int k=0;k<MAX_COORD;k++)
	{
		arr[k] = -1;
	}

	for(int k = 0; k<MAX_COORD;k++)
	{
		int t = atom_h[index].coord_id[k];

		if((t>-1)&&(!check_repeat(t,arr,MAX_COORD)))
		{
			double s[3],r1[3],r2[3];
			s[0]=s[1]=s[2]=r1[0]=r1[1]=r1[2]=r2[0]=r2[1]=r2[2]=0.0;

			s[0] = atom_h[index].sx - atom_h[t].sx;
			s[1] = atom_h[index].sy - atom_h[t].sy;
			s[2] = atom_h[index].sz - atom_h[t].sz;
			s[0] = s[0]-(int)(s[0]*2);
			s[1] = s[1]-(int)(s[1]*2);
			s[2] = s[2]-(int)(s[2]*2);
			V3mulM3(s,H_now,r1);
			for(int l=k+1;l<MAX_COORD;l++)
			{

				int t1 = atom_h[index].coord_id[l];
				/*
				if((atom[t1].coord !=atom[index].coord)||(atom[t].coord !=atom[index].coord))
				{
					arr[l] = t1;
					  arr[k] = t;
				}
				*/
				if((t1>-1)&&(!check_repeat(t1,arr,MAX_COORD)))
				{
					s[0] = atom_h[index].sx - atom_h[t1].sx;
					s[1] = atom_h[index].sy - atom_h[t1].sy;
					s[2] = atom_h[index].sz - atom_h[t1].sz;
					s[0] = s[0]-(int)(s[0]*2);
					s[1] = s[1]-(int)(s[1]*2);
					s[2] = s[2]-(int)(s[2]*2);
					V3mulM3(s,H_now,r2);
					double r[3]; r[0] = r[1] = r[2] = 0.0;
					for(int count=0;count<3;count++)
						{
							//cout << r1[count]<<"\t"<<r2[count]<<"\t";
							r[count] = r1[count]+r2[count];
						//	cout << r[count]<<" HELLO\n";
						}
					double dist=0;
					for(int count=0;count<3;count++) dist+=r[count]*r[count];
					if(sqrt(dist)<1) { arr[l] = t1;arr[k] = t;cout << "OK\n";}
					//cout << index<<"\t"<<t<<"\t"<<t1<<"\t"<<dist<<"\n";
				//	cout <<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\n";
				}
			}
		}
	}
	double H_now_inv[3][3];
	M3inv(H_now,H_now_inv);
	bool ret_val = false;
	int things_in = 0;
	for(int i = 0;i<MAX_COORD;i++)
	{
		int t = atom_h[index].coord_id[i];
		//if((t>-1)&&(!check_repeat(t,arr,MAX_COORD))&&(atom[t].coord==12))
		if((t>-1)&&(!check_repeat(t,arr,MAX_COORD)))
		{
			// specifically for jog thing - only for this case

			coords_to_insert[i][0] = 1;
			double s[3],r1[3],r2[3];
			s[0]=s[1]=s[2]=r1[0]=r1[1]=r1[2]=r2[0]=r2[1]=r2[2]=0.0;

			s[0] = atom_h[index].sx - atom_h[t].sx;
			s[1] = atom_h[index].sy - atom_h[t].sy;
			s[2] = atom_h[index].sz - atom_h[t].sz;
			s[0] = s[0]-(int)(s[0]*2);
			s[1] = s[1]-(int)(s[1]*2);
			s[2] = s[2]-(int)(s[2]*2);
			V3mulM3(s,H_now,r1);
			r2[0] = -1.0*r1[0];r2[1] = -1.0*r1[1]; r2[2] = -1.0*r1[2];
			V3mulM3(r2,H_now_inv,s);
			s[0] = atom_h[index].sx - s[0];
			s[1] = atom_h[index].sy - s[1];
			s[2] = atom_h[index].sz - s[2];


			if(s[0]>=1) s[0]=s[0]-1; if(s[0]<0) s[0]=1+s[0];
			if(s[1]>=1) s[1]=s[1]-1; if(s[1]<0) s[1]=1+s[1];
			if(s[2]>=1) s[2]=s[2]-1; if(s[2]<0) s[2]=1+s[2];


			//cout << r2[0]<<"\t"<<r2[1]<<"\t"<<r2[2]<<"\t"<<s[0]<<"\t"<<s[1]<<"\t"<<s[2]<<"\t"<< t <<" HU HU\n";
			coords_to_insert[i][1] = s[0];coords_to_insert[i][2]=s[1];coords_to_insert[i][3]=s[2];

			if(coords_to_insert[i][1]>=1.0) coords_to_insert[i][1] = coords_to_insert[i][1]-1;
			if(coords_to_insert[i][2]>=1.0) coords_to_insert[i][2] = coords_to_insert[i][2]-1;
			if(coords_to_insert[i][3]>=1.0) coords_to_insert[i][3] = coords_to_insert[i][3]-1;

			if(coords_to_insert[i][1]<0.0) coords_to_insert[i][1] = coords_to_insert[i][1]+1;
			if(coords_to_insert[i][2]<0.0) coords_to_insert[i][2] = coords_to_insert[i][2]+1;
			if(coords_to_insert[i][3]<0.0) coords_to_insert[i][3] = coords_to_insert[i][3]+1;
			ret_val = true;
		//	things_in++;if(things_in+atom[index].coord>=12) i =MAX_COORD;
			cout << "yes\n";
		}
	}
	free(arr);
	return ret_val;
}

int insert_atoms_recur(int coord_num_limit, atomic_dat *atom_now, int n_now, double H_now[3][3],bool consider, bool first_iter)
//Valid for fcc only - takes central symmetry into consideration
// very expensive - multiple neighbor lists are computed
{
	if(first_iter)
	{
		for(int i = 0;i<n;i++)
		{
			atom_now[i].BV = atom_now[i].coord;
		}

	}
	compute_CNA_and_others(atom_now,n_now, H_now,false);

	bool recur = false;
	for(int i =0;i<n_now;i++)
	{
		bool consider2= false;
		if(consider)
		{
			if((atom[i].coord==11)||(atom[i].coord >=13))
			{
				int deficit = 0;
				int hcp_n = 0;
				for(int k = 0; k<MAX_COORD;k++)
				{
					if(atom[i].coord_id[k]>-1)
					{
						if(atom[atom[i].coord_id[k]].coord!=12) deficit++;
						if(atom[atom[i].coord_id[k]].CNA == HCP) hcp_n++;
					}
				}
				if(((deficit>3)&&(hcp_n==0))) consider2 = true;
			}
		}
	//	if(((atom_now[i].coord==atom_now[i].BV)&&(atom_now[i].coord== coord_num_limit))||consider2)
		if((atom_now[i].coord<=coord_num_limit)||consider2)
		{
			cout << "entered an atom "<<i<<"\t"<<atom_now[i].coord<<"\n";
			for(int k=0;k<MAX_COORD;k++)
			{
				for(int l=0;l<LEN_H;l++)
				{
					coords_to_insert[k][l] = -1;
				}
			}
			//cout <<"getting in\n";
			recur = fill_coords_to_insert(atom_now,i,H_now);//,coords_to_insert);
			//cout << "got out\n";
			if(recur)
			{
				//cout << "HELLLLOOO CAMEEEEEE HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n";
				for(int k = 0; k<MAX_COORD;k++)
				{
					if(coords_to_insert[k][0] ==1)
					{
						int len = 1;
						double arr[len][4];
						arr[0][0] = coords_to_insert[k][1];
						arr[0][1] = coords_to_insert[k][2];
						arr[0][2] = coords_to_insert[k][3];
						arr[0][3] = 1;

						insert_atoms(arr,len);
						cout <<n<<"\n";

					}
				}
				i = n_now+1;
			}
		}
	}

	if(recur)
	{
		n_now = n;
		cout << "going into second loop\n";
		//tests
		//compute_CNA_and_others(atom_now,n_now, H_now,false);
		//save_cfg(10,H_now);
		insert_atoms_recur(coord_num_limit, atom, n_now, H_now,consider,false);

	}

}


void shift(atomic_dat *atom_now,int n_now,double H_now[3][3],double arr_shift[3])
{
	for(int tag =0; tag<n_now;tag++)
	{
		atom_now[tag].sx +=arr_shift[0];
		atom_now[tag].sy +=arr_shift[1];
		atom_now[tag].sz +=arr_shift[2];
		if(atom_now[tag].sx>=1) atom_now[tag].sx=atom_now[tag].sx-1;
		if(atom_now[tag].sx<0) atom_now[tag].sx=1+atom_now[tag].sx;


		if(atom[tag].sy>=1) atom[tag].sy=atom_now[tag].sy-1;
		if(atom[tag].sy<0) atom[tag].sy=1+atom_now[tag].sy;

		if(atom_now[tag].sz>=1) atom_now[tag].sz=atom_now[tag].sz-1;
		if(atom_now[tag].sz<0) atom_now[tag].sz=1+atom_now[tag].sz;
	}
}


void save_atoms(atomic_dat *atom_now, int n_now, double H_now[3][3])
{
	compute_CNA_and_others(atom_now,n_now, H_now,true);
	int add_atoms = 0;
	for(int i =0; i<n_now;i++)
	{
		if((atom_now[i].ackN ==6)&&(atom_now[i].coord<11)&&(atom_now[i].sy>0.5))
		{
			//atom_now[i].interface = 1;
			add_atoms++;
		}
	}
	cout << "number of atoms to be considered\t"<< add_atoms<<"\n";
	double arr_add[add_atoms][4];
	add_atoms = 0;
	for(int i =0; i<n_now;i++)
	{


		if((atom_now[i].ackN ==6)&&(atom_now[i].coord<11)&&(atom_now[i].sy>0.5))
		{
			//atom_now[i].interface = 1;
			arr_add[add_atoms][0] = atom_now[i].sx+(-0.588665+0.563803);
			arr_add[add_atoms][1] = atom_now[i].sy+(-0.842297+0.850116);
			arr_add[add_atoms][2] = atom_now[i].sz+(-0.34135+0.341399);
			arr_add[add_atoms][3] = 1;
			add_atoms++;
		}
	}
	cout << "start adding "<<add_atoms<< " atoms\t";
	insert_atoms(arr_add,add_atoms);
	cout << "n changed from\t"<< n_now <<" to "<< n<<"\t";
	cout << "done\n";

}
