/*
 * CNAandOthers.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar, testing something on Sept 30 2014
 */
#include <CNAandOthers.h>
bool compute_meandistance(atomic_dat *atom_now,int atom_num, double H_now[3][3], int nearest_n0[MAX_COORD],int nearest_n1[MAX_COORD])
{
	double mean_distance=0;
	double r[3],s[3];
	double rsqs[MAX_COORD][2];

	bool success = true;

	for(int i=0;i<MAX_COORD;i++)
	{
		rsqs[i][0] = 1e31-1;
		rsqs[i][1]=-1;
		nearest_n0[i]=-1;
		nearest_n1[i]=-1;
	}

	for(int j=0;j<MAX_COORD;j++)
	{
		int k = atom[atom_num].coord_id[j];
		if(k>-1)
		{
			s[0] = atom_now[atom_num].sx-atom_now[k].sx;
			s[1] = atom_now[atom_num].sy-atom_now[k].sy;
			s[2] = atom_now[atom_num].sz-atom_now[k].sz;
			s[0] = s[0]-(int)(s[0]*2);
			s[1] = s[1]-(int)(s[1]*2);
			s[2] = s[2]-(int)(s[2]*2);

			V3mulM3(s,H_now,r);
			rsqs[j][0] = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
			rsqs[j][1] = k;
		}
	}


	for(int i = 0;i< (MAX_COORD-1);i++)
	{
		for(int j= i;j<MAX_COORD;j++)
		{
			if(rsqs[i][0]> rsqs[j][0])
			{
				double a = rsqs[i][0];
				double b = rsqs[i][1];
				rsqs[i][0] = rsqs[j][0];
				rsqs[i][1] = rsqs[j][1];
				rsqs[j][0] = a;
				rsqs[j][1] = b;
			}
		}
    }

	mean_distance = 0.0;
	for(int i=0;i<6;i++)
	{
		if(rsqs[i][1]<0)
		{
			success = false; return(success);
		}else
		{
			mean_distance = mean_distance+rsqs[i][0];
		}
	}

	mean_distance = mean_distance/6.0;

	for(int j=0;j<MAX_COORD;j++)
	{
		if(rsqs[j][1]>-1)
		{
			if(rsqs[j][0]<1.45*mean_distance) nearest_n0[j]=rsqs[j][1];
			if(rsqs[j][0]<1.55*mean_distance) nearest_n1[j]=rsqs[j][1];
		}
	}

	return(success);
}
void compute_acklandnotation( atomic_dat *atom_now, int n_now, double H_now[3][3])
{
	int chi[8];
	int jnab,jbeg,jend,i,j;
	double rxij,ryij,rzij,rijsq;
//	double *distsq;
	int nearest_n0[MAX_COORD],nearest_n1[MAX_COORD];
	int n0, n1;

	cout << "entering compute_ackN\t"<<n_now<<"\n";

	for(int i=0;i<n_now;i++)
	{
		//cout <<"\t"<<i<<"\n";
		bool compute_ok = compute_meandistance(atom_now,i, H_now,nearest_n0,nearest_n1);
		n0=0;n1=0;
		for(int j=0;j<MAX_COORD;j++)
		{
			if(nearest_n0[j]>-1) n0++;
			if(nearest_n1[j]>-1) n1++;
		}

		chi[0] = chi[1] = chi[2] = chi[3] = chi[4] = chi[5] = chi[6] = chi[7] = 0;
		//cout <<"\t"<<i<<"\t"<<n0<<"\t"<<n1<<"\n";
		for(int j=0;j<n0;j++)
		{
			double sij[3],rij[3];
			sij[0] = atom_now[i].sx-atom_now[nearest_n0[j]].sx;
			sij[1] = atom_now[i].sy-atom_now[nearest_n0[j]].sy;
			sij[2] = atom_now[i].sz-atom_now[nearest_n0[j]].sz;
			sij[0] = sij[0]-(int)(sij[0]*2);
			sij[1] = sij[1]-(int)(sij[1]*2);
			sij[2] = sij[2]-(int)(sij[2]*2);

			V3mulM3(sij,H_now,rij);
			double normij = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
			for(int k=j+1;k<n0;k++)
			{
				double sik[3],rik[3];
				sik[0] = atom_now[i].sx-atom_now[nearest_n0[k]].sx;
				sik[1] = atom_now[i].sy-atom_now[nearest_n0[k]].sy;
				sik[2] = atom_now[i].sz-atom_now[nearest_n0[k]].sz;
				sik[0] = sik[0]-(int)(sik[0]*2);
				sik[1] = sik[1]-(int)(sik[1]*2);
				sik[2] = sik[2]-(int)(sik[2]*2);

				V3mulM3(sik,H_now,rik);
				double normik = sqrt(rik[0]*rik[0]+rik[1]*rik[1]+rik[2]*rik[2]);

				double bond_angle = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) / (normij*normik);

				// Histogram for identifying the relevant peaks
				if ((-1.00000001 <= bond_angle) && (bond_angle < -0.945)) { chi[0]++; }
				else if ((-0.945 <= bond_angle) && (bond_angle < -0.915)) { chi[1]++; }
				else if ((-0.915 <= bond_angle) && (bond_angle < -0.755)) { chi[2]++; }
				else if ((-0.755 <= bond_angle) && (bond_angle < -0.195)) { chi[3]++; }
				else if ((-0.195 <= bond_angle) && (bond_angle < 0.195)) { chi[4]++; }
				else if ((0.195 <= bond_angle) && (bond_angle < 0.245)) { chi[5]++; }
				else if ((0.245 <= bond_angle) && (bond_angle < 0.795)) { chi[6]++; }
				else if ((0.795 <= bond_angle) && (bond_angle < 1.0000001)) { chi[7]++; }
				else { cout << "bond angle not in the range\t"<< i<<"\t"<<nearest_n0[j]<<"\t"<<nearest_n0[k]<<"\t"<<bond_angle<<"\n";}
			}
		}


		double delta_bcc = 0.35*chi[4]/(double)(chi[5]+chi[6]-chi[4]);
		double delta_cp = fabs(1.-(double)chi[6]/24.);
		double delta_fcc = 0.61*(fabs((double)(chi[0]+chi[1]-6.))+(double)chi[2])/6.;
		double delta_hcp = (fabs((double)chi[0]-3.)+fabs((double)chi[0]+(double)chi[1]+(double)chi[2]+(double)chi[3]-9.))/12.;
		// Identification of the local structure according to the reference
		if (chi[0] == 7)       { delta_bcc = 0.; }
		else if (chi[0] == 6)  { delta_fcc = 0.; }
		else if (chi[0] <= 3)  { delta_hcp = 0.; }

		if (chi[7] > 0.00000)
			atom_now[i].ackN = UNKNOWN;
		else
			if (chi[4] < 3.00000)
			{
				if (n1 > 13 || n1 < 11)
					atom_now[i].ackN = UNKNOWN;
				else
					atom_now[i].ackN = ICO;
			} else
				if (delta_bcc <= delta_cp)
				{
					if (n1 < 11)
						atom_now[i].ackN = UNKNOWN;
					else
						atom_now[i].ackN = BCC;
				} else
					if (n1 > 12 || n1 < 11)
						atom_now[i].ackN = UNKNOWN;
		else
			if (delta_fcc < delta_hcp)
				atom_now[i].ackN = FCC;
	    else
			atom_now[i].ackN = HCP;
		/*
		 if((i==22537)||(i==22501))
		 {
			 cout << "***************************\t"<<i<<"\t"<<atom_now[i].ackN<<"\t"<<atom_now[i].CNA<<"\n";
		 }
		 */
		if(((atom_now[i].CNA==FCC_UNKNOWN)||(atom_now[i].CNA==BCC_UNKNOWN))) atom_now[i].CNA = atom_now[i].ackN;
		//if(atom_now[i].interface==0) atom_now[i].ackN=-1;

		/*
			if((i==22537)||(i==22501))
		 {
				cout << "***************************\t"<<i<<"\t"<<atom_now[i].ackN<<"\t"<<atom_now[i].CNA<<"\n";
		 }

		 */
	}



}

int common_neighbors(int i,int j, int atom_type,int *return_value,double len_z)
{

	//This also determines the number of common neighbors in the same plane for ring analysis -- the last two input parameters are for that

	//determine if i and j are in the same plane. If they are not assign neigh=-1  else perform the analysis
	bool inplane = false;
	double del_inplane = 1.1;
	int limit = 3*MAX_COORD;
	if(((fabs(atom[i].sz-atom[j].sz)*len_z)<del_inplane)||((atom[i].type==atom[j].type)&&(atom[i].interface==atom[j].interface)&&(atom[i].interface>0)))
	{

		inplane=true;
		(*return_value) = 0;
	}

	//common neighbors between two atoms
	int first=0, second=0, third=0;
	int store_cmn_neighs[limit];
	for(int x=0;x<10;x++) store_cmn_neighs[x]=-1;
	for(int k=0;k<MAX_COORD;k++)
	{
		if(atom[i].coord_id[k]>-1)
		{
			for(int l=0;l<MAX_COORD;l++)
			{
				if((atom[j].coord_id[l]>-1)&&(atom[i].coord_id[k]==atom[j].coord_id[l]))
				{
					store_cmn_neighs[first] = atom[i].coord_id[k];
					int rr = atom[i].coord_id[k];
					first++;
					if(inplane)
					{
						if(((fabs(atom[i].sz-atom[rr].sz)*len_z)<del_inplane)||((atom[i].type==atom[rr].type)&&(atom[i].interface==atom[rr].interface)&&(atom[i].interface>0)))
							(*return_value)++; //should change to 0.75
					}
				}
			}
		}
	}

	int store_cmn_neighs2[limit];
	for(int x=0;x<10;x++) store_cmn_neighs2[x]=-1;
	int counter = 0;
	for(int p=0;p<first;p++)
	{
		int r = store_cmn_neighs[p];
		for(int k=0;k<MAX_COORD;k++)
		{
			if(atom[r].coord_id[k]>-1)
			{
				for(int q=p+1;q<first;q++)
				{
					int s = store_cmn_neighs[q];
					if(atom[r].coord_id[k]==s)
					{
						if(!check_repeat(r,store_cmn_neighs2,limit)) {store_cmn_neighs2[counter] = r;counter++;}
						if(!check_repeat(s,store_cmn_neighs2,limit)) {store_cmn_neighs2[counter+1] = s;counter++;}
						second++;
					}
				}
			}
		}
	}

	if(atom_type==1)
	{
		if(counter==4)
		{
			counter=1;
		}else
		{
			if(counter==3) counter = 2;
		}

	}else
	{
		counter=1;
	}

	return (first*100+second*10+counter);


}


void compute_CNA( atomic_dat *atom_now, int n_now, double H_now[3][3])
{

	// For computing the CNA or any similar method, we need to consider the number of neighbors for BCC to be 14 and not 8
	// Without assuming anything, we can safely reset the neighbor lists to calculate the neighbors for CNA with new coordinates
	// But we should set the neighbor lists back to the original ones immediately after

	double len_z = sqrt(H_now[2][0]*H_now[2][0]+H_now[2][1]*H_now[2][1]+H_now[2][2]*H_now[2][2]);
	cout << "finished entering compute_CNA\t"<<n_now<<"\n";

	for(int i=0;i<n_now;i++)
	{
		int pair_data[6];
		for(int j=0;j<6;j++)
		{
			pair_data[j]=0;
		}
	//		cout << atom_now[13346].type<<"\t"<<atom_now[34691].type<<" ************\n";
	//	cout << i<<"\n";
		if((atom_now[i].coord<MAX_COORD))
		{
		//	cout << "entered first here\t"<<i<<"\n";

			int cmn_nbs_num = 0,cmn_nbs_bonds = 0, longest_nb_bonds = 0;
			for(int j=0;j<MAX_COORD;j++)
			{
				int neigh=-1;
				atom_now[i].neigh_bonds[j]= neigh;
				if(atom_now[i].coord_id[j]>-1)
				{
			//				cout <<i<<"\t" << "herere and\t";
					int struct_set = 1000;

					//This also determines the number of common neighbors in the same plane for ring analysis -- the last two input parameters are for that


					struct_set = struct_set+common_neighbors(i,atom_now[i].coord_id[j],atom_now[i].type,&neigh,len_z);
					atom_now[i].neigh_bonds[j]= neigh;

					pair_data[5] = pair_data[5]+1;
					if(struct_set == 1421) {
						pair_data[0] = pair_data[0]+1;
					}else{
						if(struct_set == 1422) {
							pair_data[1] = pair_data[1]+1;
						}else{
							if(struct_set == 1441) {
								pair_data[2] = pair_data[2]+1;
							}else{
								if(struct_set== 1661) {
									pair_data[3] = pair_data[3]+1;
								}else{
									pair_data[4] = pair_data[4]+1;
								}
							}
						}
					}

				}

			}
		}

		if(pair_data[0]==12)
		{
			atom_now[i].CNA = FCC;
		}else
		{
			if((pair_data[0]==6)&&(pair_data[1]==6))
			{
				atom_now[i].CNA = HCP;
			}else
			{
				if((pair_data[2]==6)&&(pair_data[3]==8))
				{
					atom_now[i].CNA = BCC;
				}else
				{

					if(atom_now[i].coord==12)
					{
						atom_now[i].CNA = FCC_UNKNOWN;
					}else
					{
						if(atom_now[i].coord==14)
						{
							atom_now[i].CNA = BCC_UNKNOWN;
						}else
						{
							atom_now[i].CNA = UNKNOWN;
						}
					}

				}
			}
		}

		//atom_now[i].CNA = pair_data[0]*10000+pair_data[1]*1000+pair_data[2]*100+pair_data[3]*10+pair_data[4];

	//	if(atom_now[i].interface==0) atom_now[i].CNA=-1;
	}

}



void coord_number(double H1[3][3])
{

	int jbeg,jend,jnab;
	int *coord;
	coord = (int *) malloc((n)*sizeof(int));
	double sxi,syi,szi,sij[3],rij[3],rijsq;

	int total_close = 1000;
	int close_atoms[total_close];
	int counter_close=0;

	cout << n_types << " neighbor distances are \n";
		for(int j = 0; j < n_types; j++)
		{
			cout <<j << " "<<sqrt(rcoordsq[6][j])<<"\n";
		}

	for(int i=0;i<total_close;i++)
	{
		close_atoms[i]=-1;
	}

	//cout <<"in\n";
	for(int i=0;i<n;i++)
	{
		atom[i].coord = 0;
		atom[i].interface=0;
		atom[i].neigh_pe=0.0;
		for(int j=0;j<MAX_ELEMS;j++)
		{
			atom[i].coord_agg[j] = 0;
			for(int k=0; k<SINGLE_MAX_COORD; k++)
			{
//				atom[i].coord_info[j][k] = -1;

				atom[i].coord_id[ j*SINGLE_MAX_COORD+k ] = -1;
			}
		}

	}


	for(int i=0;i<n;i++)
	{
		jbeg = nbr_ptr[i];
		jend = nbr_ptr1[i];

		sxi = atom[i].sx;
		syi = atom[i].sy;
		szi = atom[i].sz;

	  for (jnab = jbeg; jnab <= jend; jnab++)
		{
			int j = nbr_lst[jnab];
			sij[0] = sxi - atom[j].sx;
			sij[1] = syi - atom[j].sy;
			sij[2] = szi - atom[j].sz;
			sij[0] = sij[0]-(int)(sij[0]*2);
			sij[1] = sij[1]-(int)(sij[1]*2);
			sij[2] = sij[2]-(int)(sij[2]*2);

			V3mulM3(sij,H1,rij);
			rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];


			if(rijsq<=0.5)
			{
				/*
				cout << "*************************\n";
				cout << i<< " "<<atom[i].type<<" "<< atom[i].rx << " "<< atom[i].ry<< " "<< atom[i].rz<<"\n";
				cout << j<< " "<<atom[j].type<<" "<< atom[j].rx << " "<< atom[j].ry<< " "<< atom[j].rz <<"\n";
				cout << rij[0]<<" "<<rij[1]<<" "<<rij[2]<<"\n";
				cout << "*************************\n";
				*/	
				if(counter_close <total_close)
				{

					if(!check_repeat(i,close_atoms,total_close))
					{
						close_atoms[counter_close]=i;counter_close++;
					}

					if(!check_repeat(j,close_atoms,total_close))
					{
						close_atoms[counter_close]=j;counter_close++;
					}
				}else
				{
					cout << "More than " << counter_close << " atoms too close.. skipping the rest\n";
 				}
			}

			if (rijsq <= rcoordsq[atom[i].type-1][atom[j].type-1])
			{

				int i_type = atom[i].type-1;
				int j_type = atom[j].type-1;

				int i_c = atom[i].coord; int j_c=atom[j].coord;

				if(!check_repeat(j,atom[i].coord_id,MAX_COORD))
				{
					atom[i].coord_id[i_c] = j;
					atom[i].coord = atom[i].coord + 1;
					atom[i].coord_agg[j_type]++;
					atom[i].neigh_pe += atom[j].pe;
				}

				if(!check_repeat(i,atom[j].coord_id,MAX_COORD))
				{
					atom[j].coord_id[j_c] = i;
					atom[j].coord = atom[j].coord + 1;
					atom[j].coord_agg[i_type]++;
					atom[j].neigh_pe += atom[i].pe;
				}

				if(latt_cutoff[atom[i].type-1][1]!=latt_cutoff[atom[j].type-1][1])
				{
					atom[i].interface = 1;
					atom[j].interface = 1;
					if(atom[i].type==1) interface_Cu_n++;
					if(atom[j].type==1) interface_Cu_n++;
				}
				//					if((i==23232)||(i==3423)||(j==23232)||(j==3423)) cout << i<<"\t"<<atom[i].type<<"\t"<<atom[i].coord<<"\n";
			}
		}
	}

	for(int i=0;i<n;i++)
	{
		if (atom[i].interface==1)
		{
			for(int j=0;j<MAX_COORD;j++)
			{
				int r = atom[i].coord_id[j];

				if(r>=0)
				{
					//	if(i==21142) cout << r<<" before\t"<<atom[r].interface<<"\n";
					if(atom[r].interface!=1)
					{
						atom[r].interface=2;
						//						if(atom[r].type==1) interface_Cu_n--;
					}
					//if(i==21142) cout << r<<" after\t"<<atom[r].interface<<"\n";
				}
			}
		}
		if((atom[i].type==7) || (atom[i].type==11))
		{
			//2, 4, 7, 9, 11
			cout << "atom id "<< i <<" coord num is "<< atom[i].coord << " and ";
			cout << atom[i].coord_agg[1] << " "<< atom[i].coord_agg[3]<< " ";
			cout << atom[i].coord_agg[6]<< " "<< atom[i].coord_agg[8]<<" ";
			cout << atom[i].coord_agg[10]<<" " << atom[i].pe << " "<<atom[i].type<<"\n";
			/*
			for(int j=0;j<11; j++)
			{
				cout << j+1 << " "<< atom[i].coord_agg[j]<<" ";
			}
			cout <<"\n";
			*/
		}
	}

	if(counter_close>0)
	cout << "total close atoms are\t"<<counter_close<<"**************\n";
	free(coord);
}

void fill_bonds_etc()
{
	cout << "filling bonds information\n";
	cout << "also clears previously existing bonds";
	cout << "this reads the data from `atom`'s coord function'\n";
	cout << "please ensure that method is called before this function is called";
	cout << "\nThere are "<< element_bonds.size()<< "varieties of bonds";
	cout << "they are:\n";

	bonds.clear();
	angles.clear();

	for(auto it = element_bonds.begin(); it !=element_bonds.end(); it=it+3)
	{
		cout << *it<< "\t"<<*(it+1)<<"\t"<<*(it+2)<<"\n";
	}
	cout <<"\nlist of available angles:\n";
	for(auto it = element_angles.begin(); it !=element_angles.end(); it=it+4)
	{
		cout << *it<< "\t"<<*(it+1)<<"\t"<<*(it+2)<<"\t"<<*(it+3)<<"\n";
	}
	for(int i=0;i<n;i++)
	{

		for(int jj=0;jj<atom[i].coord; jj++)
		{
			int j = atom[i].coord_id[jj];
			//if(i==760) cout << "760 coords are "<< atom[i].coord<<"\n";
			//if(i==760)
			//{
			//	cout << " neighbor for 760 is "<< i<<"\t"<< atom[j].type<<"\n";
			//}


			for(auto it = element_bonds.begin(); it!=element_bonds.end(); it=it+3)
			{
					//since coord values are duplicated in the structure
					//we don't have to check i,j and j,i combinations - just i,j is sufficient
				 if( (atom[i].type==*(it+1)) && (atom[j].type==*(it+2)) )
				{
					//print out the bond length
					double rij[3],sij[3];
					sij[0] = atom[i].sx - atom[j].sx;
					sij[1] = atom[i].sy - atom[j].sy;
					sij[2] = atom[i].sz - atom[j].sz;
					sij[0] = sij[0]-(int)(sij[0]*2);
					sij[1] = sij[1]-(int)(sij[1]*2);
					sij[2] = sij[2]-(int)(sij[2]*2);

					V3mulM3(sij,Hcry,rij);
					double rijsq_here = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
					//cout << atom[i].type <<"\t"<<atom[j].type<< "\t"<<sqrt(rijsq_here)<<" bond length\n";

					bonds.push_back(*it);
					bonds.push_back(i);
					bonds.push_back(j);
				}
			}
			// for angles
			for(int kk=jj+1; kk< atom[i].coord; kk++)
			{

				int k = atom[i].coord_id[kk];
				//if(i==760)	cout << i << " "<< j <<" "<< k << " inside angles\n";
				for(auto it = element_angles.begin(); it!=element_angles.end(); it=it+4)
				{
					if( (atom[i].type == *(it+2)) &&
							(
							( (atom[j].type == *(it+1)) && (atom[k].type == *(it+3)) ) ||
							( (atom[k].type == *(it+1)) && (atom[j].type == *(it+3)) )
							)
						)
					{

						//printout angles
						double rij1[3],rij2[3],sij[3];
						sij[0] = atom[i].sx - atom[j].sx;
						sij[1] = atom[i].sy - atom[j].sy;
						sij[2] = atom[i].sz - atom[j].sz;
						sij[0] = sij[0]-(int)(sij[0]*2);
						sij[1] = sij[1]-(int)(sij[1]*2);
						sij[2] = sij[2]-(int)(sij[2]*2);

						V3mulM3(sij,Hcry,rij1);
						//double rij1_sq = rij1[0]*rij1[0]+rij1[1]*rij1[1]+rij1[2]*rij1[2];

						sij[0] = atom[i].sx - atom[k].sx;
						sij[1] = atom[i].sy - atom[k].sy;
						sij[2] = atom[i].sz - atom[k].sz;
						sij[0] = sij[0]-(int)(sij[0]*2);
						sij[1] = sij[1]-(int)(sij[1]*2);
						sij[2] = sij[2]-(int)(sij[2]*2);

						V3mulM3(sij,Hcry,rij2);
						//double rij2_sq = rij2[0]*rij2[0]+rij2[1]*rij2[1]+rij2[2]*rij2[2];

						//cout << atom[j].type << " "<< atom[i].type << " "<< atom[k].type <<" "<< acos(V3dot(rij1, rij2))*180/PI<< " Angles\n";

						angles.push_back(*it);
						angles.push_back(j);
						angles.push_back(i);
						angles.push_back(k);
					}

				}
			}

		}
	}

	cout << "number of bonds are: "<< bonds.size()*1.0/3<<"\n";
	cout << "number of angles are: "<< angles.size()*1.0/4<<"\n";

}
