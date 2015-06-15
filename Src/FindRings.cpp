/*
 * FindRings.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#include <FindRings.h>
bool RING_DATA_COM = true;
int RING_DATA_RING_SIZE = 5;
int follow = -100;
int root_n;
bool is_distance_acceptable(atomic_dat *atom_now,int i,int ring_seq[MAX_LOOP_SIZE],int len, double H_now[3][3])
{
	double inplane_cutoff_sq = 9.0*9.0;
	bool is_neighbor_of_first_atom=false;
	/*
	 int common_nbr = -1;
	 int cntr = 0;
	 while((common_nbr<0)&&(cntr<(MAX_COORD-1)))
	 {
		 if(atom[i].coord_id[cntr]>-1)
		 {
			 int f = ring_seq[len-1];
			 if(check_repeat(atom[i].coord_id[cntr],atom[f].coord_id,MAX_COORD)) common_nbr = atom[i].coord_id[cntr];
			 cntr++;
		 }
	 }
	 */

	for(int j=0;j<len;j++)
	{
		double s[3],r[3];
		s[0]=s[1]=s[2]=r[0]=r[1]=r[2]=0.0;

		int k = ring_seq[j];
		s[0] = atom_now[i].sx-atom_now[k].sx;
		s[1] = atom_now[i].sy-atom_now[k].sy;
		s[2] = atom_now[i].sz-atom_now[k].sz;

		s[0] = s[0]-(int)(s[0]*2);
		s[1] = s[1]-(int)(s[1]*2);
		s[2] = s[2]-(int)(s[2]*2);

		V3mulM3(s,H_now,r);
		//		cout << "D 2nd\t"<<i<<"\t"<<k<<"\t"<<sqrt((r[0]*r[0]+r[1]*r[1]+r[2]*r[2]))<<"\n";
		double dis = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];

		if(dis>inplane_cutoff_sq)
			{
				if(i==10770)
						{
							cout << k<<"\t"<<sqrt(dis)<< "\t"<< sqrt(inplane_cutoff_sq)<< "\t distance checks\n";
						}
				return(false);
			}

		if((len>1)&&(j>0)&&(j!=(len-1)))
		{

			//			cout <<"entered check multi bond \n";
			if(dis<(rcoordsq[atom[i].type-1][atom[k].type-1])) return(false);
			//			cout <<"passed\n";
		}

		if((j==0)&&(dis<(rcoordsq[atom[i].type-1][atom[k].type-1]))) is_neighbor_of_first_atom = true;
		//Check for first and last atoms in the array ring_seq are not neighbors with each other and with the new atom i

	}
	if(len>1)
	{
		double s[3],r[3];
		s[0]=s[1]=s[2]=r[0]=r[1]=r[2]=0.0;

		int k = ring_seq[(len-1)];
		int l = ring_seq[0];
		s[0] = atom_now[l].sx-atom_now[k].sx;
		s[1] = atom_now[l].sy-atom_now[k].sy;
		s[2] = atom_now[l].sz-atom_now[k].sz;

		s[0] = s[0]-(int)(s[0]*2);
		s[1] = s[1]-(int)(s[1]*2);
		s[2] = s[2]-(int)(s[2]*2);

		V3mulM3(s,H_now,r);
		double dis = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
		if((dis<(rcoordsq[atom[l].type-1][atom[k].type-1]))&&(is_neighbor_of_first_atom)) return(false);
	}


	if(len>1)
	{
		for(int j=0;j<len;j++)
		{
			double s[3],r[3];
			s[0]=s[1]=s[2]=r[0]=r[1]=r[2]=0.0;

			int k = ring_seq[j];

			int tot_neighs=0;
			for(int p = j+1;p<=len;p++)
			{
				int q;
				if (p< len)
				{
					q = ring_seq[p];
				}else
				{
					q = i;
				}


				s[0] = atom_now[q].sx-atom_now[k].sx;
				s[1] = atom_now[q].sy-atom_now[k].sy;
				s[2] = atom_now[q].sz-atom_now[k].sz;

				s[0] = s[0]-(int)(s[0]*2);
				s[1] = s[1]-(int)(s[1]*2);
				s[2] = s[2]-(int)(s[2]*2);

				V3mulM3(s,H_now,r);
				double dis = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
				if(dis<(rcoordsq[atom[k].type-1][atom[q].type-1]))
				{
					tot_neighs++;
					if(tot_neighs>2) return(false);
				}else
				{
					/*
					 if((q==i)&&(j!=(len-1))&&(common_nbr>-1)){

						 if(!check_repeat(common_nbr,ring_seq,len)&&(common_nbr!=i))
						 {
							 if(check_repeat(common_nbr,atom_now[k].coord_id,MAX_COORD))
							 {

								 return(false);
							 }
						 }

					 }*/

				}
			}

		}

	}

	return(true);
}

int complete_loop(int i, int j,int *counter, int ring_seq[MAX_LOOP_SIZE], atomic_dat *atom_now, double H_now[3][3])
{
	//bool check_repeat(int j,int *arr,int len)

	if ((i == follow)||(j==follow))
	{

		root_n = follow;
	}

	int len = MAX_LOOP_SIZE;

	int is_complete;
	bool distance_acceptable = false;
	//cout <<"0st\t" <<i<<"\t"<<atom_now[i].coord_id[j]<<"\t"<<atom_now[i].neigh_bonds[j]<<"\n";

	int n_tag = atom_now[i].coord_id[j];

	if((n_tag==ring_seq[0])&&(*counter>2)) return(0);

	if(check_repeat(n_tag,ring_seq,*counter)) return(1);


	if ((root_n == follow)||(i == follow)||(atom_now[i].coord_id[j]==follow))
		cout << "1st\t "<<i<<" \t "<<atom_now[i].coord_id[j]<<" \t " <<atom_now[i].neigh_bonds[j] <<"\n";

	if((atom_now[i].neigh_bonds[j]>-1)&&(atom_now[i].neigh_bonds[j]<2))
	{
		if ((root_n == follow)||(i == follow)||(atom_now[i].coord_id[j]==follow)) cout << "2nd\t "<<i<<" \t "<<atom_now[i].coord_id[j]<<"\n";
		if(is_distance_acceptable(atom_now, atom_now[i].coord_id[j],ring_seq,*counter, H_now))
		{
			if ((root_n == follow)||(i == follow)||(atom_now[i].coord_id[j]==follow)) cout <<"entered  loop for\t "<<atom_now[i].coord_id[j]<<"\n";

			ring_seq[(*counter)]=n_tag;
			if ((root_n == follow)||(i == follow)||(atom_now[i].coord_id[j]==follow)) cout << *counter << " " <<n_tag << " inside before looping multiple\n";

			(*counter)++;
			if ((root_n == follow)||(i == follow)||(atom_now[i].coord_id[j]==follow)) if((*counter>9)) cout << i<<" \t "<<n_tag<<" WHHHHHHHHHHHYYYYYYYYYYYY\n";



			for(int k=0;k<MAX_COORD;k++)
			{
				if(atom[n_tag].coord_id[k]>-1)
				{
					is_complete = complete_loop(n_tag,k,counter,ring_seq,atom_now,H_now);
					if(is_complete==0)
					{
						k=MAX_COORD;
						return(0);
					}
				}

				if(k==(MAX_COORD-1))
				{
					ring_seq[(*counter)] = -1;
					(*counter)--;
				}

				//				cout << *counter <<"\t is the counter value\n";
				/*
				 if(is_complete==-1)
				 {
					 k=MAX_COORD;
					 return(-1);
				 }
				 */

			}
		}else
		{
			return(1);
		}
	}

	return(1);
}

void compute_rings(atomic_dat *atom_now, int n_now, double H_now[3][3])
{
	cout <<"came inside the compute_rings function\n";

	/*
	 for(int i=0;i<n;i++)
	 {
		 atom[i].CNA = -1;
	 }
	 */
	int ting_tong = 0;

	for(int i=0;i<n_now;i++)
	{

//		int if ((i == 2474)||(j==2474));
//		int coord_limit=12;
//		if(atom_now[i].type==2) coord_limit=14;

		if(i==follow) cout << "interested atom is here but is it in the interface?\t" << atom_now[i].interface<<"\n";

		if((atom_now[i].interface>0)&&(atom_now[i].CNA<=6))
		{
			if(i==follow) cout << "interested atom is in the interface \n";
			int loops=0;

			for(int j=0;j<MAX_COORD;j++)
			{
				int ring_seq[MAX_LOOP_SIZE];for(int r=0;r<MAX_LOOP_SIZE;r++) ring_seq[r]=-1;
				ring_seq[0]=i;
				int counter=1;

				if ((i == follow)||(atom_now[i].coord_id[j]==follow)) cout << "neighbors are \t"<<i<<"\t"<<atom_now[i].coord_id[j]<<"\n";

				if(atom_now[i].coord_id[j]>-1)
				{
					if ((i == follow)||(atom_now[i].coord_id[j]==follow))
						{
							cout << i<<"\t"<<atom_now[i].coord_id[j]<<"\t"<<atom_now[i].neigh_bonds[j]<<"\t"<<atom_now[i].CNA<<"\n";
							//root_n = follow;
						}
					int is_complete = complete_loop(i,j,&counter,ring_seq,atom_now, H_now);
					if(is_complete==0)
					{
						//j = MAX_COORD;
						if ((i == follow)||(atom_now[i].coord_id[j]==follow)) cout <<"found one loop\t total loop count and current loop length is\t";
						loops++;
						if ((i == follow)||(atom_now[i].coord_id[j]==follow)) cout <<loops<<"\t"<<counter<<"\n";
						if(counter>MAX_LOOP_VALUE) MAX_LOOP_VALUE = counter;

						bool is_new = false;
						if((counter>2)&&(counter<16))
						{
							if(counter==3)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING3) { is_new = true; atom[ring_seq[k]].CNA = UNKNOWN_RING3;}
							}
							} else if(counter==4)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING4) { is_new = true; atom[ring_seq[k]].CNA = UNKNOWN_RING4;}
							}
							} else if(counter==5)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING5) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING5;}
							}
							} else if(counter==6)
							{	for(int k=0;k<counter;k++)
							{
							//	cout << ring_seq[k] <<"\t"<<atom[ring_seq[k]].CNA<< "\t"<< UNKNOWN_RING6 << "\t" << UNKNOWN_RING5<<"\t"<<UNKNOWN_RING4<<" in counter 6\n";
								if(atom_now[ring_seq[k]].CNA < UNKNOWN_RING6) { is_new = true;atom_now[ring_seq[k]].CNA = UNKNOWN_RING6;}
							}
							} else if(counter==7)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING7) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING7;}
							}
							} else if(counter==8)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING8) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING8;}
							}
							} else if(counter==9)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING9) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING9;}
							}
							} else if(counter==10)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING10) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING10;}
							}
							} else if(counter==11)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING11) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING11;}
							}
							} else if(counter==12)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING12) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING12;}
							}
							} else if(counter==13)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING13) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING13;}
							}
							} else if(counter==14)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING14) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING14;}
							}
							} else if(counter==15)
							{	for(int k=0;k<counter;k++)
							{
								if(atom[ring_seq[k]].CNA < UNKNOWN_RING15) { is_new = true;atom[ring_seq[k]].CNA = UNKNOWN_RING15;}
							}
							}

							if((RING_DATA_COM)&&(counter>=RING_DATA_RING_SIZE)&&(is_new))
							{
								double x_c = atom_now[ring_seq[0]].sx;
								double y_c = atom_now[ring_seq[0]].sy;
								double z_c = atom_now[ring_seq[0]].sz;
								double cm[3];
								cm[0] = x_c*counter; cm[1] = y_c*counter; cm[2] = z_c*counter;
								for(int k=1;k<counter;k++)
								{
									double s_x = atom_now[ring_seq[k]].sx;
									double s_y = atom_now[ring_seq[k]].sy;
									double s_z = atom_now[ring_seq[k]].sz;
									double ss_x = s_x-x_c;
									double ss_y = s_y-y_c;
									double ss_z = s_z-z_c;

									if(ss_x >= 0.5) ss_x= 1-ss_x;
									if(ss_y >=0.5) ss_y= 1-ss_y;
									if(ss_z >=0.5) ss_z= 1-ss_z;

									if(ss_x<= -0.5) ss_x= 1+ss_x;
									if(ss_y<= -0.5) ss_y= 1+ss_y;
									if(ss_z<= -0.5) ss_z= 1+ss_z;

									cm[0] = cm[0]+ss_x;
									cm[1] = cm[1]+ss_y;
									cm[2] = cm[2]+ss_z;
								}
								double r[3];
								cm[0] = cm[0]/counter;
								cm[1] = cm[1]/counter;
								cm[2] = cm[2]/counter;
								V3mulM3(cm,H_now,r);
								cout << "RING COM\t"<< r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<< counter<<"\t"<<ting_tong<<"\t"<<i<<"\n";
								ting_tong++;
							}

						}


						//for(int k=0;k<counter;k++) {if(atom[ring_seq[k]].CNA<counter) atom[ring_seq[k]].CNA = counter;}
					}else
					{
						//cout <<"there is no loop yet for this set\n";

					}
					//cout << is_complete<<"$$\t"<< atom_now[i].CNA<<"\n";;
				}
			}
		}
	}

	cout << "MAX LOOP VALUE ISSSSSSSSSSS\t"<< MAX_LOOP_VALUE<<"\n";
}
