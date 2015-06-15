/*
 * ChangeAtoms.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#include <ChangeAtoms.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>

int random_seed = -1;

std::set<int> find_duplicate_atoms(double cutoff_dist = 0.09)
{
	std::set<int> duplicate_atoms;
	std::set<int> keep_atoms;
	int jbeg,jend,jnab;
	int *coord;
	//bool start = false;

	double sxi,syi,szi,sij[3],rij[3],rijsq;

	//int total_close = 1000;
	//int close_atoms[total_close];
	int counter_close=0;

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

			V3mulM3(sij,Hcry,rij);
			rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];

			if(rijsq<=cutoff_dist)
			{
				//if(!start) { duplicate_atoms.insert(j); start =false;}

				//just checking to see if two different types of atoms are overlapping
				//should not happen in this case
				if(strcmp(atom[i].elem, atom[j].elem)!=0)
				{
					cout << i << " "<< j<< " "<< atom[i].elem;
					cout << " "<< atom[j].elem << " ******* WHAT THE HOLY hell!\n";
				}

				//check if j is already in the list of keep_atoms and if not,
				//mark it for deletion and mark i for keeping
				auto find_keep_j = keep_atoms.find(j);
				auto find_atom = find_keep_j;
				auto find_dup_j = duplicate_atoms.find(j);
				auto find_keep_i = keep_atoms.find(i);
				auto find_dup_i = duplicate_atoms.find(i);

				if(find_keep_j == keep_atoms.end())
				{
					if(find_dup_j == duplicate_atoms.end())
					{
						if(find_dup_i == duplicate_atoms.end())
						{
							if(find_keep_i == keep_atoms.end())
							{
								keep_atoms.insert(j);
								duplicate_atoms.insert(i);
							}else
							{
								duplicate_atoms.insert(j);
							}
						}
						{
							keep_atoms.insert(j);
						}
					}else{
						if(find_dup_i == duplicate_atoms.end())
						{
							if(find_keep_i == keep_atoms.end())
							{
								keep_atoms.insert(i);
							}else
							{

							}
						}else
						{
							//keep_atoms.insert(i);
							//duplicate_atoms.erase(i);
						}
					}
				}else
				{
					if(find_dup_i == duplicate_atoms.end())
					{
						if(find_keep_i == keep_atoms.end())
						{
							duplicate_atoms.insert(i);
						}else
						{
							duplicate_atoms.insert(i);
							keep_atoms.erase(i);
						}
					}
				}
			}
		}
	}
	cout << "keeping atoms\n";
	for(auto it = keep_atoms.begin(); it != keep_atoms.end(); it++)
	{
		cout << (*it) << " "<<atom[*it].elem<<"\n";
	}
	cout << "removing atoms\n";

	for(auto it = duplicate_atoms.begin(); it != duplicate_atoms.end(); it++)
	{
		cout << (*it) << " "<<atom[*it].elem<<"\n";
	}

	return duplicate_atoms;
}

void delete_atoms(int *arr, int len)
{

	int n_new = n-len;
	atomic_dat *atom_new;
	atom_new = (struct atomic_dat *) malloc((n_new+3)*sizeof(struct atomic_dat));
	int j = 0;
	for(int i=0;i<n;i++)
	{
		if(!check_repeat(i,arr,len))
		{
			copy_atom_indv(atom,atom_new, i,j);
			j++;
		}
	}
	n = n_new;
	free(atom);
	atom = atom_new;

}

void delete_atom_single(int tag)
{
	int n_new = n-1;
	atomic_dat *atom_new;
	atom_new = (struct atomic_dat *) malloc((n_new+3)*sizeof(struct atomic_dat));
	int j = 0;
	for(int i=0;i<n;i++)
	{
		if(i!=tag)
		{
			copy_atom_indv(atom,atom_new, i,j);
			j++;
		}
	}
	n = n_new;
	free(atom);
	atom = atom_new;
}

void insert_atoms(double (*arr)[4], int len)
{

	int n_new = n+len;
	atomic_dat *atom_new;
	atom_new = (struct atomic_dat *) malloc((n_new+3)*sizeof(struct atomic_dat));
	int j = 0;
	for(int i=0;i<n;i++)
	{

		copy_atom_indv(atom,atom_new, i,j);
		j++;
	}

	for(int i =n;i<n_new;i++)
	{
		atom_new[i].sx = arr[(i-n)][0];
		atom_new[i].sy = arr[(i-n)][1];
		atom_new[i].sz = arr[(i-n)][2];
		atom_new[i].type = arr[(i-n)][3];
	}
	n = n_new;
	free(atom);
	atom = atom_new;

}

void insert_atoms_specific(double (*arr)[5], int len)
{

	int n_new = n+len;
	atomic_dat *atom_new;
	atom_new = (struct atomic_dat *) malloc((n_new+3)*sizeof(struct atomic_dat));
int j = 0;
	for(int i=0;i<n;i++)
	{

		for(int k=0;k<len;k++)
		{
			if(arr[k][4]==i)
			{
				atom_new[j].sx = arr[k][0];
				atom_new[j].sy = arr[k][1];
				atom_new[j].sz = arr[k][2];
				atom_new[j].type = arr[k][3];
				j++;
			}
		}
		copy_atom_indv(atom,atom_new, i,j);
		j++;
	}
	n = n_new;
	free(atom);
	atom = atom_new;
}




void delete_duplicates()
{
	auto get_values = find_duplicate_atoms();
	std::vector<int> dup_atoms(get_values.begin(), get_values.end());

	delete_atoms(dup_atoms.data(), dup_atoms.size());
}

double sq_atomic_dist(int i, int j)
{
	double sij[3], rij[3];
	sij[0] = atom[i].sx - atom[j].sx;
	sij[1] = atom[i].sy - atom[j].sy;
	sij[2] = atom[i].sz - atom[j].sz;
	sij[0] = sij[0]-(int)(sij[0]*2);
	sij[1] = sij[1]-(int)(sij[1]*2);
	sij[2] = sij[2]-(int)(sij[2]*2);
	V3mulM3(sij,Hcry,rij);
	//cout <<"\n"<< i << " "<< j<< " "<<atom[i].type << " "<< atom[j].type<<"\n";
	//cout << i << " " << j<<" "<< rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]<<"\n";
	return rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
}

std::vector<int> find_atoms_to_change_randomly(int type_now, double fraction_or_natoms, double min_dist)
{
	cout << "finding atoms randomly that fit the criteria***\n";
	std::vector<int> type_info;
	std::vector<int> type_info2;
	std::vector< vector<int> > coords_change_type;

	//store all the candidate atoms into `type_info`
	for (int i=0;i<n;i++)
	{
		if(atom[i].type==type_now)
		{
			type_info.push_back(i);
		}
	}
	int rand_num;
	// determine the number of atoms
	int atoms_to_be_replaced=0;
	if (fraction_or_natoms < 1.0)
		atoms_to_be_replaced = fraction_or_natoms*type_info.size() + 0.5;
	else
			atoms_to_be_replaced = (int) fraction_or_natoms;

	cout << "number of atoms that fit the criteria are\t" << atoms_to_be_replaced<<"\n";
	int atoms_will_be_replaced = 0;
	//initialize a random seed uses time so we can use it many times without a problem
	if(random_seed == -1)
	{
		std::cout << "no random seed given; taking a pseudo-random value\n";
		std::srand(time(NULL));
	}else
	{
		std::cout << "random seed is " << random_seed<<"\n";
		std::srand(random_seed);
	}

//

	// do a while loop to find the number of
	bool reached_end = false;
	if(atoms_to_be_replaced == 0) reached_end = true;
	cout << "original size of possible atoms "<< type_info.size()<<"\n";
	while (!reached_end)
	{
		// create a random number between 0 number of eligible replaceble atoms
		rand_num = rand() % type_info.size();
		// find the atom position in type_info for that random number
		int tag = type_info[rand_num];
		// cout << type_info.size()<< " " << rand_num <<" "<<tag<<" are the size and random integer generated ";

		//check if that atom is close to any of the atoms already marked for change
		auto it = type_info2.begin();
		bool change_atom = true;
		while (it != type_info2.end())
		{
			if(min_dist*min_dist > sq_atomic_dist(tag, *it))
			{
				change_atom = false;
				//cout << tag << " " << *it << " are close so ";
				it = type_info2.end();

			}else
			{
				it++;
			}
		}
		//if the atom is not close to any other
		//add it to vector of atoms that can be changed
		if(change_atom)
		{
			type_info2.push_back(tag);
			atoms_will_be_replaced++;
			//cout << " added\n";
		}
		type_info.erase(type_info.begin()+rand_num);

		if((atoms_will_be_replaced >= atoms_to_be_replaced)||(type_info.size()==0))
			reached_end = true;

	}
	//return std::tuple<int, vector<int> > change_info(atoms_will_be_replaced, type_info2);
	cout <<" criteria is "<< atoms_will_be_replaced << " "<< atoms_to_be_replaced << " " <<type_info.size()<<"\n";
	return type_info2;

}



void delete_atoms_randomly(int type, double fraction_or_natoms, double min_dist)
{

	std::vector<int> selected_atoms = find_atoms_to_change_randomly(type, fraction_or_natoms, min_dist);
	delete_atoms(selected_atoms.data(),selected_atoms.size());
	int number_remaining  = 0;
	for(int i=0;i<n;i++)
	{
		if(atom[i].type==type)
		{
			number_remaining++;
		}
	}
	cout << "number of atoms of type "<< type << " after delete are "<< number_remaining<<"\n";

}


int change_atom_types_randomly(int type_now, int type_needed, double fraction_or_natoms, double min_dist)
{
	cout << " inside change atom type randomly function\n";
	std::vector<int> selected_atoms = find_atoms_to_change_randomly(type_now, fraction_or_natoms, min_dist);
	cout << "size of the selected_atoms vector is\t"<< selected_atoms.size()<<"\n";
	for (auto it = selected_atoms.begin(); it != selected_atoms.end(); it++)
	{
		atom[*it].type = type_needed;
		strncpy(atom[*it].elem, element_names[type_needed-1].c_str(), sizeof(atom[*it].elem));
		atom[*it].ma = element_weights[type_needed-1];
		atom[*it].charge = element_charges[type_needed-1];
	}
	return selected_atoms.size();

}

/*
int change_atom_types_by_sub_Al(int type_now, int type_needed, double fraction_or_natoms, double min_dist)
{
	cout << " inside change atom type by a set of rules\n";
	cout << "these rules are base on placing cations next to substitued Al"

	//	Rules are as follows
	//	First remove (or substitute those ions that have no Al neighbors
	//	Second remove (or substitute) ions that have 1 Al neighbor
	//	go on until no more ions need to be removed


	std::vector<int> sub0;
	std::vector<int> sub1;
	std::vector<int> sub2;
	std::vector<int> sub3;

	for(int i =0; i <n, i++)
	{
		if(atom[i].type == type_now)
		{
			``
		}
	}

}
*/

int change_bridging_oxygen(int root_type, int type_of_neighbor, int type_needed)
{
	int number_changed = 0;
	int number_of_al =0;
	int number_of_bridge_ox = 0;
	for (int i=0;i<n;i++)
	{
		if(atom[i].type==root_type)
		{
			number_of_al++;
			for(int k=0;k<atom[i].coord; k++)
			{
				int tag = atom[i].coord_id[k];
				if (atom[tag].type==type_of_neighbor)
				{
					number_of_bridge_ox++;
					atom[tag].type = type_needed;
					strncpy(atom[tag].elem, element_names[type_needed-1].c_str(), sizeof(atom[tag].elem));
					atom[tag].ma = element_weights[type_needed-1];
					atom[tag].charge = element_charges[type_needed-1];
					number_changed++;
				}
			}
		}
	}
	cout << "changed is "<<number_changed<<"\n";
	cout << "number Al is "<< number_of_al<<"\n";
	cout << "number of bridge Ox is "<< number_of_bridge_ox<<"\n";
}
