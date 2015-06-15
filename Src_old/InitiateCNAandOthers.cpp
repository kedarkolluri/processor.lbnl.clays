/*
 * InitiateCNAandOthers.cpp
 *
 *  Created on: Jul 20, 2009
 *      Author: kedar
 */

#include <InitiateCNAandOthers.h>

 double CNA_RATIO_1 = 0.899031812;
 double CNA_RATIO_2 = 1.20710676;
 double CNA_RATIO_12 = 1.0233219;

int prepare_nbrlist(double H1[3][3], int max_number_atoms)
{
	// Divide the entire cell into Boxes assuming H matrix exists
	mx = (int) (sqrt(H1[0][0]*H1[0][0]+H1[0][1]*H1[0][1]+H1[0][2]*H1[0][2])/rlst);
	my = (int) (sqrt(H1[1][0]*H1[1][0]+H1[1][1]*H1[1][1]+H1[1][2]*H1[1][2])/rlst);
	mz = (int) (sqrt(H1[2][0]*H1[2][0]+H1[2][1]*H1[2][1]+H1[2][2]*H1[2][2])/rlst);

	if(mx==0) mx = 1;  if(my==0) my = 1;  if(mz==0) mz = 1;

	ncell = mx*my*mz;

	int maxnbr = n*max_number_atoms; // assuming 100 atoms maximum per each atom would contribute to the energy of the system

	cout << mx<<"\t"<<my<<"\t"<<mz<<"\t"<<ncell<<" cell linked list parameters\n";



//	free(nbr_ptr);free(nbr_ptr1);free(nbr_lst);free(map);free(head);free(list);

	if(nbr_ptr!=NULL) { free(nbr_ptr); nbr_ptr = NULL;}
	if(nbr_ptr1!=NULL) { free(nbr_ptr1); nbr_ptr1=NULL;}
	if(nbr_lst!=NULL) { free(nbr_lst); nbr_lst = NULL;}
	if(map!=NULL) { free(map); map = NULL;}
	if(head!= NULL) { free(head); head = NULL;	}
	if(list != NULL) { free(list) ; list = NULL;}

	nbr_ptr = (int *) malloc((n+3)*sizeof(int));
	nbr_ptr1 = (int *) malloc((n+3)*sizeof(int));
	nbr_lst = (int *) malloc((maxnbr+3)*sizeof(int));
	map = (int *) malloc((ncell*13+3)*sizeof(int));
	head = (int *) malloc((ncell+3)*sizeof(int));
	list = (int *) malloc((n+3)*sizeof(int));
	mapcells();
	make_nbr_lst(H1);

}

void set_neighbordistances(char *filename)
{
	ifstream cutoff_file;
	string tmp_line;
	string *ptr_tmp_line;
	int num_entries;
	cout << filename<<"\n";
	cutoff_file.open(filename);
	cout << "file opened\n";
	if(cutoff_file.good())
	{
		cout << "file is good\n";
		ptr_tmp_line = get_next_splits(cutoff_file, num_entries);
		if(ptr_tmp_line[0] == "types")
		{
			cout << "here in n_types\n";
			n_types = atoi(ptr_tmp_line[1].c_str());
			cout << "determined n types as\t"<< n_types<<"\n";
		}else
		{
			cout << "first line in the cutoff line should be of the form \' types num_of_types_integer \' but it is not..\n";
			cout << "..exiting\n";
			exit(1);
		}
		rcoordsq = (double**) malloc((n_types)*sizeof(double *));
		latt_cutoff = (double**) malloc((n_types)*sizeof(double *));
		element_names = (string*) malloc(n_types*sizeof(string));

		for (int i=0;i<n_types;i++)
		{
			rcoordsq[i] = (double*) malloc((n_types)*sizeof(double));
			latt_cutoff[i] = (double*) malloc((2)*sizeof(double));

		}
		cout << "finished initialized the stuf\n";
		delete [] ptr_tmp_line;
		while(!cutoff_file.eof())
		{
			ptr_tmp_line = get_next_splits(cutoff_file, num_entries);
			cout << "got the values\n";
			if(ptr_tmp_line[0]=="lattice_param_block")
			{
				cout << "inside lattice_parameters\n";
				if(num_entries != n_types*3+1)
				{
					cout << "not correct number of entries .. exiting\n"; exit(1);
				}else
				{
					for (int i =1;i <num_entries;i=i+3)
					{
						cout << "before setting the value "<< atoi(ptr_tmp_line[i].c_str())<< "\n";
						latt_cutoff[atoi(ptr_tmp_line[i].c_str())-1][0] = atof(ptr_tmp_line[i+1].c_str());
						latt_cutoff[atoi(ptr_tmp_line[i].c_str())-1][1] = atof(ptr_tmp_line[i+2].c_str());
						cout << "lattice parameters "<< atoi(ptr_tmp_line[i].c_str())-1<< " "<<latt_cutoff[atoi(ptr_tmp_line[i].c_str())-1][0]<<"\t";
						cout << latt_cutoff[atoi(ptr_tmp_line[i].c_str())-1][1]<<"\n";
					}
				}
			}
			if(ptr_tmp_line[0]=="cutoffs")
			{
				if(num_entries !=3*n_types*(n_types+1)/2+1)
				{
					cout << " incorrect number of entries 1 1 cutofff 1 2 cutoff 1 3 cutoff etc... is the right format..exiting\n"; exit(1);
				}else
				{
					for (int i =1;i <num_entries;i=i+3)
					{
						rcoordsq[atoi(ptr_tmp_line[i].c_str())-1][atoi(ptr_tmp_line[i+1].c_str())-1] = atof(ptr_tmp_line[i+2].c_str())*atof(ptr_tmp_line[i+2].c_str())*latt_cutoff[atoi(ptr_tmp_line[i].c_str())-1][0]*latt_cutoff[atoi(ptr_tmp_line[i+1].c_str())-1][0];
						if((atoi(ptr_tmp_line[i].c_str())) -1!= (atoi(ptr_tmp_line[i+1].c_str())-1))
						rcoordsq[atoi(ptr_tmp_line[i+1].c_str())-1][atoi(ptr_tmp_line[i].c_str())-1] = atof(ptr_tmp_line[i+2].c_str())*atof(ptr_tmp_line[i+2].c_str())*latt_cutoff[atoi(ptr_tmp_line[i].c_str())-1][0]*latt_cutoff[atoi(ptr_tmp_line[i+1].c_str())-1][0];

						cout << atoi(ptr_tmp_line[i].c_str())-1<< " "<<atoi(ptr_tmp_line[i+1].c_str())-1<<" "<<rcoordsq[atoi(ptr_tmp_line[i].c_str())-1][atoi(ptr_tmp_line[i+1].c_str())-1]<<"\t";
						cout << rcoordsq[atoi(ptr_tmp_line[i+1].c_str())-1][atoi(ptr_tmp_line[i].c_str())-1] << "\n";
					}
				}
			}

			if(ptr_tmp_line[0]=="element_names")
			{
				if(num_entries !=2*n_types+1)
				{
					cout << " incorrect number and format of element names.. exiting\n"; exit(1);
				}else
				{
					for(int i =1;i<num_entries; i=i+2)
					{
					//	cout << atoi(ptr_tmp_line[i].c_str())-1 << "\t"<<ptr_tmp_line[i+1] <<"as \t"<<element_names[atoi(ptr_tmp_line[i].c_str())-1]<<" df\n";
					//	element_names[atoi(ptr_tmp_line[i].c_str())-1].assign("ting");
					}
				}
			}
			delete [] ptr_tmp_line;
		}



	}else
	{
		cout << "cutoffs file is corrupted..terminating\n";
		exit(1);
	}
	cutoff_file.close();

}

void compute_CNA_and_others( atomic_dat *atom_now, int n_now, double H_now[3][3])
{

//rcoordsq[0][0] = a_Cu*CNA_RATIO_1*a_Cu*CNA_RATIO_1; rcoordsq[1][1] =a_Nb*CNA_RATIO_2*a_Nb*CNA_RATIO_2; rcoordsq[0][1] =a_Cu*a_Nb*CNA_RATIO_12*CNA_RATIO_12; rcoordsq[1][0] = rcoordsq[0][1];

cout << rcoordsq[0][0] <<"\t"<< rcoordsq[1][1] <<"\t"<< rcoordsq[0][1] <<"\t"<<rcoordsq[1][0] <<"\n";

cout << "INITIATING STUFF                                DKLFHS;LSKDHFLSDFHSDLF\n";
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			cout << sqrt(rcoordsq[i][j])<<"\t";
		}
	}
	cout << "\n";
	prepare_nbrlist(H_now,100);
	coord_number(H_now);


	compute_CNA(atom_now,n_now,H_now);
	compute_acklandnotation(atom_now,n_now,H_now);

	compute_rings(atom_now,n_now,H_now);
	/*
		rcoordsq[0][0] = 3.085365*3.085365; rcoordsq[1][1] =3.08*3.08; rcoordsq[0][1] =3.534891*3.534891; rcoordsq[1][0] = 3.534891*3.534891; //	3.984418*3.984418//  3.534891*3.534891

	 prepare_nbrlist(H_now,100);
	 coord_number(H_now);
	 */


}

void compute_CNA_and_others( atomic_dat *atom_now, int n_now, double H_now[3][3], bool extra)
{


	cout << rcoordsq[0][0] <<"\t"<< rcoordsq[1][1] <<"\t"<< rcoordsq[0][1] <<"\t"<<rcoordsq[1][0] <<"\n";
	cout << "INITIATING STUFF                                DKLFHS;LSKDHFLSDFHSDLF\n";

	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			cout << sqrt(rcoordsq[i][j])<<"\t";
		}
	}
	cout << "\n";
	prepare_nbrlist(H_now,100);
	coord_number(H_now);

	if(extra)
	{
	compute_CNA(atom_now,n_now,H_now);
	compute_acklandnotation(atom_now,n_now,H_now);

	compute_rings(atom_now,n_now,H_now);
	}
	/*
		rcoordsq[0][0] = 3.085365*3.085365; rcoordsq[1][1] =3.08*3.08; rcoordsq[0][1] =3.534891*3.534891; rcoordsq[1][0] = 3.534891*3.534891; //	3.984418*3.984418//  3.534891*3.534891

	 prepare_nbrlist(H_now,100);
	 coord_number(H_now);
	 */


}
