#include "main.h"
bool MAIN_SAVE_LAMMPS = false;
bool MAIN_SAVE_FULL = false;
bool SEQ_PROCESS = true;
bool WRITE_XYZ_FORMAT = false;
bool ATOMIC_RESET = false;
bool CONVERT_VESTA = false;
bool MAKE_ILLITE = false;
bool KEEP_GHOSTS = false;
double del_z = 0.0;
int random_seed = -1;
bool LAMMPS_CHARGE = false;
bool LAMMPS_MOLECULE = false;




int main (int argc, char *argv[])
{
	inputfilename.assign("dat_lammps.");

	bool cutoff_file = false;

	for(int tt = 1;tt < argc;tt++)
	{
		if (strcmp(argv[tt], "convert_VESTA")==0)
		{
			// we will assume for now that we can take in only 1 file
			if(tt+2< argc)
			{
				//inputfilename is the stub
				if(strcmp(argv[tt+1],"filename")==0)
				{
					SEQ_PROCESS = false;
					UTILITIES_ZIP = false;
					CONVERT_VESTA = true;

					inputfilename = argv[tt+2];
          std::cout << "filename is "<< inputfilename<<"\n";
					tt = tt+2;
          std::cout << argc<<" "<<tt<<"\n";

				}else
				{
					cout << "to convert, please say `filename` and then follow with name of the file\n";

					exit(1);
				}
        if(tt+1< argc)
        {
          if(strcmp(argv[tt+1],"make_illite")==0)
  				{
  					MAKE_ILLITE = true;
  					tt=tt+1;
  					if(strcmp(argv[tt+1], "keep_ghosts")==0)
  					{
  						KEEP_GHOSTS = true;
  						tt=tt+1;
  					}
  					if(strcmp(argv[tt+1], "increase_z")==0)
  					{
  						del_z = atof(argv[tt+2]);
  						cout << "********************del z value is "<<del_z<<"\n";
  						tt=tt+2;
  					}
  				}
        }

			}else
			{
				cout << "not enough information to convert, exiting\n";
				exit(1);
			}
		}else if(strcmp(argv[tt], "rand_seed")==0)
		{

			random_seed = atoi(argv[tt+1]);

		} else if(strcmp(argv[tt], "stub")==0)
		{

			inputfilename.assign(argv[tt]);

		} else if(strcmp(argv[tt],"start")==0)
		{
			filenumber_start = atoi(argv[tt+1]);
			filenumber_end = filenumber_start+1;
			filenumber_interval = 1;
		}else if(strcmp(argv[tt],"end")==0)
		{
			filenumber_end = atoi(argv[tt+1]);
		}else if(strcmp(argv[tt], "interval")==0)
		{
			filenumber_interval = atoi(argv[tt+1]);

		}else	if(strcmp(argv[tt],"zip_unzip")==0)
		{
			if(strcmp(argv[tt+1],"1")==0)
			UTILITIES_ZIP = false;
		}else if(strcmp(argv[tt],"SAVE_LAMMPS")==0)
		{
				MAIN_SAVE_LAMMPS = true;
				if(tt+1<argc)
				{
					if(strcmp(argv[tt+1],"CHARGE")==0)
					{
						LAMMPS_CHARGE = true;
						tt++;
						if(tt+1<argc)
						{
							if(strcmp(argv[tt+1],"MOLECULE")==0)
								LAMMPS_MOLECULE = true;
								tt++;
						}
					}else if(strcmp(argv[tt+1],"MOLECULE")==0)
					{
						LAMMPS_MOLECULE = true;
						tt++;
						if(tt+1<argc)
						{
							if(strcmp(argv[tt+1],"CHARGE")==0)
								LAMMPS_CHARGE = true;
								tt++;
							}
						}
				}
		}else if(strcmp(argv[tt],"SAVE_FULL")==0)
		{
			MAIN_SAVE_FULL = true;
		}else if(strcmp(argv[tt], "XYX")==0)
		{
			WRITE_XYZ_FORMAT = true;
		}else if(strcmp(argv[tt],"SEQ")==0)
		{
			if((strcmp(argv[tt+1],"NO")==0)||(strcmp(argv[tt+1],"No")==0)||(strcmp(argv[tt+1],"no")==0))
			{
				SEQ_PROCESS = false;
			}
		}else if(strcmp(argv[tt],"ATOMIC_RESET")==0)
		{
			ATOMIC_RESET = true;
			cout << "entered in ATOMIC RESET\n";
		}else if(strcmp(argv[tt],"CUTOFF_FILE")==0)
		{
			if(tt+1<argc)
			{
				//set_neighbordistances(argv[tt+1]);
				tt++;
			}else
			{
				cout << "input error; file containing cutoffs is not there... exiting";
				exit(1);
			}
		}else
		{
			cout << "the parameter "<< argv[tt]<< " does not exist\n";
		}
	}


	if(CONVERT_VESTA)
	{
    std::cout << "entered into vista stuff\n";
    simcell mdcell = read_xyz_VESTA(inputfilename.c_str());
    //std::cout << mdcell.n <<"\n";
    std::cout << "exited 2\n";
  }
  std::cout << "exited bro\n";
/*
		//read_xyz_VESTA(inputfilename.c_str());

		set_types_to_atom_from_element_name();
		if(MAKE_ILLITE)
		{
			// THIS CODE USED TO CONVERT VESTA XYZ to LAMMPS
			//***this commenting is done to use this block to process xyz files***
			prepare_nbrlist(Hcry, 100);
			coord_number(Hcry);

			int total_K = 0;
			for(int i=0;i<n;i++)
			if(atom[i].type==7) total_K++;
			cout << "before removing duplicates, number of atoms are "<< n<< "\n";
			cout << "before removing duplicates, number of K atoms are "<< total_K<<"\n";


			delete_duplicates();
			cout << "after removing duplicates, the number of atoms are "<<n<<"\n";

			save_lammps(10,Hcry);
			save_xyz_VESTA(10,Hcry);
			//specifically for K, which is 7
			// find the number of K ions
			total_K = 0;
			for(int i=0;i<n;i++)
			if(atom[i].type==7) total_K++;
			cout << "after removing duplicates, number of K atoms are "<< total_K<<"\n";

			//change few of Ge atoms to Al
			//but min distance between 2 Al atoms is more than 4Angstorms
			int changed_ge_atoms=0;

			//change 0.35 of type 10 atoms to type 4
			//minimum distance between changed atoms is 5.0 Angstorms
			//changed_ge_atoms = change_atom_types_randomly(10, 4, 0.35, 5.0);

			changed_ge_atoms = change_atom_types_randomly(10, 4, 0.35, 5.0);

			cout << total_K-changed_ge_atoms<<" Cations to remove\n";
			if((total_K-changed_ge_atoms) < 0) {
				cout << "more atoms to be removed than there are\n";
				cout << "something is wrong\n";
				cout <<"please check input data file...\n";
				exit(1);
			}

			// delete `total_K-changed_ge_atoms` atoms (>1 means number) of type 7
			// the minimum distance btween removed should be 5.0
			// delete_atoms_randomly(7, total_K-changed_ge_atoms, 5.0);
			if( !KEEP_GHOSTS )
			{
				delete_atoms_randomly(7, total_K-changed_ge_atoms, 0.0);

			}else
			{
				change_atom_types_randomly(7, 11, total_K-changed_ge_atoms);
			}


			// change the remaining Ge atoms to Si
			// Change 1.0 fraction of type 10 atoms to type 2 atoms
			// change_atom_types_randomly(10, 2, 1.0);
			change_atom_types_randomly(10, 2, 0.999999999);

			prepare_nbrlist(Hcry, 100);

			coord_number(Hcry);
			//calling above method is necessary for below method to work
			//Find atoms of type 5 bonded to atoms of type 4 and change their type to 9
			change_bridging_oxygen(4, 5, 9) ;

		}

		prepare_nbrlist(Hcry, 100);
		coord_number(Hcry);
		int number_of_K_atoms = 0;

		int counter_K = 0;
		for(int tag=0;tag<n;tag++)
		{
			if(atom[tag].type==7)
			{
				counter_K++;
			}
		}
		cout << counter_K <<" the number of K+ atoms currently available\n";

		//do some `filling` of bonds -- add bonds
	fill_bonds_etc();
	if(MAIN_SAVE_LAMMPS) save_lammps(20,Hcry);
  save_xyz_VESTA(20,Hcry);

	}
  */
	if(SEQ_PROCESS)
	{
		std::cout << "Enetered into SEQ Process\n\n";
	//	seq_process_lammps_new();
	}
	return (0);
}
