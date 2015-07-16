#include <simulationparams.h>


atom_types_info set_atom_types_info(char *filename)
{
	ifstream cutoff_file;
	string tmp_line;
	string *ptr_tmp_line;
	int num_entries;
  atom_types_info atominfo;
	cutoff_file.open(filename);
	if(cutoff_file.good())
	{
		ptr_tmp_line = get_next_splits(cutoff_file, num_entries);
		if(ptr_tmp_line[0] == "types")
		{
			cout << "here in n_types\n";
			atominfo.n_types = atoi(ptr_tmp_line[1].c_str());
			cout << "determined n types as\t"<< atominfo.n_types<<"\n";
		}else
		{
			cout << "first line in the cutoff line should be of the form \' types num_of_types_integer \' but it is not..\n";
			cout << "..exiting\n";
			exit(1);
		}

		cout << "slurped types values\n";

		delete [] ptr_tmp_line;
    int cutoff_entries = 0;
		while(!cutoff_file.eof())
		{
			ptr_tmp_line = get_next_splits(cutoff_file, num_entries);
			if(ptr_tmp_line[0]=="lattice_param_block")
			{
				cout << "inside lattice_parameters\n";
				if(num_entries != atominfo.n_types*3+1)
				{
					cout << "not correct number of entries .. exiting\n"; exit(1);
				}else
				{
					for (int i =1;i <num_entries;i=i+3)
					{
            int local_type = atoi(ptr_tmp_line[i].c_str()) ;
            double local_cutoff = atof(ptr_tmp_line[i+1].c_str());
            //not using this for now but we could use it later if we want
            int local_block = atof(ptr_tmp_line[i+2].c_str());
            atominfo.latt_cutoff[local_type] = local_cutoff;
					}
				}
			}

			if(ptr_tmp_line[0]=="cutoffs")
			{
        cutoff_entries+=num_entries-1;

				for (int i =1;i <num_entries;i=i+3)
				{
          int localtype1 = atoi(ptr_tmp_line[i].c_str());
          int localtype2 = atoi(ptr_tmp_line[i+1].c_str());
          double localratio = atof(ptr_tmp_line[i+2].c_str());

					double localcutoff = localratio * localratio * atominfo.latt_cutoff[localtype1] * atominfo.latt_cutoff[localtype2];

          atominfo.rcoordsq[std::make_tuple(localtype1, localtype2)] = localcutoff;
          if(localtype1 != localtype2)
          atominfo.rcoordsq[std::make_tuple(localtype2, localtype1)] = localcutoff;
        }
      }

			if(ptr_tmp_line[0]=="element_names")
			{
        cout << "number of entries in element_names line are "<< num_entries<<"\n";
				if(num_entries >2*atominfo.n_types+1)
				{
					cout << "number and format of element names is more than types given in the file ";
          cout << "this is not accpetable.. exiting\n"; exit(1);
				}else
				{
					for(int i =1;i<num_entries; i=i+2)
					{

					  atominfo.element_names[atoi(ptr_tmp_line[i].c_str())].assign(ptr_tmp_line[i+1]);
					}
				}
			}

      if(ptr_tmp_line[0]=="element_weights")
      {
        cout << "\n\nentered element weights " << num_entries<<" "<< 2*atominfo.n_types+1<<"\n";
        if(num_entries >2*atominfo.n_types+1)
        {
          cout << "number and format of element weights is more than types given in the file";
          cout << "this is not accpetable.. exiting\n"; exit(1);
        }else
        {
          for(int i =1;i<num_entries; i=i+2)
          {
            atominfo.element_weights[atoi(ptr_tmp_line[i].c_str())] = atof(ptr_tmp_line[i+1].c_str());
          }
        }
      }

      if(ptr_tmp_line[0]=="element_charges")
      {
        cout << "\n\nentered element charges\n\n";
        if(num_entries >2*atominfo.n_types+1)
        {
          cout << "number and format of element charges is more than types given in the file";
          cout << "this is not accpetable.. exiting\n"; exit(1);
        }else
        {
          for(int i =1;i<num_entries; i=i+2)
          {
            atominfo.element_charges[atoi(ptr_tmp_line[i].c_str())] = atof(ptr_tmp_line[i+1].c_str());
          }
        }
      }

      if(ptr_tmp_line[0]=="element_bonds")
      {
        cout << "element bonds\n";

        if((num_entries-1) % 3 != 0)
        {
          cout << "wrong information - not exact number of bond parameters (3 req)\n";
          cout <<"exiting...\n";
          exit(1);
        }else
        {
          for (int i =1;i <num_entries;i=i+3)
          {
            int local1 = atoi(ptr_tmp_line[i].c_str());
            int local2 = atoi(ptr_tmp_line[i+1].c_str());
            int local3 = atoi(ptr_tmp_line[i+2].c_str());
            auto bondtuple= std::make_tuple(local1, local2, local3);
            atominfo.element_bonds.push_back(bondtuple);
          }
        }
      }

      if(ptr_tmp_line[0]=="element_angles")
      {
        cout << "element angles\n";
        if((num_entries-1) % 4 != 0)
        {
          cout << "wrong information - not exact number of angle parameters (4 req)\n";
          cout <<"exiting...\n";
          exit(1);
        }else
        {
          for (int i =1;i <num_entries;i=i+4)
          {
            int local1 = atoi(ptr_tmp_line[i].c_str());
            int local2 = atoi(ptr_tmp_line[i+1].c_str());
            int local3 = atoi(ptr_tmp_line[i+2].c_str());
            int local4 = atoi(ptr_tmp_line[i+3].c_str());
            auto angletuple = std::make_tuple(local1, local2, local3, local4);
            atominfo.element_angles.push_back(angletuple);
          }
        }
      }

			delete [] ptr_tmp_line;
		}
    if(cutoff_entries !=3*atominfo.n_types*(atominfo.n_types+1)/2)
    {
      cout << " incorrect number of entries 1 1 cutofff 1 2 cutoff 1 3 cutoff etc...exiting\n";
      cout << "expected "<<3*atominfo.n_types*(atominfo.n_types+1)/2<< " but got "<< cutoff_entries <<"\n";
      exit(1);
    }
	}else
	{
		cout << "cutoffs file is corrupted..terminating\n";
		exit(1);
	}
	cutoff_file.close();
  return atominfo;
}
