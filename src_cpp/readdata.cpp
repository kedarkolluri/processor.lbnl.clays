/*
 * ReadData.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */
#include <readdata.h>

simcell read_xyz_VESTA(const char *filename1)
{
	FILE *fptr;
	simcell mdcell;
	mdcell.n = -1;
	int num_val;

	ifstream inputfile;
	char filename[80]="";
	strcat(filename, filename1);
	string tmp_line;
	std::cout << "filename "<< filename<<"\n";
	unzip(filename);
	inputfile.open(filename);
	if(!(inputfile.good()))
	{
		cout << filename << " file open failed\n";
		return(mdcell);
	}else
	{
		bool read_coords = false;
		std::string *ptr_tmp_line1;
		ptr_tmp_line1 = get_next_splits(inputfile, num_val);
		mdcell.n = atoi(ptr_tmp_line1[0].c_str());
		delete [] ptr_tmp_line1;
		ptr_tmp_line1 = get_next_splits(inputfile, num_val);

		double crystal_data[6];
		crystal_data[0] = atof(ptr_tmp_line1[0].c_str());
		crystal_data[1] = atof(ptr_tmp_line1[1].c_str());
		crystal_data[2] = atof(ptr_tmp_line1[2].c_str());
		crystal_data[3] = atof(ptr_tmp_line1[3].c_str())/180*PI;
		crystal_data[4] = atof(ptr_tmp_line1[4].c_str())/180*PI;
		crystal_data[5] = atof(ptr_tmp_line1[5].c_str())/180*PI;

		delete [] ptr_tmp_line1;

		for(int i=0;i<6; i++)
		{
			mdcell.crystal[i] = crystal_data[i];
		}
		crystal_H(mdcell.crystal,mdcell.Hcry);
		M3inv(mdcell.Hcry,mdcell.Hinv);

		cout << "number of atoms are\t"<< mdcell.n <<"\n";
		int tag = 0;
		cout << "initalized the structure\n";
		while(!inputfile.eof())
		{
			std::string *ptr_tmp_line;
			ptr_tmp_line = get_next_splits(inputfile, num_val);

			if(num_val>0)
			{
				atomic_dat atom;

				for(int ih = 0; ih<3;ih++)
				{
					atom.r[ih] = atof(ptr_tmp_line[ih+1].c_str());
					atom.u[ih] = atom.r[ih];
				}

				atom.type = 0;
				atom.molID = 1;

				strncpy(atom.elem, ptr_tmp_line[0].c_str(), sizeof(atom.elem));

				V3mulM3(atom.r,mdcell.Hinv,atom.s);
				for(int ih=0;ih<3;ih++)
				{
					if(atom.s[ih] >= 1.0 ) atom.s[ih] = atom.s[ih]-1.0;
					if(atom.s[ih] <0.0 ) atom.s[ih] = atom.s[ih]+1.0;
				}
				mdcell.atomdata[tag] = atom;
				tag++;
			}
			delete [] ptr_tmp_line;

		}
	}
	inputfile.close();
	return mdcell;
}
/*
int read_lammps(char *filename)
{
	FILE *fptr;
	int num_val;
	ifstream inputfile;
	string tmp_line;

	static bool g_prev_success = true;
	if(g_prev_success)
	{
		if(atom != NULL)
		free(atom);

	}
	unzip(filename);
	cout <<"filename is \t"<<filename<<" in here \n";
	inputfile.open(filename);
	if(!(inputfile.good()))
	{
		cout << filename<<"\t:file open failed\n";
		g_prev_success = false;
		return(1);
	}
	else{
		g_prev_success = true;
		bool read_coords = false;
		int n_counter = -1;
		int id_counter = -1;
		string *ptr_tmp_line1;
		ptr_tmp_line1 = get_next_splits(inputfile, num_val);
		delete[] ptr_tmp_line1;
		ptr_tmp_line1 = get_next_splits(inputfile, num_val);

		delete[] ptr_tmp_line1;
		ptr_tmp_line1 = get_next_splits(inputfile, num_val);
		delete[] ptr_tmp_line1;
		ptr_tmp_line1 = get_next_splits(inputfile, num_val);
		n = atoi(ptr_tmp_line1[0].c_str());
		delete[] ptr_tmp_line1;
		cout << "number of atoms are\t"<<n<<"\n";;
		string *ref_line;
		int iii =0;
		atom = (struct atomic_dat *) malloc((n+5)*sizeof(struct atomic_dat));
		while(!inputfile.eof())
		{
//			cout << "am here at the top\n";
			string *ptr_tmp_line;
//			cout << "going to get values\n";
			ptr_tmp_line = get_next_splits(inputfile, num_val);
//			cout << "got values ..\n";
			if((ptr_tmp_line[0]=="ITEM:")&&(ptr_tmp_line[1]=="BOX")&&(ptr_tmp_line[2]=="BOUNDS"))
			{
				int prev_num_val  = num_val;
				ptr_tmp_line = get_next_splits(inputfile, num_val);
				xlo = atof(ptr_tmp_line[0].c_str()); xhi = atof(ptr_tmp_line[1].c_str()); if(num_val==3) xy = atof(ptr_tmp_line[2].c_str());
				ptr_tmp_line = get_next_splits(inputfile, num_val);
				ylo =atof(ptr_tmp_line[0].c_str()); yhi = atof(ptr_tmp_line[1].c_str()); if(num_val==3) xz = atof(ptr_tmp_line[2].c_str());
				ptr_tmp_line = get_next_splits(inputfile, num_val);
				zlo = atof(ptr_tmp_line[0].c_str()); zhi = atof(ptr_tmp_line[1].c_str()); if(num_val==3) yz = atof(ptr_tmp_line[2].c_str());

				lz = (zhi-zlo);
				xlo = max(xlo,xlo-xy);
				xlo = max(xlo,xlo-xz);
				xhi = min(xhi,xhi-xy);
				xhi = min(xhi,xhi-xz);
				ylo = max(ylo,ylo-yz);
				yhi = min(yhi,yhi-yz);
				lx = xhi-xlo;
				ly = yhi-ylo;
				cout << lx <<"\t"<<ly<<"\t"<<lz<<"\t"<<xy <<"\t"<<xz <<"\t"<<yz<<"\n";

				H[0][0] = lx;H[0][1]=0.0;H[0][2]=0.0;H[1][0]=xy;H[1][1]=ly;H[1][2]=0.0;H[2][0]=xz;H[2][1]=yz;H[2][2]=lz;
				for(int a=0;a<3;a++)
				{
					for(int b=0;b<3;b++)
					{
						Hcry[a][b] = H[a][b];
				//		cout << Hcry[a][b]<<"\t";
					}
				//	cout << "\n";
				}
				H_crystal(H,crystal0);
				crystal_H(crystal0,Hcry);
				//double Hcry_inv1[3][3];
				//M3inv(Hcry,Hcry_inv1);
				M3inv(Hcry,Hcry_inv);

			}
			if((ptr_tmp_line[0]=="ITEM:")&&(ptr_tmp_line[1]=="ATOMS"))
			{
				//ref_line = ptr_tmp_line;
//				cout << "am here  in the refline thing\n";
				ref_line = new string[num_val];
				for (int i =0;i<num_val;i++)
				{
					//ref_line[i] = ptr_tmp_line[i];
					ref_line[i].append(ptr_tmp_line[i]);
//					cout << ref_line[i]<<"\n";
//					strncpy(ref_line[i], ptr_tmp_line[i]);
//					istringstream iss2(ptr_tmp_line[i]);
//					cout << iss2 <<"\n";
//					 iss2 >> ref_line[i];
;				}

				int i =0;
				while (id_counter<0) {if(ref_line[i+2]=="id") id_counter = i;};
				read_coords = true;
				delete [] ptr_tmp_line;
				ptr_tmp_line = get_next_splits(inputfile, num_val);
			}
			//cout << "BEFORE READ COORDS ***********************************\n";
			if(read_coords)
			{
				n_counter++;
//				cout << n_counter << "\n";
				int tag = atoi(ptr_tmp_line[id_counter].c_str());
				tag = tag-1;
				//cout << n_counter <<"\t"<<tag<<" \n";
				for (int i =0;i<num_val;i++)
				{
					//cout << tag<< " " <<i<<"\t"<<atof(ptr_tmp_line[i].c_str())<<"\t"<<atoi(ptr_tmp_line[i].c_str())<<"\n";
					if(ref_line[i+2]=="xs")
					{
						atom[tag].sx =atof(ptr_tmp_line[i].c_str());
						if(atom[tag].sx>=1) atom[tag].sx = atom[tag].sx-1;
						if(atom[tag].sx<0) atom[tag].sx = atom[tag].sx+1;
					}
					if(ref_line[i+2]=="ys")
					{
						atom[tag].sy =atof(ptr_tmp_line[i].c_str());
						if(atom[tag].sy>=1) atom[tag].sy = atom[tag].sy-1;
						if(atom[tag].sy<0) atom[tag].sy = atom[tag].sy+1;
					}
					if(ref_line[i+2]=="zs")
					{
						atom[tag].sz =atof(ptr_tmp_line[i].c_str());
						if(atom[tag].sz>=1) atom[tag].sz = atom[tag].sz-1;
						if(atom[tag].sz<0) atom[tag].sz = atom[tag].sz+1;
					}

					if(ref_line[i+2]=="vx") {atom[tag].vx =atof(ptr_tmp_line[i].c_str());}
					if(ref_line[i+2]=="vy") {atom[tag].vy =atof(ptr_tmp_line[i].c_str());}
					if(ref_line[i+2]=="vz"){atom[tag].vz =atof(ptr_tmp_line[i].c_str());}

					if(ref_line[i+2]=="type") {atom[tag].type =atoi(ptr_tmp_line[i].c_str());}
					if(ref_line[i+2]=="mass") {atom[tag].ma =atof(ptr_tmp_line[i].c_str());}

					if(ref_line[i+2]=="c_energy") {atom[tag].pe =atof(ptr_tmp_line[i].c_str());}
					if(ref_line[i+2]=="q") {atom[tag].charge =atof(ptr_tmp_line[i].c_str());}
					if(ref_line[i+2]=="element") { strncpy(atom[tag].elem,ptr_tmp_line[i].c_str(),sizeof(atom[tag].elem));}
				}



				if(n_counter>=n-1) read_coords = false;

			}
//			cout << " am here\n";
			delete [] ptr_tmp_line;
//			cout << "done it\n";

		}
		delete [] ref_line;

	}

	inputfile.close();
	zip(filename);
	//Adjustments
	// 1. if element is not there, 1 is Cu, 2 is Nb and others are ZZ
	// 2. set the number of types
	for (int i =0;i <n;i++)
	{	//if(i==0) {cout << strcmp(atom[i].elem,"")<<"\n";cout << strcmp("abcd","abcd")<<"\n";cout << strcmp("efgh","abcd")<<"\n";}
		//cout << sizeof(atom[i].elem) <<" "<<i<<"\n";
		if(strcmp(atom[i].elem,"")==0)
		{
			//cout << i<< " hi\n";
			char abcd[80]="";
			if(atom[i].type == 1)
			{
				strcat(abcd,"Mg");
			}else if(atom[i].type == 2)
			{
				strcat(abcd,"O");
			}else
			{
				strcat(abcd,"Z");
			}
			strncpy(atom[i].elem, abcd,sizeof(atom[i].elem));
			//cout << abcd <<" "<<atom[i].elem<<"\n";
			//strncpy(atom[i].elem, element_names[atom[i].type-1].c_str(),sizeof(atom[i].elem));
			//cout << abcd <<" "<<atom[i].elem<<"\n";
			//cout << abcd <<"\n";
		}


		if(n_types< atom[i].type) n_types = atom[i].type;

	}

	cout << "total types of atoms are\t"<< n_types <<"\n";
	//cout << "exiting here\n";

	return(0);
}
*/
