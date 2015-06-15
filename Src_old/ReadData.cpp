/*
 * ReadData.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */
#include <ReadData.h>
//#include <F_read_lammps_dump.h>

int read_A(void)
{
	double dummy1;
	int dummy2;
	char *s;
	cout << inputfilename<<"\n";
	fptr = fopen(inputfilename,"r");
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&n);
	cout << n<<"\n";
	cout << dummy1<<"\n";

	fscanf(fptr,"%lf %lf %lf",&H0[0][0],&H0[0][1],&H0[0][2]);
	fscanf(fptr,"%lf %lf %lf",&H0[1][0],&H0[1][1],&H0[1][2]);
	fscanf(fptr,"%lf %lf %lf",&H0[2][0],&H0[2][1],&H0[2][2]);

	// FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2
	for (int i=0;i<3;i++) H0[2][i] = 2*H0[2][i];
	// FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2

	M3inv(H0, H0_inv);



	for (int i=0;i<=2;i++)
	{for (int j=0;j<=2;j++)
	{cout << H0[i][j]<<"\t";}
		cout <<"\n";
	}
	cout <<"\n\n";
	for (int i=0;i<=2;i++)
	{for (int j=0;j<=2;j++)
	{cout << H0_inv[i][j]<<"\t";}
		cout <<"\n";
	}

	cout << "before collecting atoms\n";
	free(atom);
	atom = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));


	for(int i=0;i<n;i++)
	{
		fscanf(fptr,"%d %lf %lf %lf",&atom[i].type,&atom[i].rx,&atom[i].ry,&atom[i].rz);
	}
	fclose(fptr);
	cout << atom[5].rx<<"\t"<<atom[5].ry<<"\t"<<atom[5].rz<<"\n";
	cout << "after collecting atoms\n";

	for(int i=0;i<n;i++)
	{
		double r[3],s[3];
		r[0]=atom[i].rx;r[1]=atom[i].ry;r[2]=atom[i].rz;
		atom[i].ux=atom[i].rx;atom[i].uy=atom[i].ry;atom[i].uz=atom[i].rz;

		V3mulM3(r,H0_inv,s);
		atom[i].sx = s[0];atom[i].sy=s[1];atom[i].sz=s[2];


		if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
		if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
		if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;

		if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
		if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
		if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;


		/*
		 if((atom[i].sz >z_start)&&(atom[i].sz<z_end))
		 {
			 interface_atoms++;

		 }
		 */

		if(atom[i].type==1)
		{
			atom[i].ma = 63.546;
		}else
		{
			if(atom[i].type==2)
			{
				atom[i].ma = 92.90638;
			}
		}

	}

	cout << atom[5].ma<<"\t"<<atom[5].sx<<"\t"<<atom[5].sy<<"\t"<<atom[5].sz<<" 88888888888888888888888\n";


	//identify atoms close to the interface - adhoc procedure only - select atoms with sz between 0.55 and 0.61
	/*
	 tag_array_interface_atoms = (int *) malloc((interface_atoms+3)*sizeof(int));
	 interface_atoms = 0;

	 for (int i = 0;i<n;i++)
	 {
		 if((atom[i].sz >z_start)&&(atom[i].sz<z_end))
		 {
			 tag_array_interface_atoms[interface_atoms] = i;
			 interface_atoms++;

		 }
	 }

	 sort(interface_atoms, tag_array_interface_atoms);
	 */
	cout <<"\n"<<interface_atoms<<" total atoms\n";

	//identify atoms close to the interface - adhoc procedure only - select atoms with sz between 0.55 and 0.61 - END




	H_crystal(H0,crystal0);
	crystal_H(crystal0,Hcry);

	for (int i=0;i<=2;i++)
	{for (int j=0;j<=2;j++)
	{cout << Hcry[i][j]<<"\t";}
		cout <<"\n";
	}


	for (int i=0;i<6;i++)
	{
		cout << crystal0[i]<<"\t";
	}
	cout <<"\n";


	return(0);
}

int read_A_exact(char *filename)
{
	double dummy1;
	int dummy2;
	char *s;
	cout << filename<<"\n";
	fptr = fopen(filename,"r");
	if(fptr!=NULL)
	{
	fscanf(fptr,"%d",&dummy2);
	cout << dummy2<<"\n";

	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&dummy2);
	fscanf(fptr,"%d",&n);
	cout << n<<"\n";
	cout << dummy1<<"\n";

	fscanf(fptr,"%lf %lf %lf",&H0[0][0],&H0[0][1],&H0[0][2]);
	fscanf(fptr,"%lf %lf %lf",&H0[1][0],&H0[1][1],&H0[1][2]);
	fscanf(fptr,"%lf %lf %lf",&H0[2][0],&H0[2][1],&H0[2][2]);

	// FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2
	for (int i=0;i<3;i++) H0[2][i] = H0[2][i];
	// FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2

	M3inv(H0, H0_inv);



	for (int i=0;i<=2;i++)
	{for (int j=0;j<=2;j++)
	{cout << H0[i][j]<<"\t";}
		cout <<"\n";
	}
	cout <<"\n\n";
	for (int i=0;i<=2;i++)
	{for (int j=0;j<=2;j++)
	{cout << H0_inv[i][j]<<"\t";}
		cout <<"\n";
	}

	cout << "before collecting atoms\n";
	free(atom);

	atom = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

	cout << "before collecting atoms - allocated\n";

	for(int i=0;i<n;i++)
	{
		fscanf(fptr,"%d %lf %lf %lf",&atom[i].type,&atom[i].rx,&atom[i].ry,&atom[i].rz);
//		cout << i<<"\t"<<n<<"\n";
	}
	fclose(fptr);
	cout << atom[5].rx<<"\t"<<atom[5].ry<<"\t"<<atom[5].rz<<"\n";
	cout << "after collecting atoms\n";

	for(int i=0;i<n;i++)
	{
		double r[3],s[3];
		r[0]=atom[i].rx;r[1]=atom[i].ry;r[2]=atom[i].rz;
		atom[i].ux=atom[i].rx;atom[i].uy=atom[i].ry;atom[i].uz=atom[i].rz;

		V3mulM3(r,H0_inv,s);
		atom[i].sx = s[0];atom[i].sy=s[1];atom[i].sz=s[2];


		if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
		if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
		if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;

		if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
		if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
		if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;


		/*
		 if((atom[i].sz >z_start)&&(atom[i].sz<z_end))
		 {
			 interface_atoms++;

		 }
		 */

		if(atom[i].type==1)
		{
			atom[i].ma = 63.546;
		}else
		{
			if(atom[i].type==2)
			{
				atom[i].ma = 92.90638;
			}
		}

	}

	cout << atom[5].ma<<"\t"<<atom[5].sx<<"\t"<<atom[5].sy<<"\t"<<atom[5].sz<<" 88888888888888888888888\n";


	//identify atoms close to the interface - adhoc procedure only - select atoms with sz between 0.55 and 0.61
	/*
	 tag_array_interface_atoms = (int *) malloc((interface_atoms+3)*sizeof(int));
	 interface_atoms = 0;

	 for (int i = 0;i<n;i++)
	 {
		 if((atom[i].sz >z_start)&&(atom[i].sz<z_end))
		 {
			 tag_array_interface_atoms[interface_atoms] = i;
			 interface_atoms++;

		 }
	 }

	 sort(interface_atoms, tag_array_interface_atoms);
	 */

	//identify atoms close to the interface - adhoc procedure only - select atoms with sz between 0.55 and 0.61 - END




	H_crystal(H0,crystal0);
	crystal_H(crystal0,Hcry);

	for (int i=0;i<=2;i++)
	{for (int j=0;j<=2;j++)
	{cout << Hcry[i][j]<<"\t";}
		cout <<"\n";
	}


	for (int i=0;i<6;i++)
	{
		cout << crystal0[i]<<"\t";
	}
	cout <<"\n";

	}else
	{
		return -1;
	}
	return(0);
}

int read_xian_ming_blas(char *filename)
{
  double dummy1;
  int dummy2;
  char s[80];
  cout << filename<<"\n";
  fptr = fopen(filename,"r");
  if(fptr!=NULL)
  {
  fscanf(fptr,"%d",&n);
  fscanf(fptr,"%d %d",&dummy1,&dummy2);
  cout << dummy2<<"\t"<< dummy1<<"\n";
  cout << n<<"\n";
  fscanf(fptr,"%s",s);
  cout <<s<<" out of reading the stuff\n";
  fscanf(fptr,"%s",s);
  cout <<s<<" out of reading the stuff\n";
  fscanf(fptr,"%s",s);
  cout <<s<<" out of reading the stuff\n";
  fscanf(fptr,"%lf %lf %lf",&H0[0][0],&H0[1][1],&H0[2][2]);
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) if (i!=j) H0[i][j] = 0;
  // FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2
  fscanf(fptr,"%s",s);
  M3inv(H0, H0_inv);



  for (int i=0;i<=2;i++)
  {for (int j=0;j<=2;j++)
  {cout << H0[i][j]<<"\t";}
          cout <<"\n";
  }
  cout <<"\n\n";
  for (int i=0;i<=2;i++)
  {for (int j=0;j<=2;j++)
  {cout << H0_inv[i][j]<<"\t";}
          cout <<"\n";
  }

  cout << "before collecting atoms\n";
  free(atom);

  atom = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

  cout << "before collecting atoms - allocated\n";

  for(int i=0;i<n;i++)
  {
          fscanf(fptr,"%lf %lf %lf %d",&atom[i].rx,&atom[i].ry,&atom[i].rz,&atom[i].type);
  }
  fclose(fptr);
  cout << atom[5].rx<<"\t"<<atom[5].ry<<"\t"<<atom[5].rz<<"\n";
  cout << "after collecting atoms\n";

  for(int i=0;i<n;i++)
  {
          double r[3],s[3];
          r[0]=atom[i].rx;r[1]=atom[i].ry;r[2]=atom[i].rz;
          atom[i].ux=atom[i].rx;atom[i].uy=atom[i].ry;atom[i].uz=atom[i].rz;

          V3mulM3(r,H0_inv,s);
          atom[i].sx = s[0];atom[i].sy=s[1];atom[i].sz=s[2];


          if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
          if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
          if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;

          if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
          if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
          if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;

          if(atom[i].type==1)
          {
                  atom[i].ma = 63.546;
          }else
          {
                  if(atom[i].type==2)
                  {
                          atom[i].ma = 92.90638;
                  }
          }

  }

  cout << atom[5].ma<<"\t"<<atom[5].sx<<"\t"<<atom[5].sy<<"\t"<<atom[5].sz<<" 111\n";




  H_crystal(H0,crystal0);
  crystal_H(crystal0,Hcry);

  for (int i=0;i<=2;i++)
  {for (int j=0;j<=2;j++)
  {cout << Hcry[i][j]<<"\t";}
          cout <<"\n";
  }


  for (int i=0;i<6;i++)
  {
          cout << crystal0[i]<<"\t";
  }
  cout <<"\n";

  }else
  {
          return -1;
  }
  return(0);
}


int read_xian_ming_blas2(char *filename)
{
  double dummy1,dummy12,dummy13;
  int dummy2;
  char s[80];
  cout << filename<<"\n";
  fptr = fopen(filename,"r");
  if(fptr!=NULL)
  {
  fscanf(fptr,"%d",&n);
  fscanf(fptr,"%d %d",&dummy1, &dummy2);
  cout << dummy2<<"\t"<< dummy1<<"\n";
  cout << n<<"\n";
  fscanf(fptr,"%lf %lf %lf",&H0[0][0],&H0[1][1],&H0[2][2]);
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) if (i!=j) H0[i][j] = 0;
//  H0[2][2]=3*H0[2][2];
  // FOR LAMMPS - create vaccume on top and bottom by multiplying H[2][*] by 2
  M3inv(H0, H0_inv);



  for (int i=0;i<=2;i++)
  {for (int j=0;j<=2;j++)
  {cout << H0[i][j]<<"\t";}
          cout <<"\n";
  }
  cout <<"\n\n";
  for (int i=0;i<=2;i++)
  {for (int j=0;j<=2;j++)
  {cout << H0_inv[i][j]<<"\t";}
          cout <<"\n";
  }

  cout << "before collecting atoms\n";
  free(atom);

  atom = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

  cout << "before collecting atoms - allocated\n";

  for(int i=0;i<n;i++)
  {
          fscanf(fptr,"%lf %lf %lf %d %lf %lf %lf",&atom[i].rx,&atom[i].ry,&atom[i].rz,&atom[i].type,&dummy1,&dummy12,&dummy13);
  }
  fclose(fptr);
  cout << atom[5].rx<<"\t"<<atom[5].ry<<"\t"<<atom[5].rz<<"\n";
  cout << "after collecting atoms\n";

  for(int i=0;i<n;i++)
  {
          double r[3],s[3];

        // atom[i].rz=atom[i].rz-1.4;

          r[0]=atom[i].rx;r[1]=atom[i].ry;r[2]=atom[i].rz;
          atom[i].ux=atom[i].rx;atom[i].uy=atom[i].ry;atom[i].uz=atom[i].rz;

          V3mulM3(r,H0_inv,s);
//          atom[i].sx = s[0];atom[i].sy=s[1];atom[i].sz=s[2]-0.5;
          atom[i].sx = s[0];atom[i].sy=s[1];atom[i].sz=s[2];

          if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
          if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
          if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;

          if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
          if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
          if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;

//          if(atom[i].sz>0.50) atom[i].type=2;

          if(atom[i].type==1)
          {
                  atom[i].ma = 63.546;
          }else
          {
                  if(atom[i].type==2)
                  {
                          atom[i].ma = 92.90638;
                  }
          }

  }

  cout << atom[5].ma<<"\t"<<atom[5].sx<<"\t"<<atom[5].sy<<"\t"<<atom[5].sz<<" 111\n";




  H_crystal(H0,crystal0);
  crystal_H(crystal0,Hcry);

  for (int i=0;i<=2;i++)
  {for (int j=0;j<=2;j++)
  {cout << Hcry[i][j]<<"\t";}
          cout <<"\n";
  }


  for (int i=0;i<6;i++)
  {
          cout << crystal0[i]<<"\t";
  }
  cout <<"\n";

  }else
  {
          return -1;
  }
  return(0);
}



int read_cfg(char *filename, atomic_dat *atom_fill)
{

	//char filename[80]="";
	char getstring[STRING_LIMIT];
	ifstream inputfile;
	char str[80];
	//strcat(filename,Filename_prefix);
	//sprintf(str,"%d",i);
	//strcat(filename,str);strcat(filename,Filename_sufix);
	cout <<"filename is \t"<<filename<<" in here \n";
	inputfile.open(filename,ifstream::in);
	double atom_mass;
	int it=0;
	char *str_seq;
	int data_seq=3;
	double ma_l; int type;
	int i =0;
	int initial_setup = true;
	while((inputfile.good())&&(!inputfile.eof()))
	{
		inputfile.getline(getstring,STRING_LIMIT);
		//cout << getstring<<"\n";
		it++;
		char *pch;
		pch = strtok(getstring," =\t");
		int wait = -1;
		int current = 0;

		char first_string[80]="";

		while(pch!=NULL)
		{
			current++;
			if(!(pch[0]=='#'))
			{
				// RULES FOR CFG FORMAT ARE ENTERED HERE
				if(initial_setup)
				{
					//cout << pch<<"\t";
					if((current==1)&&(strcmp(pch,"Number")==0)) {wait = 4;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"A")==0)) {wait = -1;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(1,1)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(1,2)")==0)){wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(1,3)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(2,1)")==0)){wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(2,2)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(2,3)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(3,1)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(3,2)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"H0(3,3)")==0)) {wait = 2;strcat(first_string,pch);}


					if((current==1)&&(strcmp(pch,"Transform(1,1)")==0)){wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"Transform(1,2)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"Transform(1,3)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"Transform(2,1)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"Transform(2,2)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"Transform(2,3)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"Transform(3,1)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"Transform(3,2)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"Transform(3,3)")==0)) {wait = 2;strcat(first_string,pch);}

					if((current==1)&&(strcmp(pch,"eta(1,1)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"eta(1,2)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"eta(1,3)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"eta(2,1)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"eta(2,2)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"eta(2,3)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"eta(3,1)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"eta(3,2)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"eta(3,3)")==0)) {wait = 2;strcat(first_string,pch);}
					if((current==1)&&(strcmp(pch,"entry_count")==0)) {wait = 2;strcat(first_string,pch);}



					if((current==wait)&&(strcmp(first_string,"Number")==0))
					{
						n = atoi(pch);
						cout << n<<"\n";
						atom = (struct atomic_dat *) malloc((n+5)*sizeof(struct atomic_dat));
					}
					if((current==wait)&&(strcmp(first_string,"H0(1,1)")==0)) {H[0][0] = strtod(pch,NULL);}
					if((current==wait)&&(strcmp(first_string,"H0(1,2)")==0)){H[0][1] = strtod(pch,NULL);}
					if((current==wait)&&(strcmp(first_string,"H0(1,3)")==0)){H[0][2] = strtod(pch,NULL);}
					if((current==wait)&&(strcmp(first_string,"H0(2,1)")==0)) {H[1][0] = strtod(pch,NULL);}
					if((current==wait)&&(strcmp(first_string,"H0(2,2)")==0)) {H[1][1] = strtod(pch,NULL);}
					if((current==wait)&&(strcmp(first_string,"H0(2,3)")==0)) {H[1][2] = strtod(pch,NULL);}
					if((current==wait)&&(strcmp(first_string,"H0(3,1)")==0)) {H[2][0] = strtod(pch,NULL);}
					if((current==wait)&&(strcmp(first_string,"H0(3,2)")==0)) {H[2][1] = strtod(pch,NULL);}
					if((current==wait)&&(strcmp(first_string,"H0(3,3)")==0)) {H[2][2] = strtod(pch,NULL);}

					/*
					 if((current==wait)&&(strcmp(first_string,"Transform(1,1)")==0)){wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"Transform(1,2)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"Transform(1,3)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"Transform(2,1)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"Transform(2,2)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"Transform(2,3)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"Transform(3,1)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"Transform(3,2)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"Transform(3,3)")==0)) {wait = 2;strcat(first_string,pch);}

					 if((current==wait)&&(strcmp(first_string,"eta(1,1)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"eta(1,2)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"eta(1,3)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"eta(2,1)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"eta(2,2)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"eta(2,3)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"eta(3,1)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"eta(3,2)")==0)) {wait = 2;strcat(first_string,pch);}
					 if((current==wait)&&(strcmp(first_string,"eta(3,3)")==0)) {wait = 2;strcat(first_string,pch);}
					 */
					//cout << H[0][0]<<"\t"<<H[1][1]<<"\t"<<H[1][2]<<" HA HAH AHHAHA\n";
					//cout << "first string value is \t"<<first_string<<"\n";
					if((current==wait)&&(strcmp(first_string,"entry_count")==0)) {data_seq=atoi(pch); initial_setup = false;cout <<"going to get out\n";}
				} else

				{
					//cout << "out in else\t"<<pch<<"\t"<< current<<"\n";
					if((strcmp(pch,"63.546")==0)&&(current==1))
					{
						type = 1;
						ma_l = 63.546;

					}else if((strcmp(pch,"Cu")==0)&&(current==2))
					{
						type = 1;
						ma_l = 63.546;
					}else if((strcmp(pch,"92.90638")==0)&&(current==1))
					{
						type = 2;
						ma_l = 92.90638;
					}else if((strcmp(pch,"Nb")==0)&&(current==2))
					{
						type = 2;
						ma_l = 92.90638;
					}else if((strcmp(pch,"196.97")==0)&&(current==1))
					{
						type = 1;
						ma_l = 196.97;
					}else if((strcmp(pch,"Au")==0)&&(current==2))
					{
						type = 1;
						ma_l = 196.97;
					} else
					{
						if(current==3)
						{

							atom[i].sx = strtod(pch,NULL);
							atom[i].type = type;
							atom[i].ma = ma_l;
								//cout << i<<"\t"<<atom[i].type <<"\t"<<atom[i].sx<<"  inside data stuff\n";
						}else if(current==4)
						{
							atom[i].sy = strtod(pch,NULL);
							//	cout <<atom[i].sy<<"\t";
						}else if(current==5)
						{
							atom[i].sz = strtod(pch,NULL);
							//	cout<<atom[i].sz<<"\n";
							i++;
						}

					}

				}
				//cout <<wait<<endl;
				pch = strtok(NULL," =\t");


			}else
			{
				pch = NULL;
			}

		}

	}
	inputfile.close();
	cout << i<<"\n";
	for(int j=0;j<3;j++)
	{
		for(int k=0;k<3;k++)
		{
			Hcry[j][k] = H[j][k];
		}
	}

	H_crystal(H,crystal0);
		crystal_H(crystal0,Hcry);
	return(0);
}

int read_lammps_specific(char *filename, atomic_dat *atom_fill, bool create_H, bool new_format, double H_here[3][3])
{
	FILE *fptr;
	static bool prev_success = true;
	if(prev_success)
		{
			if(atom_fill != NULL)
			free(atom_fill);

		}
	unzip(filename);
	//printf("opening data file %s ...\n",filename);
	fptr = fopen(filename,"r");
	if(fptr==NULL)
	{
		cout << filename<<"\t:file open failed\n";
		prev_success = false;
		return(1);

	}else
	{
		prev_success = true;
		char str3[80];
		char str32[80];
		char str33[80];

		double time_step,dummy2,dummy3;


		fscanf(fptr,"%s %s",str32,str33);
	 //   cout << str32<<"\t"<<str33<<" 1 \n";
	    fscanf(fptr,"%d",&time_step);
	    fscanf(fptr,"%s %s",str32,str33);
	   // cout << str32<<"\t"<<str33<<" 2 \n";
	    fscanf(fptr,"%s %s",str32,str33);
	   // cout << str32<<"\t"<<str33<<" 3 \n";
	    fscanf(fptr,"%d",&n);
	    fscanf(fptr,"%s %s",str32,str33);
	   // cout << str32<<"\t"<<str33<<" 4 \n";
	    fscanf(fptr,"%s",str32);
	    //cout << str32<<"5 \n";
		if(new_format)
		{
			char str_a[80],str_b[80],str_c[80];
			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
			//cout << str_a<<"\t"<<str_c<<" format changed\n";
		}
		xy = 0.0;xz=0.0;yz=0.0;
		if(new_format)
		{

			fscanf(fptr,"%lf %lf %lf",&xlo,&xhi,&xy);
			//cout << xlo <<"\tdfd\t"<<xhi<<"\t"<<xy<<"\n";
			//lx = (xhi - xlo);
			fscanf(fptr,"%lf %lf %lf",&ylo,&yhi,&xz);
			//ly = (yhi - ylo);
//			cout << ylo<<"\t ly " << yhi<<"\t"<<xz<<"\n";;

			double dummy_z;
			fscanf(fptr,"%lf",&zlo);
			//	    zlo = dummy3;
			fscanf(fptr,"%lf",&zhi);
			fscanf(fptr,"%lf",&yz);
	//		cout << yz<<" yz is \n";
			lz = (zhi-zlo);

			xlo = max(xlo,xlo-xy);
			xlo = max(xlo,xlo-xz);

			xhi = min(xhi,xhi-xy);
			xhi = min(xhi,xhi-xz);

			ylo = max(ylo,ylo-yz);
			yhi = min(yhi,yhi-yz);


			lx = xhi-xlo;
			ly = yhi-ylo;

		//	cout << zlo<<"\t z "<< zhi<<"\n";
			fscanf(fptr,"%s %s",str32,str33);
			//cout << str32<<"\t"<<str33<<"\n";
		//	cout << lx<<"\t"<<ly<<"\t"<<lz<<"\n";

			char str_a[80],str_b[80],str_c[80];

			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
		//	cout << str_a<<"\t"<<str_c<<" format changed\n";
			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
		//	cout << str_a<<"\t"<<str_c<<" format changed\n";
			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
		//	cout << str_a<<"\t"<<str_c<<" format changed\n";
			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
		//	cout << str_a<<"\t"<<str_c<<" format changed\n";
			/*
			 if((xy*xz)>0)
			 {
			 }else
			 {
			 }
			 */
		}else
		{
			fscanf(fptr,"%lf %lf",&xlo,&xhi);
		//	cout << xlo <<"\tdfd\t"<<xhi<<"\n";
			lx = (xhi - xlo);
		//	cout << xlo<<"\t"<<lx<<"\n";
			fscanf(fptr,"%lf %lf",&ylo,&yhi);
			//	    fscanf(fptr,"%f",&dummy2);
			//	    fscanf(fptr,"%f",&dummy2);
			ly = (yhi - ylo);
		//	cout << ylo<<"\t ly " << yhi<<"\n";;

			double dummy_z;
			fscanf(fptr,"%lf",&zlo);
			//	    zlo = dummy3;
			fscanf(fptr,"%lf",&zhi);
			lz = (zhi-zlo);
			cout << zlo<<"\t z "<< zhi<<"\n";
			fscanf(fptr,"%s %s",str32,str33);
//			cout << str32<<"\t"<<str33<<"\n";
//			cout << lx<<"\t"<<ly<<"\t"<<lz<<"\n";
		}

		if(create_H)
		{
			H_here[0][0] = lx;H_here[0][1]=0.0;H_here[0][2]=0.0;
			H_here[1][0]=xy;H_here[1][1]=ly;H_here[1][2]=0.0;H_here[2][0]=xz;H_here[2][1]=yz;H_here[2][2]=lz;

		}

		double H_herecry_inv1[3][3];
		M3inv(H_here,H_herecry_inv1);

		atom_fill = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

		double rx,ry,rz,vx,vy,vz,fx,fy,fz,sx,sy,sz,ux,uy,uz,ke,pe;
		int tag,type;
		double test1,test2;

		for (int i = 0; i < n; i++)
		{


			//fscanf(fptr,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tag, &type, &sx, &sy, &sz,&ux,&uy,&uz,&vx,&vy,&vz,&pe);
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &test1, &test2, &sx, &sy, &sz,&ux,&uy,&uz,&vx,&vy,&vz,&pe);
			tag = (int) test1;
			type = (int) test2;
// -0.3 for moving in and +0.3 for moving back out
			atom_fill[tag-1].sx = sx;//+0.3;
			atom_fill[tag-1].sy = sy;
			atom_fill[tag-1].sz = sz;

			if(atom_fill[tag-1].sx>=1) atom_fill[tag-1].sx=atom_fill[tag-1].sx-1;
			if(atom_fill[tag-1].sx<0) atom_fill[tag-1].sx=1+atom_fill[tag-1].sx;


			if(atom_fill[tag-1].sy>=1) atom_fill[tag-1].sy=atom_fill[tag-1].sy-1;
			if(atom_fill[tag-1].sy<0) atom_fill[tag-1].sy=1+atom_fill[tag-1].sy;

			if(atom_fill[tag-1].sz>=1) atom_fill[tag-1].sz=atom_fill[tag-1].sz-1;
			if(atom_fill[tag-1].sz<0) atom_fill[tag-1].sz=1+atom_fill[tag-1].sz;


			atom_fill[tag-1].pe = pe;
			atom_fill[tag-1].type = type;

			if(atom_fill[tag-1].type==1)
			{
				atom_fill[tag-1].ma = 63.546;
				atom_fill[tag-1].pe = atom_fill[tag-1].pe; //subtract cohesive energy from the energy of the atom
			}else
			{
				if(atom_fill[tag-1].type==2)
				{
					atom_fill[tag-1].ma = 92.90638;
					atom_fill[tag-1].pe = atom_fill[tag-1].pe;
				}
			}

			atom_fill[tag-1].vx = vx;
			atom_fill[tag-1].vy = vy;
			atom_fill[tag-1].vz = vz;

			atom_fill[tag-1].ux = ux;
			atom_fill[tag-1].uy = uy;
			atom_fill[tag-1].uz = uz;

			atom_fill[tag-1].fx = ux;
			atom_fill[tag-1].fy = uy;
			atom_fill[tag-1].fz = uz;

		}

		fclose(fptr);
		zip(filename);

	}

	//cout << "out of reading lammps\n";
	return(0);
}


int read_lammps(char *filename, atomic_dat *atom_fill, bool create_H, bool new_format)
{
	cout << "READING LAMMPS IN SPECIFIC FORMAT OF Cu-Nb RELATED... CHECK FOR THIS IF ERRORS ARE OCCURING\n";
	FILE *fptr;
	static bool prev_success = true;
	if(prev_success)
		{
			if(atom != NULL)
			free(atom);

		}
	unzip(filename);
	//printf("opening data file %s ...\n",filename);
	fptr = fopen(filename,"r");
	if(fptr==NULL)
	{
		cout << filename<<"\t:file open failed\n";
		prev_success = false;
		return(1);

	}else
	{
		prev_success = true;
		char str3[80];
		char str32[80];
		char str33[80];

		double time_step,dummy2,dummy3;


		fscanf(fptr,"%s %s",str32,str33);
	 //   cout << str32<<"\t"<<str33<<" 1 \n";
	    fscanf(fptr,"%d",&time_step);
	    fscanf(fptr,"%s %s",str32,str33);
	   // cout << str32<<"\t"<<str33<<" 2 \n";
	    fscanf(fptr,"%s %s",str32,str33);
	   // cout << str32<<"\t"<<str33<<" 3 \n";
	    fscanf(fptr,"%d",&n);
	    fscanf(fptr,"%s %s",str32,str33);
	   // cout << str32<<"\t"<<str33<<" 4 \n";
	    fscanf(fptr,"%s",str32);
	    //cout << str32<<"5 \n";
		if(new_format)
		{
			char str_a[80],str_b[80],str_c[80];
			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
			//cout << str_a<<"\t"<<str_c<<" format changed\n";
		}
		xy = 0.0;xz=0.0;yz=0.0;
		if(new_format)
		{

			fscanf(fptr,"%lf %lf %lf",&xlo,&xhi,&xy);
			//cout << xlo <<"\tdfd\t"<<xhi<<"\t"<<xy<<"\n";
			//lx = (xhi - xlo);
			fscanf(fptr,"%lf %lf %lf",&ylo,&yhi,&xz);
			//ly = (yhi - ylo);
//			cout << ylo<<"\t ly " << yhi<<"\t"<<xz<<"\n";;

			double dummy_z;
			fscanf(fptr,"%lf",&zlo);
			//	    zlo = dummy3;
			fscanf(fptr,"%lf",&zhi);
			fscanf(fptr,"%lf",&yz);
	//		cout << yz<<" yz is \n";
			lz = (zhi-zlo);

			xlo = max(xlo,xlo-xy);
			xlo = max(xlo,xlo-xz);

			xhi = min(xhi,xhi-xy);
			xhi = min(xhi,xhi-xz);

			ylo = max(ylo,ylo-yz);
			yhi = min(yhi,yhi-yz);


			lx = xhi-xlo;
			ly = yhi-ylo;

		//	cout << zlo<<"\t z "<< zhi<<"\n";
			fscanf(fptr,"%s %s",str32,str33);
			//cout << str32<<"\t"<<str33<<"\n";
		//	cout << lx<<"\t"<<ly<<"\t"<<lz<<"\n";

			char str_a[80],str_b[80],str_c[80];

			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
		//	cout << str_a<<"\t"<<str_c<<" format changed\n";
			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
		//	cout << str_a<<"\t"<<str_c<<" format changed\n";
			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
		//	cout << str_a<<"\t"<<str_c<<" format changed\n";
			fscanf(fptr,"%s %s %s",str_a,str_b,str_c);
		//	cout << str_a<<"\t"<<str_c<<" format changed\n";
			/*
			 if((xy*xz)>0)
			 {
			 }else
			 {
			 }
			 */
		}else
		{
			fscanf(fptr,"%lf %lf",&xlo,&xhi);
		//	cout << xlo <<"\tdfd\t"<<xhi<<"\n";
			lx = (xhi - xlo);
		//	cout << xlo<<"\t"<<lx<<"\n";
			fscanf(fptr,"%lf %lf",&ylo,&yhi);
			//	    fscanf(fptr,"%f",&dummy2);
			//	    fscanf(fptr,"%f",&dummy2);
			ly = (yhi - ylo);
		//	cout << ylo<<"\t ly " << yhi<<"\n";;

			double dummy_z;
			fscanf(fptr,"%lf",&zlo);
			//	    zlo = dummy3;
			fscanf(fptr,"%lf",&zhi);
			lz = (zhi-zlo);
			cout << zlo<<"\t z "<< zhi<<"\n";
			fscanf(fptr,"%s %s",str32,str33);
//			cout << str32<<"\t"<<str33<<"\n";
//			cout << lx<<"\t"<<ly<<"\t"<<lz<<"\n";
		}

		if(create_H)
		{
			H[0][0] = lx;H[0][1]=0.0;H[0][2]=0.0;H[1][0]=xy;H[1][1]=ly;H[1][2]=0.0;H[2][0]=xz;H[2][1]=yz;H[2][2]=lz;
			for(int a=0;a<3;a++)
			{
				for(int b=0;b<3;b++)
				{
					Hcry[a][b] = H[a][b];
					cout << Hcry[a][b]<<"\t";
				}
				cout << "\n";
			}
			H_crystal(H,crystal0);
			crystal_H(crystal0,Hcry);
		}

		double Hcry_inv1[3][3];
		M3inv(Hcry,Hcry_inv1);
		M3inv(Hcry,Hcry_inv);

		atom = (struct atomic_dat *) malloc((n+5)*sizeof(struct atomic_dat));

		double rx,ry,rz,vx,vy,vz,fx,fy,fz,sx,sy,sz,ux,uy,uz,ke,pe;
		int tag,type;
		double test1,test2;

		for (int i = 0; i < n; i++)
		{

			//fscanf(fptr,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			//  &tag,&type &rx, &ry, &rz,&vx,&vy,&vz,&fx,&fy,&fz);
			/*
			 fscanf(fptr,"%d %d %lf %lf %lf %lf %lf %lf", &tag, &type, &rx, &ry, &rz,&vx,&vy,&vz);

			 atom[tag-1].rx = rx;
			 atom[tag-1].ry = ry;
			 atom[tag-1].rz = rz;

			 atom[tag-1].type = type;

			 atom[tag-1].vx = vx;
			 atom[tag-1].vy = vy;
			 atom[tag-1].vz = vz;
			 */

			//fscanf(fptr,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tag, &type, &sx, &sy, &sz,&ux,&uy,&uz,&vx,&vy,&vz,&pe);
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &test1, &test2, &sx, &sy, &sz,&ux,&uy,&uz,&vx,&vy,&vz,&pe);
			tag = (int) test1;
			type = (int) test2;
// -0.3 for moving in and +0.3 for moving back out
			atom[tag-1].sx = sx;//-0.3;
			atom[tag-1].sy = sy;
			atom[tag-1].sz = sz;


			/*
			 atom[tag-1].rx = sx*(xhi-xlo);
			 atom[tag-1].ry = sy*(yhi-ylo);
			 atom[tag-1].rz = sz*(zhi-zlo);
			 double r[3],s[3];
			 r[0] = atom[tag-1].rx;
			 r[1] = atom[tag-1].ry;
			 r[2] = atom[tag-1].rz;

			 V3mulM3(r,Hcry_inv1,s);

			 atom[tag-1].sx = s[0];
			 atom[tag-1].sy = s[1];
			 atom[tag-1].sz = s[2];
			 if(tag-1==0)
			 {
				 cout << atom[tag-1].rx<<"\t"<<atom[tag-1].ry<<"\t"<<atom[tag-1].rz<<"\t"<<(ux/sx)<<"\t"<<uy/sy<<"\t"<<uz/sz<<" ++++++++\n";
			 }
			 */

			if(atom[tag-1].sx>=1) atom[tag-1].sx=atom[tag-1].sx-1;
			if(atom[tag-1].sx<0) atom[tag-1].sx=1+atom[tag-1].sx;


			if(atom[tag-1].sy>=1) atom[tag-1].sy=atom[tag-1].sy-1;
			if(atom[tag-1].sy<0) atom[tag-1].sy=1+atom[tag-1].sy;

			if(atom[tag-1].sz>=1) atom[tag-1].sz=atom[tag-1].sz-1;
			if(atom[tag-1].sz<0) atom[tag-1].sz=1+atom[tag-1].sz;


			atom[tag-1].pe = pe;
			atom[tag-1].type = type;

			if(atom[tag-1].type==1)
			{
				atom[tag-1].ma = 63.546;
				atom[tag-1].pe = atom[tag-1].pe; //subtract cohesive energy from the energy of the atom
				strncpy(atom[tag-1].elem,"Cu",sizeof(atom[i].elem));

			}else
			{
				if(atom[tag-1].type==2)
				{
					atom[tag-1].ma = 92.90638;
					atom[tag-1].pe = atom[tag-1].pe;
					strncpy(atom[tag-1].elem ,"Nb",sizeof(atom[i].elem));
				}else
				{
					atom[tag-1].ma = 63.546;
					atom[tag-1].pe = atom[tag-1].pe;
					strncpy(atom[tag-1].elem,"Z",sizeof(atom[i].elem));
				}
			}

			atom[tag-1].vx = vx;
			atom[tag-1].vy = vy;
			atom[tag-1].vz = vz;

			atom[tag-1].ux = ux;
			atom[tag-1].uy = uy;
			atom[tag-1].uz = uz;

			atom[tag-1].fx = ux;
			atom[tag-1].fy = uy;
			atom[tag-1].fz = uz;

		}

		fclose(fptr);
		zip(filename);
		/*
			double H_c_inv[3][3];
		 M3inv(Hcry, H_c_inv);

		 for(int i=0;i<n;i++)
		 {
			 double r[3],s[3];
			 r[0]=atom[i].rx;r[1]=atom[i].ry;r[2]=atom[i].rz;
			 V3mulM3(r,H_c_inv,s);
			 atom[i].sx = s[0];atom[i].sy=s[1];atom[i].sz=s[2];
			 if(atom[i].type==1)
			 {
				 atom[i].ma = 63.546;
			 }else
			 {
				 if(atom[i].type==2)
				 {
					 atom[i].ma = 92.90638;
				 }
			 }
		 }
		 */

	}

	//cout << "out of reading lammps\n";
	return(0);
}

int read_lammps_general(char *filename)
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

	}else
	{
		g_prev_success = true;
		bool read_coords = false;
		int n_counter = -1;
		int id_counter = -1;
		string *ptr_tmp_line1;
		ptr_tmp_line1 = get_next_splits(inputfile, num_val);
		delete[] ptr_tmp_line1;
//		cout << "hi there\n";
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



