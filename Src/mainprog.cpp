/*
 *  fileport.cpp
 *
 *
 *  Created by kedar on 2/27/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#include "mainprog.h"

//double a_Cu=3.615;
bool MAIN_SAVE_LAMMPS = false;
bool MAIN_SAVE_FULL = false;
bool CHAIN_MIN = false;
bool CHAIN_NEB = false;
bool SEQ_PROCESS = true;
int CHAIN_MIN_NUMBER = 50;
int CHAIN_NEB_NUMBER = 18;
int CHAIN_START_NUMBER = 1000;
char REF_STRING[512] = "";
bool SAVE_GULP = false;
char GULP_SHELL_ATOM[80]="";
char GULP_SHELL_TYPE[80]="";
double GULP_CORE_CHARGE = 0.0;
double GULP_SHELL_CHARGE = 0.0;
char STRING1[512] = "";
char STRING2[512]="";
bool GEOMETRY= false;
bool MINIMIZE_CHAIN_SEQUENCE = false;
int MINIMIZE_CHAIN_START = 1;
int MINIMIZE_CHAIN_END = 6;
bool WRITE_FORMAT_A = false;
bool WRITE_XYZ_FORMAT = false;
bool PERFORM_DISRIGISTRY = false;
bool ATOMIC_RESET = false;
bool CONVERT_VESTA = false;
bool MAKE_ILLITE = false;
bool KEEP_GHOSTS = false;
double del_z = 0.0;

// this the disrigistry type; for now 1 means 2 files comparison
// 2 means 1 inputfile and one reference vectors between 2 planes
// 3 means 1 inputfile and 1 reference vector all round

int DISRIGISTRY_TYPE =1;

int collate_volumes2()
{

	int all_cu_nb[7][19];
	for(int i =0;i<7;i++)
	{
		for(int j=0;j<19;j++)
		{
			all_cu_nb[i][j] = 0;
		}
	}
	for(int i=0;i<7;i++)
	{
//			cout <<i<<" inner i\n";
		int Cu_[5];
		int Cu__[7];
		int Nb_[7];
		if(i ==0)
		{
			Cu_[0] = 25163; Cu_[1] = 25078;Cu_[2] = 26714;Cu_[3] = 26717;Cu_[4] = 25077;
			Nb_[0]=51836;Nb_[1]=51033;;Nb_[2]=50187;;Nb_[3]=50203;;Nb_[4]=50204;;Nb_[5]=51048;;Nb_[6]=51047;
			Cu__[0]=25075;Cu__[1]=25072;Cu__[2]=25074;Cu__[3]=26713;Cu__[4]=26626;Cu__[5]=24990;Cu__[6]=25073;

		} else if(i ==1)
		{
			Cu_[0] = 26717; Cu_[1] = 26714;Cu_[2] = 26715;Cu_[3] = 26632;Cu_[4] = 26630;
			Nb_[0]=51048;Nb_[1]=51047;;Nb_[2]=50203;;Nb_[3]=49305;;Nb_[4]=50216;;Nb_[5]=51058;;Nb_[6]=50204;
			Cu__[0]=24990;Cu__[1]=25073;Cu__[2]=26713;Cu__[3]=26628;Cu__[4]=26627;Cu__[5]=26629;Cu__[6]=26626;

		}else if(i ==2)
		{
			Cu_[0] = 26632; Cu_[1] = 26715;Cu_[2] = 28268;Cu_[3] = 28271;Cu_[4] = 26631;
			Nb_[0]=51058;Nb_[1]=50204;;Nb_[2]=49305;;Nb_[3]=49319;;Nb_[4]=49320;;Nb_[5]=50217;;Nb_[6]=50216;
			Cu__[0]=26629;Cu__[1]=26626;Cu__[2]=26628;Cu__[3]=28267;Cu__[4]=28180;Cu__[5]=26545;Cu__[6]=26627;

		}else if(i ==3)
		{
			Cu_[0] = 28271; Cu_[1] = 28268;Cu_[2] = 28269;Cu_[3] = 28186;Cu_[4] = 28184;
			Nb_[0]=50217;Nb_[1]=50216;;Nb_[2]=49319;;Nb_[3]=48395;;Nb_[4]=49330;;Nb_[5]=50225;;Nb_[6]=49320;
			Cu__[0]=26545;Cu__[1]=26627;Cu__[2]=28267;Cu__[3]=28182;Cu__[4]=28181;Cu__[5]=28183;Cu__[6]=28180;


			//Nb_[0]=50217;Nb_[1]=50216;;Nb_[2]=50225;;Nb_[3]=49320;;Nb_[4]=51065;;Nb_[5]=51058;;Nb_[6]=51059;
		}else if(i ==4)
		{
			Cu_[0] = 28186; Cu_[1] = 28269;Cu_[2] = 29775;Cu_[3] = 29778;Cu_[4] = 28185;
			Nb_[0]=51065;Nb_[1]=50217;;Nb_[2]=49320;;Nb_[3]=49330;;Nb_[4]=49331;;Nb_[5]=50226;;Nb_[6]=50225;
			Cu__[0]=28183;Cu__[1]=28180;Cu__[2]=28182;Cu__[3]=29774;Cu__[4]=29692;Cu__[5]=28102;Cu__[6]=28181;


			//Nb_[0]=50225;Nb_[1]=49320;;Nb_[2]=48395;;Nb_[3]=48407;;Nb_[4]=48408;;Nb_[5]=49331;;Nb_[6]=49330;
		}else if(i ==5)
		{
			Cu_[0] = 29778; Cu_[1] = 29775;Cu_[2] = 29776;Cu_[3] = 29698;Cu_[4] = 29696;
			Nb_[0]=50226;Nb_[1]=50225;;Nb_[2]=49330;;Nb_[3]=48408;;Nb_[4]=49337;;Nb_[5]=50230;;Nb_[6]=49331;
			Cu__[0]=28102;Cu__[1]=28181;Cu__[2]=29774;Cu__[3]=29694;Cu__[4]=29693;Cu__[5]=29695;Cu__[6]=29692;


		}else if(i ==6)
		{
			Cu_[0] = 29698; Cu_[1] = 29776;Cu_[2] = 31132;Cu_[3] = 31135;Cu_[4] = 29697;
			Nb_[0]=50230;Nb_[1]=49331;;Nb_[2]=48408;;Nb_[3]=48416;;Nb_[4]=48417;;Nb_[5]=49338;;Nb_[6]=49337;
			Cu__[0]=29695;Cu__[1]=29692;Cu__[2]=29694;Cu__[3]=31131;Cu__[4]=31064;Cu__[5]=29619;Cu__[6]=29693;


		}
		for(int j=0;j<19;j++)
		{
			if(j<5){ all_cu_nb[i][j] = Cu_[j];}else if(j<12) {all_cu_nb[i][j]=Nb_[j-5];}else {all_cu_nb[i][j] = Cu__[j-12];};
		}
	}

	for(int j=1;j<7;j++)
	{
		int n_val = 57143;
		double init_cu_volume = 0.0;
		double init_nb_volume = 0.0;
		double init_cu_volume__ = 0.0;
		double init_total_volume = 0.0;
		double init_volume_array[n_val];
		cout <<"\n";

		for(int i=0;i<7;i++)
		{

		//cout <<"\n";
		//	cout <<j<<"\t"<<i<<"\t";
			FILE *fptr_what_values;
			char filename[80]="pvp.",str[20];
			sprintf(str,"%d",j);
			strcat(filename,str);
			strcat(filename,".");
			sprintf(str,"%d",i);
			strcat(filename,str);
			fptr_what_values = fopen(filename,"r");


			double volumes[n_val];
			double volume_total=0.0;
			for(int ii=0;ii<n_val;ii++)
			{
				fscanf(fptr_what_values,"%lf",&volumes[ii]);
				volume_total+=volumes[ii];
				if(i==0) init_volume_array[ii] = volumes[ii];
			}
			fclose(fptr_what_values);

			for(int k =0;k<19;k++)
			{

			//		cout << volumes[all_cu_nb[i][k]]<<"\t";
			}
			//cout <<"\n";
			double volume_cu = 0;
			double volume_nb = 0;
			double volume_diff_cu = 0;
			double volume_diff_nb = 0;
			double volume_diff_cu__ = 0;
			double volume_cu__ = 0;
			for(int k=0;k<5;k++)
			{
				volume_cu +=volumes[all_cu_nb[i][k]];
				volume_diff_cu += volumes[all_cu_nb[i][k]]-init_volume_array[all_cu_nb[i][k]];

			}
			for(int k=5;k<12;k++)
			{
				volume_nb +=volumes[all_cu_nb[i][k]];
				volume_diff_nb += volumes[all_cu_nb[i][k]]-init_volume_array[all_cu_nb[i][k]];

			}
			for(int k=12;k<19;k++)
			{
				volume_cu__ +=volumes[all_cu_nb[i][k]];
				volume_diff_cu__ += volumes[all_cu_nb[i][k]]-init_volume_array[all_cu_nb[i][k]];

			}
			if(i==0)
			{
				init_cu_volume = volume_cu;
				init_nb_volume = volume_nb;
				init_cu_volume__ = volume_cu__;
				init_total_volume = volume_total;
			}
			printf("%d %d %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf \n",j,i,volume_cu,volume_nb,volume_cu__,volume_total,volume_diff_cu, volume_diff_nb,volume_diff_cu__,volume_cu-init_cu_volume,volume_nb-init_nb_volume,volume_cu__-init_cu_volume__,volume_total-init_total_volume);

		}
	}
	return 0;
}


int collate_volumes()
{

	int all_cu_nb[7][12];
	for(int i =0;i<7;i++)
	{
		for(int j=0;j<12;j++)
		{
			all_cu_nb[i][j] = 0;
		}
	}
	for(int i=0;i<7;i++)
	{
//			cout <<i<<" inner i\n";
		int Cu_[5];
		int Nb_[7];
		if(i ==0)
		{
			Cu_[0] = 25163; Cu_[1] = 25078;Cu_[2] = 26714;Cu_[3] = 26717;Cu_[4] = 25077;
			Nb_[0]=51836;Nb_[1]=51033;;Nb_[2]=50187;;Nb_[3]=50203;;Nb_[4]=50204;;Nb_[5]=51048;;Nb_[6]=51047;
		} else if(i ==1)
		{
			Cu_[0] = 26717; Cu_[1] = 26714;Cu_[2] = 26715;Cu_[3] = 26632;Cu_[4] = 26630;
			Nb_[0]=51048;Nb_[1]=51047;;Nb_[2]=50203;;Nb_[3]=49305;;Nb_[4]=50216;;Nb_[5]=51058;;Nb_[6]=50204;
		}else if(i ==2)
		{
			Cu_[0] = 26632; Cu_[1] = 26715;Cu_[2] = 28268;Cu_[3] = 28271;Cu_[4] = 26631;
			Nb_[0]=51058;Nb_[1]=50204;;Nb_[2]=49305;;Nb_[3]=49319;;Nb_[4]=49320;;Nb_[5]=50217;;Nb_[6]=50216;
		}else if(i ==3)
		{
			Cu_[0] = 28271; Cu_[1] = 28268;Cu_[2] = 28269;Cu_[3] = 28186;Cu_[4] = 28184;
			Nb_[0]=50217;Nb_[1]=50216;;Nb_[2]=49319;;Nb_[3]=48395;;Nb_[4]=49330;;Nb_[5]=50225;;Nb_[6]=49320;

			//Nb_[0]=50217;Nb_[1]=50216;;Nb_[2]=50225;;Nb_[3]=49320;;Nb_[4]=51065;;Nb_[5]=51058;;Nb_[6]=51059;
		}else if(i ==4)
		{
			Cu_[0] = 28186; Cu_[1] = 28269;Cu_[2] = 29775;Cu_[3] = 29778;Cu_[4] = 28185;
			Nb_[0]=51065;Nb_[1]=50217;;Nb_[2]=49320;;Nb_[3]=49330;;Nb_[4]=49331;;Nb_[5]=50226;;Nb_[6]=50225;

			//Nb_[0]=50225;Nb_[1]=49320;;Nb_[2]=48395;;Nb_[3]=48407;;Nb_[4]=48408;;Nb_[5]=49331;;Nb_[6]=49330;
		}else if(i ==5)
		{
			Cu_[0] = 29778; Cu_[1] = 29775;Cu_[2] = 29776;Cu_[3] = 29698;Cu_[4] = 29696;
			Nb_[0]=50226;Nb_[1]=50225;;Nb_[2]=49330;;Nb_[3]=48408;;Nb_[4]=49337;;Nb_[5]=50230;;Nb_[6]=49331;

		}else if(i ==6)
		{
			Cu_[0] = 29698; Cu_[1] = 29776;Cu_[2] = 31132;Cu_[3] = 31135;Cu_[4] = 29697;
			Nb_[0]=50230;Nb_[1]=49331;;Nb_[2]=48408;;Nb_[3]=48416;;Nb_[4]=48417;;Nb_[5]=49338;;Nb_[6]=49337;

		}
		for(int j=0;j<12;j++)
		{
			if(j<5){ all_cu_nb[i][j] = Cu_[j];}else {all_cu_nb[i][j]=Nb_[j-5];}
		}
	}






//	cout << "entered here\n";
	for(int j=1;j<7;j++)
	{
//		cout <<j<<" outer j\n";
		int n_val = 57143;
		double init_cu_volume = 0.0;
		double init_nb_volume = 0.0;
		double init_total_volume = 0.0;
		double init_volume_array[n_val];
		cout <<"\n";

		for(int i=0;i<7;i++)
		{

			cout <<"\n";
			cout <<j<<"\t"<<i<<"\t";
			FILE *fptr_what_values;
			char filename[80]="pvp.",str[20];
			sprintf(str,"%d",j);
			strcat(filename,str);
			strcat(filename,".");
			sprintf(str,"%d",i);
			strcat(filename,str);
//			cout << "filename is \t"<<filename<<"\n";
			fptr_what_values = fopen(filename,"r");


			double volumes[n_val];
			double volume_total=0.0;
			for(int ii=0;ii<n_val;ii++)
			{
				fscanf(fptr_what_values,"%lf",&volumes[ii]);
				volume_total+=volumes[ii];
				if(i==0) init_volume_array[ii] = volumes[ii];
//				cout << volumes[ii]<<"\t"<<volume_total<<"\t"<<ii<<"\n";
			}

			fclose(fptr_what_values);
//			cout <<"finished values\n";

			for(int k =0;k<12;k++)
			{

					cout << volumes[all_cu_nb[i][k]]<<"\t";
			}
			//cout <<"\n";
			double volume_cu = 0;
			double volume_nb = 0;
			double volume_diff_cu = 0;
			double volume_diff_nb = 0;
			for(int k=0;k<5;k++)
			{
				volume_cu +=volumes[all_cu_nb[i][k]];
				volume_diff_cu += volumes[all_cu_nb[i][k]]-init_volume_array[all_cu_nb[i][k]];
			//	if(i==5)
			//	{
			//		cout << "********\t"<<all_cu_nb[i][k]<<"\t"<<volumes[all_cu_nb[i][k]]<<"\t"<<init_volume_array[all_cu_nb[i][k]]<<"\t"<<volume_cu<<"\t"<<volume_diff_cu<<"\n";
			//	}
			}
			for(int k=5;k<12;k++)
			{
				volume_nb +=volumes[all_cu_nb[i][k]];
				volume_diff_nb += volumes[all_cu_nb[i][k]]-init_volume_array[all_cu_nb[i][k]];

			}
			if(i==0)
			{
				init_cu_volume = volume_cu;
				init_nb_volume = volume_nb;
				init_total_volume = volume_total;
			}
//			printf("%d %d %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf \n",j,i,volume_cu,volume_nb,volume_total,volume_diff_cu, volume_diff_nb,volume_cu-init_cu_volume,volume_nb-init_nb_volume,volume_total-init_total_volume);

		}
	}

	return 0;

}

int find_data_for_rings2(int pp)
{
	int all_cu_nb[7][19];
		for(int i =0;i<7;i++)
		{
			for(int j=0;j<19;j++)
			{
				all_cu_nb[i][j] = 0;
			}
		}
		for(int i=0;i<7;i++)
		{
	//			cout <<i<<" inner i\n";
			int Cu_[5];
			int Cu__[7];
			int Nb_[7];
			if(i ==0)
			{
				Cu_[0] = 25163; Cu_[1] = 25078;Cu_[2] = 26714;Cu_[3] = 26717;Cu_[4] = 25077;
				Nb_[0]=51836;Nb_[1]=51033;;Nb_[2]=50187;;Nb_[3]=50203;;Nb_[4]=50204;;Nb_[5]=51048;;Nb_[6]=51047;
				Cu__[0]=25075;Cu__[1]=25072;Cu__[2]=25074;Cu__[3]=26713;Cu__[4]=26626;Cu__[5]=24990;Cu__[6]=25073;

			} else if(i ==1)
			{
				Cu_[0] = 26717; Cu_[1] = 26714;Cu_[2] = 26715;Cu_[3] = 26632;Cu_[4] = 26630;
				Nb_[0]=51048;Nb_[1]=51047;;Nb_[2]=50203;;Nb_[3]=49305;;Nb_[4]=50216;;Nb_[5]=51058;;Nb_[6]=50204;
				Cu__[0]=24990;Cu__[1]=25073;Cu__[2]=26713;Cu__[3]=26628;Cu__[4]=26627;Cu__[5]=26629;Cu__[6]=26626;

			}else if(i ==2)
			{
				Cu_[0] = 26632; Cu_[1] = 26715;Cu_[2] = 28268;Cu_[3] = 28271;Cu_[4] = 26631;
				Nb_[0]=51058;Nb_[1]=50204;;Nb_[2]=49305;;Nb_[3]=49319;;Nb_[4]=49320;;Nb_[5]=50217;;Nb_[6]=50216;
				Cu__[0]=26629;Cu__[1]=26626;Cu__[2]=26628;Cu__[3]=28267;Cu__[4]=28180;Cu__[5]=26545;Cu__[6]=26627;

			}else if(i ==3)
			{
				Cu_[0] = 28271; Cu_[1] = 28268;Cu_[2] = 28269;Cu_[3] = 28186;Cu_[4] = 28184;
				Nb_[0]=50217;Nb_[1]=50216;;Nb_[2]=49319;;Nb_[3]=48395;;Nb_[4]=49330;;Nb_[5]=50225;;Nb_[6]=49320;
				Cu__[0]=26545;Cu__[1]=26627;Cu__[2]=28267;Cu__[3]=28182;Cu__[4]=28181;Cu__[5]=28183;Cu__[6]=28180;


				//Nb_[0]=50217;Nb_[1]=50216;;Nb_[2]=50225;;Nb_[3]=49320;;Nb_[4]=51065;;Nb_[5]=51058;;Nb_[6]=51059;
			}else if(i ==4)
			{
				Cu_[0] = 28186; Cu_[1] = 28269;Cu_[2] = 29775;Cu_[3] = 29778;Cu_[4] = 28185;
				Nb_[0]=51065;Nb_[1]=50217;;Nb_[2]=49320;;Nb_[3]=49330;;Nb_[4]=49331;;Nb_[5]=50226;;Nb_[6]=50225;
				Cu__[0]=28183;Cu__[1]=28180;Cu__[2]=28182;Cu__[3]=29774;Cu__[4]=29692;Cu__[5]=28102;Cu__[6]=28181;


				//Nb_[0]=50225;Nb_[1]=49320;;Nb_[2]=48395;;Nb_[3]=48407;;Nb_[4]=48408;;Nb_[5]=49331;;Nb_[6]=49330;
			}else if(i ==5)
			{
				Cu_[0] = 29778; Cu_[1] = 29775;Cu_[2] = 29776;Cu_[3] = 29698;Cu_[4] = 29696;
				Nb_[0]=50226;Nb_[1]=50225;;Nb_[2]=49330;;Nb_[3]=48408;;Nb_[4]=49337;;Nb_[5]=50230;;Nb_[6]=49331;
				Cu__[0]=28102;Cu__[1]=28181;Cu__[2]=29774;Cu__[3]=29694;Cu__[4]=29693;Cu__[5]=29695;Cu__[6]=29692;


			}else if(i ==6)
			{
				Cu_[0] = 29698; Cu_[1] = 29776;Cu_[2] = 31132;Cu_[3] = 31135;Cu_[4] = 29697;
				Nb_[0]=50230;Nb_[1]=49331;;Nb_[2]=48408;;Nb_[3]=48416;;Nb_[4]=48417;;Nb_[5]=49338;;Nb_[6]=49337;
				Cu__[0]=29695;Cu__[1]=29692;Cu__[2]=29694;Cu__[3]=31131;Cu__[4]=31064;Cu__[5]=29619;Cu__[6]=29693;


			}
			for(int j=0;j<19;j++)
			{
				if(j<5){ all_cu_nb[i][j] = Cu_[j];}else if(j<12) {all_cu_nb[i][j]=Nb_[j-5];}else {all_cu_nb[i][j] = Cu__[j-12];};
			}
		}
		double pe_Cu = 0.0;
		double pe_Cu__ = 0.0;
		double pe_Nb = 0.0;
		double pe = 0.0;
		for(int i =0;i<19;i++)
		{
			if(i<5) {pe_Cu += atom[all_cu_nb[pp][i]].pe;atom[all_cu_nb[pp][i]].disrigistry[2] = 2;}
			if((i>4)&&(i<12)) {pe_Nb += atom[all_cu_nb[pp][i]].pe;atom[all_cu_nb[pp][i]].disrigistry[2] = 3;}
			if((i>11)) {pe_Cu__ += atom[all_cu_nb[pp][i]].pe;atom[all_cu_nb[pp][i]].disrigistry[2] = 1;}
			pe += atom[all_cu_nb[pp][i]].pe;


		}
		cout << "energy of the 19 member ring is\t";
		cout << pe_Cu<<"\t"<<pe_Nb<<"\t"<<pe_Cu__<<"\t"<<pe<<"\n";
		cout << "*******************************\n";
		/*
		FILE *fptr_patch;
		char filename[80]="color_patch.",str[20];
		if(pp==-1) pp=0;
		sprintf(str,"%06d",pp);
		strcat(filename,str);
		strcat(filename,".usr");
		fptr_patch=fopen(filename,"w");
		for(int kk=0;kk<n;kk++)
		{
			double a1=0;
			double a2=0;
			double a3=0.5;
			if(atom[kk].disrigistry[2]==2)
			{
				a1 = 0.9;a2=1.0;a3=0.0;
			}else if(atom[kk].disrigistry[2] == 3)
			{
				a1 = 1.0;a2=0.0;a3=0.0;
			}else if(atom[kk].disrigistry[2] == 1)
			{
				a1 = 0.0;a2=1.0;a3=0.0;
			}

			fprintf(fptr_patch,"%lf %lf %lf\n",a1,a2,a3);
		}
		fclose(fptr_patch);
		*/
	double sij[3],rij[3],rij_actual;
	double sum_rij;
	int counter = 0;
	double s[3];
	double Cu_r[3];
	double Nb_r[3];
	double Cu_r__[3];
	for(int i=0;i<3;i++)
	{
		s[i] = Cu_r[i] =Nb_r[i]=Cu_r__[i] = 0.0;
	}

	for(int i =0;i<5;i++)
	{
		int j = all_cu_nb[pp][i];
		s[0] +=atom[j].sx;s[1]+=atom[j].sy;s[2]+=atom[j].sz;


	}
	s[0] = s[0]/5.0;s[1]=s[1]/5; s[2] = s[2]/5;
	V3mulM3(s,Hcry,Cu_r);
	s[0] = s[1] = s[2] = 0.0;
	for(int i =5;i<12;i++)
	{
		int j = all_cu_nb[pp][i];
		s[0] +=atom[j].sx;s[1]+=atom[j].sy;s[2]+=atom[j].sz;


	}
	s[0] = s[0]/7.0;s[1]=s[1]/7; s[2] = s[2]/7;
	V3mulM3(s,Hcry,Nb_r);
	s[0] = s[1] = s[2] = 0.0;
	for(int i =12;i<19;i++)
	{
		int j = all_cu_nb[pp][i];
		s[0] +=atom[j].sx;s[1]+=atom[j].sy;s[2]+=atom[j].sz;


	}
	s[0] = s[0]/7.0;s[1]=s[1]/7; s[2] = s[2]/7;
	V3mulM3(s,Hcry,Cu_r__);
	cout <<"************************************************\n";
	for(int i =0;i<3;i++)
	{
		cout << Cu_r[i]<<"\t";
	}
	for(int i =0;i<3;i++)
	{
		cout << Nb_r[i]<<"\t";
	}
	for(int i =0;i<3;i++)
	{
		cout << Cu_r__[i]<<"\t";
	}
	cout <<"\n**************************************************\n";
	/*
	for(int i =0;i<5;i++)
	{
		for(int j=5;j<12;j++)
		{
			sij[0] = atom[all_cu_nb[pp][i]].sx - atom[all_cu_nb[pp][j]].sx;
			sij[1] = atom[all_cu_nb[pp][i]].sy - atom[all_cu_nb[pp][j]].sy;
			sij[2] = atom[all_cu_nb[pp][i]].sz - atom[all_cu_nb[pp][j]].sz;
			sij[0] = sij[0]-(int)(sij[0]*2);
			sij[1] = sij[1]-(int)(sij[1]*2);
			sij[2] = sij[2]-(int)(sij[2]*2);

			V3mulM3(sij,Hcry,rij);
			rij_actual = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
			sum_rij+=rij_actual;
			counter++;
		//	cout << pp<<"\t"<< counter<<"\t"<<rij_actual<<"\n";


		}
	}
	*/

	return 0;
}

int find_data_for_rings(int pp)
{
	int Cu_[5];
	int Nb_[7];
	if(pp ==0)
	{
		Cu_[0] = 25163; Cu_[1] = 25078;Cu_[2] = 26714;Cu_[3] = 26717;Cu_[4] = 25077;
		Nb_[0]=51836;Nb_[1]=51033;;Nb_[2]=50187;;Nb_[3]=50203;;Nb_[4]=50204;;Nb_[5]=51048;;Nb_[6]=51047;
	} else if(pp ==1)
	{
		Cu_[0] = 26717; Cu_[1] = 26714;Cu_[2] = 26715;Cu_[3] = 26632;Cu_[4] = 26630;
		Nb_[0]=51048;Nb_[1]=51047;;Nb_[2]=50203;;Nb_[3]=49305;;Nb_[4]=50216;;Nb_[5]=51058;;Nb_[6]=50204;
	}else if(pp ==2)
	{
		Cu_[0] = 26632; Cu_[1] = 26715;Cu_[2] = 28268;Cu_[3] = 28271;Cu_[4] = 26631;
		Nb_[0]=51058;Nb_[1]=50204;;Nb_[2]=49305;;Nb_[3]=49319;;Nb_[4]=49320;;Nb_[5]=50217;;Nb_[6]=50216;
	}else if(pp ==3)
	{
		Cu_[0] = 28271; Cu_[1] = 28268;Cu_[2] = 28269;Cu_[3] = 28186;Cu_[4] = 28184;
		//Nb_[0]=50217;Nb_[1]=50216;;Nb_[2]=49319;;Nb_[3]=48395;;Nb_[4]=49330;;Nb_[5]=50225;;Nb_[6]=49320;

		Nb_[0]=50217;Nb_[1]=50216;;Nb_[2]=50225;;Nb_[3]=49320;;Nb_[4]=51065;;Nb_[5]=51058;;Nb_[6]=51059;
	}else if(pp ==4)
	{
		Cu_[0] = 28186; Cu_[1] = 28269;Cu_[2] = 29775;Cu_[3] = 29778;Cu_[4] = 28185;
		//Nb_[0]=50217;Nb_[1]=51065;;Nb_[2]=50226;;Nb_[3]=49331;;Nb_[4]=49330;;Nb_[5]=49320;;Nb_[6]=50225;

		Nb_[0]=50225;Nb_[1]=49320;;Nb_[2]=48395;;Nb_[3]=48407;;Nb_[4]=48408;;Nb_[5]=49331;;Nb_[6]=49330;
	}else if(pp ==5)
	{
		Cu_[0] = 29778; Cu_[1] = 29775;Cu_[2] = 29776;Cu_[3] = 29698;Cu_[4] = 29696;
		Nb_[0]=50225;Nb_[1]=49330;;Nb_[2]=48408;;Nb_[3]=50230;;Nb_[4]=49337;;Nb_[5]=50226;;Nb_[6]=49331;

	}else if(pp ==6)
	{
		Cu_[0] = 29698; Cu_[1] = 29776;Cu_[2] = 31132;Cu_[3] = 31135;Cu_[4] = 29697;
		Nb_[0]=50230;Nb_[1]=49331;;Nb_[2]=49338;;Nb_[3]=48416;;Nb_[4]=48408;;Nb_[5]=48417;;Nb_[6]=49337;

	}

//	sij[0] = sxi - atom[j].sx;
//	sij[1] = syi - atom[j].sy;
//	sij[2] = szi - atom[j].sz;
//	sij[0] = sij[0]-(int)(sij[0]*2);
//	sij[1] = sij[1]-(int)(sij[1]*2);
//	sij[2] = sij[2]-(int)(sij[2]*2);

//	V3mulM3(sij,H1,rij);
//	rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
//	cout << i<<"\t"<<j<<"\t"<<atom[i].type<<"\t"<<atom[j].type<<"\t"<<rijsq<<"\n";

	double Cu_s[3];
	double Cu_r[3];
	double Nb_r[3];
	for(int i=0;i<3;i++)
	{
		Cu_s[i] = Cu_r[i] =Nb_r[i]=0.0;
	}

	for(int i =0;i<5;i++)
	{
		int j = Cu_[i];
		Cu_s[0] +=atom[j].sx;Cu_s[1]+=atom[j].sy;Cu_s[2]+=atom[j].sz;


	}
	Cu_s[0] = Cu_s[0]/5.0;Cu_s[1]=Cu_s[1]/5; Cu_s[2] = Cu_s[2]/5;
	V3mulM3(Cu_s,Hcry,Cu_r);
	cout <<"\n*************************************************\n";
	cout << pp<< "\t"<<Cu_r[0]<<"\t"<<Cu_r[1]<<"\t"<<Cu_r[2]<<"\t";
	for(int i=0;i<3;i++)
	{
			Cu_s[i] =0.0;
	}
	for(int i=0;i<7;i++)
	{
		int j = Nb_[i];
		Cu_s[0] +=atom[j].sx;Cu_s[1]+=atom[j].sy;Cu_s[2]+=atom[j].sz;

	}
	Cu_s[0] = Cu_s[0]/7.0;Cu_s[1]=Cu_s[1]/7.; Cu_s[2] = Cu_s[2]/7.;
	V3mulM3(Cu_s,Hcry,Nb_r);
	cout <<Nb_r[0]<<"\t"<<Nb_r[1]<<"\t"<<Nb_r[2]<<"\t";
	cout <<sqrt((Cu_r[0]-Nb_r[0])*(Cu_r[0]-Nb_r[0])+(Cu_r[1]-Nb_r[1])*(Cu_r[1]-Nb_r[1]))<<"\n";
	cout <<"\n*************************************************\n";

	return 0;
}

int temp_cout_neighbors(double H1[3][3])
{
	int array_n[5];
	int jbeg,jend,jnab;
	int *coord;
	coord = (int *) malloc((n)*sizeof(int));
	double sxi,syi,szi,sij[3],rij[3],rijsq;

	int total_close = 1000;
	int close_atoms[total_close];
	int counter_close=0;
/*
	array_n[0] = 28268;
	array_n[1] = 28269;
	array_n[2] = 28186;
	array_n[3] = 28184;
	array_n[4] = 28271;
*/
	array_n[0] = 31268;
	array_n[1] = 31268;
	array_n[2] = 31268;
	array_n[3] = 31268;
	array_n[4] = 31268;

	cout <<"\n\n\n\n\n\n";
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

				if(((i==array_n[0])||(i==array_n[1])||(i==array_n[2])||(i==array_n[3])||(i==array_n[4])||(j==array_n[0])||(j==array_n[1])||(j==array_n[2])||(j==array_n[3])||(j==array_n[4]))&&(atom[i].type!=atom[j].type))
				{
				sij[0] = sxi - atom[j].sx;
				sij[1] = syi - atom[j].sy;
				sij[2] = szi - atom[j].sz;
				sij[0] = sij[0]-(int)(sij[0]*2);
				sij[1] = sij[1]-(int)(sij[1]*2);
				sij[2] = sij[2]-(int)(sij[2]*2);

				V3mulM3(sij,H1,rij);
				rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
				cout << i<<"\t"<<j<<"\t"<<atom[i].type<<"\t"<<atom[j].type<<"\t"<<rijsq<<"\n";

				}
			}
		}
	cout <<"\n\n\n\n\n\n";
	return 0;
}


double compute_interface_energy(atomic_dat *atom_now, int n_now, double interface_pe[5][2])
{
	cout << "entered interface energy\n";
	for(int i=0;i<6;i++)
	{
		for(int j=0;j<2;j++)
		{
			interface_pe[i][j]=0.0;
		}
	}

	double ring_pe = 0;
	double ring_cn = 0;
	int local_c2=0;
	double interface_energy=0.0;

	for(int i=0;i<n;i++)
	{
		interface_energy+=atom[i].pe;
		if((atom_now[i].interface==1)&&(atom_now[i].type==1))
		{
			interface_pe[0][0] += atom_now[i].pe;
			interface_pe[0][1]++;
		}else if((atom_now[i].interface==1)&&(atom_now[i].type==2))
		{
			interface_pe[1][0] += atom_now[i].pe;
			interface_pe[1][1]++;
		}else if((atom_now[i].interface==2)&&(atom_now[i].type==1))
		{
			interface_pe[2][0] += atom_now[i].pe;
			interface_pe[2][1]++;

		}else if((atom_now[i].interface==2)&&(atom_now[i].type==2))
		{

			interface_pe[3][0] += atom_now[i].pe;
			interface_pe[3][1]++;

		}else
		{
			bool collect =false;
			int interface_collect=0;
			int atom_type = -1;
			for(int c=0;c<MAX_COORD;c++)
			{
				int j = atom_now[i].coord_id[c];
				if(j>-1)
				{
					if((atom_now[j].interface!=0)&&(atom_now[i].interface==0))
					{
						collect = true;
						interface_collect = atom_now[j].interface;
						c = MAX_COORD;
						atom_type = atom_now[j].type;
					}
				}
			}

			if((collect)&&(interface_collect!=0))
			{
				if(atom_type==2)
				{
					interface_pe[3][0]+=atom_now[i].pe;
					interface_pe[3][1]++;
				}else if(atom_type==1)
				{
					interface_pe[2][0]+=atom_now[i].pe;
					interface_pe[2][1]++;
				}
			}
		}


	}


	for(int i=0;i<4;i++)
	{
		cout <<" interface energy:\t";
		for(int j=0;j<2;j++)
		{
			cout << interface_pe[i][j]<<"\t";
		}

		cout << interface_pe[i][0]/(interface_pe[i][1]*1.0)<<"\n";
	}

	cout <<" interface energy total:\t"<< interface_pe[0][0]+interface_pe[1][0]+interface_pe[2][0]+interface_pe[3][0]<<"\n";

	cout << "RING NEIGHBORHOOD ENERGY\t"<<ring_cn<<"\t"<<ring_pe<<"\n";
	return interface_energy;
}


int read_file()
{
	if (format==1)
	{
		read_A();
	}
	return(0);
}


double compute_center_of_mass( double *x_cm,double *y_cm,double *z_cm, int atom_type)
{
	*x_cm = 0.0;*y_cm = 0.0; *z_cm=0.0;
	double tot_mass = 0.0;
	double massf=Cu_massf;
	if(atom_type==1) massf = 1.0;
	if(atom_type==2) massf = 0.0;
	for(int i=0;i<n;i++)
	{
		double r[3];
		double s[3];
		double ratio=1;
		s[0] = atom[i].sx; s[1]=atom[i].sy;s[2]=atom[i].sz;
		V3mulM3(s,Hcry,r);

		if(atom[i].type==1)
		{
			ratio = massf;

		}else
		{
			if(atom[i].type==2)
			{
				ratio = 1- massf;
			}
		}

		tot_mass = tot_mass+ratio;

		*x_cm=*x_cm+r[0]*ratio; *y_cm = *y_cm+r[1]*ratio;*z_cm = *z_cm+r[2]*ratio;
	}


	*x_cm = *x_cm/tot_mass;*y_cm=*y_cm/tot_mass;*z_cm=*z_cm/tot_mass;
	return 0.0;
}




int compute_rms_displacement()
{
	int dt_len = (total_number_of_files*1.0)/2.0;
	cout <<"\t"<<dt_len<<" IS THE THING\n";
	double rms_combined[dt_len],norm_combined[dt_len];
	double rms_combined_z[dt_len];


	for(int i=0;i<dt_len;i++)
	{
		rms_combined[i]=0.0;
		rms_combined_z[i]=0.0;
		norm_combined[i]=0.0;

	}
	cout << total_number_of_files<<"\t"<<dt_len<<"\n";
	// Test for the first atom in the "store_snapshots" variable



	FILE *fptr;

	FILE *fptr2;
	FILE *fptr3;
	FILE *fptr4;
	FILE *fptr5;
	FILE *fptr6;
	FILE *fptr7;
	FILE *fptr8;
	FILE *fptr9;
	FILE *fptr10;
	FILE *fptr11;
	FILE *fptr12;
	FILE *fptr13;
	FILE *fptr14;
	FILE *fptr15;

	fptr2 = fopen("all_rms_data","w");
	fptr3 = fopen("test_data_1","w");
	fptr4 = fopen("test_data_2","w");
	fptr5 = fopen("test_data_3","w");
	fptr6 = fopen("test_data_4","w");
	fptr7 = fopen("test_data_5","w");
	fptr8 = fopen("test_data_6","w");

	fptr9 = fopen("test_data_7","w");
	fptr10 = fopen("test_data_8","w");
	fptr11 = fopen("test_data_9","w");
	fptr12 = fopen("test_data_10","w");
	fptr13 = fopen("test_data_11","w");
	fptr14 = fopen("test_data_12","w");
	fptr15 = fopen("test_data_13","w");


	for (int atom_number = 0;atom_number<interface_atoms;atom_number++)
	{
		double rms[dt_len], rms_z[dt_len],norm[dt_len];

		for(int i=0;i<dt_len;i++)
		{
			rms[i]=0.0;
			rms_z[i]=0.0;
			norm[i]=0.0;
		}

		for(int i = 0;i<total_number_of_files;i++)
		{
			double x = store_snapshots[i*(interface_atoms*4)+(atom_number*4)+1];
			double y = store_snapshots[i*(interface_atoms*4)+(atom_number*4)+2];
			double z = store_snapshots[i*(interface_atoms*4)+(atom_number*4)+3];

			if (tag_array_interface_atoms[atom_number]==24393)
			{
				fprintf(fptr3,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==24627)
			{
				fprintf(fptr4,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==24314)
			{
				fprintf(fptr5,"%lf %lf %lf\n",x,y,z);
			}
			if (tag_array_interface_atoms[atom_number]==24047)
			{
				fprintf(fptr6,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==26242)
			{
				fprintf(fptr7,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==22947)
			{
				fprintf(fptr8,"%lf %lf %lf\n",x,y,z);
			}


			if (tag_array_interface_atoms[atom_number]==22613)
			{
				fprintf(fptr9,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==23123)
			{
				fprintf(fptr10,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==22903)
			{
				fprintf(fptr11,"%lf %lf %lf\n",x,y,z);
			}
			if (tag_array_interface_atoms[atom_number]==21743)
			{
				fprintf(fptr12,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==21568)
			{
				fprintf(fptr13,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==21884)
			{
				fprintf(fptr14,"%lf %lf %lf\n",x,y,z);
			}
			if(tag_array_interface_atoms[atom_number]==20854)
			{
				fprintf(fptr15,"%lf %lf %lf\n",x,y,z);
			}





			int final_file_number = min(total_number_of_files,i+dt_len);

			for(int j=i+1;j<final_file_number;j++)
			{
				double x1 = store_snapshots[j*(interface_atoms*4)+(atom_number*4)+1];
				double y1 = store_snapshots[j*(interface_atoms*4)+(atom_number*4)+2];
				double z1 = store_snapshots[j*(interface_atoms*4)+(atom_number*4)+3];
				int interval = j-i;
				rms[interval-1] = rms[interval-1]+(x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);
				rms_z[interval-1] = rms_z[interval-1]+(z-z1)*(z-z1);
				norm[interval-1] = norm[interval-1]+1;

			}
		}

		for(int i = 0;i<dt_len;i++)
		{
			if(norm[i]>0)
			{
				rms_combined[i] = rms_combined[i]+rms[i];
				rms_combined_z[i] = rms_combined_z[i]+rms_z[i];
				norm_combined[i] = norm_combined[i]+norm[i];
				rms[i] = rms[i]/norm[i];
				rms_z[i] = rms_z[i]/norm[i];
				fprintf(fptr2,"%d %d %d %lf %lf\n",atom_number,tag_array_interface_atoms[atom_number],(i+1),rms[i],norm[i]);
			}
		}

		if((tag_array_interface_atoms[atom_number]==24393)||(tag_array_interface_atoms[atom_number]==24627)||(tag_array_interface_atoms[atom_number]==24314)
		   ||(tag_array_interface_atoms[atom_number]==24047)||(tag_array_interface_atoms[atom_number]==26242)||(tag_array_interface_atoms[atom_number]==22947))
		{

			char filename[80]="rms.",str[80];
			sprintf(str,"%d",atom_number);
			strcat(filename,str);
			fptr = fopen(filename,"w");

			for(int i = 0;i<dt_len;i++)
			{
				if(norm[i]>0)
				{
					fprintf(fptr,"%d %lf %lf %lf\n",i,rms[i], rms_z[i],norm[i]);
				}
			}

			fclose(fptr);

		}


		if((tag_array_interface_atoms[atom_number]==22613)||(tag_array_interface_atoms[atom_number]==23123)||(tag_array_interface_atoms[atom_number]==22903)
		   ||(tag_array_interface_atoms[atom_number]==21743)||(tag_array_interface_atoms[atom_number]==21568)||(tag_array_interface_atoms[atom_number]==21884)||(tag_array_interface_atoms[atom_number]==20854))

		{

			char filename[80]="rms.",str[80];
			sprintf(str,"%d",atom_number);
			strcat(filename,str);
			fptr = fopen(filename,"w");

			for(int i = 0;i<dt_len;i++)
			{
				if(norm[i]>0)
				{
					fprintf(fptr,"%d %lf %lf %lf\n",i,rms[i], rms_z[i],norm[i]);
				}
			}

			fclose(fptr);

		}


	}
	fclose(fptr2);
	fclose(fptr3);
	fclose(fptr4);
	fclose(fptr5);
	fclose(fptr6);
	fclose(fptr7);
	fclose(fptr8);
	fclose(fptr9);
	fclose(fptr10);
	fclose(fptr11);
	fclose(fptr12);
	fclose(fptr13);
	fclose(fptr14);
	fclose(fptr15);

	char filename1[80]="combined.rms";
	fptr = fopen(filename1,"w");

	for(int i = 0;i<dt_len;i++)
	{
		if(norm_combined[i]>0)
		{

			rms_combined[i] = rms_combined[i]/norm_combined[i];
			rms_combined_z[i] = rms_combined_z[i]/norm_combined[i];

			fprintf(fptr,"%d %lf %lf %lf\n",i,rms_combined[i],rms_combined_z[i],norm_combined[i]);
		}
	}
	fclose(fptr);
	return(0);
}


// In the attempt to reduce the number of particle trajectory
// 1. atoms in the center were selected from ".structure.out" file from sz coordinate (between z_start and z_end) (the assumption is the atoms at the\
//    interface do not drastically move away from their initial positions).
// 2.  The atom list was then sorted for easier navigation through data files
int msd_interface_atoms()
{
	if (filenumber_start<0) filenumber_start = 0;
    if (filenumber_end<0) filenumber_end = 99999999;
    if(filenumber_interval<=0) filenumber_interval = 10;

	double xu_1,yu_1,zu_1,xu_2,yu_2,zu_2;
	int n_check;
	bool start = 0;
	int counter=0;

	int array_len = int((filenumber_end-filenumber_start)*1.0/(filenumber_interval*1.0)+0.5)+5;

	cout << array_len<<" array length \n";

	cout << "allocating memory\t"<<sizeof(double)*(interface_atoms*4*array_len+4)*1.0/1048576.0<<" MB\n";
	store_snapshots = (double *) malloc((interface_atoms*4*array_len+4)*sizeof(double));
	int counter_snapshots=0;
	FILE *fptr;
	fptr = fopen("all_data","w");
	FILE *fptr_cm;
	fptr_cm = fopen("center_of_mass","w");

	FILE *fptr_interface_pe;

	fptr_interface_pe = fopen("interface_pe","w");

	bool first_snapshot = true;
	double x_cm_first=0.0,y_cm_first = 0.0,z_cm_first = 0.0;
	for(int pp=filenumber_start;pp<=filenumber_end;pp+=filenumber_interval)
	{
		total_number_of_files++;
		char filename[80]="dat.",str[80];
		sprintf(str,"%d",pp);
		strcat(filename,str);
		int is_read = read_lammps(filename, atom,false,false);
		if(is_read==0)
		{
			prepare_nbrlist(Hcry,100);
			coord_number(Hcry);
			compute_CNA_and_others(atom,n, Hcry);

			if(first_snapshot)
			{
				compute_center_of_mass(&x_cm_first,&y_cm_first,&z_cm_first,0);
				first_snapshot = false;
				cout << "EENNNNNNNNNTTTTTTTEEEEEEEEERRRRRRRREEEEEEEEDDDDDDD HHHHHHHHEEEEEEEEEERRRRRRREEEEEEEEEE\n";
			}
			double x_cm_Cu,y_cm_Cu,z_cm_Cu;
			double x_cm_Nb,y_cm_Nb,z_cm_Nb;
			double x_cm_all,y_cm_all,z_cm_all;
			double x_corr=0.0,y_corr=0.0,z_corr=0.0;

			compute_center_of_mass(&x_cm_Cu,&y_cm_Cu,&z_cm_Cu,1);
			compute_center_of_mass(&x_cm_Nb,&y_cm_Nb,&z_cm_Nb,2);
			compute_center_of_mass(&x_cm_all,&y_cm_all,&z_cm_all,0);

			fprintf(fptr_cm,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",pp,x_cm_Cu,y_cm_Cu,z_cm_Cu,x_cm_Nb,y_cm_Nb,z_cm_Nb,x_cm_all,y_cm_all,z_cm_all);

			//save_cfg(pp,Hcry);
			//save_lammps(pp,Hcry);

			int counter_tmp = 0;
			double interface_pe1 = 0.0;
			double interface_pe2 = 0.0;
			for(int i = 0;i<n;i++)
			{
				if(i == tag_array_interface_atoms[counter_tmp])
				{

					if(atom[i].interface==1) interface_pe1 = interface_pe1+atom[i].pe;
					if(atom[i].interface==2) interface_pe2 = interface_pe2+atom[i].pe;

					counter_tmp++;

					store_snapshots[counter_snapshots]   = pp;
					if(atom[i].type==1)
					{
						store_snapshots[counter_snapshots+1] = atom[i].ux-x_cm_all;
						store_snapshots[counter_snapshots+2] = atom[i].uy-x_cm_all;
						store_snapshots[counter_snapshots+3] = atom[i].uz-x_cm_all;
					}else
					{
						if(atom[i].type==2)
						{
							store_snapshots[counter_snapshots+1] = atom[i].ux-x_cm_all;
							store_snapshots[counter_snapshots+2] = atom[i].uy-x_cm_all;
							store_snapshots[counter_snapshots+3] = atom[i].uz-x_cm_all;
						}
					}

					fprintf(fptr,"%d %d %d %lf %lf %lf %d\n",atom[i].interface,pp,i,store_snapshots[counter_snapshots+1],store_snapshots[counter_snapshots+2],store_snapshots[counter_snapshots+3],atom[i].CNA);

					counter_snapshots+=4;
				}
			}
			cout <<" out of looping\n";
			fprintf(fptr_interface_pe,"%d %lf %lf\n",pp,interface_pe1,interface_pe2);
		}
	}
	fclose(fptr);
	fclose(fptr_cm);
	fclose(fptr_interface_pe);

	//compute_rms_displacement();

	return(0);
}


int seq_process_lammps()
{
	if (filenumber_start<0) filenumber_start = 0;
    if (filenumber_end<0) filenumber_end = 99999999;
    if(filenumber_interval<=0) filenumber_interval = 10;

	double xu_1,yu_1,zu_1,xu_2,yu_2,zu_2;
	int n_check;
	bool start = 0;
	int counter=0;

	int array_len = int((filenumber_end-filenumber_start)*1.0/(filenumber_interval*1.0)+0.5)+5;

	cout << array_len<<" array length \n";

	cout << "allocating memory\t"<<sizeof(double)*(interface_atoms*4*array_len+4)*1.0/1048576.0<<" MB\n";
	FILE *fptr;
	FILE *fptr_cm;
	for(int pp=filenumber_start;pp<=filenumber_end;pp+=filenumber_interval)
	{
		total_number_of_files++;
		char filename[80]="dat.",str[80];
		sprintf(str,"%d",pp);
		strcat(filename,str);
		int is_read = read_lammps(filename, atom,false,false);
		if(is_read==0)
		{
			prepare_nbrlist(Hcry,100);
			coord_number(Hcry);
			compute_CNA_and_others(atom,n, Hcry);
			save_cfg(pp,Hcry);
			save_lammps(pp,Hcry);

		}
	}


	//compute_rms_displacement();

	return(0);
}

void save_additional_color(int pp, atomic_dat *atom_now, int n_now)
{
		FILE *fptr;
	char filename[80]="addcolor_spl.",str[20];
	if(pp==-1) pp=0;
	sprintf(str,"%06d",pp);
	strcat(filename,str);
	fptr=fopen(filename,"w");

	for(int i=0;i<n_now;i++)
	{
		if(atom[i].interface==0)
		{
			fprintf(fptr,"%lf %lf %lf\n",-1.0,-1.0,-1.0);
		}else
		{
			if(atom[i].type==2)
			{
				fprintf(fptr,"%lf %lf %lf\n",-1.0,-1.0,-1.0);
			}else
			{
				if(atom[i].interface==1)
				{
					if(atom[i].CNA==6)
					{
						fprintf(fptr,"%lf %lf %lf\n",0.0,0.0,0.50);
					}else if((atom[i].CNA==1)||(atom[i].CNA==2))
					{
						fprintf(fptr,"%lf %lf %lf\n",0.0,0.250,0.250);
					}else if(atom[i].CNA==8)
					{
						fprintf(fptr,"%lf %lf %lf\n",0.0,0.750,1.00);
					}else if(atom[i].CNA==9)
					{
						fprintf(fptr,"%lf %lf %lf\n",0.0,0.0,0.0);
					}else
					{
						fprintf(fptr,"%lf %lf %lf\n",1.0,1.0,1.0);
					}

				}else
				{
					bool paint = false;
					if(atom[i].CNA==2)
					{
						//fprintf(fptr,"%lf %lf %lf\n",0.75,.50,0.0);
						paint = true;
					}else
					{

						int unknown_neighbors=0;
						if(atom[i].CNA!=0)
						{

							for(int j=0;j<atom[i].coord;j++)
							{
								if(atom[i].coord_id[j]>=0)
								{
										if((atom[atom[i].coord_id[j]].CNA==2)&&(atom[atom[i].coord_id[j]].interface==2))
										{
											paint = true;
										}else if((atom[atom[i].coord_id[j]].CNA==6)&&(atom[atom[i].coord_id[j]].interface==2))
										{
											unknown_neighbors++;
										}

								}
							}

						}
						if((unknown_neighbors>2)&&(!paint)) paint = true;

					}


					if(paint)
					{
						fprintf(fptr,"%lf %lf %lf\n",0.75,.50,0.0);
						//fprintf(fptr,"%lf %lf %lf\n",-1.00,-1.0,-1.0);
					}else
					{
						fprintf(fptr,"%lf %lf %lf\n",-1.00,-1.0,-1.0);
					}

				}
			}
		}
	}

	fclose(fptr);

}



void determine_cm_rings(int pp, atomic_dat *atom_now, int n_now, double Hcry_now[3][3])
{
	FILE *fptr_ringscm;
	double cm_s[3],cm_r[3];
	static double rxfirst=0.0;
	static double ryfirst=0.0;
	static double rzfirst=0.0;
	static double cxfirst = 0.0;
	static double cyfirst=0.0;
	static double czfirst=0.0;
	static bool start_accepted = false;
	static double rxprev=0.0;
	static double ryprev=0.0;
	static double rzprev=0.0;
	static double cxprev=0.0;
	static double cyprev=0.0;
	static double czprev=0.0;
	static double rmax = 0.0;
	static double rr=0.0;
	static double rr_prev=0.0;
	double jumps=0.0;

	bool first =true;
	int count = 0;
	double first_x=0.0,first_y=0.0,first_z=0.0;
	for(int i=0;i<3;i++)
	{
		cm_s[i] =0.0;cm_r[i]=0.0;
	}

	for(int i=0;i<n_now;i++)
	{
		if((atom_now[i].CNA>8))
		{

			if(first)
			{
				first = false;
				first_x =atom_now[i].sx; first_y=atom_now[i].sy;first_z = atom_now[i].sz;
				count++;
			}else
			{
				double delnow_x = atom_now[i].sx-first_x;
				double delnow_y = atom_now[i].sy-first_y;
				double delnow_z = atom_now[i].sz-first_z;
				delnow_x = delnow_x-(int)(delnow_x*2);
				delnow_y = delnow_y-(int)(delnow_y*2);
				delnow_z = delnow_z-(int)(delnow_z*2);
				cm_s[0] +=delnow_x;
				cm_s[1] +=delnow_y;
				cm_s[2] +=delnow_z;
				count++;

			}
		}
	}

	if(count >0)
	{
		cm_s[0] = cm_s[0]/count;
		cm_s[1] = cm_s[1]/count;
		cm_s[2] = cm_s[2]/count;
		cm_s[0] += first_x;
		cm_s[1] += first_y;
		cm_s[2] += first_z;
		cm_s[0] = cm_s[0]-(int)(cm_s[0]*2);
		cm_s[1] = cm_s[1]-(int)(cm_s[1]*2);
		cm_s[2] = cm_s[2]-(int)(cm_s[2]*2);
		if(start_accepted)
		{
			cm_s[0] = cm_s[0]-(int)((cm_s[0]-cxfirst)*2);
			cm_s[1] = cm_s[1]-(int)((cm_s[1]-cyfirst)*2);
			cm_s[2] = cm_s[2]-(int)((cm_s[2]-czfirst)*2);
		}


	V3mulM3(cm_s,Hcry_now,cm_r);

	if(!start_accepted)
	{
		if((int((count+4)/5))<4)
		{
		rxfirst = cm_r[0];ryfirst = cm_r[1];rzfirst = cm_r[2];
		cxfirst = cm_s[0];cyfirst = cm_s[1];czfirst = cm_s[2];
		rxprev = cm_r[0];ryprev = cm_r[1];rzprev = cm_r[2];
		cxprev = cm_s[0];cyprev = cm_s[1];czprev = cm_s[2];
		start_accepted = true;
		}

	}else
	{
		if((int((count+4)/5))<4)
		{
			double ds[3];
			double dr[3];
			ds[0] = cm_s[0]-cxfirst;
			ds[1] = cm_s[1]-cyfirst;
			ds[2] = cm_s[2]-czfirst;
			ds[0] = ds[0]-(int)(ds[0]*2);
			ds[1] = ds[1]-(int)(ds[1]*2);
			ds[2] = ds[2]-(int)(ds[2]*2);
			V3mulM3(ds,Hcry_now,dr);
			rr= sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			rxprev = cm_r[0];ryprev = cm_r[1];rzprev = cm_r[2];
			cxprev = cm_s[0];cyprev = cm_s[1];czprev = cm_s[2];
		}


	}

	if(start_accepted)
	{
		jumps+=fabs(rr-rr_prev);
		rr_prev=rr;
		fptr_ringscm = fopen("rings_cm","a");
		fprintf(fptr_ringscm,"%d %24.16lf %24.16lf %24.16lf %24.16lf %d %lf\n",pp,cm_r[0],cm_r[1],cm_r[2],rr,(int(0.5+(int((count+4)/5)))/2),jumps);
		fclose(fptr_ringscm);
	}
	else
	{
		fptr_ringscm = fopen("rings_cm","w");
		fclose(fptr_ringscm);
	}
	}

//	char command[512]="echo ",str[80];
//	sprintf(str,"%d",(int)(rmax*10));
//	strcat(command,str);
//	strcat(command, " > ../max_movement");
//	execute_system_command(command);



}
/*
void create_jogpair()
{
	double r[3]; for(int i=0;i<3;i++) r[i]=0.0;

	//20418

}
*/
int seq_just_atomeye()
{
	cout << "entered seq process 2\n";
	//exit(1);
	if (filenumber_start<0) filenumber_start = 0;
    if (filenumber_end<0) filenumber_end = 99999999;
    if(filenumber_interval<=0) filenumber_interval = 10;

	int n_check;
	bool start = 0;
	int counter=0;
	atomic_dat *atom_ref;
	double H_ref[3][3];

	FILE *fptr;
	FILE *fptr_interface;
	fptr_interface = fopen("interface.energy","w");
	fclose(fptr_interface);
	for(int pp=filenumber_start;pp<=filenumber_end;pp+=filenumber_interval)
	{
		total_number_of_files++;
		char filename[80]="dat.",str[80];
		sprintf(str,"%d",pp);
		strcat(filename,str);
		int is_read = read_lammps(filename, atom,true,true);
		if(is_read==0)
		{
			save_cfg(pp,Hcry);
//			save_cfg_interface(pp,Hcry,1,1);
//			save_cfg_interface(pp,Hcry,1,0);
//			save_lammps(pp,Hcry);

		}
	}
	return 0;
}

void save_atom_indices(int num)
{
	FILE *fptr1;

	char filename[80]="atom_indices_new.",str[80];
	sprintf(str,"%d",num);
	strcat(filename,str);
	fptr1 = fopen(filename,"w");

	for(int i =0;i<n;i++)
	{
		/*
		if(atom[i].interface!=0)
		{
			fprintf(fptr1,"%d\n",1);
		}else
		{
			bool set = false;
			for(int j=0;j<MAX_COORD;j++)
			{
				int k = atom[i].coord_id[j];
				if((j>-1)&&(atom[k].interface !=0))
				{
					set = true;
					j = MAX_COORD;
				}

			}
			if(set)
			{
				fprintf(fptr1,"%d\n",1);
			}else
			{
				fprintf(fptr1,"%d\n",0);
			}
		}
		*/
		bool wrote = false;
		if(((atom[i].type==1)&&(atom[i].coord==9))||((atom[i].type==2)&&(atom[i].coord==10)))
		{
			bool not_surface = false;
			int a1,a2;

			for(int j=0;j<MAX_COORD;j++)
			{
				if(atom[i].coord_id[j]>-1)
				{
					if(atom[atom[i].coord_id[j]].type!=atom[i].type)
					{
						not_surface = true;
						j=MAX_COORD;
					}
				}
			}
			if(!not_surface) {wrote = true; fprintf(fptr1,"%d\n",0);}
		}
		if(!wrote)
		fprintf(fptr1,"%d\n",1);
	}
	fclose(fptr1);
}

void atomic_reset3(atomic_dat *atom_now,atomic_dat *atom_ol,double H_now[3][3])
{
	static int abcd=0;
	abcd++;
	int arr_size=0;
	for(int i=0;i<n;i++)
	{
		if((atom_now[i].interface==1)&&(atom_now[i].type==1))
		{
			arr_size++;
		}
	}
	double now_int_atoms[arr_size][6];
	double ol_int_atoms[arr_size][6];
	int arr_size_now=0,arr_size_ol=0;

	for(int i=0;i<n;i++)
	{
		if((atom_now[i].interface==1)&&(atom_now[i].type==1))
		{
			now_int_atoms[arr_size_now][0]=i;
			now_int_atoms[arr_size_now][1]=-1;
			now_int_atoms[arr_size_now][2]=atom_now[i].sx;
			now_int_atoms[arr_size_now][3]=atom_now[i].sy;
			now_int_atoms[arr_size_now][4]=atom_now[i].sz;
			now_int_atoms[arr_size_now][5]=-1;
			arr_size_now++;
		}
		if((atom_ol[i].interface==1)&&(atom_ol[i].type==1))
		{
			ol_int_atoms[arr_size_ol][0]=i;
			ol_int_atoms[arr_size_ol][1]=-1;
			ol_int_atoms[arr_size_ol][2]=atom_ol[i].sx;
			ol_int_atoms[arr_size_ol][3]=atom_ol[i].sy;
			ol_int_atoms[arr_size_ol][4]=atom_ol[i].sz;
			ol_int_atoms[arr_size_ol][5]=-1;
			arr_size_ol++;
		}
	}
	cout << arr_size<< " is the size of the array here *****\n";
	for (int i=0;i<arr_size;i++)
	{

		double r[3];
		double s[3];
		double s1[3];
		s1[0] =ol_int_atoms[i][2];s1[1] = ol_int_atoms[i][3]; s1[2]=ol_int_atoms[i][4];
		double rijsq=1e6;
		int curr_j=ol_int_atoms[i][0];
		int curr_jj=i;
		if(i%1000==0)			cout << i<<" before ";
		int p=0;
		//if(abcd==0) p =0; else p=i;
		for (int j=p;j<arr_size;j++)
		{
			s[0] = s1[0] - now_int_atoms[j][2];
			s[1] = s1[1] - now_int_atoms[j][3];
			s[2] = 0;

			s[0] = s[0]-(int)(s[0]*2);
			s[1] = s[1]-(int)(s[1]*2);
			s[2] = s[2]-(int)(s[2]*2);

			V3mulM3(s,H_now,r);
			double rijsq_tmp = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
		//	if(ol_int_atoms[i][0]==9161) cout << "aaa " << rijsq_tmp<<" "<<j<<" "<<rijsq<<" "<<curr_j<<"\n";
			if(((abcd>1)&&(now_int_atoms[j][1]!=0))||(abcd==1))
			{
				if(rijsq_tmp<rijsq)
				{
					rijsq= rijsq_tmp;
					curr_j=j;
				}
			}

		}
		ol_int_atoms[i][1]=curr_j;
		ol_int_atoms[i][5]=sqrt(rijsq);
		now_int_atoms[curr_j][1]=0;
	}
	int count_nms[arr_size];
	for(int i=0;i<arr_size;i++) count_nms[i]=0;
	for(int i=0;i<arr_size;i++)
	{
		int i_val =ol_int_atoms[i][0];
		int j_val_t = ol_int_atoms[i][1];
		int j_val = now_int_atoms[j_val_t][0];
		double dist = ol_int_atoms[i][5];
		count_nms[j_val_t]++;
		if(abcd!=0) cout << "bbb "<<dist<<"\t"<<i_val<<"\t"<<j_val<<"\n";
	}
	for(int i=0;i<arr_size;i++)
	{
		int i_val =ol_int_atoms[i][0];
		int j_val_t = ol_int_atoms[i][1];
		int j_val = now_int_atoms[j_val_t][0];
		if(count_nms[i]==0) now_int_atoms[i][1]=-1;
		if(count_nms[j_val_t]!=1) ol_int_atoms[i][1]=-1;
	}

	//tests
	for(int i=0;i<arr_size;i++) count_nms[i]=0;
	for(int i=0;i<arr_size;i++)
	{
			int i_val =ol_int_atoms[i][0];
			int j_val_t = ol_int_atoms[i][1];
			int j_val = now_int_atoms[j_val_t][0];
			double dist = ol_int_atoms[i][5];
			count_nms[j_val_t]++;
		//	cout << "bbb "<<dist<<"\t"<<i_val<<"\t"<<j_val<<"\n";
	}
	for(int i=0;i<arr_size;i++)
	{
		if(count_nms[i]!=1) cout << count_nms[i]<<" count ss is\n";
	}
	//tests - end
	int second_level_olcount=0;
	int second_level_nowcount=0;
	int here_count = 0;
	for(int i=0;i<arr_size;i++)
	{
		int i_val =ol_int_atoms[i][0];
		int j_val_t = ol_int_atoms[i][1];
		int j_val = now_int_atoms[j_val_t][0];
		double dist = ol_int_atoms[i][5];
		if(now_int_atoms[i][1]==-1)
		{
			second_level_nowcount++;
			//cout << now_int_atoms[i][0]<<" un accounted now \n";
		}
		if(ol_int_atoms[i][1]==-1)
		{
			//	cout << i_val << " un accounted old\n";
				second_level_olcount++;
		}else
		{
			swap_atom(atom_now,j_val,i_val);
			atom_now[i_val].interface=2;
			atom_now[j_val].interface=2;
			atom_ol[i_val].interface=2;
			here_count++;

		}
	}
	int abcd1=0,abcd2=0;
	for(int i=0;i<n;i++)
	{
		if((atom_now[i].interface==1)&&(atom_now[i].type==1)) abcd1++;
		if((atom_ol[i].interface==1)&&(atom_ol[i].type==1)) abcd2++;
	}
	cout << "after swapping, the interface atoms in now and old are "<<abcd1<<"\t"<<abcd2<<"\n";
	cout << "matching is "<< second_level_olcount<<"\t"<<second_level_nowcount<<" "<<here_count<<"\n";
	//if((second_level_count_check==0)&&(second_level_count>0))
	//	{
	//		cout << "*************GOING IN AGAIN************ "<< second_level_count<<" "<<second_level_count_check<<"\n";
			if(abcd<1) atomic_reset3(atom_now,atom_ol,H_now);
	//	}

}

void atomic_reset2_(atomic_dat *atom_now,atomic_dat *atom_ol,double H_now[3][3])
{
		double max_dist=-100;
		cout << "entered inside atomic reset\n";
		double **arr_list;
		arr_list = (double **) malloc((n+3)*(n+3)*sizeof(double));
		for(int i=0;i<(n+3)*(n+3);i++)
		{
			arr_list[i] = (double *) malloc(3*sizeof(double));
			for(int k=0;k<3;k++)
			{
				arr_list[i][k]=-1;
			}
		}
		int numerouno=0;
		for (int i=n-1;i>=0;i=i-1)
		{
			double r[3];
			double s[3];
			double s1[3];
			s1[0] =atom_ol[i].sx;s1[1] = atom_ol[i].sy; s1[2]=atom_ol[i].sz;
			double rijsq=1e6;
			int curr_j=i;
			if(i%1000==0)			cout << i<<" before ";
			for(int j=n-1;j>=0;j=j-1)
			{
				if((atom_now[j].type==atom_ol[i].type)&&(atom_now[j].interface==atom_ol[i].interface))
				{
					s[0] = s1[0] - atom_now[j].sx;
					s[1] = s1[1] - atom_now[j].sy;
					s[2] = s1[2] - atom_now[j].sz;

					s[0] = s[0]-(int)(s[0]*2);
					s[1] = s[1]-(int)(s[1]*2);
					s[2] = s[2]-(int)(s[2]*2);

					V3mulM3(s,H_now,r);
					double rijsq_tmp = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
					if(rijsq_tmp<rijsq)
					{
						rijsq= rijsq_tmp;
						curr_j=j;
					}
				}

			}
			arr_list[numerouno][0]=i;
			arr_list[numerouno][1]=curr_j;
			arr_list[numerouno][2] =rijsq;

			if(rijsq>max_dist) max_dist = rijsq;

			if(i%1000==0) cout <<"after and max is " <<sqrt(max_dist) <<" "<<sqrt(rijsq);

			if(i%1000==0) cout <<" swapped\n";

		}
		int *arr_list2, *arr_list3;
		arr_list2 = (int *) malloc((n+3)*sizeof(int));
		arr_list3 = (int *) malloc((n+3)*sizeof(int));
		for(int i=0;i<n;i++) { arr_list2[i]=arr_list3[i]=0;};
		for(int i =0;i<(n+3)*(n+3);i++)
		{
			if(arr_list[i][0]>-1)
			{
			arr_list2[((int) arr_list[i][0])]++;
			arr_list3[((int) arr_list[i][1])]++;
			}
		}
		for(int i =0;i<n;i++)
		{
			if(arr_list2[i]!=1)
			{
				cout << "first list "<<i<< " "<<arr_list2[i]<<"\n";
			}
			if(arr_list3[i]!=1)
			{
				cout << "first list "<<i<< " "<<arr_list3[i]<<"\n";
			}
		}

}

void atomic_reset(atomic_dat *atom_now,atomic_dat *atom_ol,double H_now[3][3])
{
	// this is for atoms where I think stuff is not happening for now everywhere but at the interface
		double max_dist=-100;
		int num_extra=0;
		cout << "entered inside atomic reset\n";
			for (int i=n-1;i>=0;i=i-1)
			{
			if(!((atom_ol[i].type==1)&&(atom_ol[i].interface==1)))
				//if(true)
				{
					double r[3];
					double s[3];
					double s1[3];
					s1[0] =atom_ol[i].sx;s1[1] = atom_ol[i].sy; s1[2]=atom_ol[i].sz;
					double rijsq=1e6;
					int curr_j=-1;
					if(i%1000==0)			cout << i<<" aaa before ";
					for(int j=i;j>=0;j=j-1)
					{
						if((atom_now[j].type==atom_ol[i].type)&&(atom_now[j].interface==atom_ol[i].interface))
						{
							if(true)//if(!((atom_now[j].type==1)&&(atom_now[j].interface==1)))
							{
								s[0] = s1[0] - atom_now[j].sx;
								s[1] = s1[1] - atom_now[j].sy;
								s[2] = s1[2] - atom_now[j].sz;

								s[0] = s[0]-(int)(s[0]*2);
								s[1] = s[1]-(int)(s[1]*2);
								s[2] = s[2]-(int)(s[2]*2);

								V3mulM3(s,H_now,r);
								double rijsq_tmp = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
								if(rijsq_tmp<rijsq)
								{
									rijsq= rijsq_tmp;
									curr_j=j;
								}
							}
						}
					}

					if(curr_j!=-1)
					{
						if(rijsq>6) num_extra++;
						if(rijsq>max_dist) max_dist = rijsq;
						if(i%1000==0) cout <<"after and max is " <<sqrt(max_dist) <<" "<<sqrt(rijsq);
						swap_atom(atom_now,curr_j,i);
						if(i%1000==0) cout <<" swapped ";
						if(i%1000==0) cout <<num_extra<<"\n";
					}

				}
			}

	//end- everywhere but at the interface

	//begin for interface
	atomic_reset3(atom_now,atom_ol,H_now);
	//atomic_reset2(atom,atom_ref,Hcry);


}

int seq_process_lammps2()
{
	cout << "entered seq process 2\n";
	//exit(1);
	if (filenumber_start<0) filenumber_start = 0;
    if (filenumber_end<0) filenumber_end = 99999999;
    if(filenumber_interval<=0) filenumber_interval = 10;

	int n_check;
	bool start = 0;
	int counter=0;
	atomic_dat *atom_ref;
	double H_ref[3][3];

	FILE *fptr;
	FILE *fptr_interface;
	fptr_interface = fopen("interface.energy","w");
	fclose(fptr_interface);
	for(int pp=filenumber_start;pp<=filenumber_end;pp+=filenumber_interval)
	{
		total_number_of_files++;
		char filename[80]="dat.",str[80];
		sprintf(str,"%d",pp);
		strcat(filename,str);
//		int is_read = read_lammps(filename, atom,true,true);
		int is_read = read_lammps_general(filename);
		if(is_read==0)
		{
			double interface_pe[6][2];
			compute_CNA_and_others(atom,n, Hcry);
			determine_cm_rings(pp,atom,n,Hcry);
		//	find_data_for_rings2(pp);
			double total_energy = compute_interface_energy(atom,n,interface_pe);

			//save_atom_indices(pp);

			double energy_difference = 0.0;
			double energy_difference_full = 0.0;
			double rdist = 0.0;
			int number_1=0;
			int number_2=0;
			//temp_cout_neighbors(Hcry);
		//	 find_data_for_rings( pp);

			if((pp-filenumber_start)!=0)
			{
				energy_difference = compute_energy_difference( atom_ref, atom,n, 0.05, false,&number_1);
				energy_difference_full = compute_energy_difference( atom_ref, atom,n, 0.0, true,&number_2);
//				rdist = compute_distance_scalar(atom_ref, atom, n, Hcry, Hcry);
			}
			fptr_interface = fopen("interface.energy","a");
			//fprintf(fptr_interface,"%d %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %12.8lf %d %d\n",pp,interface_pe[0][0],interface_pe[1][0],interface_pe[0][1],interface_pe[1][1],total_energy,energy_difference,energy_difference_full, rdist,number_1,number_2);
			fprintf(fptr_interface,"%d %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %24.16lf %12.8lf %d %d %d\n",pp,interface_pe[0][0],interface_pe[1][0],interface_pe[2][0],interface_pe[0][1],interface_pe[1][1],interface_pe[2][1],interface_pe[3][0],interface_pe[3][1],interface_pe[4][0],interface_pe[4][1],interface_pe[5][0],interface_pe[5][1],total_energy,energy_difference,energy_difference_full, rdist,number_1,number_2,n);
			fclose(fptr_interface);
			FILE *fptr_ring;
			char filename_ring[80] = "ring.";
			strcat(filename_ring,str);
			fptr_ring = fopen(filename_ring,"w");
			for(int i=0;i<n;i++)
			{
				if((atom[i].CNA>6)&&(atom[i].interface==1)&&(atom[i].type==1))
				{
					double r[3];
					double s[3];
					s[0] = atom[i].sx; s[1] = atom[i].sy;s[2]=atom[i].sz;
					s[0] = s[0]-0.3;
					if(s[0]>=1) s[0]=s[0]-1.0;
					if(s[0]<0) s[0] = s[0]+1;
					V3mulM3(s,Hcry,r);
					fprintf(fptr_ring,"%d %d %d %lf %lf %lf\n",i, atom[i].CNA,atom[i].type,r[0],r[1],r[2]);

				}


			}

			fclose(fptr_ring);

///			save_lammps(222,Hcry);


/*
			double H_t_now[3][3];
			M3mul(Hcry,KS1_trans_mas,H_t_now);
			for(int j=0;j<3;j++)
			{
				cout << H_t_now[j][0]<<"\t"<<H_t_now[j][1]<<"\t"<<H_t_now[j][2]<<"\n";
			}
//			write_A(inputfilename, pp, atom, H_t_now);
			write_A(inputfilename, pp, atom, Hcry);
*/
			if(WRITE_FORMAT_A)
			{
				write_A(inputfilename.c_str(), pp, atom, Hcry);
			}

//			write_A(inputfilename, pp, atom, Hcry);



			cout << "entered copy\n";
			cout << "TESTING COPY\n" << atom[15].sx<<"\t"<< pp<<"\t number is\t";

			if((pp-filenumber_start)==0)
			{
				for(int i=0;i<3;i++)
				{
					for(int j=0;j<3;j++)
					{
						H_ref[i][j] = Hcry[i][j];
						cout << H_ref[i][j]<<"\t";
					}
					cout <<" Reference H*****\n";
				}
				atom_ref = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

				copy_atomstruct(atom, atom_ref,n);
				compute_absolute_displacement_interface(atom,atom_ref,Hcry,H_ref);
				char mv_displace_0[80]="";
				strcat(mv_displace_0,"mv displacement_data");
				strcat(mv_displace_0," displacement_data_interface.");
				strcat(mv_displace_0,str);
				execute_system_command(mv_displace_0);
				compute_absolute_displacement(atom,atom_ref,Hcry,H_ref,false);
				char mv_displace[80]="";
				strcat(mv_displace,"mv displacement_data");
				strcat(mv_displace," displacement_data.");
				strcat(mv_displace,str);
				execute_system_command(mv_displace);
                compute_absolute_displacement(atom,atom_ref,Hcry,H_ref,false,2,2);
				char mv_displace1[80]="";
				strcat(mv_displace1,"mv displacement_data");
				strcat(mv_displace1," displacement_data_2.");
				strcat(mv_displace1,str);
				execute_system_command(mv_displace1);

				compute_absolute_displacement(atom,atom_ref,Hcry,H_ref,false,2,2);


				char mv_displace2[80]="";
				strcat(mv_displace2,"mv displacement_data");
				strcat(mv_displace2," displacement_data_3.");
				strcat(mv_displace2,str);
				execute_system_command(mv_displace2);

				if(MAIN_SAVE_FULL) save_cfg(pp,Hcry);
				//save_additional_color(pp,atom,n);
				save_cfg_interface(pp,Hcry,1,1);
//				save_cfg_interface(pp,Hcry,1,0);
//				save_cfg_test(pp,Hcry);
//				save_cfg_interface(pp,Hcry,2,1);
				if(MAIN_SAVE_LAMMPS) save_lammps(pp,Hcry);
				if(SAVE_GULP) save_gulp(pp, Hcry);
			//	save_lammps_specific("dat_lammps_all",pp,atom, Hcry);

			}else
			{
				if(ATOMIC_RESET)
				{
					atomic_reset(atom,atom_ref,Hcry);
					//atomic_reset3(atom,atom_ref,Hcry);
				}
				cout << "\n********************************************************************\n";
				double en_val = 0;
				int count_n_1 = 0;
				int count_n_2 = 0;
				double a_en=0;
				for(int i=0;i<n;i++)
				{

					//if(atom[i].type ==1) atom[i].pe = (atom[i].pe+3.54);
					//if(atom[i].type ==2) atom[i].pe = (atom[i].pe+7.47);
					//atom[i].pe = (atom[i].pe-atom_ref[i].pe);

					if((atom[i].CNA==8)||(atom[i].CNA==9))
					{
						a_en +=atom[i].pe;
					}

					if(fabs(atom[i].pe)>0.08)
					{
						en_val+=atom[i].pe;
						if(atom[i].type==1)
						{
						count_n_1++;
						}else
						{
							count_n_2++;
						}

					//	cout << atom[i].interface<< "\t"<<atom[i].type<<"\t"<<i<<"\t"<<atom[i].coord<<"\t"<<atom[i].CNA<<"\t"<<atom[i].coord-atom_ref[i].coord<<"\t"<<atom[i].pe<<"\n";
					}
					//atom[i].delr[3] = fabs(atom[i].sz-atom_ref[i].sz);
					//if(i==23287) cout << i<<"change in z direction is\t"<< atom[i].delr[3]<<" abcd cd \n";

				}
				cout << en_val<<"\t"<<count_n_1<<"\t"<<count_n_2<<"\t"<<a_en<<"\n";
				cout << "********************************************************************\n";

				compute_absolute_displacement_interface(atom,atom_ref,Hcry,H_ref);

				char mv_displace_0[80]="";
				strcat(mv_displace_0,"mv displacement_data");
				strcat(mv_displace_0," displacement_data_interface.");
				strcat(mv_displace_0,str);
				execute_system_command(mv_displace_0);




				compute_absolute_displacement(atom,atom_ref,Hcry,H_ref,false);


				char mv_displace[80]="";
				strcat(mv_displace,"mv displacement_data");
				strcat(mv_displace," displacement_data.");
				strcat(mv_displace,str);
				execute_system_command(mv_displace);

				compute_absolute_displacement(atom,atom_ref,Hcry,H_ref,false,1,2);


				char mv_displace1[80]="";
				strcat(mv_displace1,"mv displacement_data");
				strcat(mv_displace1," displacement_data_2.");
				strcat(mv_displace1,str);
				execute_system_command(mv_displace1);

				compute_absolute_displacement(atom,atom_ref,Hcry,H_ref,false,2,2);


				char mv_displace2[80]="";
				strcat(mv_displace2,"mv displacement_data");
				strcat(mv_displace2," displacement_data_3.");
				strcat(mv_displace2,str);
				execute_system_command(mv_displace2);


/*
				FILE *fptr_test;
				fptr_test = fopen("test_movement","w");
				for(int i=0;i<n;i++)
				{
					if(atom[i].CNA>6)
					{
					 double r1[3],s1[3],r2[3],s2[3];
					 s1[0] = atom_ref[i].sx;s1[1] = atom_ref[i].sy;s1[2] =atom_ref[i].sz;
					 s2[0] = atom[i].sx;s2[1] = atom[i].sy;s2[2] =atom[i].sz;
					 V3mulM3(s1,H_ref,r1);
					 V3mulM3(s2,Hcry,r2);
					 double dee = (r2[0]-r1[0])*(r2[0]-r1[0])+(r2[1]-r1[1])*(r2[1]-r1[1])+(r2[2]-r1[2])*(r2[2]-r1[2]);
					 fprintf(fptr_test,"%d %lf %lf %lf %lf %lf %lf %lf %d\n",i,r1[0],r1[1],r1[2],r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2],sqrt(dee),atom[i].CNA);
					 }

				}
				fclose(fptr_test);
				*/


				//int cc=pp;
				//for(int cc=0;cc<3;cc++)
				//{
				//	determine_mobility(atom_ref,atom,n,Hcry,Hcry,"mobility.out.data", filenumber_start,pp, true);

				if(MAIN_SAVE_FULL)	save_cfg(pp,Hcry);
				if(SAVE_GULP) save_gulp(pp, Hcry);
				//	save_additional_color(pp,atom,n);

				save_cfg_interface(pp,Hcry,1,1);
//					save_cfg_interface(pp,Hcry,1,0);
//					save_cfg_interface(pp,Hcry,2,1);
				//	save_lammps_specific("dat_lammps_all",pp,atom, Hcry);


				if(MAIN_SAVE_LAMMPS) save_lammps(pp,Hcry);
					//}
//					double closure_vector[3];
//					char *filename_atoms = "atoms1.list";
//					FindClosure (atom_ref, atom, n, H_ref, Hcry, true, filename_atoms, number_of_atoms, closure_vector);


			}


		}
	}

	return(0);
}

void delete_ghosts()
{
	cout << "**********************************\n";
	cout << "deleting ghost atoms....\n";
	cout << "**********************************\n";
	int delete_type = 11;

	std::vector<int> ghosts;
	for(int i = 0; i < n; i++)
	{
		if(atom[i].type==delete_type)
		{
			ghosts.push_back(i);
		}
	}
	delete_atoms(ghosts.data(), ghosts.size());

}

std::vector<double> find_distance(int i, int j, double H1[3][3])
{
	double sij[3], rij[3];
	sij[0] = atom[j].sx - atom[i].sx;
	sij[1] = atom[j].sy - atom[i].sy;
	sij[2] = atom[j].sz - atom[i].sz;
	sij[0] = sij[0]-(int)(sij[0]*2);
	sij[1] = sij[1]-(int)(sij[1]*2);
	sij[2] = sij[2]-(int)(sij[2]*2);

	V3mulM3(sij,H1,rij);
	std::vector<double> retval(rij, rij+3);
	return retval;

}

void print_no_neighbors()
{
//	int atomarray [] = {321, 305, 288};
//	for (int i = 0; i <sizeof(atomarray)/sizeof(int); i++ )
	double defaultvector[6][3] =
															{
																{2.613, 4.50915, 0.0},
																{-2.613, -4.50915, 0.0},
																{5.226, 0, 0},
																{-5.226, 0, 0},
																{2.613, -4.51, 0.0},
																{-2.613, 4.51, 0.0}
															};

	std::cout << "NeighborAnalyses Start\n";
	int check_id_val = 288;
	for(int i = 0; i < n; i++)
	{
		bool isempty[6] = {true, true, true, true, true, true};
		int tag = i;
		if(atom[i].type==7)
		{
			for (int j = 0; j < atom[tag].coord; j++)
			{
				int tag2 = atom[tag].coord_id[j];

				if(atom[tag2].type==7)
				{
					auto rvec = find_distance(tag, tag2, Hcry);
				//	if(tag==check_id_val) cout << "rvec is "<< tag+1 <<" "<<tag2+1<< " "<< rvec[0] << " "<< rvec[1]<< " "<<rvec[2]<<"\n";
					double angle_tmp = 180;
					int kk_lowest = -1;
					for(int kk=0;kk<6;kk++)
					{

						double angle1 = angle(defaultvector[kk], rvec.data());
						if(fabs(angle1)< angle_tmp) {
							kk_lowest = kk;
							angle_tmp = fabs(angle1);
						}
					}
					isempty[kk_lowest] = false;
				//	if(tag==check_id_val) cout << "mapped to " << defaultvector[kk_lowest][0]<<" "<<defaultvector[kk_lowest][1]<< " " << defaultvector[kk_lowest][2]<< " " << kk_lowest << " " <<angle_tmp<<"\n";
				}
			}
			for(int i = 0; i <sizeof(isempty)/sizeof(bool); i++)
			{
				if(isempty[i])
				{
					std::cout << tag+1<<" "<< defaultvector[i][0]<<" "<<defaultvector[i][1]<<
										" "<< defaultvector[i][2]<<"\n";
				}
			}
		}

	}
	std::cout << "NeighborAnalyses End\n";
}

int seq_process_lammps_new()
{
	cout << "entered seq process new\n";
	//exit(1);
	if (filenumber_start<0) filenumber_start = 0;
	if (filenumber_end<0) filenumber_end = 99999999;
	if(filenumber_interval<=0) filenumber_interval = 10;

	int n_check;
	bool start = 0;
	int counter=0;
	atomic_dat *atom_ref;
	double H_ref[3][3];

	FILE *fptr;
	FILE *fptr_interface;
	fptr_interface = fopen("interface.energy","w");
	fclose(fptr_interface);
	for(int pp=filenumber_start;pp<=filenumber_end;pp+=filenumber_interval)
	{
		total_number_of_files++;
		char filename[80]="dat.",str[80];
		sprintf(str,"%d",pp);
		strcat(filename,str);
		int is_read = read_lammps_general(filename);
		if(is_read==0)
		{

			cout <<  crystal0[0]<<" "<< crystal0[1]<<" ";
			cout<< crystal0[2]<<" "<<crystal0[3]*180/PI<<" ";
			cout << crystal0[4]*180/PI<<" "<<crystal0[5]*180/PI<<" pdbformat\n";

			double interface_pe[6][2];
			cout << Hcry[0][0]<<" "<< Hcry[0][1] << " "<<Hcry[0][2]<<" matrix\n";
			cout << Hcry[1][0]<< " "<< Hcry[1][1]<< " "<<Hcry[1][2]<<" matrix\n";
			cout << Hcry[2][0]<< " "<<Hcry[2][1]<< " "<< Hcry[2][2]<<" matrix\n";

			// even though we are not doing CNA, neighbor lists and stuff is built here
			cout << "processing file number "<< pp << "\n";
			delete_ghosts();

			compute_CNA_and_others(atom,n, Hcry, false);

			fill_bonds_etc();

			print_no_neighbors();


			double total_energy = compute_interface_energy(atom,n,interface_pe);
			cout << "total energy is: " << total_energy<<"\n";
			double energy_difference = 0.0;
			double energy_difference_full = 0.0;
			double rdist = 0.0;
			int number_1=0;
			int number_2=0;
			//temp_cout_neighbors(Hcry);
			//	 find_data_for_rings( pp);

			if((pp-filenumber_start)!=0)
			{
				energy_difference = compute_energy_difference( atom_ref, atom,n, 0.05, false,&number_1);
				energy_difference_full = compute_energy_difference( atom_ref, atom,n, 0.0, true,&number_2);
			}

			if((pp-filenumber_start)==0)
			{
				atom_ref = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
				copy_atomstruct(atom, atom_ref,n);

			}else
			{
				if(ATOMIC_RESET)
				{
					atomic_reset(atom,atom_ref,Hcry);
				}
			}

			if(MAIN_SAVE_FULL) save_cfg(pp,Hcry);
			if(MAIN_SAVE_LAMMPS) save_lammps(pp,Hcry);
			if(SAVE_GULP) save_gulp(pp, Hcry);


		}
	}

	return(0);
}

void analyze_neighbor_configs(int pp,atomic_dat *atom_now,atomic_dat *atom_ref,double H_now[3][3],double H_ref[3][3],int n_here)
{
	FILE *fptr,*fptr1,*fptr2,*fptr3;

	fptr = fopen("neigh_configs11.log","a");
	fptr1 = fopen("neigh_configs12.log","a");
	fptr2 = fopen("neigh_configs21.log","a");
	fptr3 = fopen("some_data","a");
	int tot_count[5];
	double avg[5];
	for(int i=0;i<5;i++)
	{
		tot_count[i]=0;
		avg[i]=0.0;
	}
	int max_count = 0;
	for(int i=0;i<n;i++)
	{
		atom_now[i].neigh_config = atom_ref[i].neigh_config;
		for(int k=0;k<4;k++)
		{
			atom_now[i].delr[k] = atom_ref[i].delr[k];
		}
	}
	for(int i=0;i<n;i++)
	{
		if(atom_now[i].interface!=0)
		{
			int jbeg = nbr_ptr[i];
			int jend = nbr_ptr1[i];
			int disc_count=0;
			int *coord_arr1;
			int *coord_arr2;

			double s_now[3],r_now[3],s_ref[3],r_ref[3],rijsq=0.0;
			s_now[0] = atom_now[i].sx;s_now[1]=atom_now[i].sy;s_now[2]=atom_now[i].sz;
			s_ref[0] = atom_ref[i].sx;s_ref[1]=atom_ref[i].sy;s_ref[2]=atom_ref[i].sz;
			V3mulM3(s_now,H_now,r_now);
			V3mulM3(s_ref,H_ref,r_ref);

			if(sqrt((r_now[0]-r_ref[0])*(r_now[0]-r_ref[0]))> 0.5*sqrt((H_now[0][0]*H_now[0][0]+H_now[0][1]*H_now[0][1]+H_now[0][2]*H_now[0][2])))
			{
				if(s_now[0]<0.5)
				{
					s_now[0]=1.0+s_now[0];
				}else
				{
					s_now[0]=s_now[0]-1.0;
				}
			}
			if(sqrt((r_now[1]-r_ref[1])*(r_now[1]-r_ref[1]))> 0.5*sqrt((H_now[1][0]*H_now[1][0]+H_now[1][1]*H_now[1][1]+H_now[1][2]*H_now[1][2])))
			{
				if(s_now[1]<0.5)
				{
					s_now[1]=1.0+s_now[1];
				}else
				{
					s_now[1]=s_now[1]-1.0;
				}
			}

			V3mulM3(s_now,H_now,r_now);

			rijsq = (r_now[0]-r_ref[0])*(r_now[0]-r_ref[0])+(r_now[1]-r_ref[1])*(r_now[1]-r_ref[1])+(r_now[2]-r_ref[2])*(r_now[2]-r_ref[2]);

			atom_now[i].delr[0] = atom_now[i].delr[0]+r_now[0]-r_ref[0];
			atom_now[i].delr[1] = atom_now[i].delr[1]+r_now[1]-r_ref[1];
			atom_now[i].delr[2] = atom_now[i].delr[2]+r_now[2]-r_ref[2];
			atom_now[i].delr[3] = atom_now[i].delr[3]+ rijsq;


			if(atom[i].coord>atom_ref[i].coord)
			{
				coord_arr1 = atom[i].coord_id;
				coord_arr2 = atom_ref[i].coord_id;
			}else
			{
				coord_arr2 = atom[i].coord_id;
				coord_arr1 = atom_ref[i].coord_id;
			}

			for (int jnab = 0; jnab <MAX_COORD; jnab++)
			{
				int j = coord_arr2[jnab];
				if(j>-1)
				{
					if(!check_repeat(j,coord_arr1,MAX_COORD))
					{
						disc_count++;
					}
				}
			}

			/*
			 for (int jnab = 0; jnab <MAX_COORD; jnab++)
			 {
				 int j = coord_arr1[jnab];
				 if(j>-1)
				 {
					 if(!check_repeat(j,coord_arr2,MAX_COORD))
					 {
						 disc_count++;
					 }
				 }
			 }
			 */


			if(disc_count>0) fprintf(fptr3,"%d %lf %lf\n",disc_count,sqrt(rijsq),sqrt((r_now[0]-r_ref[0])*(r_now[0]-r_ref[0])+(r_now[1]-r_ref[1])*(r_now[1]-r_ref[1])));


			if((disc_count>0)&&(rijsq>0.5))
			{
				double rrij = sqrt(atom_now[i].delr[0]*atom_now[i].delr[0]+atom_now[i].delr[1]*atom_now[i].delr[1]+atom_now[i].delr[2]*atom_now[i].delr[2]);
				double rrij_ref = sqrt(atom_ref[i].delr[0]*atom_ref[i].delr[0]+atom_ref[i].delr[1]*atom_ref[i].delr[1]+atom_ref[i].delr[2]*atom_ref[i].delr[2]);
				if(disc_count>1)
				{
					if(rrij>rrij_ref)
						atom_now[i].neigh_config = atom_now[i].neigh_config+disc_count;

					/* if(disc_count>0) {tot_count[0]++; avg[0] = avg[0]+disc_count;}
					if(disc_count>1) {tot_count[1]++; avg[1] = avg[1]+disc_count;}
					if(disc_count>2) {tot_count[2]++; avg[2] = avg[2]+disc_count;}
					if(disc_count>3) {tot_count[3]++; avg[3] = avg[3]+disc_count;}
					if(disc_count>4) {tot_count[4]++; avg[4] = avg[4]+disc_count;} */
				}

				double rij = sqrt(atom_now[i].delr[3]);

				//(rijsq<1000)&&
				if((atom[i].interface==1)&&atom[i].type==1)
					fprintf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",pp,i,atom_now[i].neigh_config,r_now[0],r_now[1],r_now[2],(r_now[0]-r_ref[0]),(r_now[1]-r_ref[1]),(r_now[2]-r_ref[2]),rij,rrij,atom_now[i].ux,atom_now[i].uy,atom_now[i].uz,(atom_now[i].ux-atom_ref[i].ux),(atom_now[i].uy-atom_ref[i].uy),(atom_now[i].uz-atom_ref[i].uz));

				if((atom[i].interface==2)&&atom[i].type==1)
					fprintf(fptr1,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",pp,i,atom_now[i].neigh_config,r_now[0],r_now[1],r_now[2],(r_now[0]-r_ref[0]),(r_now[1]-r_ref[1]),(r_now[2]-r_ref[2]),rij,rrij,atom_now[i].ux,atom_now[i].uy,atom_now[i].uz,(atom_now[i].ux-atom_ref[i].ux),(atom_now[i].uy-atom_ref[i].uy),(atom_now[i].uz-atom_ref[i].uz));

				if((atom[i].interface==1)&&atom[i].type==2)
					fprintf(fptr2,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",pp,i,atom_now[i].neigh_config,r_now[0],r_now[1],r_now[2],(r_now[0]-r_ref[0]),(r_now[1]-r_ref[1]),(r_now[2]-r_ref[2]),rij,rrij,atom_now[i].ux,atom_now[i].uy,atom_now[i].uz,(atom_now[i].ux-atom_ref[i].ux),(atom_now[i].uy-atom_ref[i].uy),(atom_now[i].uz-atom_ref[i].uz));


			}
			if(atom_now[i].neigh_config>max_count) max_count = atom_now[i].neigh_config;
		}

	}
	fclose(fptr);
	fclose(fptr1);
	fclose(fptr2);
	fclose(fptr3);
	//cout <<"TOTAL COUNT IS\t"<<tot_count<<"\n";
}

void seq_process_and_more()
{
	atomic_dat *atom_ref;
	double H_ref[3][3];
	bool initial=true;

	for(int pp=filenumber_start;pp<=filenumber_end;pp+=filenumber_interval)
	{
		total_number_of_files++;
		char filename[80]="dat.",str[80];
		sprintf(str,"%d",pp);
		strcat(filename,str);
		if(!initial)
		{
			cout << atom_ref[23310].neigh_config<<" just before ********************************++++++++++**************************************************************************\n";

		}
		int is_read = read_lammps(filename, atom,false,false);
		if(is_read==0)
		{
			prepare_nbrlist(Hcry,100);
			coord_number(Hcry);
			if(!initial)
			{
				analyze_neighbor_configs(pp,atom,atom_ref,Hcry,H_ref,n);
				//free(atom_ref);
				//atom_ref = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
			}else
			{
				atom_ref = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
				initial = false;
				cout <<"ENTERED HERERE ***************************%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n";
				for(int i=0;i<n;i++)
				{
					atom[i].delr[0]=atom[i].delr[1]=atom[i].delr[2]=atom[i].delr[3]=0.0;
				}
			}

			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					H_ref[i][j] = Hcry[i][j];
				}
			}

			cout << atom[23310].neigh_config<<"\t"<<atom[23310].interface<<"\t just before **********************************************************************************************************\n";

			copy_atomstruct(atom, atom_ref,n);

			cout << atom_ref[23310].neigh_config<<" just after **********************************************************************************************************\n";
			//	save_cfg(pp,Hcry);
			//	save_lammps(pp,Hcry);

		}
	}
	free(atom_ref);

	FILE *fptr;
	fptr = fopen("aux_neighconfigs_delr","w");
	double aux_prop[n][5];
	int counter1=0;
	int counter2 = 1000;
	for(int i=0;i<n;i++)
	{
		double delrr = sqrt(atom[i].delr[0]*atom[i].delr[0]+atom[i].delr[1]*atom[i].delr[1]+atom[i].delr[2]*atom[i].delr[2]);
		fprintf(fptr,"%d %d %lf\n",i,atom[i].neigh_config,delrr);
		aux_prop[i][0] = atom[i].neigh_config;
		aux_prop[i][1] = atom[i].delr[0];
		aux_prop[i][2] = atom[i].delr[1];
		aux_prop[i][3] = atom[i].delr[2];
		if((delrr>1.0)||(atom[i].neigh_config>8))
		{
			counter1++;
			aux_prop[i][4] = counter1;
		}else
		{
			if(atom[i].interface!=0)
			{
				counter2++;
				aux_prop[i][4] = counter2;
			}else
			{
				aux_prop[i][4]=5000;
			}
		}
	}
	fclose(fptr);



	for(int pp=filenumber_start;pp<=filenumber_end;pp+=filenumber_interval)
	{
		total_number_of_files++;
		char filename[80]="dat.",str[80];
		sprintf(str,"%d",pp);
		strcat(filename,str);

		int is_read = read_lammps(filename, atom,false,false);
		if(is_read==0)
		{
			compute_CNA_and_others(atom,n, Hcry);
			for(int i=0;i<n;i++)
			{
				atom[i].neigh_config=aux_prop[i][0];
				atom[i].delr[0] = aux_prop[i][1];
				atom[i].delr[1] = aux_prop[i][2];
				atom[i].delr[2] = aux_prop[i][3];
				atom[i].delr[3] = aux_prop[i][4];
			}
			save_cfg(pp,Hcry);
			//	save_lammps(pp,Hcry);

		}
	}




}


int trigger_MSD()
{

	read_file();
	prepare_nbrlist(Hcry,100);
	coord_number(Hcry);
	//save_cfg(12,Hcry);
	interface_atoms = 0.0;
	for(int i=0;i<n;i++)
	{
		if(atom[i].interface>0)
		{
			interface_atoms++;
		}
	}

	tag_array_interface_atoms = (int *) malloc((interface_atoms+3)*sizeof(int));
	interface_atoms = 0;

	for (int i = 0;i<n;i++)
	{
		if(atom[i].interface>0)
		{
			tag_array_interface_atoms[interface_atoms] = i;
			interface_atoms++;

		}
	}

	cout << "interface atoms are\t"<<interface_atoms<<"\n";

	msd_interface_atoms();
	return 0;
}



void temp_delete(char *filename)
{
	 int is_read = read_lammps(filename, atom,true,true);
		 if(is_read==0)
		 {

			 //compute_CNA_and_others(atom, n, Hcry);
			 //save_cfg(1,Hcry);

			 int len = 1;
			 int arr[len];
			 //	arr[0] = 24443;
			 //	arr[1] = 31541;
			 //	arr[2] = 35095;
			 //	arr[0] = 23859;
			 arr[0] = 516;
			 //arr[1] = 41754;

			 //		arr[2] = 18947;
			 delete_atoms(arr,len);
			 compute_CNA_and_others(atom, n, Hcry);
			 save_cfg(65,Hcry);
			 save_lammps(65,Hcry);
		 }
}


int seq_convert_A_lammps(char *prefix)
{
	if (filenumber_start<0) filenumber_start = 0;
    if (filenumber_end<0) filenumber_end = 99999999;
    if(filenumber_interval<=0) filenumber_interval = 10;

	int n_check;
	bool start = 0;
	int counter=0;
	atomic_dat *atom_ref;
	double H_ref[3][3];

	for(int pp=filenumber_start;pp<=filenumber_end;pp+=filenumber_interval)
	{
		char filename[80]="",str[80];
		strcat(filename,prefix);
		strcat(filename,".");
		sprintf(str,"%d",pp);
		strcat(filename,str);
		int is_read = read_A_exact(filename);
		if(is_read==0)
		{
			save_lammps(pp,Hcry);
		}
	}

	return 0;
}

int test_distances(int argc, char *argv[])
{
	// TESTING END
//	cout << "number of params is "<< argc<<"\n";
//	for (int i =0; i<argc;i++)
//	{
//		cout << argv[i]<< "\n";
//	}
//	cout << "\n";
	if(argc!=13) {cout << "not enough params; exiting\n"; exit(1);}
	//cout << argv[0]<<"\n";
	double rij[3],sij[3],htest[3][3],htestinv[3][3],testcrystal[6];
	double r1[3],r2[3],s1[3],s2[3];

	for(int i =0; i<3;i++)
		for (int j = 0; j<3; j++)
		{htest[i][j]=0.0; htestinv[i][j]=0.0;}

	testcrystal[0]=atof(argv[1]); testcrystal[1]=atof(argv[2]); testcrystal[2]=atof(argv[3]);
	testcrystal[3]=atof(argv[4])/180*PI; testcrystal[4]=atof(argv[5])/180*PI;
	testcrystal[5]=atof(argv[6])/180*PI;
	crystal_H(testcrystal,htest);
	M3inv(htest,htestinv);

//	r1[0] = 8.82604; r1[1] =13.6657; r1[2] = 40.0879;
//	r2[0] = 12.8078; r2[1] = 13.6714; r2[2] = 0.00917919;
	r1[0] = atof(argv[7]); r1[1] =atof(argv[8]); r1[2] = atof(argv[9]);
	r2[0] = atof(argv[10]); r2[1] = atof(argv[11]); r2[2] = atof(argv[12]);



	V3mulM3(r1,htestinv,s1);
	V3mulM3(r2,htestinv,s2);
	//cout << s1[0]<< " "<<s1[1]<<" "<<s1[2]<<"\n";
	//cout << s2[0]<< " "<<s2[1]<<" "<<s2[2]<<"\n";

	sij[0] = s1[0] - s2[0];
	sij[1] = s1[1] - s2[1];
	sij[2] = s1[2] - s2[2];
	sij[0] = sij[0]-(int)(sij[0]*2);
	sij[1] = sij[1]-(int)(sij[1]*2);
	sij[2] = sij[2]-(int)(sij[2]*2);

	V3mulM3(sij,htest,rij);
	double rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
	cout << sqrt(rijsq)<< " "<<rij[0]<<" "<<rij[1]<<" "<<rij[2]<< "\n";
	exit(1);


	// TESTING BEGIN
}


int main (int argc, char *argv[])
{
	//test_distances( argc,  argv);
/*
First change - input format needs to be changed
*/
	inputfilename.assign("dat_lammps.");
/*
	if(argc<3)
	{
		cout << "incomplete information\n";
		cout << "at minimum first file, last file, and increment is needed - sorry";
		cout << "Give dummy values 1 10 for last file and increment if only 1 "
		"file is being processed";
		return(0);
	}
*/
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

					tt = tt+2;

				}else
				{
					cout << "to convert, please say `filename` and then follow with name of the file\n";

					exit(1);
				}

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
		}else if(strcmp(argv[tt],"CNA_RATIO_1")==0)
		{
			CNA_RATIO_1 = atof(argv[tt+1]);
		}else if(strcmp(argv[tt],"CNA_RATIO_2")==0)
		{
			CNA_RATIO_2 = atof(argv[tt+1]);
		}else if(strcmp(argv[tt],"CNA_RATIO_12")==0)
		{
			CNA_RATIO_12 = atof(argv[tt+1]);
		}else if(strcmp(argv[tt],"a_1")==0)
		{
			a_Cu = atof(argv[tt+1]);
		}else if(strcmp(argv[tt],"a_2")==0)
		{
			a_Nb = atof(argv[tt+1]);
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
		}else if(strcmp(argv[tt],"SAVE_GULP")==0)
		{
			SAVE_GULP = true;
			if(tt+4 < argc){
			strcat(GULP_SHELL_ATOM, argv[tt+1]);
			strcat(GULP_SHELL_TYPE, argv[tt+2]);
			GULP_CORE_CHARGE = atof(argv[tt+3]);
			GULP_SHELL_CHARGE = atof(argv[tt+4]);
			tt=tt+4;
			}else
			{
				cout << "************************************************************\n";
				cout << "no extra parameters, so assuming only CORES and no SHELLS!\n";
				cout << "************************************************************\n";

			}
		}else if(strcmp(argv[tt], "XYX")==0)
		{
			WRITE_XYZ_FORMAT = true;

		}else if(strcmp(argv[tt],"WRITE_FORMAT_A")==0)
		{
			WRITE_FORMAT_A = true;

		}else if(strcmp(argv[tt],"neb_save_info")==0)
		{
			NEB_SAVE_INFO = atoi(argv[tt+1]);
		}else if(strcmp(argv[tt],"neb_climb")==0)
		{
			if(atoi(argv[tt+1])== 0)
				NEB_CLIMB = false;
			else
				NEB_CLIMB = true;
		}else if(strcmp(argv[tt],"neb_climb_after")==0)
		{
			NEB_CLIMB = true;
			NEB_CLIMB_AFTER = atoi(argv[tt+1]);
		}else if(strcmp(argv[tt],"neb_climb_tol")==0)
		{
			NEB_CLIMB = true;
			NEB_CLIMB_AFTER = atof(argv[tt+1]);

		}else if(strcmp(argv[tt],"neb_climb_mode")==0)
		{
			NEB_CLIMB = true;
			if(strcmp(argv[tt+1],"one")==0)
			{
				NEB_CLIMB_CONDITION = 1;
				tt++;
			}else if (strcmp(argv[tt+1],"multiple")==0)
			{
				NEB_CLIMB_CONDITION = 2;
				tt++;
			}

		}else if(strcmp(argv[tt],"neb_max_iter")==0)
		{

			NEB_MAX_ITER = atoi(argv[tt+1]);
		}else if(strcmp(argv[tt],"neb_force_tol")==0)
		{

			NEB_FORCE_TOL = atof(argv[tt+1]);
		}else if(strcmp(argv[tt],"neb_modify_linesearch")==0)
		{


			if(strcmp(argv[tt+1],"yes")==0)
			{
				NEB_MODIFY_LINE_SEARCH = true;
			}else if (strcmp(argv[tt+1],"no")==0)
			{
				NEB_MODIFY_LINE_SEARCH = false;
			}
		}else if(strcmp(argv[tt],"static_mode")==0)
		{
			if(strcmp(argv[tt+1],"none")== 0)
			{
				OPTOPTIONS_NEB_STATIC_MODE = NONE;
			}else if(strcmp(argv[tt+1],"energy")== 0)
			{
				OPTOPTIONS_NEB_STATIC_MODE = ENERGY;
				OPTOPTIONS_NEB_LIMIT = atof(argv[tt+2]);
			}else if(strcmp(argv[tt+1],"displacement")== 0)
			{
				OPTOPTIONS_NEB_STATIC_MODE = DISPLACEMENT;
				OPTOPTIONS_NEB_LIMIT = atof(argv[tt+2]);
			}
		}else if(strcmp(argv[tt],"neb_save_eachstep")==0)
		{
			if(strcmp(argv[tt+1],"yes")==0)
			{
				NEB_SAVE_EACHSTEP =true;
			}
		}else if(strcmp(argv[tt],"neb_dt")==0)
		{
			NEB_dt = atof(argv[tt+1]);
		}else if(strcmp(argv[tt],"neb_spring")==0)
		{
			NEB_springK = atof(argv[tt+1]);
		}else if(strcmp(argv[tt],"neb_slots")==0)
		{
			NEB_SLOT = atoi(argv[tt+1]);
		}else if(strcmp(argv[tt],"neb_backtrack")==0)
		{
			if(strcmp(argv[tt+1],"yes")==0)
			{
				NEB_BACKTRACK_M =true;
			}
		}else if (strcmp(argv[tt],"create_system")==0)
		{
			char *input1;char *input2;
			bool keep_xy = true; bool keepcms_together = true; bool keep_max = false;
			input1 = argv[tt+1];
			input2 = argv[tt+2];
			if(tt+3>argc)
			{
			if (strcmp(argv[tt+3],"false")==0) keep_xy = false;
			if (strcmp(argv[tt+4],"false")==0) keepcms_together = false;
			if (strcmp(argv[tt+5],"true")==0) keep_max = true;
			}
			cout << keep_xy<<"\t"<<keepcms_together<<"\t"<<keep_max<<"\n";
			create_system(input1, input2, keep_xy,keepcms_together,keep_max);

		}else if(strcmp(argv[tt],"CHAIN_SEQ_MIN")==0)
		{
			MINIMIZE_CHAIN_SEQUENCE = true;
			SEQ_PROCESS = false;
			if(strcmp(argv[tt+1],"seq_additional")==0)
			{
				ADDITIONAL_ATOMS = true;
				strcat(ADDITIONAL_ATOMS_FILE,argv[tt+2]);
				if(tt+3<argc)
				{
				if(strcmp(argv[tt+3],"no_rings")==0)
				{
					RINGS_NO = true;

				}
				}
			}

		}else if(strcmp(argv[tt],"CHAIN_MIN")==0)
		{
			if((tt+1)<argc)
			{
			if(atoi(argv[tt+1])!=0)
			{
				CHAIN_MIN_NUMBER = atoi(argv[tt+1]);
				if(atoi(argv[tt+2])!=0)
				{
					CHAIN_START_NUMBER = atoi(argv[tt+2]);
				}
			}
			}
			CHAIN_MIN = true;
			SEQ_PROCESS = false;
		}else if(strcmp(argv[tt],"CHAIN_NEB")==0)
		{
			int k = tt+1;
			if((k)<argc)
			{
			if(atoi(argv[tt+1])!=0)
			{
				CHAIN_NEB_NUMBER = atoi(argv[k]);
				if(atoi(argv[k+1])!=0)
				{
					CHAIN_START_NUMBER = atoi(argv[k+1]);
				}
			}
			}
			CHAIN_NEB = true;
			SEQ_PROCESS = false;
		}else if(strcmp(argv[tt],"SEQ")==0)
		{
			if((strcmp(argv[tt+1],"NO")==0)||(strcmp(argv[tt+1],"No")==0)||(strcmp(argv[tt+1],"no")==0))
			{
				SEQ_PROCESS = false;
			}
		}else if(strcmp(argv[tt],"GEOMETRY")==0)
		{
			if(tt+3<argc)
			{
				strcat(REF_STRING,argv[tt+1]);
				strcat(STRING1,argv[tt+2]);
				strcat(STRING2,argv[tt+3]);
			}
			GEOMETRY = true;
			cout << REF_STRING<<"\n";
			cout << STRING1<<"\n";
			cout << STRING2<<"\n";
		}else if(strcmp(argv[tt],"DISRIGISTRY")==0)
		{
			int s_atom_type = 2;
			int s_atom_interface=1;
			if(tt+3<argc)
			{
				try{
				DISRIGISTRY_TYPE=atoi(argv[tt+3]);
				s_atom_type = atoi(argv[tt+4]);
				s_atom_interface=atoi(argv[tt+5]);
				}catch(int ExNum){
				DISRIGISTRY_TYPE =1;
				}
			}
			SEQ_PROCESS=false;
			perform_disrigistry(argv[tt+1],argv[tt+2],DISRIGISTRY_TYPE,2,s_atom_type,s_atom_interface);
		}else if(strcmp(argv[tt],"ATOMIC_RESET")==0)
		{
			ATOMIC_RESET = true;
			cout << "entered in ATOMIC RESET\n";
		}else if(strcmp(argv[tt],"CUTOFF_FILE")==0)
		{
			if(tt+1<argc)
			{
				set_neighbordistances(argv[tt+1]);
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
	cout << "exited options\n";


	if(CONVERT_VESTA)
	{
		//COMMAND TO RUN
		//lbnl_processor_exec.out convert_VESTA filename small_3x_structure-Si2Al.xyz [make_illite] [keep_ghosts] CUTOFF_FILE cutoff_file.illite SAVE_LAMMPS CHARGE MOLECULE

		read_xyz_VESTA(inputfilename.c_str());
		set_types_to_atom_from_element_name();
		//prepare_nbrlist(Hcry, 100);
		//coord_number(Hcry);

		//delete_duplicates();

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

			/*
			if(del_z != 0)
			{
				cout << "we are going to change coordinates\n";
				cout << "we will keep all but K and Gh at the sam site\n";
				cout << "we will move K and Gh half del_z\n";

				double crystal2[6];
				for (int j=0; j<6;j++)
				{
					crystal2[j] = crystal0[j];
				}
				crystal2[2] = crystal2[2]+ del_z;
				double H_new[3][3], H_new_inv[3][3];

				crystal_H(crystal2,H_new);
				M3inv(H_new,H_new_inv);

				for (int i = 0; i < n; i++)
				{
					double r[3],s[3];

					if(atom[i].sx>1.0) atom[i].sx = atom[i].sx-1;
					if(atom[i].sy>1.0) atom[i].sy = atom[i].sy-1;
					if(atom[i].sz>1.0) atom[i].sz = atom[i].sz-1;

					if(atom[i].sx<0.0) atom[i].sx = 1.0+atom[i].sx;
					if(atom[i].sy<0.0) atom[i].sy = 1.0+atom[i].sy;
					if(atom[i].sz<0.0) atom[i].sz = 1.0+atom[i].sz;


					s[0] = atom[i].sx;s[1] = atom[i].sy;s[2] = atom[i].sz;
					V3mulM3(s,Hcry,r);

					if( (atom[i].type == 7) || (atom[i].type == 11))
					{

					}else
					{
						V3mulM3(r,H_new_inv,s);
					}
				}
			}
			*/

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
/*
				double rij[3], sij[3];
				if(atom[tag].sx>1.0) atom[tag].sx = atom[tag].sx-1;
				if(atom[tag].sy>1.0) atom[tag].sy = atom[tag].sy-1;
				if(atom[tag].sz>1.0) atom[tag].sz = atom[tag].sz-1;

				if(atom[tag].sx<0.0) atom[tag].sx = 1.0+atom[tag].sx;
				if(atom[tag].sy<0.0) atom[tag].sy = 1.0+atom[tag].sy;
				if(atom[tag].sz<0.0) atom[tag].sz = 1.0+atom[tag].sz;

				sij[0] = atom[tag].sx;
				sij[1] = atom[tag].sy;
				sij[2] = atom[tag].sz;
				V3mulM3(sij,Hcry,rij);
				cout << tag << " "<< rij[0]<<" "<<rij[1]<<" "<<rij[2]<<";\n";
*/
			}
		}
		cout << counter_K <<" the number of K+ atoms currently available\n";
/*
		int arry [3] = {291, 305, 309};
		for(int i =0; i<3; i++)
		{
			int j = arry[i];
			double rij[3], sij[3];

			for(int k =0;k < atom[j].coord; k++ )
			{
				int l = atom[j].coord_id[k];
				sij[0] = atom[l].sx - atom[j].sx;
				sij[1] = atom[l].sy - atom[j].sy;
				sij[2] = atom[l].sz - atom[j].sz;
				sij[0] = sij[0]-(int)(sij[0]*2);
				sij[1] = sij[1]-(int)(sij[1]*2);
				sij[2] = sij[2]-(int)(sij[2]*2);

				V3mulM3(sij,Hcry,rij);
				double rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
				cout << j << "\t"<< rij[0]<<"\t"<<rij[1]<<"\t"<<rij[2]<<" ***\n";
			}
		}
*/

		//do some `filling` of bonds -- add bonds
		fill_bonds_etc();

		// compute_neighbor_info(atom, n, Hcry);
		// default to no zip
//		save_cfg(20,Hcry);
	 if(MAIN_SAVE_LAMMPS)
		save_lammps(20,Hcry);

		save_xyz_VESTA(20,Hcry);

		//NOTES FOR DOING MORE STUFF
		//Angles between Si-Oa(6)-H and Al(3, 4)-Oa(6)-H
		//Bonds between O-H

	}

	if(GEOMETRY)
	{
		load_reference_vectors(REF_STRING);
		init_geometry(STRING1,STRING2);
	}

	if(MINIMIZE_CHAIN_SEQUENCE)
	{
		if (filenumber_start<0) filenumber_start = 0;
		if (filenumber_end<0) filenumber_end = 99999999;
		if (filenumber_interval<=0) filenumber_interval = 10;
		//create_recurring_chain(filenumber_start, filenumber_end, filenumber_interval, 1000, 30);

		int end_number;
		int start_number=format;
		int start_number_main = start_number;
		for(int ijk=MINIMIZE_CHAIN_START;ijk<MINIMIZE_CHAIN_END;ijk++)
		{


		//	reactant_product_config_changes("dat.48", "dat.50", start_number,false,true, true, true, ijk);
		//	end_number=start_number+1;
			cout << "entered prepare for chain\n";
			prepare_for_chain( filenumber_start, filenumber_end, filenumber_interval, start_number, &end_number, true, "minimize_chain.bash",true, true,ijk);
			int filenumber_start1 = start_number;
			int filenumber_end1 = end_number;
			int filenumber_interval1 = 1;
			cout << ijk <<"\t"<<"done doing whatever it is\n";
		//	create_recurring_chain(filenumber_start1, filenumber_end1, filenumber_interval1, ijk*1000, 3);
			start_number = end_number+1;
		}
	}

	if(SEQ_PROCESS)
	{
		std::cout << "Enetered into SEQ Process\n\n";
		seq_process_lammps_new();
	}


	if(CHAIN_MIN)
	{
		create_recurring_chain(filenumber_start, filenumber_end, filenumber_interval, CHAIN_START_NUMBER, CHAIN_MIN_NUMBER);
	}
	if(CHAIN_NEB)
	{
		create_recurring_chain_equidistant(filenumber_start, filenumber_end, filenumber_interval, CHAIN_START_NUMBER, CHAIN_NEB_NUMBER);
	}


	return (0);
}
