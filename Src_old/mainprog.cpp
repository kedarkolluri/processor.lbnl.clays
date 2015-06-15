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
bool PERFORM_DISRIGISTRY = false;
bool ATOMIC_RESET = false;
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
				write_A(inputfilename, pp, atom, Hcry);
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



int main (int argc, char *argv[])
{

	/* BEGIN - Tests for some crystal orientation stuff */
	double Htest[3][3];
	double Htest2[3][3];
	double crystaltest[6];
	double crystaltest2[6];
	crystaltest2[0] = 5.226; crystaltest2[1] = 9.0183; crystaltest2[2]=20.413;
	crystaltest2[3] = PI/2; crystaltest2[4] = 95.665/180*PI;crystaltest2[5] = PI/2;
	Htest[0][0] = 109.7; Htest[0][1]=0.0; Htest[0][2] = 0.0;
	Htest[1][0] = 0.0; Htest[1][1]=47.0034; Htest[1][2] = 0.0;
	Htest[2][0] = -4.639794; Htest[2][1]= 0.0; Htest[2][2] = 84.5;
	H_crystal(Htest,crystaltest);
	crystal_H(crystaltest, Htest2);

	crystaltest2[0] = 5.226; crystaltest2[1] = 9.0183; crystaltest2[2]=20.143;
	crystaltest2[3] = PI/2; crystaltest2[4] = 95.665/180*PI;crystaltest2[5] = PI/2;
	crystal_H(crystaltest2, Htest);
	H_crystal(Htest,crystaltest);
	double s_now[3];
	double r_now[3];
	s_now[0] = 0.4080; s_now[1]=0.5671; s_now[2]=0.0454;
	V3mulM3(s_now,Htest,r_now);
	cout << "abcd efgh\n";
	cout << r_now[0]<<"\t"<<r_now[1]<<"\t"<<r_now[2]<<"\n";
	r_now[2] = r_now[2]+1.0;
	double Htest_inv[3][3];
	double s_new[3];
	M3inv(Htest,Htest_inv);
	V3mulM3(r_now, Htest_inv, s_new);
	cout << "second stuff\n";
	cout << s_new[0]<<"\t"<<s_new[1]<<"\t"<<s_new[2]<<"\n";
	exit(1);
	/* END - Tests for some crystall orientation stuff */

	if(argc<5)

	{
		cout << "incomplete information\n";
		cout << "Format is : input_file_name input_file_format output_filename output_file_format( 0 for out, 1 for lammps, 2 for cfg)\n";
		return(0);
	}

	bool cutoff_file = false;
	format = atoi(argv[2]);
	target_format = atoi(argv[4]);

	inputfilename = argv[1];
	outputfilename = argv[3];
	cout << format<<"\t"<<target_format<<"\n";
	cout << inputfilename<<"\t"<< outputfilename<<"\n";

	filenumber_start = atoi(argv[5]);
	filenumber_end = atoi(argv[6]);
	filenumber_interval = atoi(argv[7]);
	//temp_delete("dat.0");
	for(int tt = 0;tt < argc;tt++)
	{
		if(strcmp(argv[tt],"zip_unzip")==0)
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
			}else
			{
				cout << "************************************************************\n";
				cout << "no extra parameters, so assuming only CORES and no SHELLS!\n";
				cout << "************************************************************\n";

			}

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
			}else if (strcmp(argv[tt+1],"multiple")==0)
			{
				NEB_CLIMB_CONDITION = 2;
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
			}else
			{
				cout << "input error; file containing cutoffs is not there... exiting";
				exit(1);
			}
		}


	}


//	int is_read = read_lammps_general("dat.0");
//	cout << "is read\t "<< is_read <<"\n";
//	compute_CNA_and_others(atom,n, Hcry);
//
//    save_cfg(0,Hcry);
//
//	//collate_volumes2();
//	exit(1);

// This is for KS

	Cu_alpha_matrix[0][0] = (3.0-sqrt(6.0))*(a_Cu+a_Nb)/a_Cu;
	Cu_alpha_matrix[0][1] = (sqrt(3.0)-sqrt(2.0))*(a_Cu-a_Nb)/a_Cu;
	Cu_alpha_matrix[0][2] = 0.0;
	Cu_alpha_matrix[1][0] = (6-2*sqrt(6.0))/sqrt(2)+3*(2-sqrt(6.0))/sqrt(2.0)*a_Nb/a_Cu;
	Cu_alpha_matrix[1][1] = (sqrt(6.0)-2)-(sqrt(6.0)-3)*a_Nb/a_Cu;
	Cu_alpha_matrix[1][2] = 0.0;
	Cu_alpha_matrix[2][0] = 0.0;
	Cu_alpha_matrix[2][1] = 0.0;
	Cu_alpha_matrix[2][2] = 1.0;

	Nb_alpha_matrix[0][0] = (sqrt(6)-2)*(a_Cu+a_Nb)/a_Nb;
	Nb_alpha_matrix[0][1] = (sqrt(3.0)-sqrt(2.0))*(a_Cu-2*a_Nb)/(2*a_Nb);
	Nb_alpha_matrix[0][2] = 0.0;
	Nb_alpha_matrix[1][0] = (6-3*sqrt(6.0))/sqrt(3)+2*(3-sqrt(6.0))/sqrt(3.0)*a_Cu/a_Nb;
	Nb_alpha_matrix[1][1] = (3-sqrt(6.0))+(0.5*sqrt(6.0)-1)*a_Cu/a_Nb;
	Nb_alpha_matrix[1][2] = 0.0;
	Nb_alpha_matrix[2][0] = 0.0;
	Nb_alpha_matrix[2][1] = 0.0;
	Nb_alpha_matrix[2][2] = 1.0;



	// This is for NW

		Cu_alpha_matrix[0][0] = 1.0;//0.977138138945160;
		Cu_alpha_matrix[0][1] = 0.0;//0.013199301634188;
		Cu_alpha_matrix[0][2] = 0.0;
		Cu_alpha_matrix[1][0] = 0;
		Cu_alpha_matrix[1][1] =  1.054339014407257;
		Cu_alpha_matrix[1][2] = 0.0;
		Cu_alpha_matrix[2][0] = 0.0;
		Cu_alpha_matrix[2][1] = 0.0;
		Cu_alpha_matrix[2][2] = 1.0;

		Nb_alpha_matrix[0][0] = 0.756711018600919;
		Nb_alpha_matrix[0][1] = 1.082670015191280;
		Nb_alpha_matrix[0][2] = 0.0;
		Nb_alpha_matrix[1][0] = 0;
		Nb_alpha_matrix[1][1] = 1;
		Nb_alpha_matrix[1][2] = 0.0;
		Nb_alpha_matrix[2][0] = 0.0;
		Nb_alpha_matrix[2][1] = 0.0;
		Nb_alpha_matrix[2][2] = 1.0;


		if(GEOMETRY)
		{
			load_reference_vectors(REF_STRING);
			init_geometry(STRING1,STRING2);
		}

	//transformation matrix for self created KS1


	KS1_trans_mas[0][0]=9.9999869941852970e-01;		KS1_trans_mas[0][1]=-1.5839654368795880e-03;	KS1_trans_mas[0][2]=3.0367409767026839e-04;
	KS1_trans_mas[1][0]=1.5841996408845647e-03;		KS1_trans_mas[1][1]=9.9999844461906551e-01;		KS1_trans_mas[1][2]=-7.7528462723164379e-04;
	KS1_trans_mas[2][0]=-3.0244560326537980e-04;	KS1_trans_mas[2][1]=7.7576469762323703e-04;	    KS1_trans_mas[2][2]=9.9999965335798380e-01;


	/*
	// TRANSFORMATION MATRIX FOR KS1 GIVEN BY MICHAEL - BEGIN
	double t_H[3][3],H_new_new[3][3], H_new_new_inv[3][3], crystal[6];

	 t_H[0][0]=0.9712353759E+02; t_H[0][1]=-.1538852993E+00; t_H[0][2]=-.7539158004E-01;
	 t_H[1][0]=-.1538852193E+00; t_H[1][1] = 0.9715180368E+02; t_H[1][2] = 0.2950252469E-01;
	 t_H[2][0] = -.7539158004E-01*2; t_H[2][1]= 0.2950252469E-01*2; t_H[2][2]= 0.1742400233E+03*2;

	 cout << " ACHTUNG - THIS KS1_TRANS_MAS IS FOR THE ONE PROVIDED BY MICHAEL AND IS NOT VALID FOR THE STRUCTURE MADE BY KEDAR !!!\n";
	 cout << " ACHTUNG - THIS KS1_TRANS_MAS IS FOR THE ONE PROVIDED BY MICHAEL AND IS NOT VALID FOR THE STRUCTURE MADE BY KEDAR !!!\n";
	 cout << " ACHTUNG - THIS KS1_TRANS_MAS IS FOR THE ONE PROVIDED BY MICHAEL AND IS NOT VALID FOR THE STRUCTURE MADE BY KEDAR !!!\n";
	 cout << " ACHTUNG - THIS KS1_TRANS_MAS IS FOR THE ONE PROVIDED BY MICHAEL AND IS NOT VALID FOR THE STRUCTURE MADE BY KEDAR !!!\n";

	 H_crystal(t_H,crystal);
	 crystal[2] = crystal[2];
	 for(int i=0;i<6;i++)
	 {
		cout << crystal[i]<<"\t";
	 }
	 cout <<"\n";
	 crystal_H(crystal,H_new_new);
	 for(int i=0;i<3;i++)
	 {
		for(int j=0;j<3;j++)
		{
			cout << H_new_new[i][j]<<"\t";
		}
		cout <<"\n";
	 }

	 double t_H_inv[3][3];
	 M3inv(H_new_new,H_new_new_inv);
	// M3mul(t_H_inv,H_new_new,KS1_trans_mas);
	  M3mul(H_new_new_inv,t_H,KS1_trans_mas);
	  for(int i=0;i<3;i++)
	 {
		for(int j=0;j<3;j++)
		{
			cout <<KS1_trans_mas[i][j]<<" HURRT\t";
		}
		cout <<"\n";
	 }
	 // TRANSFORMATION MATRIX FOR KS1 GIVEN BY MICHAEL - END
	 */

	//	9.9999869941852970e-01		-1.5839654368795880e-03	3.0367409767026839e-04
	//	1.5841996408845647e-03		9.9999844461906551e-01		-7.7528462723164379e-04
	//	-3.0244560326537980e-04	7.7576469762323703e-04	    9.9999965335798380e-01



	cout << "output Cu alpha matrix\n";
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout << Cu_alpha_matrix[i][j]<<"\t";
		}
		cout <<"\n";
	}
	cout << "output Cu alpha matrix - End\n";

	cout << "output Nb alpha matrix\n";
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout << Nb_alpha_matrix[i][j]<<"\t";
		}
		cout <<"\n";
	}
	cout << "output Nb alpha matrix - End\n";





	/*
	 // Seq processing for identifying mobility - begin
	 read_file();
	 seq_process_and_more();

	 // Seq processing for identifying mobility - end
	 */

/*
	 // Delete Atoms for KS1 -Begin
	 int is_read = read_lammps("dat.0", atom,true,true);
	 if(is_read==0)
	 {

		 //compute_CNA_and_others(atom, n, Hcry);
		 //save_cfg(1,Hcry);

		 int len = 2;
		 int arr[len];
		 //	arr[0] = 24443;
		 //	arr[1] = 31541;
		 //	arr[2] = 35095;
		 //	arr[0] = 23859;
		 arr[0] = 19357;
		 arr[1] = 19358;
		 //		arr[2] = 18947;
		 delete_atoms(arr,len);
		 compute_CNA_and_others(atom, n, Hcry);
		 save_cfg(65,Hcry);
		 save_lammps(65,Hcry);
	 }

	// Delete Atoms for KS1 -End
*/
//	cout << "entered here to read Blas' files";
//	for (int i=1;i<4;i++)
//	  {
//	    char filename[80]="mgo.desired_format.",str[80];
//	    char file_dd[80]=".dat";
//	    sprintf(str,"%d",i);
//	    strcat(filename,str);
//	    //strcat(filename,file_dd);
//	    read_xian_ming_blas2(filename);
//        save_cfg(i,Hcry);
//        save_lammps(i,Hcry);
//	  }
//exit(1);


/*
	// Fill deficient atoms
//	read_A_exact("ting.100");
//	int is_read=0;
	int is_read = read_lammps("dat.71", atom,true,true);
	if(is_read==0)
	{
		double arr_shift[3]; arr_shift[0] = 0.0;arr_shift[1]=0.4;arr_shift[2]=0.3;
//		shift(atom,n,Hcry,arr_shift);
		compute_CNA_and_others(atom,n, Hcry,true);
		save_cfg(72,Hcry);
		save_lammps(72,Hcry);
		insert_atoms_recur(10, atom, n, Hcry,false,true);
		compute_CNA_and_others(atom,n, Hcry,false);
		int remove = 0;
	//	delete_atom_single(785);
	//	insert_atoms_recur(10, atom, n, Hcry,false);
		compute_CNA_and_others(atom,n, Hcry,true);
		save_cfg(73,Hcry);
		save_lammps(73,Hcry);

		int adv_remove = 0;

		for(int i =0; i <n;i++)
		{
			if((atom[i].coord==11)||(atom[i].coord >13))
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
				if(((deficit>3)&&(hcp_n==0))||(deficit>=5)) adv_remove++;
			}
		}
		int adv_remove_arr[adv_remove];
		adv_remove = 0;
		for(int i =0; i <n;i++)
		{
			if((atom[i].coord==11)||(atom[i].coord >13))
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
				if(((deficit>3)&&(hcp_n==0))||(deficit>=5)) {adv_remove_arr[adv_remove] = i;adv_remove++;}
			}
		}

		delete_atoms(adv_remove_arr, adv_remove);

//		for(int i =0;i<n;i++)
//		{
//			if(atom[i].coord>=14) remove++;
//		}
//		int ting[remove];
//		remove =0;
//		for(int i=0;i<n;i++)
//		{
//			if((atom[i].coord>=14)||(atom[i].coord<=10))
//				{
//					ting[remove] = i;
//					remove++;
//				}
//		}
//		delete_atoms(ting, remove);

		cout << "before second\n";
		save_cfg(74,Hcry);
		save_lammps(74,Hcry);
		insert_atoms_recur(10, atom, n, Hcry,false,true);
		compute_CNA_and_others(atom,n, Hcry,false);
		save_cfg(75,Hcry);
		save_lammps(75,Hcry);
		remove = 0;
		for(int i =0;i<n;i++)
		{
			if(atom[i].coord>14) remove++;
		}
		int ting1[remove];
		remove =0;
		for(int i=0;i<n;i++)
		{
			if(atom[i].coord>14)
				{
					ting1[remove] = i;
					remove++;
				}
		}

		delete_atoms(ting1, remove);
		save_cfg(76,Hcry);
		save_lammps(76,Hcry);
		insert_atoms_recur(10, atom, n, Hcry,false,true);
		save_cfg(77,Hcry);
		save_lammps(77,Hcry);

	}

*/
/*
	 //Insert atoms KS1 -Begin
	 int is_read = read_lammps("dat.71", atom,true,true);
	 if(is_read==0)
	 {

		 //compute_CNA_and_others(atom, n, Hcry);
		 //save_cfg(1,Hcry);

//		 int len = 1;
//		 double arr[len][4];
		 //	arr[0] = 24443;
		 //	arr[1] = 31541;
		 //	arr[2] = 35095;

		 //0.630994 0.35502
		 //0.606581 0.351303
		 //0.606745333
		 //0.350902667
//		 arr[0][0] = (0.630994+0.606581)/2.0;
//		 arr[0][1] = (0.35502+0.351303)/2.0;
//		 arr[0][2] =  0.397149667-(0.404739-0.397149667);
//		 arr[0][3] = 1;
//		 insert_atoms(arr,len);
//		 compute_CNA_and_others(atom, n, Hcry);
//		 save_cfg(71,Hcry);
//		 save_lammps(71,Hcry);



		 int len = 2;
		 double arr[len][5];
//		 28777 1 0.606581 0.351303 0.394425 58.8549 34.0367 68.7246 0 0 0 -2.79544

		 arr[0][0] = 0.486399;
		 arr[0][1] = 0.603359;
		 arr[0][2] =  0.394991;
		 arr[0][3] = 1;
		 arr[0][4] = 19358;

		 arr[1][0] = 0.471781;
		 arr[1][1] = 0.622889;
		 arr[1][2] =  0.3981;
		 arr[1][3] = 1;
		 arr[1][4] = 19357;

		 insert_atoms_specific(arr,len);
		 compute_CNA_and_others(atom, n, Hcry);
		 save_cfg(91,Hcry);
		 save_lammps(91,Hcry);
	 }

	 // Insert atoms KS1 - End
*/

	/*
	 int	is_read = read_lammps("dat.insertion_betweenCuCu_vac", atom,true,true);
	 if(is_read==0)
	 {
		 compute_CNA_and_others(atom, n, Hcry);
		 save_cfg(710,Hcry);
		 save_lammps(710,Hcry);
	 }

	 */
	//	seq_process_lammps2();
	/*
	 read_A();
	 compute_CNA_and_others(atom,n,Hcry);
	 save_cfg(1,Hcry);
	 cout <<"herer\n";
	 int is_read = read_lammps("dat.152000",atom,true,true);
	 compute_CNA_and_others(atom,n,Hcry);
	 save_cfg(152000,Hcry);
	 */
/*
	read_A_exact("xyz_mformat");
//		 compute_CNA_and_others(atom, n, Hcry);
		 save_cfg(84,Hcry);
		 save_lammps(84,Hcry);
		 exit(1);

*/

/*
	int is_read = read_cfg("dat.0.cfg",atom);
		 if(is_read==0)
		 {
			 prepare_nbrlist(Hcry,100);
			 coord_number(Hcry);
			 compute_CNA_and_others(atom,n, Hcry);
			 save_cfg(34,Hcry);
			 save_lammps(34,Hcry);
			 for(int i=0;i<n;i++)
			 {
				 if((atom[i].CNA>6)&&(atom[i].interface>0))
				 {
					 double r[3];
					 double s[3];
					 s[0] = atom[i].sx; s[1] = atom[i].sy;s[2]=atom[i].sz;
					 V3mulM3(s,Hcry,r);

				 }
			 }
		 }

exit(1);
*/
	/*
	 // Delete atoms from KSmin - Begin
	 read_A_exact();
	 compute_CNA_and_others(atom, n, Hcry);
	 save_cfg(84,Hcry);
	 save_lammps(84,Hcry);

	 int len = 4;
	 int arr[len];
	 arr[0] = 2746;
	 arr[1] = 30155;
	 arr[2] = 30614;
	 arr[3] = 15315;
	 delete_atoms(arr,len);
	 compute_CNA_and_others(atom, n, Hcry);
	 save_cfg(85,Hcry);
	 save_lammps(85,Hcry);

	 //Delete atoms from KSmin -end
	 */

	/*
	 Hcry[0][0] = 97.151930;Hcry[0][1]=0.0;Hcry[0][2]=0.0;
	 Hcry[1][0] =-0.117442; Hcry[1][1] = 97.123618;Hcry[1][2]=0.0;
	 Hcry[2][0] = -0.552101; Hcry[2][1]=  0.081867; Hcry[2][2] = 174.239148;

	 int is_read = read_lammps("dat_minimized_KS1_nonortho.542", atom,false);
	 if(is_read==0)
	 {
		 prepare_nbrlist(Hcry,100);
		 coord_number(Hcry);
		 //compute_disrigistry(atom,atom_ref,Hcry,Hcry);
		 compute_CNA_and_others(atom,n, Hcry);
		 compute_disrigistry(atom,n,Hcry);

		 save_cfg(2000,Hcry);
		 save_lammps(2000,Hcry);
	 }

	 */

	//	save_cfg(1,Hcry);
	//	seq_process_lammps();
	//trigger_MSD();

/*

	 // Create disrigistry MAPS
//	 int is_read = read_lammps("/Users/kedar/Documents/lib/ref_NbNb_Cuallalpha.lammps", atom,true,true);
//	 int is_read = read_lammps("/Users/kedar/Documents/lib/ref_Nbalphaall_CuCu.lammps", atom,true,true);

	 int is_read = read_lammps("/Users/kedar/Documents/lib/ref_Nb_Nb_alpha_Cu_Cu.lammps", atom,true,true);

	 if(is_read==0)
	 {
		 //		prepare_nbrlist(Hcry,100);
		 //		coord_number(Hcry);
		 //compute_disrigistry(atom,atom_ref,Hcry,Hcry);
		 //		compute_CNA_and_others(atom,n, Hcry);
		 //		compute_disrigistry(atom,n,Hcry);

		 		save_cfg(2000,Hcry);
		 //		save_lammps(2000,Hcry);
	 }else
	 {
		 cout << "file not found\n";
	 }
	 atomic_dat *atom_ref;
	 cout << "entered copy\n";
	 cout << "TESTING COPY\n" << atom[15].sx<<"\t";
	 double H_ref[3][3];
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

	 is_read = read_lammps("dat_minimized_KS1_nonortho.561", atom,true,true);
	 if(is_read==0)
	 {
		 prepare_nbrlist(Hcry,100);
		 coord_number(Hcry);

		 compute_CNA_and_others(atom,n, Hcry);
		 //		compute_disrigistry(atom,n,Hcry);
		 compute_disrigistry(atom,atom_ref,Hcry,H_ref,true);
		 save_cfg(2020,Hcry);
		save_cfg_interface(2020,Hcry,1,1);
		save_cfg_interface(2020,Hcry,2,1);
		 //		save_lammps(2000,Hcry);
	 }

	 //Creat disrigistry maps -end

	*/


	/*

	 // begin - THIS IS FOR TILTING THE KS1 STRUCTURE BY ADOPTING THE VALUES FROM KS1.STRUCTURE.OUT
	 // THE X AND Y DIRECTIONS IN THIS FILE AND THAT IN KS1.STRUCTURE.OUT ARE NOT THE SAME, THEY ARE INTERCHANGED AND HENCE THIS INTERCHANGE THAT IS ADOPTED HERE
	 double t_H[3][3],H_new_new[3][3], H_new_new_inv[3][3];

	 t_H[0][0]=0.9712353759E+02; t_H[0][1]=-.1538852993E+00; t_H[0][2]=-.7539158004E-01;
	 t_H[1][0]=-.1538852193E+00; t_H[1][1] = 0.9715180368E+02; t_H[1][2] = 0.2950252469E-01;
	 t_H[2][0] = -.7539158004E-01; t_H[2][1]= 0.2950252469E-01; t_H[2][2]= 0.1742400233E+03;

	 H_crystal(t_H,crystal0);


	 double tmp;
	 tmp = crystal0[0];
	 crystal0[0] = crystal0[1];
	 crystal0[1] = tmp;
	 tmp = crystal0[3];
	 crystal0[3] = crystal0[4];
	 crystal0[4] = tmp;

	 crystal_H(crystal0,H_new_new);

	 for(int j=0;j<3;j++)
	 {
		 for(int k=0;k<3;k++)
		 {
			 cout << H_new_new[j][k]<<" *****************\t";
		 }
		 cout <<"\n";
	 }

	 M3inv(H_new_new,H_new_new_inv);


	 int is_read = read_lammps("dat_b.1", atom,true,false);
	 if(is_read==0)
	 {
		 for(int j=0;j<3;j++)
		 {
			 for(int k=0;k<3;k++)
			 {
				 cout << Hcry[j][k]<<" **********CRY*******\t";
			 }
			 cout <<"\n";

			 prepare_nbrlist(Hcry,100);
			 coord_number(Hcry);
			 //compute_disrigistry(atom,atom_ref,Hcry,Hcry);
			 compute_CNA_and_others(atom,n, Hcry);

			 save_cfg(999,Hcry);
			 save_lammps(999,Hcry);


		 }
		 for(int i=0;i<n;i++)
		 {
			 double r1[3],s1[3],r[3],s[3];
			 s[0]= atom[i].sx;
			 if(s[0]>=1) s[0]=s[0]-1;
			 if(s[0]<0) s[0]=1+s[0];

			 s[1] = atom[i].sy;
			 if(s[1]>=1) s[1]=s[1]-1;
			 if(s[1]<0) s[1]=1+s[1];

			 s[2] = atom[i].sz;
			 if(s[2]>=1) s[2]=s[2]-1;
			 if(s[2]<0) s[2]=1+s[2];

			 V3mulM3(s,Hcry,r);
			 V3mulM3(r,H_new_new_inv,s1); //without any transformation
			 atom[i].sx = s1[0];
			 atom[i].sy = s1[1];
			 atom[i].sz = s1[2];

		 }

		 prepare_nbrlist(H_new_new,100);
		 coord_number(H_new_new);
		 //compute_disrigistry(atom,atom_ref,Hcry,Hcry);
		 compute_CNA_and_others(atom,n, H_new_new);

		 save_cfg(1000,H_new_new);
		 save_lammps(1000,H_new_new);

	 }else
	 {
		 cout << "unable to read file\n";
	 }
	 */

	/*

	 is_read = read_lammps("dat_minimized_KS1.485", atom,false);
	 if(is_read==0)
	 {
		 prepare_nbrlist(Hcry,100);
		 coord_number(Hcry);
		 //compute_disrigistry(atom,atom_ref,Hcry,Hcry);
		 compute_CNA_and_others(atom,n, Hcry);

		 save_cfg(2000,Hcry);
		 save_lammps(2000,Hcry);
	 }
	 */

	// END - THIS IS FOR TILTING THE KS1 STRUCTURE BY ADOPTING THE VALUES FROM KS1.STRUCTURE.OUT
	// THE X AND Y DIRECTIONS IN THIS FILE AND THAT IN KS1.STRUCTURE.OUT ARE NOT THE SAME, THEY ARE INTERCHANGED AND HENCE THIS INTERCHANGE THAT IS ADOPTED HERE




	//trigger_MSD();

	//	create_system("Cu","Nb");

	/*
	 cout << "entered here\n";
	int is_read = read_lammps("dat_b.1", atom,true,true);
	 if(is_read==0)
	 {
		// prepare_nbrlist(Hcry,100);
		 //coord_number(Hcry);
		 //compute_disrigistry(atom,atom_ref,Hcry,Hcry);
		 //	compute_CNA_and_others(atom,n, Hcry);
		 		create_alpha_system(5,true);
		 save_cfg(1210,Hcry);
		 save_lammps(1210,Hcry);
	 }
	 */

	/*

	 double interface_pe[4][2];
	 FILE *fptr;

	 int start_num =900000;
	 int end_num = 1220000;
	 cout <<"entered here\n";

	 char start_str[20], end_str[20],start_str_e[20],end_str_e[20];


	 sprintf(start_str,"%d",start_num);
	 sprintf(end_str,"%d",end_num);

	 sprintf(start_str_e,"%6d",start_num);
	 sprintf(end_str_e,"%6d",end_num);

	 char start_file[80]="",end_file[80]="";
	 strcat(start_file,"dat.");
	 strcat(start_file,start_str);

	 strcat(end_file,"dat.");
	 strcat(end_file,end_str);

	 char file_ring_ref[80]= "";
	 strcat(file_ring_ref,"rings_data_ref.");
	 strcat(file_ring_ref,start_str);
	 strcat(file_ring_ref,".");
	 strcat(file_ring_ref,end_str);

	 char file_ring_final[80]= "";
	 strcat(file_ring_final,"rings_data_final.");
	 strcat(file_ring_final,start_str);
	 strcat(file_ring_final,".");
	 strcat(file_ring_final,end_str);



	 fptr = fopen(file_ring_ref,"w");

	 //int is_read = read_lammps(start_file, atom,true,true);
	 int is_read = read_cfg("dat.0.cfg",atom);
	 if(is_read==0)
	 {
		 prepare_nbrlist(Hcry,100);
		 coord_number(Hcry);
		 compute_CNA_and_others(atom,n, Hcry);

		 for(int i=0;i<n;i++)
		 {
			 if((atom[i].CNA>6)&&(atom[i].interface>0))
			 {
				 double r[3];
				 double s[3];
				 s[0] = atom[i].sx; s[1] = atom[i].sy;s[2]=atom[i].sz;
				 V3mulM3(s,Hcry,r);
				 fprintf(fptr,"%d %d %d %lf %lf %lf\n",i, atom[i].CNA,atom[i].type,r[0],r[1],r[2]);
			 }
		 }

		 fclose(fptr);


		 compute_interface_energy(atom,n,interface_pe);
		 //save_cfg_interface(start_num, Hcry, 1,1);

		 save_cfg(1000,Hcry);
	 }else
	 {
		 cout << "file not found\n";
	 }
	 atomic_dat *atom_ref;
	 cout << "entered copy\n";
	 cout << "TESTING COPY\n" << atom[15].sx<<"\t";


	 double H_ref[3][3];
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


	 fptr = fopen(file_ring_final,"w");
	 //is_read = read_lammps(end_file, atom,true,true);
	 is_read = read_cfg("dat.24.cfg",atom);
	 if(is_read==0)
	 {
		 prepare_nbrlist(Hcry,100);
		 coord_number(Hcry);

		 compute_CNA_and_others(atom,n, Hcry);
		 compute_disrigistry(atom,atom_ref,Hcry,H_ref,false);

		 char mv_disrig[80]="";
		 strcat(mv_disrig,"mv disrigistry_data");
		 strcat(mv_disrig," disrig.");
		 strcat(mv_disrig,start_str);
		 strcat(mv_disrig,"_");
		 strcat(mv_disrig,end_str);
		 execute_system_command(mv_disrig);

		 compute_slipvector(atom,atom_ref,Hcry,H_ref,false);

		 char mv_slip[80]="";
		 strcat(mv_slip,"mv slip_data");
		 strcat(mv_slip," slip_data.");
		 strcat(mv_slip,start_str);
		 strcat(mv_slip,"_");
		 strcat(mv_slip,end_str);
		 execute_system_command(mv_slip);



		 for(int i=0;i<n;i++)
		 {
			 if((atom[i].CNA>6)&&(atom[i].interface>0))
			 {
				 double r[3];
				 double s[3];
				 s[0] = atom[i].sx; s[1] = atom[i].sy;s[2]=atom[i].sz;
				 V3mulM3(s,Hcry,r);
				 fprintf(fptr,"%d %d %d %lf %lf %lf\n",i, atom[i].CNA,atom[i].type,r[0],r[1],r[2]);
			 }
		 }

		 fclose(fptr);

		 compute_interface_energy(atom,n,interface_pe);

		 save_cfg_interface(end_num,Hcry,1,1);
		 save_cfg(end_num,Hcry);
		 save_lammps(2000,Hcry);
	 }




	 //	make_chain();
	 */
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

/*
 for(int i =1;i<2;i++)
 {
	 char filename[80]="fs.",str_h[20];
	sprintf(str_h,"%d",i);
	strcat(filename,str_h);
	read_A_exact(filename);
	save_lammps(i,Hcry);
 }
*/
/*
read_A_exact(inputfilename);
compute_CNA_and_others(atom,n,Hcry);
save_atoms(atom,n, Hcry);
cout <<"here\n";
save_cfg(format,Hcry);
save_lammps(format,Hcry);
//save_cfg_interface(120,Hcry,1,1);
//save_cfg_interface(120,Hcry,1,0);
//save_cfg_interface(120,Hcry,2,1);

//write_A("ting", 100, atom, Hcry);
exit(1);
*/

/*
int is_read = read_lammps("dat.200",atom,true,true);
compute_CNA_and_others(atom,n,Hcry,true);
write_A("ting", 100, atom, Hcry);
exit(1);
*/
/*
int is_read = read_lammps("dat.71", atom,true,true);
compute_CNA_and_others(atom,n, Hcry,true);
save_atoms(atom,n, Hcry);
cout << "came out safely\n";
compute_CNA_and_others(atom,n, Hcry,true);
int len_hh = 4;
int tagging[len_hh];
tagging[0]= 90848;
tagging[1]=136174;
tagging[2]=129847;
tagging[3]=136946;
for(int c=0;c<len_hh;c++)
{
	atom[tagging[c]].interface = 1;
	for(int i =0;i<MAX_COORD;i++)
	{
		int k = atom[tagging[c]].coord_id[i];
		if(k>-1) atom[k].interface = 1;
	}
}
save_cfg_interface(710,Hcry,1,1);
save_cfg(710,Hcry);
save_lammps(710,Hcry);
*/

//create_system("dat.Cu15", "dat.Nb100", true,true,false);
//create_alpha_system(0, false);
//save_cfg(1210,Hcry);
/*
int is_read = read_lammps("dat.34", atom,true,true);
 if(is_read==0)
 {
	 prepare_nbrlist(Hcry,100);
	 coord_number(Hcry);
	 //compute_disrigistry(atom,atom_ref,Hcry,Hcry);
	 compute_CNA_and_others(atom,n, Hcry);
			create_alpha_system(0,true);
			prepare_nbrlist(Hcry,100);
			coord_number(Hcry);
			compute_CNA_and_others(atom,n, Hcry);	 save_cfg(35,Hcry);
			save_lammps(35,Hcry);
 }
*/
if(SEQ_PROCESS)
{
	seq_process_lammps2();
}


if(CHAIN_MIN)
{
	create_recurring_chain(filenumber_start, filenumber_end, filenumber_interval, CHAIN_START_NUMBER, CHAIN_MIN_NUMBER);
}
if(CHAIN_NEB)
{
	create_recurring_chain_equidistant(filenumber_start, filenumber_end, filenumber_interval, CHAIN_START_NUMBER, CHAIN_NEB_NUMBER);
}
//read_A_exact("rel.structure.out");
//compute_CNA_and_others(atom,n,Hcry);
//save_cfg(0,Hcry);
/*
 read_A_exact("KSmin.structure.out");
 compute_CNA_and_others(atom,n,Hcry);
 save_cfg(0,Hcry);
 save_lammps(0,Hcry);
*/
 /*
read_A_exact("KS1_38_41.a");
compute_CNA_and_others(atom,n,Hcry);
save_cfg(0,Hcry);
save_lammps(0,Hcry);
read_A_exact("KS1_38_41.b");
compute_CNA_and_others(atom,n,Hcry);
save_cfg(20,Hcry);
save_lammps(20,Hcry);
read_A_exact("KS1_38_41.t");
compute_CNA_and_others(atom,n,Hcry);
save_cfg(40,Hcry);
save_lammps(40,Hcry);
*/
//seq_convert_A_lammps(inputfilename);

// start NEB
/*
for(int tt = 0;tt < argc;tt++)
{
	if(strcmp(argv[tt],"neb_save_info")==0)
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
		}else if (strcmp(argv[tt+1],"multiple")==0)
		{
			NEB_CLIMB_CONDITION = 2;
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
		if (strcmp(argv[tt+3],"false")==0) keep_xy = false;
		if (strcmp(argv[tt+4],"false")==0) keepcms_together = false;
		if (strcmp(argv[tt+5],"true")==0) keep_max = true;

		cout << keep_xy<<"\t"<<keepcms_together<<"\t"<<keep_max<<"\n";
		create_system(input1, input2, keep_xy,keepcms_together,keep_max);

	}


}
*/

/*
 int chain_length = 0;
for(int i = filenumber_start;i<=filenumber_end;i=i+filenumber_interval)
{
	chain_length++;
}
int array[chain_length];
int counter=0;
for(int i = filenumber_start;i<=filenumber_end;i=i+filenumber_interval)
{
	array[counter]=i;
	counter++;
}
NEB(array,chain_length);

//end NEB
*/


	return (0);
}

