/*
 * Utilities.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */
#include "Utilities.h"
#include "global.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>


bool UTILITIES_ZIP = true;

int min(int a, int b)
{
	//find minimum of 2 integers
	if(a>b)
	{
		return b;
	}else
	{
		return a;
	}
}


double min(double a, double b)
{
	//find minimum of 2 doubles
	if(a>b)
	{
		return b;
	}else
	{
		return a;
	}
}

int max(int a, int b)
{
	if(a>b)
	{
		return a;
	}else
	{
		return b;
	}
}


double max(double a, double b)
{
	if(a>b)
	{
		return a;
	}else
	{
		return b;
	}
}
bool check_repeat(int j,int *arr,int len)
{
	/*
	check if an item already exists in the array
	to make an array work like a set
	*/
	for(int k =0;k <len;k++)
    {
		if(arr[k]==j)
		{
			return true;
		}
    }

	return false;
}



int run_command(const char *strCommand, bool vfork_form)
{
  int iForkId, iStatus;
  if(vfork_form)
    {
      iForkId = vfork();
    }else
    {
      iForkId = fork();
    }
  if (iForkId == 0)// This is the child
    {
      iStatus = execl("/bin/sh","sh","-c", strCommand, (char*) NULL);
      exit(iStatus);// We must exit here,
      // or we will have multiple
      // mainlines running...
    }
  else if (iForkId > 0)// Parent, no error
    {
      iStatus = 0;
    }
  else// Parent, with error (iForkId == -1)
    {
      iStatus = -1;
    }
  return(iStatus);
}

int execute_command_2( char *command)
{
  int iNumProc = 0, iChildiStatus = 0, iStatus = 0, iDeadId = 0;
  int iExitFlag = 0;

  iStatus = run_command(command,true);
  if (!iStatus)
    iNumProc++;


  // Wait till the commands complete
  while (iNumProc && !iExitFlag)
    {
      iDeadId = waitpid(-1, &iChildiStatus, WNOHANG);
      if (iDeadId < 0)
        {
          // Wait id error - exit the loop
          iExitFlag = 1;
          cout << "Error ****\n";
          cout << command <<"\n";
          cout << "command could not be executed .. error occured\n";
          cout << "Error ****\n";
        }
      else if (iDeadId > 0)
        {
          iNumProc--;
        }
      else  // iDeadId == 0, no processes died
        {
	  //          sleep(3);// give them time to die
        }
    }
  return 0;
}

void execute_system_command( char *command)
{
    execute_command_2(command);
}


void zip (char *filename)
{

   if(UTILITIES_ZIP)
   {
	char command[80]= " ";
    strcat(command,"gzip -vf ");
    strcat(command,filename);
    execute_command_2(command);
   }
}


FILE *pipestream( char *command)
{
	if( system(command)==-1)
	{
		cout << "failed execute of the command\n";
		cout << command<<"\n";
	}

}

void unzip (char *filename)
{
 if(UTILITIES_ZIP)
 {
	 int return_val = 0;
	 char command[80]= " ";
	 strcat(command,"gunzip -vf ");
	 strcat(command,filename);
	 strcat(command,".gz");
	 execute_command_2(command);
 }
}

//sort ascending order
void sort (int count,int *data)
{
	for(int i = 0;i< (count-1);i++)
	{
		for(int j= i;j<count;j++)
		{
			if(data[i]> data[j])
			{
				int a = data[i];
				data[i] = data[j];
				data[j] = a;
			}
		}
    }
}

void V3cross (double a[3], double b[3], double c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

void V3norm(double *B)
{
	double a = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
	B[0] /=a;  B[1] /=a;  B[2] /=a;

	return;
}

double V3dot (double *a, double *b)
{
	double dot_p;
	double an = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	double bn = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);

	dot_p = (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/an/bn;
	return dot_p;
}

double V3mulM3 (double a[3], double B[3][3], double c[3])
{
	c[0] = a[0]*B[0][0] + a[1]*B[1][0] + a[2]*B[2][0];
	c[1] = a[0]*B[0][1] + a[1]*B[1][1] + a[2]*B[2][1];
	c[2] = a[0]*B[0][2] + a[1]*B[1][2] + a[2]*B[2][2];
}


double M3mulV3(double A[3][3], double b[3], double c[3])
{
	c[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
	c[1] = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
	c[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];

}

double M3inv (double A[3][3], double B[3][3])
{
  double determinant;
  B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1];
  B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2];
  B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0];
  B[1][0] = A[1][2]*A[2][0]-A[1][0]*A[2][2];
  B[2][1] = A[2][0]*A[0][1]-A[2][1]*A[0][0];
  B[0][2] = A[0][1]*A[1][2]-A[0][2]*A[1][1];
  B[2][0] = A[1][0]*A[2][1]-A[2][0]*A[1][1];
  B[0][1] = A[2][1]*A[0][2]-A[0][1]*A[2][2];
  B[1][2] = A[0][2]*A[1][0]-A[1][2]*A[0][0];
  determinant = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
  if (determinant == 0.)
  {

	printf ("error: M3inv: determinant = %e\n"
              "matrix is singular\n", determinant);
	cout<<"\t"<<A[0][0]<<"\t"<<A[0][1]<<"\t"<<A[0][2]<<"\t"<<A[1][0]<<"\t"<<A[1][1]<<"\t"<<A[1][2]<<"\t"<<A[2][0]<<"\t"<<A[2][1]<<"\t"<<A[2][2]<<"\n";
      exit(1);
  }
  B[0][0] /= determinant;
  B[1][1] /= determinant;
  B[2][2] /= determinant;
  B[1][0] /= determinant;
  B[2][1] /= determinant;
  B[0][2] /= determinant;
  B[2][0] /= determinant;
  B[0][1] /= determinant;
  B[1][2] /= determinant;
  return (determinant);
} /* end M3inv() */

void M3mul (double A[3][3], double B[3][3], double C[3][3])
{
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
  C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
  C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];
  C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
  C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
  C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
  return;
} /* end M3mul() */

void M2mul (double A[2][2], double B[2][2], double C[2][2])
{
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
  return;
}

void M2nms (double A[2][2],double C[2][2], double *trace, double *vmises)
{

	C[0][0] = (A[0][0]*A[0][0] + A[0][1]*A[1][0])/2.0-0.5;
  C[0][1] = (A[0][0]*A[0][1] + A[0][1]*A[1][1])/2.0;
  C[1][0] = (A[1][0]*A[0][0] + A[1][1]*A[1][0])/2.0;
  C[1][1] = (A[1][0]*A[0][1] + A[1][1]*A[1][1])/2.0-0.5;

  *trace = (C[0][0]+C[1][1])/2.0;
  *vmises =sqrt(C[1][0]*C[1][0]+C[0][1]*C[0][1]+((C[0][0]-C[1][1])*(C[0][0]-C[1][1]))/6);
  return;
} /* end M3mul() */


void H_crystal(double H[3][3],double crystal[6])
{
	for(int i=0;i<3;i++) crystal[i] = sqrt(H[i][0]*H[i][0]+H[i][1]*H[i][1]+H[i][2]*H[i][2]);
	crystal[5] = acos((H[0][0]*H[1][0]+H[0][1]*H[1][1]+H[0][2]*H[1][2])/crystal[0]/crystal[1]);

	crystal[3] = acos((H[2][0]*H[1][0]+H[2][1]*H[1][1]+H[2][2]*H[1][2])/crystal[2]/crystal[1]);

	crystal[4] = acos((H[0][0]*H[2][0]+H[0][1]*H[2][1]+H[0][2]*H[2][2])/crystal[0]/crystal[2]);

}

void crystal_H(double crystal[6],double H[3][3])
{
	H[0][0] = crystal[0];H[0][1] =0.0;H[0][2]=0.0;

	H[1][0] = cos(crystal[5])*crystal[1];
	H[1][1] = sin(crystal[5])*crystal[1];
	H[1][2] = 0.0;

	H[2][0] = crystal[2]*cos(crystal[4]);
	H[2][1] = (crystal[2]*crystal[1]*cos(crystal[3])-H[2][0]*H[1][0])/H[1][1];
	H[2][2] = sqrt(crystal[2]*crystal[2]-H[2][0]*H[2][0]-H[2][1]*H[2][1]);
}

bool check_file_exists(char *filename)
{
 FILE *fptr = NULL;
 char filename2[80]="";
 strcat(filename2,filename);
 fptr = fopen(filename2,"r");
 if(fptr==NULL)
 {
	 if(UTILITIES_ZIP) strcat(filename2,".gz");
	fptr=fopen(filename2,"r");
	if(fptr==NULL)
		{
			return false;
		}else
		{
			fclose(fptr);
			return true;
		}
 }else
 {
	 fclose(fptr);
	 return true;
 }
}

bool check_file_exists_bash(char *filename)
{
 FILE *fptr = NULL;
 char filename2[80]="";
 strcat(filename2,filename);
 fptr = fopen(filename2,"r");
 if(fptr==NULL)
 {
	if( UTILITIES_ZIP) strcat(filename2,".gz");
	fptr=fopen(filename2,"r");
	if(fptr==NULL)
		{
		return false;
		}else
		{
			fclose(fptr);
			return true;
		}
 }else
 {
	 fclose(fptr);
	 return true;
 }
}

string* split_string(string temp_str, int &num_sub_strings)
{
	int ct = 0;
	istringstream iss1(temp_str);
	string sub;
	do
	{
		iss1 >> sub;
		ct++;
	} while (iss1);

	string* ptr_str = new string[ct];

	ct=0;
	istringstream iss2(temp_str);
	do
	{
		iss2 >> ptr_str[ct];
		ct++;
	} while (iss2);

	num_sub_strings = ct-1;
	return ptr_str;
}

string* get_next_splits( ifstream &inputfile_h, int &numval)
{
	string tmp_line_h;
	getline(inputfile_h,tmp_line_h);
	return split_string(tmp_line_h, numval);
}
