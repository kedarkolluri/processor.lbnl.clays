//	=================================================================	//
// Name        : F_system_funcs.h
// Author      : Richard E. Baumer
// Copyright   : Copyright 2009 Richard E. Baumer. All rights Reserved
// Description : System command functions
// Date		   : September 17, 2010
//	=================================================================	//

#ifndef F_SYSTEM_FUNCS_H
#define F_SYSTEM_FUNCS_H

//Included libraries and files | Namespaced Declaration
//	=================================================================	//
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

//	FUNCTION DECLARATION - Zip a file
//	=================================================================	//
void zip(string filename	//	Filename of file to ZIP (resulting name has ".gz" appended
		 );

//	FUNCTION DECLARATION - Unzip a file
//	=================================================================	//
void unzip(string filename	//	Filename of file to UNZIP (resulting name has ".gz" removed
		 );

//	FUNCTION DECLARATION - Split a string into an array of substrings, dimensions num_sub_strings
//	=================================================================	//
string* F_split_string(string temp_str, int &num_sub_strings);


//	FUNCTION DEFINITION - Zip a file
//	=================================================================	//
void zip (string filename)
{
	//	Declare global variables	
	ostringstream o;	//	String object	
	string zip_cmd = "gzip -vf ";
	
	//	Zip the file 	
	o.str("");
	o << zip_cmd << filename;
	system(o.str().c_str());
}

//	FUNCTION DEFINITION - Unzip a file
//	=================================================================	//
void unzip (string filename)
{
	//	Declare global variables	
	ostringstream o;	//	String object	
	string unzip_cmd = "gunzip -vf ";
	
	//	Zip the file 	
	o.str("");
	o << unzip_cmd << filename;
	system(o.str().c_str());
}

//	FUNCTION DEFINITION - Split a string
//	=================================================================	//
string* F_split_string(string temp_str, int &num_sub_strings)
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

#endif
