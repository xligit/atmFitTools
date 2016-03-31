#ifndef __KEYREAD_H__ 
#define __KEYREAD_H__

//class to read inputs from card files

#include<iostream>
#include<fstream>
#include<string>
#include<TString.h>
#include<map>

//USAGE//
//CREATE WITH keyread("<cardname.txt>")
//ACCESS DATA FROM CARD FILE WITH I.E. .getKeyD("<variable_name>")
//CARD PARAMETER FORMAT (WHITESPACE SENSITIVE)
//<TYPE>,<NAME>=<VALUE>;
//EX:
//i,par1=20;
//f,afloatpar=1.2;
//s,astringpar=characters;
//$ this particular line is a comment line
////////

using namespace std;

#define LINESIZEMAX 1000

class keyread{
public:
 keyread(const char* afile); //initializer
 const char* fname;  //name of card file
 map<string,int> imap; //maps variable name to integer data
 map<string,double> dmap; //maps variable name to double data
 map<string,float> fmap;//maps variable name to float data
 map<string,TString> smap; //maps variable name to TString 
 void readFile(); //fills the maps
 void processLine(TString sline); //processes each line of card file
 double getKeyD(string key);  //call these to return value of variable key
 float  getKeyF(string key);
 int    getKeyI(string key);
 TString getKeyS(string key);
};

#endif
