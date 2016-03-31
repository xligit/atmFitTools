#include "keyread.h"

keyread::keyread(const char* afile){
  fname = afile;
}

double keyread::getKeyD(string key){
  return dmap[key];
}

float keyread::getKeyF(string key){
  return fmap[key];
}

int keyread::getKeyI(string key){
  return imap[key];
}

TString keyread::getKeyS(string key){
  return smap[key];
}

void keyread::readFile(){
  ifstream file(fname);
  if (!file.is_open()){
    cout<<fname<<" is not a good file "<<endl;
    return;
  }
  char line[LINESIZEMAX];
  int  nread;
  TString sline;
  int nline=0;
  while(!file.eof()){
   // cout<<"reading line: "<<nline<<endl;
    nline++;
    file.getline(line,LINESIZEMAX);
    for (int i=0;i<file.gcount();i++){
      //cout<<"appending "<<line[i]<<" to tstring"<<endl;
      sline.Append(line[i]);
    }
    //cout<<"processing line: "<<endl;
    processLine(sline);
    sline.Clear();
  }
}

void keyread::processLine(TString sline){
  TString stype;
  TString skey;
  TString sval;
  int ival;
  float fval;
  double dval;
  ///get key type
  if (sline(0)=='i') stype="i";
  else if (sline(0)=='f') stype="f";
  else if (sline(0)=='d') stype="d";
  else if (sline(0)=='s') stype="s";
  else if (sline(0)=='$'){
    return;
  }
  else if (sline(0)==' ' || sline(0)=='\n'){
    return;
  }
  else {
    cout<<sline(0)<<" is not a valid type!!"<<endl;
    return;
  }
  ///get key name
  int ichar=2;
  while (sline(ichar)!='='){
    skey.Append(sline(ichar));
    ichar++;
  }
  ///get key value
  int jchar=ichar+1;
  while (sline(jchar)!=';'){
    sval.Append(sline(jchar));
    jchar++;
  }
  //place key value in appropriate map
  if (stype(0)=='i'){
    ival = sval.Atoi();
    cout<<"mapping "<<skey.Data()<<" to "<<ival<<endl;
    imap[skey.Data()] =  ival;
  }
  if (stype(0)=='f'){
    fval = sval.Atof(); 
    cout<<"mapping "<<skey.Data()<<" to "<<fval<<endl;
    fmap[skey.Data()]=fval;
  }
  if (stype(0)=='d'){
    dval=sval.Atof();
    cout<<"mapping "<<skey.Data()<<" to "<<dval<<endl;
    dmap[skey.Data()]=dval;
  }
  if (stype(0)=='s'){
    cout<<"mapping "<<skey.Data()<<" to "<<sval.Data()<<endl;
    smap[skey.Data()]=sval;
  }
   
}
