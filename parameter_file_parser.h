// File: parameter_file_parser.h
// Desc: Parser of the parameter input file to change default parameter values

#ifndef __PARAMETER_FILE_PARSER__
#define __PARAMETER_FILE_PARSER__

#include <iostream>
#include <map>
#include <string>
#include <algorithm>

template<class T=double>
class PARAMETER_FILE_PARSER
{
 public:
  PARAMETER_FILE_PARSER(string filename) : parameter_filename(filename)
  {param_map=NULL;}

  ~PARAMETER_FILE_PARSER() {if(param_map) delete param_map;}

  int Parse_Parameter_File();

  bool Get_Value(string key, int& value) 
  {
    if(!param_map) return false;
    map<string,string>::iterator i = (*param_map).find(key);
    if(i != (*param_map).end()){
      value = atoi(i->second.c_str());
      return true;
    }else return false;
  }

  bool Get_Value(string key, double& value) 
  {
    if(!param_map) return false;
    map<string,string>::iterator i = (*param_map).find(key);
    if(i!=(*param_map).end()){
      value = atof(i->second.c_str());
      return true;
    }else return false;
  }

  bool Get_Value(string key, string& value) 
  {
    if(!param_map) return false;
    map<string,string>::iterator i = (*param_map).find(key);
    if(i!=(*param_map).end()){
      string raw_val = i->second;
      size_t ip1,ip2;
      ip1 = raw_val.find('"');
      ip2 = raw_val.find('"',ip1+1);
      if(ip1!=string::npos && ip2!=string::npos){
	value = raw_val.substr(ip1+1,ip2-ip1-1);
	return true;
      }else return false;
    }else return false;
  }
 /*
  bool Get_Value(string key, bool& value) 
  {
    if(!param_map) return false;
    map<string,string>::iterator i = (*param_map).find(key);
    if(i!=(*param_map).end()){
      switch(i->second){ //make lowercase
      case "true": value = true;break;
      case "false": value = false;break;
      }
      return true;
      }else return false;
    }else return false;
  }
 */
  bool Processed() {return param_map!=NULL;}

 private:
  map<string, string> *param_map;
  string parameter_filename;
};
//*****************************************************************************
// Parsing parameter_filename and populating the parameters map
//*****************************************************************************
template<class T>
int PARAMETER_FILE_PARSER<T>::Parse_Parameter_File()
{
  ifstream pfile(parameter_filename.c_str());
  if (!pfile.is_open())
  {cout << "Unable to open parameter file:parameter_filename"<<endl;return 0;}

  param_map = new map<string, string>();

  string line,key, value;
  while (getline(pfile, line))
  {
    istringstream tokens(line);
    while(tokens >> ws && getline(tokens, key, '=') && 
	  tokens >> ws && getline(tokens, value)){
      key.erase(remove(key.begin(),key.end(),' '),key.end());
      size_t i = value.find("//");
      if(i!=string::npos) {
	value = value.substr(0,i-1);
	value.erase(remove(value.begin(),value.end(),' '),value.end());
      }
      (*param_map)[key] = value;
    }
  }
 
  pfile.close();
  return 1;
}
//*****************************************************************************
#endif
