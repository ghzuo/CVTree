#ifndef STRINGOPT_H
#define STRINGOPT_H

#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <string>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <sys/stat.h>

using namespace std;

string Ltrim(const string&);
string Rtrim(const string&);
string trim(const string&);

int separateWord(vector<string>&, string);

void addsuffix(string&, char);
string chgsuffix(const string&, const string&);
string getsuffix(const string& nm);

string toUpper(const string&);
string toLower(const string&);

int str2int(const string&);
float  str2float(const string&);
double str2double(const string& str);

template <class T> 
void str2number(const string& str, T& v){
    v = boost::lexical_cast<T>(trim(str));
};

long getFileSize(const string&);

template <typename T>
bool isValid(const string& str){
    bool res = true;
    try{
	T tmp = boost::lexical_cast<T>(trim(str));
    } 
    catch(boost::bad_lexical_cast &e){
	res = false;
    }
    return res;
};

template <class T> 
void readlist(const string& file, vector<T>& list){
    ifstream infile(file.c_str());
    if(!infile){
        cerr << "\nCannot found the input file "
             << file << endl;
        exit(1);
    }

    T item;
    while(infile >> item) 
	list.push_back(item);
    infile.close();
};
#endif
