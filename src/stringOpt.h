#ifndef STRINGOPT_H
#define STRINGOPT_H

#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <string>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <chrono>

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
   istringstream iss(str);
   iss >> v;
};

long getFileSize(const string&);
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

struct Timer{
    std::chrono::system_clock::time_point start;

    Timer():start(std::chrono::system_clock::now()){};

    double elapsed(){
        auto now = std::chrono::system_clock::now();
        return std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    };
};

#endif
