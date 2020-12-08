/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-11-15 20:20:23
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-08 14:57:55
 */

#ifndef STRINGOPT_H
#define STRINGOPT_H

#include <cctype>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

using namespace std;

/********************************************************************************
 * @brief string options
 * 
 * @return string 
 ********************************************************************************/
string Ltrim(const string &);
string Rtrim(const string &);
string trim(const string &);

int separateWord(vector<string> &, string, const string &sep = " ,");

template <class T, class A>
string strjoin(const A &begin, const A &end, const T &t) {
  ostringstream buf;
  A iter = begin;
  buf << *iter;
  for (++iter; iter != end; iter++) {
    buf << t;
    buf << *iter;
  }

  return buf.str();
}

string toUpper(const string &);
string toLower(const string &);

// template <typename T>
// bool isValid(const string& str){
//    bool res = true;
//    try{
//       T tmp = boost::lexical_cast<T>(trim(str));
//    }
//    catch(boost::bad_lexical_cast &e){
//       res = false;
//    }
//    return res;
// };

/********************************************************************************
 * @brief option on the suffix of file name
 * 
 ********************************************************************************/

void addsuffix(string &, char);
string addsuffix(const string&, char);
string addsuffix(const string&, const string&);
string chgsuffix(const string &, const string &);
string getsuffix(const string &);
string delsuffix(const string &);

/********************************************************************************
 * @brief option on convert string to number
 * 
 ********************************************************************************/
int str2int(const string &);
float str2float(const string &);
double str2double(const string &str);

template <class T> void str2number(const string &str, T &v) {
  istringstream iss(str);
  iss >> v;
};

/********************************************************************************
 * @brief options on read columns file
 * 
 ********************************************************************************/
int nColumns(const string&);

//read column-style file 
template <class T>
void readlist(const string &file, vector<T> &list, int ncol = 0) {
  ifstream infile(file.c_str());
  if (!infile) {
    cerr << "\nCannot found the input file " << file << endl;
    exit(1);
  }

  if (ncol == 0) {
    T item;
    while (infile >> item)
      list.push_back(item);
  } else {
    ncol --;
    for (string line; getline(infile, line);) {
      line = trim(line);
      if (!line.empty()) {
        vector<string> items;
        separateWord(items, line);

        if(items.size() > ncol){
          T item;
          istringstream iss(items[ncol]);
          iss >> item;
          list.emplace_back(item);
        }
      };
    }
  }
  infile.close();
};

// unique the vetor
template <class T> void uniqueWithOrder(vector<T> &list) {
  vector<T> tmpVector;
  set<T> tmpSet;
  for (auto item : list) {
    if (tmpSet.insert(item).second) {
      tmpVector.emplace_back(item);
    }
  }
  tmpVector.shrink_to_fit();
  list.swap(tmpVector);
};

/********************************************************************************
 * @brief function by sys stat
 ********************************************************************************/
long getFileSize(const string &);
bool fileExists(const string &);
bool isDirectory(const string &);

#endif
