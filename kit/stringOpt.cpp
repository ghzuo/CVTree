/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-09-01 13:03:04
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-04-26 19:52:20
 */

#include "stringOpt.h"

/// separate Word
int separateWord(vector<string> &w, string t, const string& sep) {

  w.clear();
  // convert all unvisible charater into space;
  for (unsigned ix = 0; ix < t.size(); ++ix) {
    char ch_t = t[ix];
    if (ch_t < 33 or ch_t > 126)
      t[ix] = 32;
  }

  vector<string> words_t;
  string::size_type pos = 0, prev_pos = 0;
  while ((pos = t.find_first_of(sep, pos)) != string::npos) {
    words_t.emplace_back(t.substr(prev_pos, pos - prev_pos));
    prev_pos = ++pos;
  }
  words_t.emplace_back(t.substr(prev_pos, pos - prev_pos));
  for (vector<string>::iterator it = words_t.begin(); it != words_t.end();
       ++it) {
    string s_t = *it;
    if (s_t.size() != 0)
      w.emplace_back(s_t);
  }
  return w.size();
}

/// trim string
string Ltrim(const string &str) {
  size_t pos = str.find_first_not_of(" \n\r\t");
  if (pos != string::npos)
    return str.substr(pos);
  else
    return "";
}

string Rtrim(const string &str) {
  size_t pos = str.find_last_not_of(" \n\r\t");
  if (pos != string::npos)
    return str.substr(0, pos + 1);
  else
    return "";
}

string trim(const string &str) { return Ltrim(Rtrim(str)); }

/// output the suffix of filename

string chgsuffix(const string &nm, const string &suf) {
  return nm.substr(0, nm.find_last_of('.') + 1) + suf;
}

string getsuffix(const string &nm) {
  return nm.substr(nm.find_last_of('.') + 1);
}

string delsuffix(const string &nm) {
  return nm.substr(0, nm.find_last_of('.'));
}

/// change string to number
int str2int(const string &str) { return stoi(trim(str)); }

float str2float(const string &str) { return stof(trim(str)); }

double str2double(const string &str) { return stod(trim(str)); }

/// upping and lower the charater
string toUpper(const string &s) {
  string ss(s);
  for (auto &ch : ss)
    ch = toupper(ch);
  return ss;
}

string toLower(const string &s) {
  string ss(s);
  for (auto &ch : ss)
    ch = tolower(ch);
  return ss;
}

void addsuffix(string &str, char c) {
  string::iterator iter = str.end();
  if (*(--iter) != c)
    str += c;
}

long getFileSize(const string &filename) {
  struct stat fileInfo;
  if (stat(filename.c_str(), &fileInfo) != 0)
    return -1;
  return fileInfo.st_size;
}

bool fileExists(const string &filename) {
  struct stat buffer;
  return (stat(filename.c_str(), &buffer) == 0);
};
