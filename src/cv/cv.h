/*
 * Copyright (c) 2025
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2025-12-06 Saturday 13:42:35
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-05 Tuesday 14:54:15
 */

#ifndef CV_H
#define CV_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#include "kstring.h"

using namespace std;

typedef pair<Kstr, double> CVdim;
typedef unordered_map<Kstr, double, Kstr_Hash> CVmap;
enum LPnorm { L0, L1, L2 };

struct CVvec {
  static LPnorm lp;
  vector<CVdim> vec;
  double norm;

  CVvec() : norm(0.0){};
  CVvec(size_t n) : vec(n), norm(0.0){};
  CVvec(const string &f) { read(f); };
  CVvec(CVmap &);
  CVvec(const vector<Kstr>&, double val=1.0);

  vector<CVdim>::iterator begin() { return vec.begin(); };
  vector<CVdim>::iterator end() { return vec.end(); };
  vector<CVdim>::const_iterator begin() const { return vec.begin(); };
  vector<CVdim>::const_iterator end() const { return vec.end(); };

  size_t size() const { return vec.size(); };
  void clear() { vec.clear(); norm = 0.0; };
  void setNormL0();
  void setNormL1();
  void setNormL2();
  void normalize();
  CVmap toMap();
  vector<Kstr> toVK();

  void write(const string &);
  string read(const string &);

  friend ostream &operator<<(ostream &, const CVvec &);
};

size_t cvsize(const string &);
#endif // CV_H