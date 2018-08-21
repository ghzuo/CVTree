/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-07-31 15:49:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-31 15:49:27
 */

#ifndef DISTANCE_H
#define DISTANCE_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "distmatrix.h"
#include "info.h"
#include "memory.h"
#include "method.h"
#include "stringOpt.h"

// select the Kstr class
#include "kstring.h"

using namespace std;

// read arguments
struct CVitem {
  size_t ndx;
  string fname;
  size_t nNAN;
  CVvec cv;

  CVitem() = default;
  CVitem(size_t i, size_t n) : ndx(i), nNAN(n){};
  void fill() { readcv(fname, cv); };
  void clear() { CVvec().swap(cv); };

  bool operator<(const CVitem &it) const { return it.nNAN < nNAN; }
};

struct IterStep {
  static size_t ndx;
  static void reIndex();

  vector<CVitem> cvlist;
  vector<CVitem *> introBlock;
  vector<CVitem *> interBlock;
  float size;

  IterStep(const Mdist &, const vector<string> &);
  void checkSize(float, const Mdist &);

  size_t length() const;
  void fillBlock(); // TODO: mergin the readcv into calcInDist function
  void calcInDist(Mdist &, Method *);
  void calcOutDist(Mdist &, Method *);
  void execute(Mdist &, Method *);
  string info();
};

void assignDM(const string &, bool, Mdist &);
#endif