/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-07-26 12:49:46
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 12:49:46
 */

#ifndef Method_H
#define Method_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <unordered_map>
#include <vector>

#include "distmatrix.h"
#include "readgenome.h"
#include "stringOpt.h"

// select the Kstr class
#include "kstring.h"

using namespace std;

// the parent class
struct Method {
  GeneType theg;
  string gsuff, cvsuff, cvdir;
  int kmin, kmax;

  Method() = default;

  // two ways for cvdir
  string (Method::*getCVpref)(const string &);
  string gcSameDir(const string &);
  string gcDiffDir(const string &);
  string getCVname(const string &, size_t);

  // from genome to cv
  void setg(const string &);
  void checkK(const vector<size_t> &);
  void setCVdir(const string &);
  void getCV(const string &, const vector<size_t> &, bool check = true);

  // basic function for the method
  static size_t count(const Genome &, size_t, CVmap &);
  static void markov(const CVmap &, const CVmap &, const CVmap &, double,
                     CVmap &);

  // the virtual function for different methods
  static Method *create(const string &, string gtype = "faa");
  virtual void cv(const Genome &, vector<pair<int, CVmap>> &) = 0;
  virtual double dist(const CVvec &, const CVvec &) = 0;
};

// son class for Hao method
struct HaoMethod : public Method {
  HaoMethod() {
    cvsuff = ".cv";
    kmin = 3;
  }

  void cv(const Genome &, vector<pair<int, CVmap>> &) override;
  double dist(const CVvec &, const CVvec &) override;
};

// son class for Li method
struct LiMethod : public Method {
  LiMethod() {
    cvsuff = ".ncv";
    kmin = 1;
  }

  void cv(const Genome &, vector<pair<int, CVmap>> &) override;
  double dist(const CVvec &, const CVvec &) override;
};

// son class for Li method
struct ZuoMethod : public Method {
  ZuoMethod() {
    cvsuff = ".ncv";
    kmin = 1;
  }
  void cv(const Genome &, vector<pair<int, CVmap>> &) override;
  double dist(const CVvec &, const CVvec &) override;
};

#endif
