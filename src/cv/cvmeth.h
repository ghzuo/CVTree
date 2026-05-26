/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-21 Thursday 15:54:55
 */

#ifndef CVMETH_H
#define CVMETH_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <unordered_map>
#include <vector>

#include "aligncv.h"
#include "cv.h"
#include "cvarray.h"
#include "cvbundle.h"
#include "kit.h"

struct CVmeth {
  string cvsuff;
  string cvdir;
  long kmin = 1;
  long kmax = 14;
  bool save = true;

  // the create function
  static CVmeth *ptr;
  static CVmeth *create(const string &, const string &, const string &,
                        const string &cgfile = "", bool save = true);
  static CVmeth *get() { return ptr; }

  // initial method
  void init(const string &, const string &);
  void setCVdir(const string &);

  // get the cvname for diffent cvdir
  function<string(const string &, size_t)> getCVname;

  // from genome to cv
  void checkK(const vector<size_t> &);
  void checkK(size_t);

  // execute the caculate
  void execute(const string &, const vector<size_t> &, bool chk = true);
  CVvec getCV(const string &, int);
  CVBundle getCVB(const string &, int);
  CVArray getCVA(const string &, int);
  size_t infoCVA(const string &, int);
  void _caclCVA(const string &, int, CVArray &);

  // basic function for the method
  size_t count(const Genome &, size_t, CVmap &);
  size_t count(const Gene &, size_t, CVmap &);

  // virtual function for different
  virtual void cv(const Genome &, vector<pair<int, CVmap>> &) = 0;
  virtual void cv(const Gene &, vector<pair<int, CVmap>> &) = 0;
};

// son class for Hao method
struct HaoMethod : public CVmeth {
  HaoMethod() {
    kmin = 3;
    cvsuff = ".Hao";
  };
  void cv(const Genome &, vector<pair<int, CVmap>> &) override;
  void cv(const Gene &, vector<pair<int, CVmap>> &) override;
  void markov(const CVmap &, const CVmap &, const CVmap &, double, CVmap &);

  template <typename T> void docv(const T &, vector<pair<int, CVmap>> &);
};

// son class for Li method
struct Counting : public CVmeth {
  Counting() { cvsuff = ".Count"; };
  void cv(const Genome &, vector<pair<int, CVmap>> &) override;
  void cv(const Gene &, vector<pair<int, CVmap>> &) override;
};

#endif