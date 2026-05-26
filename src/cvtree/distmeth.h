/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-14 Sunday 17:37:12
 */

#ifndef DISTMETH_H
#define DISTMETH_H

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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "aligncv.h"
#include "cv.h"
#include "distmatrix.h"
#include "info.h"
#include "memory.h"
#include "stringOpt.h"

using namespace std;

// read arguments
struct CVitem {
  size_t ndx;
  long nNAN;
  string fname;
  CVvec cv;

  CVitem() = default;
  CVitem(size_t i, const string &str) : ndx(i), nNAN(0), fname(str){};
  void fill() { cv.read(fname); };
  void clear() { cv.clear(); };
};

struct DistMeth {
  vector<CVitem> cvlist;
  vector<CVitem *> steplist;
  vector<CVitem *> restlist;
  float maxM;

  // the create function
  static DistMeth *create(const string &);

  // the base function
  void setMaxMem(float, size_t, size_t);
  void setflist(const vector<string> &);

  // divided task in steps
  float setStep(const Mdist &);
  string infoStep(long, float);
  size_t length() const;
  void updateStep(const Mdist &);

  // execute calculation
  void fillBlock();
  void setdist4NAN(Mdist &, const CVitem &, const CVitem &);
  void calcInDist(Mdist &);
  void calcOutDist(Mdist &);
  void execute(const vector<string> &, Mdist &);

  // the virtual function for different methods
  virtual void preprocess(CVvec &) = 0;
  virtual double dist(const CVvec &, const CVvec &) = 0;
};

// ... son class for different method
// ... distance scaling at L2
struct Cosine : public DistMeth {
  void preprocess(CVvec &) override;
  double dist(const CVvec &, const CVvec &) override;
};

struct Euclidean : public DistMeth {
  void preprocess(CVvec &) override;
  double dist(const CVvec &, const CVvec &) override;
};

// ... distance scaling at L1
struct InterList : public DistMeth {
  void preprocess(CVvec &) override;
  double dist(const CVvec &, const CVvec &) override;
};

struct Min2Max : public DistMeth {
  void preprocess(CVvec &) override;
  double dist(const CVvec &, const CVvec &) override;
};

// ... distance scaling at L0
struct InterSet : public DistMeth {
  void preprocess(CVvec &) override;
  double dist(const CVvec &, const CVvec &) override;
};

struct Dice : public DistMeth {
  void preprocess(CVvec &) override;
  double dist(const CVvec &, const CVvec &) override;
};

struct ItoU : public DistMeth {
  void preprocess(CVvec &) override;
  double dist(const CVvec &, const CVvec &) override;
};

#endif