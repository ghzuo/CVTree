/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-02-13 20:15:20
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 21:55:12
 */

#ifndef DISTMATRIX_H
#define DISTMATRIX_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <netcdfcpp.h>
#include <set>
#include <string>
#include <vector>

#include "stringOpt.h"
#include "tree.h"

using namespace std;

// ... only distance matrix
class Mdist {
  int ng;
  vector<double> dist;
  vector<string> name;

public:
  Mdist();
  void init(const vector<string> &);

  void _setdist(size_t, size_t, double);
  double _getdist(size_t, size_t) const;

  // read the matrix
  bool readmtx(const string &, const bool netcdf = false);
  void readmtxnc(const string &);
  void readmtxtxt(const string &);

  // output the matrix
  void writemtx(const string &, const bool netcdf = false);
  void writemtxnc(const string &);
  void writemtxtxt(const string &);

  // get exchange index and genome name
  bool name2ndx(const vector<string> &, vector<size_t> &) const;
  bool ndx2name(const vector<size_t> &, vector<string> &) const;

  // reduce the matrix by name or index list
  void reduce(const vector<size_t> &);
  void reduce(const vector<string> &);

  // extent the matrix
  void extend(const vector<string> &, const vector<double> &);
  void extend(const vector<string> &);
  void extend(const vector<double> &);

  // get the distance from other distance matrix
  void cleanName();
  void assignDM(const Mdist &, vector<size_t> &);
  void assignDM(const Mdist &);

  // check NAN distance
  int chkNAN(const vector<size_t> &, vector<pair<size_t, size_t>> &) const;
  int chkNAN(const vector<string> &, vector<pair<size_t, size_t>> &) const;
  int chkAllNAN(vector<pair<size_t, size_t>> &) const;
  int nNAN(size_t) const;
  int nNAN() const;
  bool isNAN(size_t, size_t) const;
  bool isNAN(size_t) const;
  bool hasNAN(size_t) const;
  bool hasNAN() const;
  string info() const;

  // get/set value of matrix
  double getdist(size_t, size_t) const;
  string getname(size_t) const;
  pair<size_t, size_t> getIndex(size_t) const;
  vector<string> getNameList() const;
  void setdist(size_t, size_t, double);
  void setname(size_t, const string &);

  // global options
  void resize(size_t);
  size_t size() const;
  size_t msize() const;
  size_t capacity() const;
};

#endif

#ifndef NCERR
#define NCERR
static const int NC_ERR = 10;
#endif
