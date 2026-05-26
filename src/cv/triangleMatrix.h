/*
 * Copyright (c) 2026
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2026-05-15 22:54:26
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-16 13:25:58
 */


#ifndef TRIANGLEMATRIX_H
#define TRIANGLEMATRIX_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "../kit/kit.h"
using namespace std;

struct TriangleHeader {
  string name;
  long nsize = 0;

  TriangleHeader() = default;
  explicit TriangleHeader(long n) : nsize(n){};
  explicit TriangleHeader(const string &);
  void read(gzFile &);
  void write(gzFile &) const;
  friend ostream &operator<<(ostream &, const TriangleHeader &);
};

struct Mtriangle {
  TriangleHeader header;
  vector<float> data;

  Mtriangle() = default;
  explicit Mtriangle(long n, float d0 = 0.0) : header(n) {
    data.resize(n * (n + 1) / 2, d0);
  };
  explicit Mtriangle(const string &fname) { read(fname); }

  // reset matrix size
  void reset(long n, float d0 = 0.0);

  // get/set/add value of matrix with signed integer indices
  void set(long, long, float);
  void add(long, long, float);
  float get(long, long) const;

  // get/set/add value of matrix with unsigned integer indices
  void set(size_t, size_t, float);
  void add(size_t, size_t, float);
  float get(size_t, size_t) const;

  // internal get/set/add without boundary check
  void _set(size_t, size_t, float);
  void _add(size_t, size_t, float);
  float _get(size_t, size_t) const;

  // index functions
  size_t index(size_t, size_t) const;
  size_t _index(size_t, size_t) const;
  string outIndex(size_t, size_t) const;
  size_t nVertex() const { return header.nsize; }

  // output info
  string info() const;
  void write(const string &, float);
  void read(const string &);

  // output stream
  friend ostream &operator<<(ostream &, const Mtriangle &);
};

#endif