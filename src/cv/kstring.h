/*
 * Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-03-17 15:39:23
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-05 Tuesday 19:36:41
 */

#ifndef KSTRING_H
#define KSTRING_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "fileOpt.h"
#include "genome.h"

using namespace std;

struct Kstr {
  unsigned long ks;

  Kstr();
  Kstr(const vector<Letter> &);
  Kstr(unsigned long);
  static long kmax();

  string decode() const;
  size_t length() const;

  void append(Letter);
  void addhead(Letter);

  void behead();
  void choptail();

  void forward(Letter);
  void backward(Letter);

  bool operator<(const Kstr &) const;
  bool operator>(const Kstr &) const;
  bool operator==(const Kstr &) const;
  bool operator!=(const Kstr &) const;
  bool operator()(const Kstr &, const Kstr &) const;

  friend ostream &operator<<(ostream &, const Kstr &);
};

struct Kstr_Hash {
  size_t operator()(const Kstr &r) const { return size_t(r.ks); };
};

string readvk(const string &, vector<Kstr> &);
void writevk(const string &, const vector<Kstr> &);

#endif
