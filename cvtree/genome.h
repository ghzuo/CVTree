/*
 * Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-04-19 11:37:42
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: Thu May 09 2024
 */

#ifndef READGENOME_H
#define READGENOME_H

#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "stringOpt.h"
#include "sampling.h"
using namespace std;

struct Letter {
  static char encMap[128];
  static string decMap;
  static long nbase;
  static vector<Letter> nonVoid;

  static void init(const string &, const string&);
  static void aainit();
  static void nainit();
  static void decMap4encMap();
  static void setLowercaseUpper();
  static void simplify(const string&);
  static void renewAuxiliary();
  static void initByStr(const string&);
  static char decode(size_t);

  char e;

  Letter() = default;
  Letter(char);
  
  friend ostream& operator<<(ostream&, const Letter);
};

struct Gene {
  string head;
  string seq;
  vector<Letter> code;

  Gene() = default;
  Gene(const string &str);
  size_t size() const;
  Letter operator[](size_t) const;
  vector<Letter> substr(size_t, size_t) const;
  void translate();
  friend ostream& operator<<(ostream&, const Gene&);
};

typedef vector<Gene> Genome;
size_t readFasta(const string &, Genome &);
void writeFasta(const string&, Genome&);
Genome sampleGenome(const Genome&, SampleMeth*);

#endif
