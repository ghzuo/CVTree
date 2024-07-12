/*
 * Copyright (c) 2024
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-05-06 13:07:26
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-11 19:03:51
 */

#ifndef SAMPLEMETH_H
#define SAMPLEMETH_H

#include <string>
#include <vector>
#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>

using namespace std;

struct SampleMeth {
  string name;
  
  // the create function
  static SampleMeth *create(const string &, double ratio = 0.8);

  // sample genome
  void sampleGenome(const string&, const vector<string>&);

  // virtual function for different
  virtual string wkdir() = 0;
  virtual vector<long> operator()(long size) = 0;
};

// bootstrap
struct Bootstrap : public SampleMeth {
  Bootstrap(){ name = "bootstrap";};
  string wkdir() override;
  vector<long> operator()(long) override;
};

// jackkinfe
struct Jackknife : public SampleMeth {
  double ratio;
  
  Jackknife(double);
  string wkdir() override;
  vector<long> operator()(long) override;
};

#endif // !SAMPLEMETH_H
