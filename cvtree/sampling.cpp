/*
 * Copyright (c) 2024
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-05-06 13:21:46
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-07 09:04:46
 */

#include "sampling.h"
/********************************************************************************
 * @brief for create a resample method
 *
 * @param methStr method name
 * @param r the ratio for jackknife method if required
 * @return SampleMeth*
 ********************************************************************************/
SampleMeth *SampleMeth::create(const string &methStr, double r) {

  SampleMeth *meth;
  if (methStr == "Bootstrap") {
    meth = new Bootstrap;
  } else if (methStr == "Jackknife") {
    meth = new Jackknife(r);
  } else {
    cerr << "Unknow Resample Method: " << methStr << endl;
    exit(3);
  }

  return meth;
}

/********************************************************************************
 * @brief the bootstrap method
 *
 * @param size the number of the samples
 * @return vector<long>
 ********************************************************************************/
vector<long> Bootstrap::operator()(long size) {
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> distrib(0, size - 1);

  vector<long> ndxlist(size);
  for (auto &ndx : ndxlist) {
    ndx = distrib(gen);
  }

  return std::move(ndxlist);
}

string Bootstrap::wkdir() { return name; };

/********************************************************************************
 * @brief the Jackknife method
 *
 * @param size  the number of the samples
 * @return vector<long>
 ********************************************************************************/
Jackknife::Jackknife(double r) : ratio(r) {
  if (r > 1.0 || r < 0.0) {
    cerr << "The value for jackknife should be in (0.0, 1.0)" << endl;
    exit(3);
  }
  name = "jackknife";
};

vector<long> Jackknife::operator()(long size) {
  // for random
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> distrib(0.0, 1.0);

  // the index list
  long nitem = size * ratio;
  vector<long> ndxlist(size);
  iota(ndxlist.begin(), ndxlist.end(), 0);

  for (auto i = 0; i < nitem; ++i) {
    long rest = size - i;
    long sele = i + rest * distrib(gen);
    long tmp = ndxlist[i];
    ndxlist[i] = ndxlist[sele];
    ndxlist[sele] = tmp;
  }

  // cut the output
  ndxlist.resize(nitem);
  return std::move(ndxlist);
}

string Jackknife::wkdir() { return name + '/' + to_string(ratio); };

/********************************************************************************
 * @brief for testing of the functions
 *
 ********************************************************************************/
#ifdef TESTING
int main(int argc, char *argv[]) {
  SampleMeth *mth;
  long n = 100;
  if (argc > 1)
    n = atoi(argv[1]);
  if (argc > 2) {
    double r = atof(argv[2]);
    mth = SampleMeth::create("Jackknife", r);
  } else {
    mth = SampleMeth::create("Bootstrap");
  }

  vector<long> xx = (*mth)(n);
  long ndx = 0;
  for (auto x : xx) {
    cout << ndx++ << "  " << x << endl;
  }
}
#endif