/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-12 16:57:14
 */

#ifndef CVTREE_H
#define CVTREE_H

#include "cvmeth.h"
#include "distmeth.h"
#include "kit.h"
#include "marktree.h"
#include "sampling.h"
#include "treemeth.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// read arguments
struct Args {
  string program, dmName, treeName;
  string refdm;
  vector<string> glist;
  vector<string> flist;
  vector<string> blist;
  vector<size_t> klist;
  float memorySize;

  CVmeth *cmeth;
  DistMeth *dmeth;
  TreeMeth *tmeth;
  SampleMeth *smeth;

  Args(int, char **);
  void usage();
};

void maintree(const Args &);
void doSampleTest(const Args &);

void initDMs(const vector<size_t> &, const vector<string> &, const string &,
             vector<pair<size_t, Mdist>>&);
void getCVs(CVmeth *, const vector<string> &,
            const vector<pair<size_t, Mdist>> &);
void getDM(DistMeth *, CVmeth *, const vector<string> &, const string &,
           pair<size_t, Mdist> &);
void getTree(TreeMeth *, const string &, pair<size_t, Mdist> &);
void bootTree(const string&, const string&, const vector<string>&);

#endif
