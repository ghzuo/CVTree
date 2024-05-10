/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-06 18:07:46
 */

#ifndef CVTREE_H
#define CVTREE_H

#include "kit.h"
#include "cvmeth.h"
#include "distmeth.h"
#include "treemeth.h"
#include "marktree.h"
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

void maintree(const Args&);
void initMainDM(const Args&, vector<pair<size_t, Mdist>> &);
void getMainCV(const Args&, const vector<pair<size_t, Mdist>> &);

void doSampleTest(const Args&);
void getSampleCV(const Args&);

void onetree(const Args&, const string&, pair<size_t, Mdist> &, Node*&);

#endif
