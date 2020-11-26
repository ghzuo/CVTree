/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-11-14 11:33:59
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 21:52:05
 */

#ifndef DM2TREE_H
#define DM2TREE_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <regex>
#include <string>
#include <vector>
//#include <omp.h>

#include "distmatrix.h"
#include "info.h"
#include "treemeth.h"
#include "stringOpt.h"
#include "tree.h"
using namespace std;

// read arguments
struct Args {
  string program, distfile, outfile;
  vector<size_t> splist;
  TreeMeth* meth;

  Args(int, char **);
  void usage();
};

#endif
