/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-04-19 11:37:42
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-19 23:17:49
 */

#ifndef CV_H
#define CV_H

#include "method.h"
#include "info.h"
struct Args {
  string program;
  vector<string> flist;
  vector<size_t> klist;
  Method* meth;

  Args(int, char **);
  void usage();
};

#endif
