/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-08-01 22:27:43
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-08-01 22:27:43
 */
#ifndef INFO_H
#define INFO_H

#include "stringOpt.h"

using namespace std;
struct Info {
  bool quiet;
  int dep;
  Timer mytimer;

  Info();
  ~Info();

  void operator()(const string&, int idep=0);
};

#endif