/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-08-27 14:15:00
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-08-27 14:15:00
 */

#ifndef TREEMETH_H
#define TREEMETH_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <string>
#include <vector>

#include "distmatrix.h"
#include "tree.h"

using namespace std;

// the starTree
typedef vector<Node *> StarTree;
struct Neighbor {
  int first, second;
  double dd;

  Neighbor() : first(0), second(0), dd(numeric_limits<double>::max()){};
};

struct TreeMeth {
  // the create function
  static TreeMeth *create(const string &);

  // virtual function for different methods
  virtual Node *tree(Mdist &) = 0;
};

// The neighbor Joint Method
struct NeighborJoint : TreeMeth {
  Node *tree(Mdist &) override;

  // the distance from the star point
  void lenStar(const StarTree &, const Mdist &);

  // reset the distance of the nearest neighbor
  void njnearest(const Mdist &, StarTree &, Neighbor &);

  // joint the two neighbors
  // reset the distance of the nearest neighbor
  void joint(Mdist &, StarTree &, Neighbor &);
};

#endif