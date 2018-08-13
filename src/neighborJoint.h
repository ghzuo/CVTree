/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-04-19 11:37:42
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 21:52:13
 */

#ifndef NEIGHBORJOINT_H
#define NEIGHBORJOINT_H

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
typedef pair<int, int> Neighbor;

// build the tree by neighbor joint algorithm
Node *neighborJoint(Mdist &);

// the distance from the star point
void lenStar(const StarTree &, const Mdist &);

// reset the distance of the nearest neighbor
void njnearest(const Mdist &, StarTree &, Neighbor &);

// joint the two neighbors
// reset the distance of the nearest neighbor
void joint(Mdist &, StarTree &, Neighbor &);

#endif
