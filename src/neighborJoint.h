#ifndef NEIGHBORJOINT_H
#define NEIGHBORJOINT_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <limits>

#include "global.h"
#include "tree.h"
#include "distmatrix.h"

using namespace std;

// the starTree
typedef vector<Node*> StarTree;

// build the tree by neighbor joint algorithm
Node* neighborJoint(Mdist&);

// the distance from the star point
void lenStar(const StarTree&, const Mdist&);

// reset the distance of the nearest neighbor
void njnearest(const Mdist&, StarTree&, StarTree::iterator&, StarTree::iterator&);

// joint the two neighbors

// reset the distance of the nearest neighbor
void joint(Mdist&, StarTree&, StarTree::iterator&, StarTree::iterator&);
void recjoint(Mdist&, StarTree&, StarTree::iterator&, StarTree::iterator&);

#endif
