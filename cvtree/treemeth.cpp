/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-11-21 23:24:12
 */

#include "treemeth.h"
/********************************************************************************
 * @brief the base class of tree method
 *
 * @param methStr
 * @return TreeMeth*
 ********************************************************************************/
TreeMeth *TreeMeth::create(const string &methStr) {
  TreeMeth *meth;
  if (methStr == "NJ") {
    meth = new NeighborJoint;
  } else if (methStr == "UPGMA") {
    meth = new UPGMA;
  } else {
    cerr << "Unknow Tree Method: " << methStr << endl;
    exit(3);
  }

  return meth;
};

void TreeMeth::startree(const Mdist &dm, Node* aTree, StarTree& nodes) {
  // the distance matrix and leafs list are copy from the main program
  // so it will be changed in the main program
  nodes.resize(dm.size());
  for (size_t i = 0; i < dm.size(); ++i) {
    nodes[i] = aTree->reproduct(i, dm.getname(i));
  }
}

/********************************************************************************
 * @brief the functions for neighbor joint method
 *
 * @param dm
 * @return Node*
 ********************************************************************************/
/// For neighborJoint
void NeighborJoint::tree(Node*& aTree, Mdist &dm) {
  // the function to neighbor joint the start tree
  StarTree nodes;
  startree(dm, aTree, nodes);
  Node *outgrp = nodes[0];

  // get the new length of the star tree
  lenStar(nodes, dm);
  int nNode = nodes.size() - 3;

  // get the nj tree
  for (int i = 0; i < nNode; ++i) {
    // get two nearest items and set length and distance matrix
    Neighbor nb;
    njnearest(dm, nodes, nb);

    // joint the two nearest neighbors
    joint(dm, nodes, nb);
  }

  // get the root of the tree which include three branches
  (*nodes[0]).length =
      0.5 * ((*nodes[0]).length - dm.getdist((*nodes[1]).id, (*nodes[2]).id));
  (*nodes[1]).length =
      0.5 * ((*nodes[1]).length - dm.getdist((*nodes[2]).id, (*nodes[0]).id));
  (*nodes[2]).length =
      0.5 * ((*nodes[2]).length - dm.getdist((*nodes[0]).id, (*nodes[1]).id));

  // add the root
  for (auto &nd : nodes)
    (*aTree).addChild(nd);

  // set the outgroup
  aTree = (*aTree).resetroot(outgrp);
}

void NeighborJoint::lenStar(const StarTree &vn, const Mdist &dm) {
  for (auto &np : vn) {
    double len = 0.0;
    size_t idp = np->id;
    for (auto &nd : vn) {
      if (np != nd)
        len += dm.getdist(idp, (*nd).id);
    }
    np->length = len;
  }
};

#pragma omp declare reduction(mindd                                            \
                              : Neighbor                                       \
                              : omp_out =                                      \
                                    omp_in.dd < omp_out.dd ? omp_in : omp_out)

void NeighborJoint::njnearest(const Mdist &dm, StarTree &nodes, Neighbor &nb) {

  // get the nearest neighbor
  size_t nNode(nodes.size());
  double m2star(nNode - 2);
  auto half = nNode / 2;
  nb.dd = numeric_limits<double>::max();
  vector<double> length(nNode);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < nNode; ++i) {
      length[i] = nodes[i]->length / m2star;
    }
#pragma omp barrier

      //... for parallel rearranged the triangle loop into retangle loop
#pragma omp for reduction(mindd : nb)
    for (size_t i = 1; i <= half; ++i) {
      //.. one half
      nb.dd += length[i];
      for (size_t j = 0; j < i; ++j) {
        double ddxy = dm._getdist(nodes[j]->id, nodes[i]->id) - length[j];
        if (ddxy < nb.dd) {
          nb.first = j;
          nb.second = i;
          nb.dd = ddxy;
        }
      }
      nb.dd -= length[i];

      //.. other half
      auto ni = nNode - i;
      nb.dd += length[ni];
      for (size_t j = 0; j < ni; ++j) {
        double ddxy = dm._getdist(nodes[j]->id, nodes[ni]->id) - length[j];
        if (ddxy < nb.dd) {
          nb.first = j;
          nb.second = ni;
          nb.dd = ddxy;
        }
      }
      nb.dd -= length[ni];
    }
  }
};

void NeighborJoint::joint(Mdist &dm, StarTree &nodes, Neighbor &nb) {

  // reset the length branch
  Node *nx = nodes[nb.first];
  Node *ny = nodes[nb.second];
  size_t idx = (*nx).id;
  size_t idy = (*ny).id;

  double mstar = double(nodes.size());
  double dx = (*nx).length;
  double dy = (*ny).length;
  double dxdy = (dx - dy) / (mstar - 2);
  double dxy = dm._getdist(idx, idy);
  (*nx).length = 0.5 * (dxy + dxdy);
  (*ny).length = 0.5 * (dxy - dxdy);

  // set the parent nodes and add the two nodes to the new node
  Node *nz = nx->reproduct(idx);
  (*nz).length = 0.5 * (dx + dy - dxy * mstar);
  (*nz).addChild(ny);
  (*nz).addChild(nx);

  // replace the first node by new nodes and delete the second node
  nodes[nb.first] = nz;
  nodes.erase(nodes.begin() + nb.second);

  // set the distrance between node z and other nodes
  // use the nx index to the index of nz in distrance matrix
  for (size_t i = 0; i < nb.first; ++i) {
    size_t idu = nodes[i]->id;
    double duxy = dm._getdist(idu, idx) + dm._getdist(idu, idy);
    double duz = 0.5 * (duxy - dxy);
    dm.setdist(idx, idu, duz);
    nodes[i]->length += (duz - duxy);
  }

  for (size_t i = nb.first + 1; i < nodes.size(); ++i) {
    size_t idu = nodes[i]->id;
    double duxy = dm._getdist(idx, idu) + dm.getdist(idy, idu);
    double duz = 0.5 * (duxy - dxy);
    dm.setdist(idx, idu, duz);
    nodes[i]->length += (duz - duxy);
  }
}

/********************************************************************************
 * @brief functions for UPGMA method
 *
 * @param dm
 * @return Node*
 ********************************************************************************/
/// For UPGMA
void UPGMA::tree(Node*& aTree, Mdist &dm) {
  // the function to UPGMA the start tree
  StarTree nodes;
  startree(dm, aTree, nodes);

  // run the task for combine into tree
  int nNode = nodes.size();
  for (int i = 1; i < nNode; ++i) {
    // get the closest node pair
    Neighbor nb;
    closest(dm, nodes, nb);

    // combine the closest node pair
    combine(dm, nodes, nb);
  }

  // reduce the branch length
  resetLength(nodes.front());
  delete aTree;
  aTree = nodes.front();
};

// get the closest node pair
void UPGMA::closest(const Mdist &dm, StarTree &nodes, Neighbor &nb) {
  int nNode = nodes.size();
  nb.dd = numeric_limits<double>::max();
  for (int i = 0; i < nNode - 1; ++i) {
    size_t idx = nodes[i]->id;
    for (int j = i + 1; j < nNode; ++j) {
      size_t idy = nodes[j]->id;
      double dist = dm.getdist(idx, idy);
      if (dist < nb.dd) {
        nb.first = i;
        nb.second = j;
        nb.dd = dist;
      }
    }
  }
}

// combine the closest node pair and update the distance
void UPGMA::combine(Mdist &dm, StarTree &nodes, Neighbor &nb) {
  // the two combining nodes
  Node *ny = nodes[nb.second];
  ny->length = 0.5 * nb.dd;
  nodes.erase(nodes.begin() + nb.second);
  Node *nx = nodes[nb.first];
  nx->length = 0.5 * nb.dd;
  nodes.erase(nodes.begin() + nb.first);

  // the new node and replace first node
  Node *nz = nx->reproduct(nx->id);
  nz->nleaf = nx->nleaf + ny->nleaf;
  nz->addChild(nx);
  nz->addChild(ny);

  // update the distance matrix
  for (auto nu : nodes) {
    double dist = dm.getdist(nu->id, nx->id) * nx->nleaf +
                  dm.getdist(nu->id, ny->id) * ny->nleaf;
    dist /= nz->nleaf;
    dm.setdist(nz->id, nu->id, dist);
  }
  nodes.push_back(nz);
}

// update the branch length
void UPGMA::resetLength(Node *nd) {
  if (!nd->isLeaf()) {
    nd->length -= nd->children.front()->length;
    for (auto &nc : nd->children) {
      resetLength(nc);
    }
  }
}