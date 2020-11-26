/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-03-17 15:39:23
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-11-26 21:15:25
 */

#ifndef TREE_H
#define TREE_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "stringOpt.h"
using namespace std;

typedef pair<string, string> str2str;

struct Node {
  typedef vector<Node *> Children;

  string name;
  size_t id;
  double length;
  double bootstrap;
  Node *parent;
  Children children;

  Node();
  Node(size_t);
  Node(size_t, const string &);
  Node(size_t, const vector<Node *> &);

  void clear();
  void addChild(Node *);
  void deleteChild(Node *);
  void getDescendants(vector<Node *> &);
  void getLeafs(vector<Node *> &);
  bool isLeaf();

  Node *resetroot(const string &);
  Node *resetroot(Node *);

  void _outnwk(ostream &);
  void outnwk(ostream &);
  void _innwk(istream &);
  void _nwkItem(const string&);
  void innwk(istream &);

  void renewId(const unordered_map<string, size_t> &);
};

#endif
