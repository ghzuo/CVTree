#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <math.h>
#include <algorithm>

#include "stringOpt.h"
#include "global.h"
using namespace std;

struct Node{
    typedef vector<Node*> Children;

    string name;
    size_t id;
    double length;
    Node* parent;
    Children children;

    Node();
    Node(size_t);
    Node(size_t, const string&);
    Node(size_t, const vector<Node*>&);

    void clear();
    void addChild(Node*);
    void deleteChild(Node*);
    void getDescendants(vector<Node*>&);
    void getLeafs(vector<Node*>&);
    bool isLeaf();    

    Node* resetroot(const string&);
    Node* resetroot(Node*);

    void _outnwk(ostream&);
    void outnwk(ostream&);
};

#endif
