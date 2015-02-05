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

typedef pair<string, string> str2str;

struct Node{
    typedef vector<Node*> Children;

    string name;
    size_t id;
    double length;
    Node* parent;
    Children children;

    size_t taxSize, nleaf, nxleaf, nUpLeaf; // nxleaf is the number of unclassfied leafs
    bool unclassified, uploaded;

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
    
    void checkUnclassified();
    void checkUploaded();
    void getDefineLeafs(vector<Node*>&);
    void getUndefineLeafs(vector<Node*>&);
    void getUndefineNames(vector<string>&);
    void chkLeafsName(vector<string>&, vector<string>&, bool);

    string getStrainName();
    void _getPrediction(string&);
    void outPrediction(ostream&);

    Node* resetroot(const string&);
    Node* resetroot(Node*);

    void _outnwk(ostream&);
    void outnwk(ostream&);
    void _innwk(istream&);
    void innwk(istream&);

    void _injson(istream&);
    void _getStr(istream&, string&);
    void _getKeyValue(string&, string&);
    void injson(istream&);
    void outjson(ostream&);
    void outjsonAbbr(ostream&);

    void renewId(const unordered_map<string,size_t>&);
    void reinitTree();
    void chgLeafName(const str2str&);
};

#endif
