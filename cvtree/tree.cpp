/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-03-17 15:39:23
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-11-26 21:14:14
 */

#include "tree.h"
const size_t N_FORKS(2);

/*******************************************************************/
/********* Member Functions For Node Class *************************/
/*******************************************************************/
Node::Node() : name(""), id(0), length(NAN), parent(NULL) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n) : name(""), id(n), length(NAN), parent(NULL) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n, const string &str)
    : name(str), id(n), length(NAN), parent(NULL) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n, const vector<Node *> &vn)
    : name(""), id(n), length(NAN), parent(NULL), children(vn) {
  children.reserve(N_FORKS);
};

bool Node::isLeaf() { return children.empty(); };

void Node::addChild(Node *nd) {
  children.emplace_back(nd);
  (*nd).parent = this;
};

void Node::deleteChild(Node *nd) {
  Children::iterator iter = find(children.begin(), children.end(), nd);
  children.erase(iter);
};

void Node::clear() {
  vector<Node *> nodes;
  getDescendants(nodes);
  vector<Node *>::iterator iter = nodes.begin();
  vector<Node *>::iterator iterEnd = nodes.end();
  for (; iter != iterEnd; ++iter)
    delete *iter;
  children.clear();
}

void Node::getDescendants(vector<Node *> &nds) {
  if (!isLeaf()) {
    for (const auto &nd : children) {
      nds.emplace_back(nd);
      (*nd).getDescendants(nds);
    }
  }
};

void Node::getLeafs(vector<Node *> &nodes) {
  if (isLeaf()) {
    nodes.emplace_back(this);
  } else {
    for (const auto &nd : children)
      (*nd).getLeafs(nodes);
  }
};

Node *Node::resetroot(const string &str) {
  vector<Node *> leafs;
  getLeafs(leafs);
  Node *og = NULL;
  for (Node *nd : leafs) {
    if ((*nd).name == str) {
      og = nd;
      break;
    }
  }

  if (og == NULL)
    return og;
  return resetroot(og);
}

Node *Node::resetroot(Node *np) {
  Node *outgrp = np;
  vector<Node *> nlist;
  while ((*np).parent != NULL) {
    nlist.emplace_back((*np).parent);
    np = (*np).parent;
  }

  for (vector<Node *>::reverse_iterator iter = nlist.rbegin();
       iter != nlist.rend() - 1; ++iter) {
    vector<Node *>::reverse_iterator next = iter + 1;
    (**iter).length = (**next).length;
    (**next).length = NAN;
    (**iter).deleteChild(*next);
    (**next).parent = NULL;
    (**next).addChild(*iter);
  }
  Node *theRoot = nlist.front();

  // set the outgroup as the last child
  Children::reverse_iterator iter = (*theRoot).children.rbegin();
  if ((*iter) != outgrp) {
    for (++iter; iter != (*theRoot).children.rend(); ++iter) {
      if ((*iter) == outgrp) {
        *iter = (*theRoot).children.back();
        (*theRoot).children.back() = outgrp;
      }
    }
  }

  return theRoot;
};

void Node::_outnwk(ostream &os) {
  if (!isLeaf()) {
    os << "(";
    Children::const_iterator iter = children.begin();
    (*(*iter))._outnwk(os);
    for (++iter; iter != children.end(); ++iter) {
      os << ",";
      (*(*iter))._outnwk(os);
    }
    os << ")";
  }

  string forename = name.substr(name.find_first_of('|') + 1);
  if (name.find(' ') == std::string::npos)
    os << forename;
  else
    os << '"' << forename << '"';

  if (!std::isnan(length))
    os << ":" << fixed << setprecision(5) << length;
};

void Node::outnwk(ostream &os) {
  _outnwk(os);
  os << ";" << endl;
};

void Node::_nwkItem(const string &str) {
  vector<string> words;
  separateWord(words, str, ":");
  if (isLeaf()) {
    name = words.front();
  } else {
    try {
      bootstrap = str2double(words.front());
    } catch (exception e) {
      name = words.front();
    }
  }
  length = str2double(words.back());
};

void Node::_innwk(istream &is) {
  Node *np = new Node;
  addChild(np);
  string brStr;

  while (is.good()) {
    char c = is.get();
    if (c == ',') {
      if (!brStr.empty()) {
        np->_nwkItem(brStr);
        brStr.clear();
      }

      np = new Node;
      addChild(np);
    } else if (c == ')') {
      if (!brStr.empty()) {
        np->_nwkItem(brStr);
        brStr.clear();
      }
      break;
    } else if (c == '(') {
      np->_innwk(is);
    } else if (c != '"' && c != '\n' && c != '\t' && c != '\'') {
      brStr.push_back(c);
    }
  }
};

void Node::innwk(istream &is) {

  char c = is.get();
  string brStr;
  while (c != ';') {
    if (c == '(') {
      _innwk(is);
    } else if (c != '"' && c != '\n' && c != '\t' && c != '\'') {
      brStr.push_back(c);
    }

    if (is.good()) {
      c = is.get();
    } else {
      cerr << "some wrong in the nwk file" << endl;
      exit(1);
    }
  }
};

void Node::renewId(const unordered_map<string, size_t> &mgi) {
  if (isLeaf()) {
    id = mgi.find(name)->second;
  } else {
    for (auto &nd : children)
      (*nd).renewId(mgi);
  }
};
