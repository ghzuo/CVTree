/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-03-17 15:39:23
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 22:03:36
 */

#include "tree.h"
const size_t N_FORKS(4);

/*******************************************************************/
/********* Member Functions For Node Class *************************/
/*******************************************************************/
Node::Node()
    : name(""), id(0), length(NAN), parent(NULL), taxSize(0), taxLevel(0),
      nleaf(0), nxleaf(0), unclassified(false), uploaded(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n)
    : name(""), id(n), length(NAN), parent(NULL), taxSize(0), taxLevel(0),
      nleaf(0), nxleaf(0), unclassified(false), uploaded(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n, const string &str)
    : name(str), id(n), length(NAN), parent(NULL), taxSize(0), taxLevel(0),
      nleaf(0), nxleaf(0), unclassified(false), uploaded(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n, const vector<Node *> &vn)
    : name(""), id(n), length(NAN), parent(NULL), children(vn), taxSize(0),
      nleaf(0), nxleaf(0), unclassified(false), uploaded(false) {
  children.reserve(2);
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
  for (Node *nd : leafs)
    if ((*nd).name == str) {
      og = nd;
      break;
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

  if (name.find(' ') == std::string::npos)
    os << name;
  else
    os << '"' << name << '"';

  if (!std::isnan(length))
    os << ":" << fixed << setprecision(5) << length;
};

void Node::outnwk(ostream &os) {
  _outnwk(os);
  os << ";" << endl;
};

void Node::_innwk(istream &is) {
  Node *np = new Node;
  addChild(np);
  string *sp = &((*np).name);
  string nstr;

  while (is.good()) {
    char c = is.get();
    if (c == ',') {
      if (!nstr.empty()) {
        (*np).length = str2double(nstr);
        nstr.clear();
      }

      np = new Node;
      sp = &((*np).name);
      addChild(np);
    } else if (c == ')') {
      if (!nstr.empty()) {
        (*np).length = str2double(nstr);
        nstr.clear();
      }
      break;
    } else if (c == '(') {
      (*np)._innwk(is);
    } else if (c == ':') {
      sp = &nstr;
    } else if (c != '"' && c != '\n' && c != '\t' && c != '\'') {
      (*sp).push_back(c);
    }
  }
};

void Node::innwk(istream &is) {

  char c = is.get();
  while (c != ';') {
    if (c == '(') {
      _innwk(is);
    } else if (c != '"' && c != '\n' && c != '\t' && c != '\'') {
      name.push_back(c);
    }

    if (is.good()) {
      c = is.get();
    } else {
      cerr << "some wrong in the nwk file" << endl;
      exit(1);
    }
  }
};

void Node::checkUnclassified() {
  if (isLeaf()) {
    if (name.find(">Unclassified<") != string::npos) {
      nxleaf = 1;
    } else {
      nleaf = 1;
    }
  } else {
    for (const auto &nd : children) {
      (*nd).checkUnclassified();
      nxleaf += (*nd).nxleaf;
      nleaf += (*nd).nleaf;
    }
  }

  if (nleaf == 0)
    unclassified = true;
};

void Node::checkUploaded() {
  if (isLeaf()) {
    if (name.find(".UPLOAD") != string::npos) {
      nUpLeaf = 1;
      uploaded = true;
    }
  } else {
    bool isUpload = true;
    for (const auto &nd : children) {
      (*nd).checkUploaded();
      nUpLeaf += (*nd).nUpLeaf;
      isUpload = isUpload && (*nd).uploaded;
    }
    uploaded = isUpload;
  }
};

void Node::getDefineLeafs(vector<Node *> &nodes) {
  if (isLeaf()) {
    if (!unclassified)
      nodes.emplace_back(this);
  } else {
    for (const auto &nd : children)
      (*nd).getDefineLeafs(nodes);
  }
};

void Node::injson(const string& file){
   ifstream ijson(file.c_str());
   if(! ijson){
      cerr << " Cannot found the input file " << file << endl;
      exit(4);
   }
   injson(ijson);
}

void Node::injson(istream &is) {
  is.get();
  _injson(is);
};

void Node::_injson(istream &is) {
  char c = is.get();
  bool isValue(false);
  string key, value;
  while (c != '}') {
    if (c == '"') {
      if (isValue)
        _getStr(is, value);
      else
        _getStr(is, key);
    } else if (c == ':') {
      isValue = true;
    } else if (c == '{') {
      Node *nd = new Node;
      addChild(nd);
      (*nd)._injson(is);
    } else if (c == ',' && isValue) {
      _getKeyValue(key, value);
      isValue = false;
    }
    c = is.get();
  }
  if (isValue)
    _getKeyValue(key, value);
};

void Node::_getKeyValue(string &key, string &value) {
  if (key.compare("name") == 0) {

    // find the arch
    size_t lpos = value.find('{');
    size_t rpos = value.find('}');
    size_t spos = value.find('/');
    size_t ppos = value.find('+');

    // get the name
    name = value.substr(0, lpos);

    // get the unclassified number
    if (ppos == string::npos) {
      nxleaf = 0;
    } else {
      nxleaf = str2int(value.substr(ppos + 1, rpos - ppos - 1));
      rpos = ppos;
    }

    // get the classified number and taxonomy size
    if (spos == string::npos) {
      taxSize = str2int(value.substr(lpos + 1, rpos - lpos - 1));
      nleaf = taxSize;
    } else {
      nleaf = str2int(value.substr(lpos + 1, spos - lpos - 1));
      taxSize = str2int(value.substr(spos + 1, rpos - spos - 1));
    }
  } else if (key.compare("length") == 0) {
    length = str2double(value);
  } else if (key.compare("upload") == 0) {
    uploaded = true;
  } else if (key.compare("unclassified")) {
    unclassified = true;
  }
};

void Node::_getStr(istream &is, string &str) {
  str.clear();
  char c = is.get();
  while (c != '"') {
    str += c;
    c = is.get();
  }
};

void Node::outjson(ostream &os) {

  os << "{";
  if (!name.empty()) {
    os << '"' << "name"
       << "\":\"" << name << "{" << nleaf;
    if (nleaf > 0 && nleaf < taxSize)
      os << "/" << taxSize;
    if (nxleaf != 0) {
      os << "+" << nxleaf;
      // if(nleaf==0 && nxleaf < taxSize)
      // 	os << "/" << taxSize;
    }
    os << "}\"";
  }

  if (uploaded)
    os << ',' << '"' << "upload"
       << "\":\""
       << "true" << '"';

  if (unclassified)
    os << ',' << '"' << "unclassified"
       << "\":\""
       << "true" << '"';

  if (nleaf == taxSize) {
    os << ',' << '"' << "ntype"
       << "\":\""
       << "Coincide" << '"';
  } else if (nleaf == 0 && nxleaf == taxSize) {
    os << ',' << '"' << "ntype"
       << "\":\""
       << "Coincide" << '"';
  }

  if (!std::isnan(length))
    os << ',' << '"' << "length"
       << "\":\"" << fixed << setprecision(5) << length << '"';

  if (!isLeaf()) {
    os << ",\"children\":[";

    Children::const_iterator iter = children.begin();
    (*(*iter)).outjson(os);
    for (++iter; iter != children.end(); ++iter) {
      os << ",";
      (*(*iter)).outjson(os);
    }
    os << "]";
  }
  os << "}";
};

void Node::outjsonAbbr(ostream &os) {

  os << "{";
  if (!name.empty()) {
    os << '"' << "n"
       << "\":\"" << name << "{" << nleaf;
    if (nleaf < taxSize)
      os << "/" << taxSize;
    if (nxleaf != 0) {
      os << "+" << nxleaf;
      // if(nleaf==0 && nxleaf < taxSize)
      // 	os << "/" << taxSize;
    }
    os << "}\"";
  }

  if (uploaded)
    os << ',' << '"' << "u"
       << "\":\""
       << "1" << '"';

  if (unclassified)
    os << ',' << '"' << "x"
       << "\":\""
       << "1" << '"';

  if (nleaf == taxSize)
    os << ',' << '"' << "t"
       << "\":\""
       << "C" << '"';
  else if (nleaf == 0 && nxleaf == taxSize) {
    os << ',' << '"' << "t"
       << "\":\""
       << "C" << '"';
  }

  // if(!isnan(length))
  // 	os << ',' << '"' << "l" << "\":\"" << fixed << setprecision(5) << length
  // <<'"';

  if (!isLeaf()) {
    os << ",\"c\":[";

    Children::const_iterator iter = children.begin();
    (*(*iter)).outjsonAbbr(os);
    for (++iter; iter != children.end(); ++iter) {
      os << ",";
      (*(*iter)).outjsonAbbr(os);
    }
    os << "]";
  }
  os << "}";
};

void Node::renewId(const unordered_map<string, size_t> &mgi) {
  if (isLeaf()) {
    id = mgi.find(name)->second;
  } else {
    for (auto &nd : children)
      (*nd).renewId(mgi);
  }
};

string Node::getStrainName() {
  return name.substr(name.find_last_of('>') + 1);
};

void Node::outPrediction(ostream &os) {
  vector<Node *> leafs;
  getLeafs(leafs);
  for (auto &nd : leafs) {
    if ((*nd).unclassified) {
      string p;
      (*nd)._getPrediction(p);
      os << (*nd).getStrainName() << "\t" << p << endl;
    }
  }
};

void Node::_getPrediction(string &p) {
  if (parent == NULL) {
    p = name;
  } else {
    if ((*parent).nleaf == (*parent).taxSize) {
      p = (*parent).name;
    } else {
      (*parent)._getPrediction(p);
    }
  }
};

void Node::chkLeafsName(vector<string> &defLeafs, vector<string> &undefLeafs,
                        bool quite) {
  vector<Node *> leafs;
  getLeafs(leafs);
  set<string> nameSet;
  for (auto &leaf : leafs) {
    if (nameSet.find((*leaf).name) == nameSet.end()) {
      nameSet.insert((*leaf).name);
    } else {
      if (!quite)
        cerr << "repeat strain: " << name << endl;
      (*leaf).name += to_string(nameSet.size());
      nameSet.insert((*leaf).name);
    }
    if ((*leaf).unclassified) {
      (*leaf).id = undefLeafs.size();
      undefLeafs.emplace_back((*leaf).name);
    } else {
      (*leaf).id = defLeafs.size();
      defLeafs.emplace_back((*leaf).name);
    }
  }
}

void Node::reinitTree() {
  if (isLeaf()) {
    name.erase(name.find('|'), 1);
  } else {
    name = "";
    for (const auto &nd : children)
      (*nd).reinitTree();
  }

  nxleaf = 0;
  nleaf = 0;
  taxSize = 0;
  unclassified = false;
}

void Node::getUndefineLeafs(vector<Node *> &nodes) {
  if (children.empty()) {
    if (unclassified)
      nodes.emplace_back(this);
  } else {
    for (auto nd : children)
      (*nd).getUndefineLeafs(nodes);
  }
};

void Node::getUndefineNames(vector<string> &names) {
  if (children.empty()) {
    if (unclassified)
      names.emplace_back(name);
  } else {
    for (auto nd : children)
      (*nd).getUndefineNames(names);
  }
};
