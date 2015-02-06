#include "tree.h"
const size_t N_FORKS(4);

/*******************************************************************/
/********* Member Functions For Node Class *************************/
/*******************************************************************/
Node::Node():name(""),id(0),length(NAN),parent(NULL){
    children.reserve(N_FORKS);
};

Node::Node(size_t n):name(""),id(n),length(NAN),parent(NULL){
    children.reserve(N_FORKS);
};

Node::Node(size_t n, const string& str):name(str),id(n),length(NAN),parent(NULL){
    children.reserve(N_FORKS);
};

Node::Node(size_t n, const vector<Node*>& vn)
    :name(""),id(n),length(NAN), parent(NULL), children(vn){
        children.reserve(N_FORKS);
};

bool Node::isLeaf(){
    return children.empty();
};

void Node::addChild(Node* nd){
    children.emplace_back(nd);
    (*nd).parent = this;
};

void Node::deleteChild(Node* nd){
    Children::iterator iter = find(children.begin(),children.end(), nd);
    children.erase(iter);
};

void Node::clear(){
    vector<Node*> nodes;
    getDescendants(nodes);
    vector<Node*>::iterator iter=nodes.begin();
    vector<Node*>::iterator iterEnd=nodes.end();
    for(;iter!=iterEnd; ++iter)
	delete *iter;
    children.clear();
}

void Node::getDescendants(vector<Node*>& nds){
    if(!isLeaf()){
	for(const auto &nd : children){
	    nds.emplace_back(nd);
	    (*nd).getDescendants(nds);
	}
    }
};

void Node::getLeafs(vector<Node*>& nodes){
    if(isLeaf()){
	nodes.emplace_back(this);
    }else{
	for(const auto &nd : children)
	    (*nd).getLeafs(nodes);
    }
};

Node* Node::resetroot(const string& str){
    vector<Node*> leafs;
    getLeafs(leafs);
    Node* og = NULL;
    for(Node* nd : leafs)
	if((*nd).name == str){
	    og = nd; break;
	}

    if(og == NULL)
	return og;
    return resetroot(og);
}

Node* Node::resetroot(Node* np){
    Node* outgrp = np;
    vector<Node*> nlist;
    while((*np).parent != NULL){
	nlist.emplace_back((*np).parent);
	np = (*np).parent;
    }

    for(vector<Node*>::reverse_iterator iter = nlist.rbegin(); 
	iter != nlist.rend()-1; ++iter){
	vector<Node*>::reverse_iterator next = iter+1;
	(**iter).length = (**next).length;
	(**next).length = NAN;
	(**iter).deleteChild(*next);
	(**next).parent = NULL;
	(**next).addChild(*iter);
    }
    Node* theRoot = nlist.front();

    // set the outgroup as the last child
    Children::reverse_iterator iter = (*theRoot).children.rbegin();
    if((*iter) != outgrp){
        for( ++iter; iter != (*theRoot).children.rend(); ++iter){
            if((*iter)==outgrp){
                *iter = (*theRoot).children.back();
                (*theRoot).children.back() = outgrp;
            }
        }
    }
    
    return theRoot;
};

void Node::_outnwk(ostream& os){
    if(!isLeaf()){
	os << "(";
	Children::const_iterator iter = children.begin();
	(*(*iter))._outnwk(os);
	for(++iter; iter!=children.end(); ++iter){
	    os << ",";
	    (*(*iter))._outnwk(os);
	}
	os << ")";
    }

    if(name.find(' ')==std::string::npos)
	os << name ;
    else
	os << '"' << name << '"';

    if(!std::isnan(length))
	os << ":" << fixed << setprecision(5) << length;
};

void Node::outnwk(ostream& os){
    _outnwk(os);
    os << ";" << endl;
};
