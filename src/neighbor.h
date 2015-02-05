#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <limits>
//#include <omp.h>

#include "stringOpt.h"
#include "global.h"
#include "tree.h"
#include "distmatrix.h"
#include "neighborJoint.h"
using namespace std;

// read arguments
struct Args{
    string program, distfile, outfile;
    vector<size_t> splist;
    bool netcdf;
    
    Args(int, char**);
    void usage();
};

// select the genomes to build the tree
void selectLeafs(const Mdist&, const string&, vector<Node*>&);


#endif
