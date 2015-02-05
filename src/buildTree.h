#ifndef BUILDTREE_H
#define BUILDTREE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <regex>
#include <unordered_map>

#include "global.h"
#include "stringOpt.h"
#include "distmatrix.h"
#include "memory.h"
#include "neighborJoint.h"

// select the Kstr class
#include "kstring.h"

using namespace std;

// read arguments
struct Args{
    string program, outfile, mtxfile, extdir, orgdir;
    vector<string> suflist, extlist, orglist;
    vector<size_t> orgIndex;
    unordered_map<string,string> taxmap;
    bool outtax, netcdf;
    float memorySize;
    
    Args(int, char**);
    void usage();
};

struct IterStep{
    size_t first, second;
    float size;

    IterStep();
    IterStep(size_t ibeg):first(ibeg),second(ibeg),size(0.0){};
};

void checkCVsize(const Args&, size_t, vector<IterStep>&);

#endif
