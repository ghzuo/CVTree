#ifndef CVTREE_H
#define CVTREE_H

#include "global.h"
#include "stringOpt.h"
#include "readgenome.h"
#include "kstring.h"

struct Args{
    string program, gtype, outdir, indir;
    vector<pair<string,size_t>> flist;
    set<int> klist;
    set<int> slist;
    int method;
    string suffix;

    Args(int, char**);
    void kcheck(int);
    void usage();
};

size_t count(const Genome&, size_t, CVmap&);
void subtract(const CVmap&, const CVmap&, const CVmap&, double, CVmap&);
bool bySecond(const pair<string,long>&, const pair<string,long>&);

#endif
