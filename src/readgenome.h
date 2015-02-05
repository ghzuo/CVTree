#ifndef READGENOME_H
#define READGENOME_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <iomanip>

#include "global.h"
#include "stringOpt.h"
using namespace std;

typedef string Gene;
typedef vector<Gene> Genome;

struct GeneType{
    char mc[128];
    vector<char> letters;
    
    GeneType(const string&);

    void aainit();
    void nainit();

    size_t readgene(string&, Genome&) const;
    void checkgene(string&) const;
};

#endif
