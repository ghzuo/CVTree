#ifndef KSTRING_H
#define KSTRING_H

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <math.h> 
#include <zlib.h>
#include <unordered_map>

#include "global.h"
#include "readgenome.h"

struct Kstr{
    static vector<char> charSet;
    unsigned long ks;

    Kstr();
    Kstr(const string&);
    Kstr(unsigned long);
    static int init(const vector<char>&);

    string decode() const;
    size_t length() const;

    void append(char);
    void addhead(char);

    void behead();
    void choptail();

    void forward(char);
    void backward(char);

    bool operator<(const Kstr&) const;
    bool operator>(const Kstr&) const;
    bool operator==(const Kstr&) const;
    bool operator!=(const Kstr&) const;
    bool operator()(const Kstr&, const string&) const;

    friend ostream& operator<<(ostream&, const Kstr&);

private:
    static char cmap[128];
    static size_t nbase;
};

struct Kstr_Hash{
    size_t operator()(const Kstr& r) const{
	return size_t(r.ks);
    };
};

typedef pair<Kstr,double> CVdim;
typedef unordered_map<Kstr,double,Kstr_Hash> CVmap;
typedef vector<CVdim> CVvec;

void writecv(const CVmap&, const string&);
void writecv(const CVvec&, const string&);
double readcv(const string&, CVvec&);
size_t cvsize(const string&);

double module(const CVvec&);
void normalize(CVvec&);
double dist(const CVvec&, const CVvec&);

void readvk(const string&, vector<Kstr>&);
void writevk(const string&, const vector<Kstr>&);

#endif
