#ifndef DISTMATRIX_H
#define DISTMATRIX_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <netcdfcpp.h>

#include "tree.h"
#include "global.h"

using namespace std;


// ... only distance matrix
class Mdist{
    int ng;
    vector<double> dist;
    vector<string> name;

    void _setdist(size_t, size_t, double);
    double _getdist(size_t, size_t) const;

public:
    Mdist();

    // read the matrix
    void readmtx(const string&);
    void readmtxnc(const string&);
    void readmtxtxt(const string&);

    // output the matrix
    void writemtx(const string&);
    void writemtxnc(const string&);
    void writemtxtxt(const string&);

    // reduce the matrix by name or index list
    void reduce(const vector<size_t>&);
    void reduce(const vector<string>&);

    // extent the matrix
    void extend(const vector<string>&, const vector<double>&);
    void extend(const vector<string>&);
    void extend(const vector<double>&);

    // get/set value of matrix
    double getdist(size_t, size_t) const;
    string getname(size_t) const;
    void setdist(size_t, size_t, double);
    void setname(size_t, const string&);

    // global options
    void resize(size_t);
    size_t size() const;
    size_t msize() const;
    size_t capacity() const;
};

#endif

#ifndef NCERR
#define NCERR
static const int NC_ERR = 10;
#endif

