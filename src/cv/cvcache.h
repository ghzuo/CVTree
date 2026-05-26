/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-11-25 11:34:53
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-16 10:05:55
 */

#ifndef CVCACHE_H
#define CVCACHE_H

#include "cvbrick.h"
#include "cvmeth.h"
#include "distmatrix.h"
#include "kit.h"

// Handle Class for write CV database
typedef map<Kstr, vector<CVKinfo>> KstrIndex;
struct CVCacheHandle {
  size_t nBrick = 0;
  string filename;
  fstream fcv;
  fstream fnm;
  vector<CVGinfo> name;
  KstrIndex kstr;

  CVCacheHandle(const string &fname) : filename(fname){};

  // get CV cache data into cache files
  bool has(const string &gn);
  void fillCV(CVmeth*, const vector<string>&, vector<CVGinfo> &, size_t, bool app = false);
  void insertCV(const CVvec &, CVGinfo, vector<pair<Kstr, CVatom>> &,
                vector<pair<CVKinfo, CVatom>> &);

  // write cv into file
  void writeOneCV(const vector<pair<Kstr, CVatom>> &,
                  const vector<pair<CVKinfo, CVatom>> &);
  void appendCVbrick(const Kstr &, const CVatom &);
  void insertCVatom(const CVKinfo &, const CVatom &);

  // for check cache
  void checkCache();
  void readName(const string &);
  void readCache(const string &);

  // output the data into CV array
  void writecva();
  size_t writeblock(const pair<Kstr, vector<CVKinfo>> &, ostream &);
};

#endif // !CVCACHE_H
