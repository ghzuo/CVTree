/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-05 Tuesday 15:11:38
 */

#include "cv.h"

LPnorm CVvec::lp = L2;

CVvec::CVvec(CVmap &cm) {
  // set vec
  vec = extractAndSort(cm);
  // set norm
  setNormL2();
};

CVvec::CVvec(const vector<Kstr> &vk, double val) {
  vec.reserve(vk.size());
  for (const Kstr &k : vk)
    vec.push_back(CVdim(k, val));
}

CVmap CVvec::toMap() {
  CVmap cm(vec.begin(), vec.end());
  return cm;
}

vector<Kstr> CVvec::toVK() {
  vector<Kstr> vk;
  vk.reserve(vec.size());
  for (auto &c : vec)
    vk.push_back(c.first);
  return vk;
}

void CVvec::normalize() {
  for (auto &cd : vec)
    cd.second /= norm;
}

void CVvec::setNormL0() { norm = double(vec.size()); }

void CVvec::setNormL1() {
  norm = 0.0;
  for (auto &item : vec)
    norm += abs(item.second);
}

void CVvec::setNormL2() {
  norm = 0.0;
  for (auto &item : vec)
    norm += item.second * item.second;
  norm = sqrt(norm);
};

// output cv vector file and gzip it
void CVvec::write(const string &file) {
  unsigned long size = vec.size();

  gzFile fp;
  if ((fp = gzopen(file.c_str(), "wb")) == NULL) {
    cerr << "Error happen on write cvfile: " << file << endl;
    exit(1);
  }

  gzwrite(fp, &norm, sizeof(double));
  gzwrite(fp, &size, sizeof(unsigned long));
  gzwrite(fp, vec.data(), size * sizeof(CVdim));
  // for(const CVdim& s : cv)
  // 	gzwrite(fp, &s, sizeof(CVdim));

  // write Letter char set
  gzputs(fp, Letter::decMap.c_str());
  gzclose(fp);
}

// Read the cv from binary gz file
// read the cv from operation on the cv vector
// used for the distance calculate
string CVvec::read(const string &file) {
  gzFile fp;
  if ((fp = gzopen(file.c_str(), "rb")) == NULL) {
    cerr << "CV file not found: \"" << file << '"' << endl;
    exit(1);
  }

  // get the head (norm & size) of the cv file
  pair<double, size_t> tmp;
  long len = gzread(fp, (char *)&tmp, sizeof(CVdim));
  if (len != sizeof(CVdim)) {
    cerr << "Failed to read cv file " << file << endl;
    exit(10);
  }
  norm = sqrt(tmp.first);
  size_t size = tmp.second;

  // read the cv
  vec.resize(size);
  gzread(fp, (char *)vec.data(), size * sizeof(CVdim));

  // read char set
  string charSet;
  gzline(fp, charSet);
  gzclose(fp);

  return charSet;
}

ostream &operator<<(ostream &os, const CVvec &cv) {
  os << "CVvec: " << cv.size() << " dims, norm: " << cv.norm << endl;
  for (const CVdim &c : cv.vec)
    os << c.first << ':' << c.second << '\n';
  return os;
}

size_t cvsize(const string &filename) {
  gzFile fp;
  if ((fp = gzopen(filename.c_str(), "rb")) == NULL) {
    cerr << "CV file not found: \"" << filename << '"' << endl;
    exit(1);
  }

  double inner;
  gzread(fp, (char *)&inner, sizeof(double));
  size_t size;
  gzread(fp, (char *)&size, sizeof(size_t));
  gzclose(fp);

  return size * sizeof(CVdim);
};
