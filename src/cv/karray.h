/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-12-03 22:14:17
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-20 Wednesday 20:47:14
 */

#ifndef KARRAY_H
#define KARRAY_H

#include "cv.h"
//#include "string.h"

using namespace std;

// basic items of CV brick
struct Kitem {
  unsigned int index = 0;
  float value = NAN;

  Kitem() = default;
  Kitem(unsigned ndx, float val) : index(ndx), value(val){};
  Kitem(pair<unsigned, float> p) : index(p.first), value(p.second){};
  bool operator<(const Kitem &rh) const { return index < rh.index; };
  bool operator<(unsigned &ndx) const { return index < ndx; };
  friend ostream &operator<<(ostream &, const Kitem &);
};

// the infomation of Kstr Dimension
struct KdimInfo {
  Kstr kstr;
  pair<size_t, size_t> index;

  KdimInfo() = default;
  KdimInfo(const KdimInfo &rhs) : kstr(rhs.kstr), index(rhs.index){};
  KdimInfo(const Kstr &ks) : kstr(ks), index(0, 0){};
  KdimInfo(const Kstr &ks, size_t i, size_t j) : kstr(ks), index(i, j){};

  bool operator==(const KdimInfo &rh) const { return kstr == rh.kstr; };
  bool operator<(const KdimInfo &rh) const { return kstr < rh.kstr; };

  friend ostream &operator<<(ostream &, const KdimInfo &);
  friend istream &operator>>(istream &, KdimInfo &);
};

// the infomation of CV dimension
struct CVdimInfo {
  double len = NAN;
  double lasso = NAN;
  double norm = NAN;

  CVdimInfo() = default;
  CVdimInfo(const CVdimInfo &rhs)
      : len(rhs.len), lasso(rhs.lasso), norm(rhs.norm){};
  CVdimInfo(double v) : len(v), lasso(v), norm(v){};
  CVdimInfo(const CVvec &);
  CVdimInfo(const CVmap &);

  friend ostream &operator<<(ostream &, const CVdimInfo &);
};

#endif // !KARRAY_H
       // end of karray.h