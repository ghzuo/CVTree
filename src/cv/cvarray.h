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
 * @Last Modified Time: 2026-05-21 Thursday 16:38:37
 */

#ifndef CVARRAY_H
#define CVARRAY_H

#include <tuple>
#include <unordered_map>
#include <vector>

#include "cvbundle.h"
#include "karray.h"
#include "kit.h"

using namespace std;

struct Kblock {
  vector<Kitem>::const_iterator _begin;
  vector<Kitem>::const_iterator _end;

  Kblock() = default;
  Kblock(vector<Kitem>::const_iterator b, vector<Kitem>::const_iterator e)
      : _begin(b), _end(e){};

  vector<Kitem>::const_iterator begin() const { return _begin; };
  vector<Kitem>::const_iterator end() const { return _end; };
  size_t size() const { return _end - _begin; }
};

struct CVAinfo {
  size_t nCV;
  size_t nKstr;
  size_t nItem;

  CVAinfo() = default;
  CVAinfo(size_t nc, size_t nk, size_t ni) : nCV(nc), nKstr(nk), nItem(ni){};
  CVAinfo(const CVAinfo &rhs)
      : nCV(rhs.nCV), nKstr(rhs.nKstr), nItem(rhs.nItem){};
  CVAinfo(const string &fn) { read(fn); };
  void read(const string &);
  friend ostream &operator<<(ostream &, const CVAinfo &);
};

struct CVArray {
  vector<KdimInfo> kdi;
  vector<CVdimInfo> cvdi;
  vector<Kitem> data;
  vector<float> norm;

  CVArray() = default;
  CVArray(const CVArray &rhs)
      : kdi(rhs.kdi), cvdi(rhs.cvdi), data(rhs.data), norm(rhs.norm){};
  CVArray(const vector<CVmap> &cvs) { set(cvs); };
  CVArray(const CVBundle &cvs) { set(cvs); };
  CVArray(const string &fname) : CVArray(fname, CVvec::lp){};
  CVArray(const string &fname, enum LPnorm normType) {
    read(fname);
    setNorm(normType);
  }

  // void set(const CVBundle &);
  // void set(const vector<CVmap> &);

  template <typename T> void set(const T &cvs) {
    // get the cvdiminfo
    cvdi.reserve(cvs.size());
    for (const auto &cv : cvs) {
      cvdi.emplace_back(CVdimInfo(cv));
    }

    // align kstr into map
    size_t nItem = 0;
    unordered_map<Kstr, vector<Kitem>, Kstr_Hash> khash;
    for (auto i = 0; i < cvs.size(); ++i) {
      for (const auto &cd : cvs[i]) {
        khash[cd.first].emplace_back(Kitem(i, cd.second));
        nItem++;
      }
    }
    vector<pair<Kstr, vector<Kitem>>> kvec = extractAndSort(khash);

    // set CVArray based on map
    kdi.reserve(kvec.size());
    data.reserve(nItem);
    for (auto &kar : kvec) {
      size_t ibeg = data.size();
      size_t iend = ibeg + kar.second.size();
      kdi.emplace_back(KdimInfo(kar.first, ibeg, iend));
      data.insert(data.end(), kar.second.begin(), kar.second.end());
    }
  }

  void setNorm();
  void setNorm(enum LPnorm);
  void calcNorm();
  void calcNorm(enum LPnorm);
  Kblock getKblock(size_t) const;

  void read(const string &);
  void write(const string &) const;
  void calcCVdimInfo();

  size_t cvSize() const { return norm.size() == 0 ? cvdi.size() : norm.size(); }
  size_t kstrSize() const { return kdi.size(); }
  size_t itemSize() const { return data.size(); }
  string info() const;

  friend ostream &operator<<(ostream &, const CVArray &);
};

template <typename T>
void alignSortVector(const vector<T> &va, const vector<T> &vb,
                     vector<pair<size_t, size_t>> &aln) {
  auto itb = vb.begin();
  for (auto ita = va.begin(); ita != va.end(); ++ita) {
    itb = lower_bound(itb, vb.end(), *ita);
    if (*ita == *itb) {
      aln.emplace_back(ita - va.begin(), itb - vb.begin());
      ++itb;
    }
  }
}
#endif // !CVARRAY_H
