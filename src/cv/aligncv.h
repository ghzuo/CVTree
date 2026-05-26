/*
 * Copyright (c) 2025
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2025-12-08 Monday 10:49:30
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-12 Friday 16:26:04
 */

#ifndef ALIGNCV_H
#define ALIGNCV_H

#include "cv.h"

// CVblock: a segment of the CVvec
typedef vector<CVdim>::const_iterator CViter;
struct CVblock {
  CViter begin, end;

  CVblock(const CViter &a, const CViter &b) : begin(a), end(b){};

  bool empty() const { return begin == end; };
  long length() const { return distance(begin, end); };
  CViter front() const { return begin; };
  CViter last() const { return end - 1; };
  CViter mid() const { return begin + length() / 2; };

  CViter getLower(const CViter &obj) const {
    return lower_bound(begin, end, *obj, [](const CVdim &a, const CVdim &b) {
      return a.first < b.first;
    });
  };

  CViter getUpper(const CViter &obj) const {
    return upper_bound(begin, end, *obj, [](const CVdim &a, const CVdim &b) {
      return a.first < b.first;
    });
  };

  bool pop() {
    ++begin;
    return empty();
  };

  bool resetBegin(CViter &obj) {
    begin = getLower(obj);
    return empty();
  };

  bool resetEnd(CViter &obj) {
    end = getUpper(obj);
    return empty();
  };
};

bool fitBoundary(CVblock &, CVblock &);
bool fitBegin(CVblock &, CVblock &);
bool fitEnd(CVblock &, CVblock &);

double align(CVblock &, CVblock &);
size_t nInterSection(CVblock &, CVblock &);
double overlap(CVblock &, CVblock &);

#endif // ALIGNCV_H
