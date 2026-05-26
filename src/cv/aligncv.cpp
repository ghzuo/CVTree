/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-12 Friday 16:11:50
 */

#include "aligncv.h"


// Function for CVblock
bool fitBegin(CVblock &block1, CVblock &block2) {

  // fit the begin terminal
  if (block1.begin->first < block2.begin->first) {
    if (block1.resetBegin(block2.begin))
      return false;
  } else if (block1.begin->first > block2.begin->first) {
    if (block2.resetBegin(block1.begin))
      return false;
  }

  return true;
}

bool fitEnd(CVblock &block1, CVblock &block2) {

  // get the upper bound
  CViter iter1Last = block1.last();
  CViter iter2Last = block2.last();
  if (iter1Last->first < iter2Last->first) {
    if (block2.resetEnd(iter1Last))
      return false;
  } else if (iter1Last->first > iter2Last->first) {
    if (block1.resetEnd(iter2Last))
      return false;
  }
  return true;
}

bool fitBoundary(CVblock &block1, CVblock &block2) {
  // get the lower bound
  if (!fitBegin(block1, block2))
    return false;

  // get the upper bound
  if (!fitEnd(block1, block2))
    return false;

  return true;
}


// get dot product by sequentially align
double align(CVblock &block1, CVblock &block2) {

  double d(0);
  for (;;) {
    if (block1.begin->first == block2.begin->first) {
      d += block1.begin->second * block2.begin->second;
      if (block1.pop())
        return d;
      if (block2.pop())
        return d;
    } else {
      if (block1.begin->first < block2.begin->first) {
        if (block1.pop())
          return d;
      } else {
        if (block2.pop())
          return d;
      }
    }
  }
};

// get the number of intersection
size_t nInterSection(CVblock &block1, CVblock &block2) {

  // the number record
  size_t nI(0);

  // get the distance after reset bound with binary search
  if (fitBegin(block1, block2)) {
    for (;;) {
      if (block1.begin->first == block2.begin->first) {
        nI++;
        if (block1.pop())
          break;
        if (block2.pop())
          break;
      } else {
        if (block1.begin->first < block2.begin->first) {
          if (block1.pop())
            break;
        } else {
          if (block2.pop())
            break;
        }
      }
    }
  }

  return nI;
};

// get the overlap of kmer
double overlap(CVblock &block1, CVblock &block2) {

  double sumMin(0);
  if (fitBegin(block1, block2)) {
    for (;;) {
      if (block1.begin->first == block2.begin->first) {
        sumMin += min(block1.begin->second, block1.begin->second);
        if (block1.pop())
          break;
        else if (block2.pop())
          break;
      } else if (block1.begin->first < block2.begin->first) {
        if (block1.pop())
          break;
      } else {
        if (block2.pop())
          break;
      }
    }
  }
  return sumMin;
};

