/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-07-31 15:50:49
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-31 15:50:50
 */

#include "distmeth.h"
extern Info theInfo;

//// set max memory
void DistMeth::setMaxMem(float ms, int ng, int nk) {
  float maxNameLen = 2048.0;
  float giga = 1073741824.0;
  float bs = maxNameLen * ng + nk * ng * (ng - 1) * sizeof(double) / 2 + giga;
  maxM = ms - bs;
};

////// for DistMeth
void DistMeth::init(const vector<string> &flist) {
  cvlist.clear();
  for (size_t i = 0; i < flist.size(); ++i) {
    cvlist.emplace_back(CVitem(i, flist[i]));
  }
};

float DistMeth::setStep(const Mdist &dm) {

  // get the cvitem sort by NAN
  vector<pair<int, CVitem *>> nanItems;
  for (size_t i = 0; i < dm.size(); ++i) {
    int n = dm.nNAN(i);
    if (n > 0) {
      nanItems.emplace_back(make_pair(n, &cvlist[i]));
    }
  }
  sort(nanItems.begin(), nanItems.end());

  // divid the nanItems into two blocks
  float size(0.0);
  interBlock.clear();
  introBlock.clear();
  for (auto iter = nanItems.rbegin(); iter != nanItems.rend(); ++iter) {
    if (size < maxM) {
      auto it = iter + 1;
      for (; it != nanItems.rend(); ++it) {
        if (dm.isNAN(iter->second->ndx, it->second->ndx)) {
          size += cvsize(iter->second->fname);
          introBlock.emplace_back(iter->second);
          break;
        }
      }
      if (it == nanItems.rend()) {
        interBlock.emplace_back(iter->second);
      }
    } else {
      interBlock.emplace_back(iter->second);
    }
  }

#if defined(_OPENMP)
  if (interBlock.size() < omp_get_max_threads()) {
    introBlock.insert(introBlock.begin(), interBlock.begin(), interBlock.end());
    vector<CVitem *>().swap(interBlock);
  }
#endif

  return size;
}

void DistMeth::cleanStep() {
  for (auto &item : introBlock) {
    item->clear();
  }
}

size_t DistMeth::length() const { return introBlock.size(); };

void DistMeth::fillBlock() {
#pragma omp parallel for
  for (auto i = 0; i < introBlock.size(); ++i) {
    introBlock[i]->fill();
  }
  return;
};

void DistMeth::calcInDist(Mdist &dm) {
// for the intro-distances between the genomes
#pragma omp parallel
#pragma omp single
  {
    for (auto i = 0; i < introBlock.size(); ++i) {
      for (auto j = i + 1; j < introBlock.size(); ++j) {
#pragma omp task
        if (dm.isNAN(introBlock[i]->ndx, introBlock[j]->ndx)) {
          dm.setdist(introBlock[i]->ndx, introBlock[j]->ndx,
                     dist(introBlock[i]->cv, introBlock[j]->cv));
        }
      }
    }
  }
  return;
};

void DistMeth::calcOutDist(Mdist &dm) {
// for the intro-distances between the genomes
#pragma omp parallel for
  for (auto i = 0; i < interBlock.size(); ++i) {
    // read an orginal cv
    interBlock[i]->fill();

    // obtain the distances between the genome and the extend genomes
    for (auto j = 0; j < introBlock.size(); ++j) {
      if (dm.isNAN(interBlock[i]->ndx, introBlock[j]->ndx)) {
        dm.setdist(interBlock[i]->ndx, introBlock[j]->ndx,
                   dist(interBlock[i]->cv, introBlock[j]->cv));
      }
    }
    interBlock[i]->clear();
  }
  return;
};

void DistMeth::execute(const vector<string> &flist, Mdist &dm) {

  init(flist);
  int ndx(0);

  do {
    // check the memory and set steps
    theInfo(infoStep(++ndx, setStep(dm)), 1);

    // read the extend cv of the step
    fillBlock();
    theInfo("Complete read CVs", 1);

    // calculate the inner distances
    calcInDist(dm);
    theInfo("Complete calculate of intro-block");

    // calculate the outer distances
    if (!interBlock.empty()) {
      calcOutDist(dm);
      theInfo("Complete calculate of inter-block", -1);
    } else {
      theInfo("No inter-block calculate are required", -1);
    }
    // clean the cvs
    cleanStep();

    // output complete infomation
    theInfo("Complete Distance Calculate of the Block", -1);

  } while (dm.hasNAN());
}

string DistMeth::infoStep(int ndx, float size) {
  string str = "Start the calculate step " + to_string(ndx);
  str += "\n" + to_string(introBlock.size()) + "/" +
         to_string(introBlock.size() + interBlock.size()) +
         " CVs will be resident in memory,\nwhich size is " +
         to_string(size / 1073741824) + "G";
  return str;
}

double Cosine::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  // get the distance after reset bound with binary search
  if (fitBegin(block1, block2)) {
    return 0.5 * (1.0 - align(block1, block2));
  } else {
    return 1.0;
  }

  // for binary
  // if(fitBoundary(block1,block2)){
  //     return 0.5*(1.0 - binaryAlign(block1, block2));
  //     //return 0.5*(1.0 - shrink(block1, block2));
  // }else{
  //     return 1.0;
  // }
};

double InterSet::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  // get the distance after reset bound with binary search
  if (fitBegin(block1, block2)) {
    double n(0);
    for (;;) {
      if (block1.begin->first == block2.begin->first) {
        n++;
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
    return 1.0 - n / sqrt(double(cv1.size()) * double(cv2.size()));

  } else {
    return 1.0;
  }
}

double InterList::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  // get the distance after reset bound with binary search
  if (fitBegin(block1, block2)) {
    double n(0);
    double sum(0);
    for (;;) {
      if (block1.begin->first == block2.begin->first) {
        n += min(block1.begin->second, block1.begin->second);
        sum += block1.begin->second;
        sum += block2.begin->second;
        if (block1.pop()) {
          while (!block2.pop()) {
            sum += block2.begin->second;
          }
          break;
        }
        if (block2.pop()) {
          while (!block1.pop()) {
            sum += block1.begin->second;
          }
          break;
        }
      } else {
        if (block1.begin->first < block2.begin->first) {
          sum += block1.begin->second;
          if (block1.pop()) {
            while (!block2.pop()) {
              sum += block2.begin->second;
            }
            break;
          }
        } else {
          sum += block2.begin->second;
          if (block2.pop()) {
            while (!block1.pop()) {
              sum += block1.begin->second;
            }
            break;
          }
        }
      }
    }
    return 1.0 - 2 * n / sum;
  } else {
    return 1.0;
  }
}