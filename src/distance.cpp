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

#include "distance.h"
extern Info theInfo;

size_t IterStep::ndx = 0;
void IterStep::reIndex() { ndx = 0; }

////// for IterStep
IterStep::IterStep(const Mdist &dm, const vector<string> &flist) : size(0.0) {
  for (size_t i = 0; i < dm.size(); ++i) {
    size_t n = dm.nNAN(i);
    if (n > 0) {
      CVitem it(i, n);
      cvlist.emplace_back(it);
    }
  }

  // set the file name
  for (auto &it : cvlist) {
    it.fname = flist[it.ndx];
  }

  ++ndx;
};

void IterStep::checkSize(float maxM, const Mdist &dm) {

  // sort by the number of NAN
  sort(cvlist.begin(), cvlist.end());

  for (size_t i = 0; i < cvlist.size(); ++i) {
    if (size < maxM) {
      size_t j = i + 1;
      for (; j < cvlist.size(); ++j) {
        if (dm.isNAN(cvlist[i].ndx, cvlist[j].ndx)) {
          size += cvsize(cvlist[i].fname);
          introBlock.emplace_back(&cvlist[i]);
          break;
        }
      }
      if (j == cvlist.size()) {
        interBlock.emplace_back(&cvlist[i]);
      }
    } else {
      interBlock.emplace_back(&cvlist[i]);
    }
  }
}

size_t IterStep::length() const { return introBlock.size(); };

void IterStep::fillBlock() {
#pragma omp parallel for
  for (auto i = 0; i < introBlock.size(); ++i) {
    introBlock[i]->fill();
  }
  return;
};

void IterStep::calcInDist(Mdist &dm, Method *meth) {
// for the intro-distances between the genomes
#pragma omp parallel
#pragma omp single
  {
    for (auto i = 0; i < introBlock.size(); ++i) {
      for (auto j = i + 1; j < introBlock.size(); ++j) {
#pragma omp task
        if (dm.isNAN(introBlock[i]->ndx, introBlock[j]->ndx)) {
          dm.setdist(introBlock[i]->ndx, introBlock[j]->ndx,
                     meth->dist(introBlock[i]->cv, introBlock[j]->cv));
        }
      }
    }
  }
  return;
};

void IterStep::calcOutDist(Mdist &dm, Method *meth) {
// for the intro-distances between the genomes
#pragma omp parallel for
  for (auto i = 0; i < interBlock.size(); ++i) {
    // read an orginal cv
    interBlock[i]->fill();

    // obtain the distances between the genome and the extend genomes
    for (auto j = 0; j < introBlock.size(); ++j) {
      if (dm.isNAN(interBlock[i]->ndx, introBlock[j]->ndx)) {
        dm.setdist(interBlock[i]->ndx, introBlock[j]->ndx,
                   meth->dist(interBlock[i]->cv, introBlock[j]->cv));
      }
    }
    interBlock[i]->clear();
  }
  return;
};

void IterStep::execute(Mdist &dm, Method *meth) {
  // output step info
  theInfo.output(info(), 1);

  // read the extend cv of the step
  fillBlock();
  theInfo.output("Complete read CVs", 1);

  // calculate the inner distances
  calcInDist(dm, meth);
  theInfo.output("Complete calculate of intro-block");

  // calculate the outer distances
  if (!interBlock.empty()) {
    calcOutDist(dm, meth);
    theInfo.output("Complete calculate of inter-block", -1);
  } else {
    theInfo.output("No inter-block calculate are required", -1);
  }

  // output complete infomation
  theInfo.output("Complete Distance Calculate of the Block", -1);
}

void assignDM(const string &refdm, bool netcdf, Mdist &dm) {
  // assign distance by reference DM
  vector<string> dmlist;
  separateWord(dmlist, refdm);
  for (auto &fnm : dmlist) {
    Mdist xdm;
    bool readxdm;
#pragma omp critical
    { readxdm = xdm.readmtx(fnm, netcdf); }
    if (readxdm) {
      xdm.cleanName();
      dm.assignDM(xdm);
      if (!dm.hasNAN())
        return;
    }
  }
};

string IterStep::info() {
  string str = "Start the calculate block " + to_string(ndx);
  str += "\n" + to_string(introBlock.size()) +
         " CVs will be loaded, which size is " + to_string(size / 1073741824) +
         "G";
  return str;
}