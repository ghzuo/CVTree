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
  sort(cvlist.begin(), cvlist.end());
  npos = cvlist.size();

  for (auto &it : cvlist) {
    it.fname = flist[it.ndx];
  }

  ++ndx;
};

void IterStep::checkSize(float maxM) {
  npos = 0;
  for (; npos < cvlist.size(); ++npos) {
    size += cvsize(cvlist[npos].fname);
    if (size > maxM) {
      ++npos;
      break;
    }
  }
}

size_t IterStep::length() const { return npos; };

void IterStep::fillBlock() {
#pragma omp parallel for
  for (auto i = 0; i < npos; ++i) {
    cvlist[i].fill();
  }
  return;
};

void IterStep::calcInDist(Mdist &dm, Method *meth) {
// for the intro-distances between the genomes
#pragma omp parallel
#pragma omp single
  {
    for (auto i = 0; i < npos; ++i) {
      for (auto j = i + 1; j < npos; ++j) {
#pragma omp task
        if (dm.isNAN(cvlist[i].ndx, cvlist[j].ndx)) {
          dm.setdist(cvlist[i].ndx, cvlist[j].ndx,
                     meth->dist(cvlist[i].cv, cvlist[j].cv));
        }
      }
    }
  }
  return;
};

void IterStep::calcOutDist(Mdist &dm, Method *meth) {
// for the intro-distances between the genomes
#pragma omp parallel for
  for (auto i = npos; i < cvlist.size(); ++i) {
    // read an orginal cv
    cvlist[i].fill();

    // obtain the distances between the genome and the extend genomes
    for (auto j = 0; j < npos; ++j) {
      if (dm.isNAN(cvlist[i].ndx, cvlist[j].ndx)) {
        dm.setdist(cvlist[i].ndx, cvlist[j].ndx,
                   meth->dist(cvlist[i].cv, cvlist[j].cv));
      }
    }
    cvlist[i].clear();
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
  if (npos < cvlist.size()) {
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
  str += "\n" + to_string(npos) + " CVs will be loaded, which size is " +
         to_string(size / 1073741824) + "G";
  return str;
}