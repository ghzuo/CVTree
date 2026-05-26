/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-14 Sunday 17:40:21
 */

#include "distmeth.h"

/// the create function
DistMeth *DistMeth::create(const string &methStr) {

  // create the distance method
  DistMeth *meth;
  if (methStr == "Cosine") {
    meth = new Cosine();
  } else if (methStr == "Euclidean") {
    meth = new Euclidean();
  } else if (methStr == "InterList") {
    meth = new InterList();
  } else if (methStr == "Min2Max") {
    meth = new Min2Max();
  } else if (methStr == "InterSet") {
    meth = new InterSet();
  } else if (methStr == "Jaccard" || methStr == "ItoU") {
    meth = new ItoU();
  } else if (methStr == "Dice") {
    meth = new Dice();
  } else {
    cerr << "Unknow Distance Method: " << methStr << endl;
    exit(3);
  }

  return meth;
}

/// set max memory
void DistMeth::setMaxMem(float ms, size_t ng, size_t nk) {
  float maxNameLen = 2048.0;
  float giga = 1073741824.0;
  float bs =
      maxNameLen * ng + float(nk) * ng * (ng - 1) * sizeof(double) / 2 + giga;
  maxM = ms - bs;
};

/// for DistMeth
void DistMeth::setflist(const vector<string> &flist) {
  cvlist.clear();
  restlist.clear();
  cvlist.reserve(flist.size());
  restlist.reserve(flist.size());
  for (size_t i = 0; i < flist.size(); ++i) {
    cvlist.emplace_back(CVitem(i, flist[i]));
    restlist.emplace_back(&cvlist[i]);
  }
};

float DistMeth::setStep(const Mdist &dm) {
  // sort restlist by number of NAN
  sort(restlist.begin(), restlist.end(),
       [](auto &a, auto &b) { return a->nNAN < b->nNAN; });

  // chop the large NAN number item to step list
  float size(0.0);
  steplist.clear();
  for (auto iter = restlist.rbegin(); iter != restlist.rend(); ++iter) {
    if (size < maxM) {
      steplist.emplace_back(*iter);
      size += cvsize((*iter)->fname);
    } else {
      break;
    }
  }
  size_t nrest = restlist.size() - steplist.size();
  restlist.resize(nrest);

#ifdef _OPENMP
  if (restlist.size() < omp_get_max_threads()) {
    steplist.insert(steplist.end(), restlist.begin(), restlist.end());
    vector<CVitem *>().swap(restlist);
  }
#endif

  return size;
}

void DistMeth::updateStep(const Mdist &dm) {
  // clean cv in the steplist
  for (auto &item : steplist) {
    item->clear();
  }

  // get the nonNAN cvitems in restlist
  for (auto i = 1; i < restlist.size(); ++i) {
    for (auto j = 0; j < i; ++j) {
      if (dm.isNAN(restlist[i]->ndx, restlist[j]->ndx)) {
        ++restlist[i]->nNAN;
        ++restlist[j]->nNAN;
      }
    }
  }

  // get the nonNAN list to restlist
  vector<CVitem *> nanlist;
  for (auto &ci : restlist)
    if (ci->nNAN > 0)
      nanlist.emplace_back(ci);
  restlist.swap(nanlist);
}

size_t DistMeth::length() const { return steplist.size(); };

void DistMeth::fillBlock() {
#pragma omp parallel for
  for (auto i = 0; i < steplist.size(); ++i) {
    steplist[i]->fill();
    preprocess(steplist[i]->cv);
  }
  return;
};

void DistMeth::setdist4NAN(Mdist &dm, const CVitem &it1, const CVitem &it2) {
  if (dm.isNAN(it1.ndx, it2.ndx)) {
    dm.setdist(it1.ndx, it2.ndx, dist(it1.cv, it2.cv));
  }
};

void DistMeth::calcInDist(Mdist &dm) {
  // for the intro-distances between the genomes
  OMP4TriAngleLoop loopOpt(steplist.size());
#pragma omp parallel for
  for (auto i = loopOpt.outBeg; i < loopOpt.outEnd; ++i) {
    for (auto j = loopOpt.inBeg; j < i; ++j) {
      setdist4NAN(dm, *steplist[i], *steplist[j]);
    }

    auto ir = loopOpt.inEnd - i;
    for (auto j = loopOpt.inBeg; j < ir; ++j) {
      setdist4NAN(dm, *steplist[ir], *steplist[j]);
    }
  }

  return;
};

void DistMeth::calcOutDist(Mdist &dm) {
// for the intro-distances between the genomes
#pragma omp parallel for
  for (auto i = 0; i < restlist.size(); ++i) {
    // read an orginal cv
    restlist[i]->fill();
    preprocess(restlist[i]->cv);

    // obtain the distances between the genome and the extend genomes
    for (auto j = 0; j < steplist.size(); ++j) {
      setdist4NAN(dm, *restlist[i], *steplist[j]);
    }
    restlist[i]->clear();
  }
  return;
};

void DistMeth::execute(const vector<string> &flist, Mdist &dm) {

  setflist(flist);
  updateStep(dm);
  long ndx(0);

  while (!restlist.empty()) {
    // check the memory and set steps
    theInfo(infoStep(++ndx, setStep(dm)), 1);

    // read the extend cv of the step
    fillBlock();
    theInfo("Complete read CVs", 1);

    // calculate the inner distances
    calcInDist(dm);
    theInfo("Complete calculate of main block");

    // calculate the outer distances
    if (!restlist.empty()) {
      calcOutDist(dm);
      theInfo("Complete calculate of rest block", -1);
    } else {
      theInfo("No rest block calculate are required", -1);
    }
    // clean the cvs
    updateStep(dm);

    // output complete infomation
    theInfo("Complete Distance Calculate of the Block", -1);
  }
}

string DistMeth::infoStep(long ndx, float size) {
  string str = "Start the calculate step " + to_string(ndx);
  str += "\n" + to_string(steplist.size()) + "/" +
         to_string(steplist.size() + restlist.size()) +
         " CVs will be resident in memory,\nwhich size is " +
         to_string(size / 1073741824) + "G";
  return str;
}

///.........................
/// Three method based on vector
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

void Cosine::preprocess(CVvec &cv) { cv.normalize(); }

double Euclidean::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  double dd = cv1.norm * cv1.norm + cv2.norm * cv2.norm;
  dd -= 2.0 * align(block1, block2);
  return sqrt(dd);
};

void Euclidean::preprocess(CVvec &cv) { cv.normalize(); }

double InterList::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  return 1.0 - 2.0 * overlap(block1, block2) / (cv1.norm + cv2.norm);
};

void InterList::preprocess(CVvec &cv) { cv.setNormL1(); }

double Min2Max::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  double qOverlap = overlap(block1, block2);
  double qTotal = cv1.norm + cv2.norm - qOverlap;
  return 1.0 - qOverlap / qTotal;
};

void Min2Max::preprocess(CVvec &cv) { cv.setNormL1(); }

///..............
/// Three method for set of kmer
double InterSet::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  return 1.0 -
         double(nInterSection(block1, block2)) / sqrt(cv1.norm * cv2.norm);
}

void InterSet::preprocess(CVvec &cv) { cv.setNormL0(); }

double Dice::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  return 1.0 -
         2.0 * double(nInterSection(block1, block2)) / (cv1.norm + cv2.norm);
}
void Dice::preprocess(CVvec &cv) { cv.setNormL0(); }

double ItoU::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  size_t nI = nInterSection(block1, block2);
  size_t nU = cv1.size() + cv2.size() - nI;
  return 1.0 - double(nI) / double(nU);
}

void ItoU::preprocess(CVvec &cv) { cv.setNormL0(); }
