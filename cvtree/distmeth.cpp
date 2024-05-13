/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-12 21:30:28
 */

#include "distmeth.h"
extern Info theInfo;

/// for CV item
void CVitem::fill() {

  //.. get the cv and its norm
  pair<double, string> meta = readcv(fname, cv);
  double norm = meta.first;

  //.. revise norm according to lp of method
  switch (DistMeth::lp) {
  case L0:
    norm = double(cv.size());
    break;
  case L1:
    norm = 0;
    for (auto &item : cv)
      norm += abs(item.second);
  case L2:
    break;
  }

  // do normlize if method required
  if (DistMeth::normalize) {
    for (auto &cvdim : cv) {
      cvdim.second /= norm;
    }
  }
}

///.. whether do normlize default
bool DistMeth::normalize = false;

///.. the type of norm in CVitem for the DistMeth
LPnorm DistMeth::lp = L2;

/// the create function
DistMeth *DistMeth::create(const string &methStr, bool normal) {

  // change whether do normalize based on parameter
  normalize = normal;

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
  float bs = maxNameLen * ng + nk * ng * (ng - 1) * sizeof(double) / 2 + giga;
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
  sort(restlist.begin(), restlist.end(), [](auto& a, auto& b){return a->nNAN < b->nNAN;});

  // chop the large NAN number item to step list
  float size(0.0);
  steplist.clear();
  for (auto iter = restlist.rbegin(); iter != restlist.rend(); ++iter) {
    if (size < maxM) {
      steplist.emplace_back(*iter);
      size += cvsize((*iter)->fname);
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

void DistMeth::updateStep(const Mdist& dm) {
  // clean cv in the steplist
  for (auto &item : steplist) {
    item->clear();
  }

  // get the nonNAN cvitems in restlist
  for(auto i=1; i<restlist.size(); ++i){
    for(auto j=0; j<i; ++j ){
      if(dm.isNAN(restlist[i]->ndx,restlist[j]->ndx)){
        ++restlist[i]->nNAN;
        ++restlist[j]->nNAN;
      }
    }
  }

  // get the nonNAN list to restlist
  vector<CVitem*> nanlist;
  for(auto& ci : restlist)
    if(ci->nNAN > 0)
      nanlist.emplace_back(ci);
  restlist.swap(nanlist);
}

size_t DistMeth::length() const { return steplist.size(); };

void DistMeth::fillBlock() {
#pragma omp parallel for
  for (auto i = 0; i < steplist.size(); ++i) {
    steplist[i]->fill();
  }
  return;
};

void DistMeth::setdist4NAN(Mdist &dm, const CVitem &cv1, const CVitem &cv2) {
  if (dm.isNAN(cv1.ndx, cv2.ndx)) {
    dm.setdist(cv1.ndx, cv2.ndx, dist(cv1, cv2));
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

  while(!restlist.empty()) {
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
double Cosine::dist(const CVitem &cv1, const CVitem &cv2) {
  CVblock block1(cv1.cv.begin(), cv1.cv.end());
  CVblock block2(cv2.cv.begin(), cv2.cv.end());

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

double InterList::dist(const CVitem &cv1, const CVitem &cv2) {
  CVblock block1(cv1.cv.begin(), cv1.cv.end());
  CVblock block2(cv2.cv.begin(), cv2.cv.end());

  return 1.0 - 2.0 * overlap(block1, block2) / (cv1.norm + cv2.norm);
};

double Min2Max::dist(const CVitem &cv1, const CVitem &cv2) {
  CVblock block1(cv1.cv.begin(), cv1.cv.end());
  CVblock block2(cv2.cv.begin(), cv2.cv.end());

  double qOverlap = overlap(block1, block2);
  double qTotal = cv1.norm + cv2.norm - qOverlap;
  return 1.0 - qOverlap / qTotal;
};

/// function to summary the rest vector
double Euclidean::ddTail(CVblock &blk) {
  double dd(0);
  while (!blk.pop())
    dd += blk.begin->second * blk.begin->second;
  return dd;
};

double Euclidean::dist(const CVitem &cv1, const CVitem &cv2) {
  CVblock block1(cv1.cv.begin(), cv1.cv.end());
  CVblock block2(cv2.cv.begin(), cv2.cv.end());

  // get the distance after reset bound with binary search
  double dd(0);
  for (;;) {
    if (block1.begin->first == block2.begin->first) {
      double dx = block1.begin->second - block2.begin->second;
      dd += dx * dx;
      if (block1.pop()) {
        dd += ddTail(block2);
        break;
      } else if (block2.pop()) {
        dd += ddTail(block1);
        break;
      }
    } else if (block1.begin->first < block2.begin->first) {
      if (block1.pop()) {
        dd += ddTail(block2);
        break;
      }
    } else {
      if (block2.pop()) {
        dd += ddTail(block1);
        break;
      }
    }
  }
  return sqrt(dd);
};

///..............
/// Three method for set of kmer
double InterSet::dist(const CVitem &cv1, const CVitem &cv2) {
  CVblock block1(cv1.cv.begin(), cv1.cv.end());
  CVblock block2(cv2.cv.begin(), cv2.cv.end());

  return 1.0 -
         double(nInterSection(block1, block2)) / sqrt(cv1.norm * cv2.norm);
}

double Dice::dist(const CVitem &cv1, const CVitem &cv2) {
  CVblock block1(cv1.cv.begin(), cv1.cv.end());
  CVblock block2(cv2.cv.begin(), cv2.cv.end());

  return 1.0 -
         2.0 * double(nInterSection(block1, block2)) / (cv1.norm + cv2.norm);
}

double ItoU::dist(const CVitem &cv1, const CVitem &cv2) {
  CVblock block1(cv1.cv.begin(), cv1.cv.end());
  CVblock block2(cv2.cv.begin(), cv2.cv.end());

  size_t nI = nInterSection(block1, block2);
  size_t nU = cv1.cv.size() + cv2.cv.size() - nI;
  return 1.0 - double(nI) / double(nU);
}
