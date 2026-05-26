/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-19 Tuesday 16:29:08
 */

#include "similarity.h"

/**************************************************************
 * the similar methods
 **************************************************************/
Similarity *Similarity::create(const string &methStr) {

  // create the distance method
  Similarity *meth;
  if (methStr == "Cosine") {
    meth = new Cosine();
    CVvec::lp = L2;
  } else if (methStr == "Euclidean") {
    meth = new Euclidean();
    CVvec::lp = L2;
  } else if (methStr == "InterList") {
    meth = new InterList();
    CVvec::lp = L1;
  } else if (methStr == "Min2Max") {
    meth = new Min2Max();
    CVvec::lp = L1;
  } else if (methStr == "InterSet") {
    meth = new InterSet();
    CVvec::lp = L0;
  } else if (methStr == "Jaccard" || methStr == "ItoU") {
    meth = new ItoU();
    CVvec::lp = L0;
  } else if (methStr == "Dice") {
    meth = new Dice();
    CVvec::lp = L0;
  } else {
    cerr << "Unknow Distance Method: " << methStr << endl;
    exit(3);
  }

  return meth;
}

void Similarity::calcSim(const CVArray &cva, const CVArray &cvb,
                         Mrectangle &mtx) {
  vector<pair<size_t, size_t>> aln;
  alignSortVector(cva.kdi, cvb.kdi, aln);
  calcSim(cva, cvb, aln, mtx);
};

void Similarity::calcSim(const CVArray &cva, const CVArray &cvb,
                         const vector<pair<size_t, size_t>> &aln,
                         Mrectangle &mtx) {
  for (auto &it : aln) {
    Kblock kba = cva.getKblock(it.first);
    Kblock kbb = cvb.getKblock(it.second);

    for (auto ka = kba.begin(); ka != kba.end(); ++ka) {
      for (auto kb = kbb.begin(); kb != kbb.end(); ++kb) {
        mtx.add(ka->index, kb->index, calcDim(ka->value, kb->value));
      }
    }
  }

  for (auto i = 0; i < cva.norm.size(); ++i) {
    for (auto j = 0; j < cvb.norm.size(); ++j) {
      mtx.set(i, j, scale(mtx.get(i, j), cva.norm[i], cvb.norm[j]));
    }
  }
};

void Similarity::calcSim(const CVArray &cva, Mtriangle &mtx) {
  for (size_t k = 0; k < cva.kdi.size(); ++k) {
    Kblock kb = cva.getKblock(k);

    for (auto ita = kb.begin(); ita != kb.end(); ++ita) {
      for (auto itb = ita+1; itb != kb.end(); ++itb) {
        mtx.add(size_t(ita->index), size_t(itb->index), calcDim(ita->value, itb->value));
      }
    }
  }

  for (size_t i = 1; i < cva.norm.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      mtx.set(i, j, scale(mtx.get(i, j), cva.norm[i], cva.norm[j]));
    }
  }
};

///.........................
/// Three method based on vector
float Cosine::calcDim(float a, float b) { return a * b; }
float Cosine::scale(float val, float aNorm, float bNorm) {
  return val / (aNorm * bNorm);
}

float Euclidean::calcDim(float a, float b) { return (a - b) * (a - b); }
float Euclidean::scale(float val, float aNorm, float bNorm) {
  // for vector a{a1, a2, 0} and b{b1, 0, b3}
  // d^2 = (a1-b1)^2 + a2^2 + b3^2 = (a1^2 + a2^2) + (b1^2 + b3^2) - 2*a1*b1
  // d = sqrt(2 - 2 * a1 *b1)
  return 1 - 0.5 * sqrt(2) * sqrt(1 - val);
}

// ... distance scaling at L1

float InterList::calcDim(float a, float b) { return min(a, b); }
float InterList::scale(float val, float aNorm, float bNorm) {
  return 2.0 * val / (aNorm + bNorm);
}

float Min2Max::calcDim(float a, float b) { return min(a, b); }
float Min2Max::scale(float val, float aNorm, float bNorm) { return val; }

// ... distance scaling at L0
float InterSet::calcDim(float a, float b) { return 1.0; }
float InterSet::scale(float val, float aNorm, float bNorm) {
  return val / sqrt(aNorm * bNorm);
}

float Dice::calcDim(float a, float b) { return 1.0; }

float Dice::scale(float val, float aNorm, float bNorm) {
  return 2.0 * val / (aNorm + bNorm);
}

float ItoU::calcDim(float a, float b) { return 1.0; }

float ItoU::scale(float val, float aNorm, float bNorm) {
  return val / (aNorm + bNorm - val);
}
