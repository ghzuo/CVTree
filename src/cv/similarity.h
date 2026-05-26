/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-19 Tuesday 16:29:29
 */

#ifndef SIMILARITY_H
#define SIMILARITY_H

#include "cvarray.h"
#include "kit.h"
#include "rectangleMatrix.h"
#include "triangleMatrix.h"

using namespace std;

struct Similarity {
  // the create function
  float mindist;
  static Similarity *create(const string &);

  // get the similarity matrix without give alignment
  void calcSim(const CVArray &, const CVArray &, Mrectangle &);
  // get the similarity matrix with give alignment
  void calcSim(const CVArray &, const CVArray &,
               const vector<pair<size_t, size_t>> &, Mrectangle &);
  // get the similarity matrix within a single CVArray
  void calcSim(const CVArray &, Mtriangle &);
  // scale the value at the end
  virtual float calcDim(float, float) = 0;
  virtual float scale(float, float, float) = 0;
};

// ... son class for different method
// ... distance scaling at L2
struct Cosine : public Similarity {
  float calcDim(float, float) override;
  float scale(float, float, float) override;
};

struct Euclidean : public Similarity {
  float calcDim(float, float) override;
  float scale(float, float, float) override;
};

// ... distance scaling at L1
struct InterList : public Similarity {
  float calcDim(float, float) override;
  float scale(float, float, float) override;
};

struct Min2Max : public Similarity {
  float calcDim(float, float) override;
  float scale(float, float, float) override;
};

// ... distance scaling at L0
struct InterSet : public Similarity {
  float calcDim(float, float) override;
  float scale(float, float, float) override;
};

struct Dice : public Similarity {
  float calcDim(float, float) override;
  float scale(float, float, float) override;
};

struct ItoU : public Similarity {
  float calcDim(float, float) override;
  float scale(float, float, float) override;
};
#endif