/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-12-03 22:15:34
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-10 Wednesday 18:36:10
 */

#include "karray.h"

ostream &operator<<(ostream &os, const Kitem &ki) {
  os << ki.index << ":" << ki.value;
  return os;
};

ostream &operator<<(ostream &os, const KdimInfo &kdi) {
  os << kdi.kstr << "\t" << kdi.index.first << "\t" << kdi.index.second;
  return os;
};

CVdimInfo::CVdimInfo(const CVvec& cv){
  len = float(cv.size());
  lasso = 0.0;
  for(const auto& cd : cv)
    lasso += cd.second;
  norm = cv.norm;
};

CVdimInfo::CVdimInfo(const CVmap& cm) {
  len = float(cm.size());
  lasso = 0.0;
  norm = 0.0;
  for(const auto& cd : cm){
    lasso += cd.second;
    norm += cd.second * cd.second;
  }
  norm = sqrt(norm);
};


ostream &operator<<(ostream &os, const CVdimInfo &cdi) {
  os << cdi.len << "\t" << cdi.lasso << "\t" << cdi.norm;
  return os;
};
