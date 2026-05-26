/*
 * Copyright (c) 2025
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2025-12-06 Saturday 18:00:39
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-10 Wednesday 14:38:22
 */

#include "cvbundle.h"

string readCVBundle(const string &file, CVBundle &cvs) {
  // open file to read
  gzFile fp;
  if ((fp = gzopen(file.c_str(), "rb")) == NULL) {
    cerr << "CV file not found: \"" << file << '"' << endl;
    exit(1);
  }

  // get the cvs size
  size_t ncv;
  gzread(fp, &ncv, sizeof(size_t));
  // read the cvs: read size and then the cvs
  for (size_t i = 0; i < ncv; i++) {
    size_t size;
    gzread(fp, (char *)&size, sizeof(size_t));
    CVvec cv(size);
    gzread(fp, (char *)cv.vec.data(), size * sizeof(CVdim));
    cvs.emplace_back(cv);
  }

  // read char set
  string charSet;
  gzline(fp, charSet);
  gzclose(fp);

  return charSet;
};

void writeCVBundle(const string &file, const CVBundle &cvs) {
  gzFile fp;
  if ((fp = gzopen(file.c_str(), "wb")) == NULL) {
    cerr << "Error happen on write cvfile: " << file << endl;
    exit(1);
  }

  size_t ncv = cvs.size();
  gzwrite(fp, &ncv, sizeof(size_t));

  for (auto &cv : cvs) {
    size_t size = cv.size();
    gzwrite(fp, &size, sizeof(size_t));
    gzwrite(fp, cv.vec.data(), size * sizeof(CVdim));
  }

  // write Letter char set
  gzputs(fp, Letter::decMap.c_str());
  gzclose(fp);
};