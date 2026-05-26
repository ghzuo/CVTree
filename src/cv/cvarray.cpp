/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-09 5:10:18
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-21 Thursday 16:40:25
 */

#include "cvarray.h"
// for CV array info/header
void CVAinfo::read(const string &fname) {
  // open and test file
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL) {
    cerr << "CV file not found: \"" << gzfile << '"' << endl;
    exit(1);
  }

  // get size of the cvarray
  gzread(fp, (char *)this, sizeof(CVAinfo));
  gzclose(fp);
};

ostream &operator<<(ostream &os, const CVAinfo &hd) {
  os << hd.nCV << "\t" << hd.nKstr << "\t" << hd.nItem;
  return os;
};

void CVArray::setNorm() {
  setNorm(CVvec::lp);
}

void CVArray::setNorm(enum LPnorm lp) {
  norm.clear();
  norm.reserve(cvdi.size());
  switch (lp) {
  case L0:
    for (auto &ci : cvdi)
      norm.emplace_back(ci.len);
    break;
  case L1:
    for (auto &ci : cvdi)
      norm.emplace_back(ci.lasso);
    break;
  case L2:
    for (auto &ci : cvdi)
      norm.emplace_back(ci.norm);
    break;
  }
  cvdi.clear();
};

Kblock CVArray::getKblock(size_t ndx) const {
  auto _beg = data.begin() + kdi[ndx].index.first;
  auto _end = data.begin() + kdi[ndx].index.second;
  Kblock kb(_beg, _end);
  return kb;
};

void CVArray::read(const string &fname) {
  try {
    // open file to read
    gzFile fp;
    string gzfile = addsuffix(fname, ".gz");
    if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL)
      throw runtime_error("Cannot open file for reading: " + gzfile);

    // get size of the cvarray
    CVAinfo hd;
    gzread(fp, (char *)&hd, sizeof(CVAinfo));

    // read the cvdiminfo
    cvdi.resize(hd.nCV);
    gzread(fp, (char *)cvdi.data(), sizeof(CVdimInfo) * hd.nCV);

    // read the kdiminfo
    kdi.resize(hd.nKstr);
    gzread(fp, (char *)kdi.data(), sizeof(KdimInfo) * hd.nKstr);

    // read the data
    data.resize(hd.nItem);
    gzread(fp, (char *)data.data(), sizeof(Kitem) * hd.nItem);

    // close file
    gzclose(fp);

  } catch (std::exception &e) {
    cerr << "Error reading file: " << fname << "\n" << e.what() << endl;
    exit(1);
  } catch (...) {
    cerr << "Error reading file: " << fname << "\n" << endl;
    exit(1);
  }
};

void CVArray::write(const string &fname) const {
  // open and test file
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "wb")) == NULL) {
    cerr << "Error happen on write cvfile: " << gzfile << endl;
    exit(1);
  }

  // write the size of CVArray
  CVAinfo hd(cvdi.size(), kdi.size(), data.size());
  gzwrite(fp, &hd, sizeof(CVAinfo));

  // write the diminfo and data
  gzwrite(fp, cvdi.data(), cvdi.size() * sizeof(CVdimInfo));
  gzwrite(fp, kdi.data(), kdi.size() * sizeof(KdimInfo));
  gzwrite(fp, data.data(), data.size() * sizeof(Kitem));

  // close file
  gzclose(fp);
};

void CVArray::calcCVdimInfo() {
  cvdi.clear();
  cvdi.resize(norm.size(), CVdimInfo(0.0));
  for(auto& d : data) {
    cvdi[d.index].len++;
    cvdi[d.index].lasso += d.value;
    cvdi[d.index].norm += d.value * d.value;
  }
  for(auto& cd : cvdi)
    cd.norm = sqrt(cd.norm);
};

void CVArray::calcNorm() {
  calcNorm(CVvec::lp);
}

void CVArray::calcNorm(enum LPnorm lp) {
  fill(norm.begin(), norm.end(), 0.0);
  if (lp == L0){
    for(auto& d : data)
      norm[d.index] ++;
  } else if(lp == L1){
    for(auto& d : data)
      norm[d.index] += abs(d.value);
  } else {
    for(auto& d : data)
      norm[d.index] += d.value * d.value;
    for(auto& d : norm)
      d = sqrt(d);
  }
};

string CVArray::info() const {
  string infoStr;
  infoStr += to_string(cvSize()) + " CVs, ";
  infoStr += to_string(kstrSize()) + " Kmers, ";
  infoStr += to_string(itemSize()) + " items";
  return infoStr; 
}



ostream &operator<<(ostream &os, const CVArray &cva) {
  // output cvdiminfo
  os << "The number of CV is " << cva.cvdi.size() << "\nThe number of K is "
     << cva.kdi.size() << "\nThe number of items is " << cva.data.size()
     << endl;

  os << "\n==== The Norms ========================\n";
  for (const auto &cd : cva.cvdi) {
    os << cd << "\n";
  }
  os << endl;

  // output data according to kdiminfo
  os << "\n==== The kmer items ====================\n";
  for (const auto &kd : cva.kdi) {
    os << kd.kstr << "\t";
    for (auto i = kd.index.first; i < kd.index.second; ++i) {
      os << cva.data[i] << " ";
    }
    os << endl;
  }

  return os;
};
