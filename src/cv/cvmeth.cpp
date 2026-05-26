/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-21 Thursday 15:55:27
 */

#include "cvmeth.h"

CVmeth *CVmeth::ptr = nullptr;
CVmeth *CVmeth::create(const string &methStr, const string &cvdir,
                       const string &gtype, const string &cgfile, bool save) {

  // Don't recreat method
  if (ptr != nullptr)
    return ptr;

  if (methStr == "Hao") {
    ptr = new HaoMethod;
  } else if (methStr == "Count") {
    ptr = new Counting;
  } else {
    cerr << "Unknow CV Method: " << methStr << endl;
    exit(3);
  }

  // set letter type for sequence
  Letter::init(gtype, cgfile);

  // set parameters for CV method
  ptr->kmax = Kstr::kmax();
  ptr->setCVdir(cvdir);

  // cache the cv file
  ptr->save = save;

  return ptr;
};

void CVmeth::checkK(const vector<size_t> &klist) {
  if (klist.front() < kmin || klist.back() > kmax) {
    cerr << "The range of k value is [" << kmin << "," << kmax << "]" << endl;
    exit(3);
  }
}

void CVmeth::checkK(size_t k) {
  if (k < kmin || k > kmax) {
    cerr << "The range of k value is [" << kmin << "," << kmax << "]" << endl;
    exit(3);
  }
}

void CVmeth::setCVdir(const string &str) {
  cvdir = str;
  if (str.empty()) {
    getCVname = [this](const string &str, size_t k) {
      return str + cvsuff + to_string(k) + ".gz";
    };
  } else {
    addsuffix(cvdir, '/');
    mkpath(cvdir);
    getCVname = [this](const string &str, size_t k) {
      return cvdir + getFileName(str) + cvsuff + to_string(k) + ".gz";
    };
  }
}

void CVmeth::execute(const string &gfile, const vector<size_t> &klist,
                     bool chk) {

  vector<pair<int, CVmap>> mcv;
  // check the existed cvfile
  if (chk) {
    for (auto k : klist) {
      string cvfile = getCVname(gfile, k);
      if (!gzvalid(cvfile)) {
        CVmap cv;
        mcv.emplace_back(make_pair(k, cv));
      }
    }
  } else {
    for (auto k : klist) {
      CVmap cv;
      mcv.emplace_back(make_pair(k, cv));
    }
  }

  // calculate the CV
  if (!mcv.empty()) {
    // read genomes
    Genome genome;
    readFasta(gfile, genome);

    // get the cv of all K of the genome
    cv(genome, mcv);

    // write down CVs
    for (auto item : mcv) {
      string outfile = getCVname(gfile, item.first);
      CVvec tmp(item.second);
      tmp.write(outfile);
    }
  }
};

CVvec CVmeth::getCV(const string &fname, int k) {
  CVvec aCV;
  string cvfile = getCVname(fname, k);
  if (gzvalid(cvfile)) {
    aCV.read(cvfile);
  } else {
    vector<pair<int, CVmap>> mcv{make_pair(k, CVmap())};
    // read genomes
    Genome genome;
    readFasta(fname, genome);
    // get the cv in map format
    cv(genome, mcv);
    // convert CV map into CV vector
    aCV = CVvec(mcv.front().second);
    // save CV into file if required
    if (save)
      aCV.write(cvfile);
  }
  return aCV;
};


size_t CVmeth::infoCVA(const string &fname, int k) {
  string cvfile = getCVname(fname, k);
  if (gzvalid(cvfile)) {
    CVAinfo cvi(cvfile);
    return cvi.nCV;
  } else {
    CVArray cva;
    _caclCVA(fname, k, cva);
    cva.write(cvfile);
    return cva.cvdi.size();
  }
};

CVArray CVmeth::getCVA(const string &fname, int k) {
  string cvfile = getCVname(fname, k);
  CVArray cva;
  if (gzvalid(cvfile)) {
    cva.read(cvfile);
  } else {
    _caclCVA(fname, k, cva);
    if (save)
      cva.write(cvfile);
  }
  cva.setNorm();
  return cva;
};

void CVmeth::_caclCVA(const string &fname, int k, CVArray &cva) {
  // read genomes
  Genome genome;
  readFasta(fname, genome);

  // get the cv in map format
  vector<CVmap> cvs;
  for (auto &gene : genome) {
    vector<pair<int, CVmap>> mcv{make_pair(k, CVmap())};
    cv(gene, mcv);
    cvs.emplace_back(std::move(mcv.front().second));
  }
  cva = CVArray(cvs);
}

CVBundle CVmeth::getCVB(const string &fname, int k) {
  CVBundle cvb;
  string cvfile = getCVname(fname, k);
  if (gzvalid(cvfile)) {
    readCVBundle(cvfile, cvb);
  } else {
    // read genomes
    Genome genome;
    readFasta(fname, genome);

    // get the cv in map format
    cvb.reserve(genome.size());
    vector<pair<int, CVmap>> mcv{make_pair(k, CVmap())};
    for (auto &gene : genome) {
      cv(gene, mcv);
      cvb.emplace_back(CVvec(mcv.front().second));
      mcv.front().second.clear();
    }

    // save the gene cv for genome into cvbundle
    if (save)
      writeCVBundle(cvfile, cvb);
  }

  // return the size of cvbundle
  return cvb;
};

// count the kmers
size_t CVmeth::count(const Genome &genome, size_t k, CVmap &cv) {
  size_t n(0);
  for (const auto &gene : genome) {
    // the number of kstring of the gene
    n += count(gene, k, cv);
  }
  return n;
};

size_t CVmeth::count(const Gene &gene, size_t k, CVmap &cv) {
  // get the first k string
  Kstr ks(gene.substr(0, k));
  CVmap::iterator iter = cv.find(ks);
  if (iter == cv.end()) {
    cv[ks] = 1.0;
  } else {
    ++(iter->second);
  }

  // get the next kstring
  for (int i = k; i < gene.size(); ++i) {
    ks.behead();
    ks.append(gene[i]);
    CVmap::iterator iter = cv.find(ks);
    if (iter == cv.end()) {
      cv[ks] = 1.0;
    } else {
      ++(iter->second);
    }
  }

  return gene.size() - k + 1;
};

// The Markov Method
void Counting::cv(const Genome &genome, vector<pair<int, CVmap>> &vcv) {
  for (auto &item : vcv) {
    count(genome, item.first, item.second);
  }
};

void Counting::cv(const Gene &gene, vector<pair<int, CVmap>> &vcv) {
  for (auto &item : vcv) {
    count(gene, item.first, item.second);
  }
};

// Hao method based the Markov Model
void HaoMethod::cv(const Genome &genome, vector<pair<int, CVmap>> &vcv) {
  docv(genome, vcv);
}

void HaoMethod::cv(const Gene &gene, vector<pair<int, CVmap>> &vcv) {
  docv(gene, vcv);
}

template <typename T>
void HaoMethod::docv(const T &seq, vector<pair<int, CVmap>> &vcv) {

  // require k to count
  set<int> slist;
  for (auto &item : vcv) {
    slist.insert(item.first);
    slist.insert(item.first - 1);
    slist.insert(item.first - 2);
  }

  // count k string
  map<int, CVmap> mvc;
  map<int, double> nstr;
  for (auto &k : slist) {
    CVmap cv;
    nstr[k] = count(seq, k, cv);
    mvc[k] = cv;
  }

  // get the cv with subtract
  for (auto &item : vcv) {
    int k = item.first;
    double factor = nstr[k] * nstr[k - 2] / (nstr[k - 1] * nstr[k - 1]);
    markov(mvc[k], mvc[k - 1], mvc[k - 2], factor, item.second);
  }
};

void HaoMethod::markov(const CVmap &mck, const CVmap &mckM1, const CVmap &mckM2,
                       double factor, CVmap &cv) {

  CVmap::const_iterator iter;
  for (const auto &cd : mckM2) {
    Kstr ksM2 = cd.first;
    double nksM2 = cd.second;
    for (auto &c : Letter::nonVoid) {
      Kstr ksM1A = ksM2;
      ksM1A.append(c);
      iter = mckM1.find(ksM1A);
      if (iter != mckM1.end()) {
        double nksM1A = iter->second;

        for (auto &d : Letter::nonVoid) {
          Kstr ksM1B = ksM2;
          ksM1B.addhead(d);
          iter = mckM1.find(ksM1B);
          if (iter != mckM1.end()) {
            double nksM1B = iter->second;
            double nks0 = factor * nksM1B * nksM1A / nksM2;

            Kstr ks = ksM1B;
            ks.append(c);
            iter = mck.find(ks);
            double nks = 0;
            if (iter != mck.end())
              nks = iter->second;
            cv[ks] = (nks - nks0) / nks0;
          }
        }
      }
    }
  }
};