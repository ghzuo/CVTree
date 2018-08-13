/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-07-26 12:49:56
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 12:49:56
 */

#include "method.h"
/***********************************************
 ******* for the parent class
 ************************************************/
void Method::setg(const string &gtype) {
  // set the genome type
  gsuff = "." + gtype;
  cvsuff = gsuff + cvsuff;

  // init genome type read file
  theg.init(gtype);

  // get the genome letters map to check the sequcne
  kmax = Kstr::init(theg.letters);
};

void Method::checkK(const vector<size_t> &klist) {
  if (*(klist.begin()) < kmin || *(--klist.end()) > kmax) {
    cerr << "The range of k value is [" << kmin << "," << kmax << "]" << endl;
    exit(3);
  }
}

string Method::gcSameDir(const string &str) { return str + cvsuff; };

string Method::gcDiffDir(const string &str) {
  return cvdir + str.substr(str.find_last_of('/') + 1) + cvsuff;
};

void Method::setCVdir(const string &str) {
  cvdir = str;
  if (str.empty()) {
    getCVpref = &Method::gcSameDir;
  } else {
    mkdir(cvdir.c_str(), 0755); 
    getCVpref = &Method::gcDiffDir;
  }
}

string Method::getCVname(const string &gname, size_t k) {
  return (this->*getCVpref)(gname) + to_string(k) + ".gz";
};

void Method::getCV(const string &gname, const vector<size_t> &klist, bool check) {

  vector<pair<int, CVmap>> mcv;
  // check the existed cvfile
  if (check) {
    for (auto k : klist) {
      string cvfile = getCVname(gname, k);
      if (!fileExists(cvfile)) {
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
    string gfile = gname + gsuff;
    Genome genome;
    theg.readgene(gfile, genome);

    // get the cv of all K of the genome
    cv(genome, mcv);

    // write down CVs
    for (auto item : mcv) {
      string outfile = getCVname(gname, item.first);
      writecv(item.second, outfile);
    }
  }
};

size_t Method::count(const Genome &genome, size_t k, CVmap &cv) {
  size_t n(0);
  for (const auto &gene : genome) {
    // the number of kstring of the gene
    n += (gene.size() - k + 1);

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
  }
  return n;
};

void Method::markov(const CVmap &mck, const CVmap &mckM1, const CVmap &mckM2,
                    double factor, CVmap &cv) {

  CVmap::const_iterator iter;
  for (const auto &cd : mckM2) {
    Kstr ksM2 = cd.first;
    double nksM2 = cd.second;
    for (auto &c : ksM2.charSet) {
      Kstr ksM1A = ksM2;
      ksM1A.append(c);
      iter = mckM1.find(ksM1A);
      if (iter != mckM1.end()) {
        double nksM1A = iter->second;

        for (auto &d : ksM2.charSet) {
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

Method *Method::create(const string &str, string gtype) {

  Method *meth;

  if (str == "Hao" || str == "CVTree") {
    meth = new HaoMethod;
  } else if (str == "Li" || str == "InterSet") {
    meth = new LiMethod;
  } else if (str == "Zuo" || str == "InterList") {
    meth = new ZuoMethod;
  } else {
    cerr << "Unknow Method: " << str << endl;
    exit(3);
  }

  (*meth).setg(gtype);

  return meth;
};

/***********************************************
 ******* for the Hao Method class
 ************************************************/
void HaoMethod::cv(const Genome &genome, vector<pair<int, CVmap>> &vcv) {

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
    nstr[k] = count(genome, k, cv);
    mvc[k] = cv;
  }

  // get the cv with subtract
  for (auto &item : vcv) {
    int k = item.first;
    double factor = nstr[k] * nstr[k - 2] / (nstr[k - 1] * nstr[k - 1]);
    markov(mvc[k], mvc[k - 1], mvc[k - 2], factor, item.second);
  }
};

double HaoMethod::dist(const CVvec &cv1, const CVvec &cv2) {
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

/***********************************************
 ******* for the Li Method class
 ************************************************/
// get cv
void LiMethod::cv(const Genome &genome, vector<pair<int, CVmap>> &vcv) {
  for (auto &item : vcv) {
    count(genome, item.first, item.second);
  }
};

// get distance
double LiMethod::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  // get the distance after reset bound with binary search
  if (fitBegin(block1, block2)) {
    double n(0);
    for (;;) {
      if (block1.begin->first == block2.begin->first) {
        n++;
        if (block1.pop())
          break;
        if (block2.pop())
          break;
      } else {
        if (block1.begin->first < block2.begin->first) {
          if (block1.pop())
            break;
        } else {
          if (block2.pop())
            break;
        }
      }
    }
    return 1.0 - n / sqrt(double(cv1.size()) * double(cv2.size()));

  } else {
    return 1.0;
  }
}
/***********************************************
 ******* for the Zuo Method class
 ************************************************/
void ZuoMethod::cv(const Genome &genome, vector<pair<int, CVmap>> &vcv) {
  for (auto &item : vcv) {
    count(genome, item.first, item.second);
  }
};

// get distance
double ZuoMethod::dist(const CVvec &cv1, const CVvec &cv2) {
  CVblock block1(cv1.begin(), cv1.end());
  CVblock block2(cv2.begin(), cv2.end());

  // get the distance after reset bound with binary search
  if (fitBegin(block1, block2)) {
    double n(0);
    double sum(0);
    for (;;) {
      if (block1.begin->first == block2.begin->first) {
        n += min(block1.begin->second, block1.begin->second);
        sum += block1.begin->second;
        sum += block2.begin->second;
        if (block1.pop()) {
          while (!block2.pop()) {
            sum += block2.begin->second;
          }
          break;
        }
        if (block2.pop()) {
          while (!block1.pop()) {
            sum += block1.begin->second;
          }
          break;
        }
      } else {
        if (block1.begin->first < block2.begin->first) {
          sum += block1.begin->second;
          if (block1.pop()) {
            while (!block2.pop()) {
              sum += block2.begin->second;
            }
            break;
          }
        } else {
          sum += block2.begin->second;
          if (block2.pop()) {
            while (!block1.pop()) {
              sum += block1.begin->second;
            }
            break;
          }
        }
      }
    }
    return 1.0 - 2 * n / sum;
  } else {
    return 1.0;
  }
}

/***********************************************
 ******* command methods
 ************************************************/
