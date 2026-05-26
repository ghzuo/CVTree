/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-11-28 22:22:12
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-17 17:39:35
 */

#include "cvcache.h"

/********************************************************************************
 * @brief Section for read the previous data
 ********************************************************************************/
// function for name list index
void CVCacheHandle::checkCache() {
  string dbcache = filename + ".data.cache";
  string gncache = filename + ".name.cache";
  if (fileExists(dbcache) && fileExists(gncache)) {
    readName(gncache);
    readCache(dbcache);
  } else {
    cerr << "The database files are incompleted, please check!" << endl;
    exit(2);
  }
}

void CVCacheHandle::readName(const string &file) {
  ifstream fgn(file.c_str());
  for (string line; getline(fgn, line);) {
    if (!line.empty()) {
      CVGinfo cgi;
      istringstream buf(line);
      buf >> cgi;
      name.emplace_back(cgi);
    }
  }
  fgn.close();
};

void CVCacheHandle::readCache(const string &file) {
  ifstream fcv(file.c_str());
  CVbrick cb;
  while (fcv.read(reinterpret_cast<char *>(&cb), sizeof(cb))) {
    auto iter = kstr.find(cb.label.kstr);
    if (iter != kstr.end()) {
      iter->second.emplace_back(CVKinfo(nBrick++, cb.getlen()));
    } else {
      kstr[cb.label.kstr].emplace_back(CVKinfo(nBrick++, cb.getlen()));
    }
  }
};

/********************************************************************************
 * @brief Section for insert data
 ********************************************************************************/
bool CVCacheHandle::has(const string &gn) {
  return binary_search(name.begin(), name.end(), gn);
};

void CVCacheHandle::fillCV(CVmeth *cmeth, const vector<string> &flist,
                           vector<CVGinfo> &cglist, size_t k, bool app) {
  // open cache file for write
  string dbcache = filename + ".data.cache";
  string gncache = filename + ".name.cache";
  if (app) {
    checkCache();
    fcv.open(dbcache.c_str(), ios::out | ios::binary | ios::app);
    fnm.open(gncache.c_str(), ios::out | ios::app);
  } else {
    mkpath(filename);
    fcv.open(dbcache.c_str(), ios::out | ios::binary);
    fnm.open(gncache.c_str(), ios::out);
  }

  // fill the missing cv item
  for (auto i = 0; i < cglist.size(); ++i) {
    auto cgi = cglist[i];
    if (isnan(cgi.norm)) {
      // get CV and CV genome info
      CVvec cv;
      cgi.norm = cmeth->getcv(flist[i], k, cv);
      cgi.len = cv.size();

      // insert CV into data index and get wirte file items
      vector<pair<Kstr, CVatom>> bricks;
      vector<pair<CVKinfo, CVatom>> atoms;
      insertCV(cv, cgi, bricks, atoms);

      // write CV into file
      writeOneCV(bricks, atoms);
    }
    theInfo("Complete the " + to_string(cgi.ndx) + " CV");
  }
  fcv.close();
  fnm.close();
}

void CVCacheHandle::insertCV(const CVvec &cv, CVGinfo cgi,
                             vector<pair<Kstr, CVatom>> &bricks,
                             vector<pair<CVKinfo, CVatom>> &atoms) {
  // add cv into index of database
  theInfo("Begin to align CV, number = " + to_string(cgi.len));
  cgi.lasso = 0;
  for (auto i = 0; i < cv.size(); ++i) {
    cgi.lasso += abs(cv[i].second);
    CVatom ci(cgi.ndx, cv[i].second);
    Kstr ks = cv[i].first;
    auto iks = kstr.lower_bound(ks);
    if (iks->first != ks) {
      iks = kstr.emplace_hint(iks, ks, vector<CVKinfo>());
      iks->second.emplace_back(nBrick++, 1);
      bricks.emplace_back(ks, ci);
    } else {
      if (iks->second.back().len < BSIZE) {
        atoms.emplace_back(iks->second.back(), ci);
        ++(iks->second.back().len);
      } else {
        iks->second.emplace_back(nBrick++, 1);
        bricks.emplace_back(ks, ci);
      }
    }
  }

  // get the cv info and index transition
  cgi.ndx = name.size();
  name.emplace_back(cgi);
  fnm << cgi << endl;
}

void CVCacheHandle::writeOneCV(const vector<pair<Kstr, CVatom>> &bricks,
                               const vector<pair<CVKinfo, CVatom>> &atoms) {
  // append a new CV bricks into file
  theInfo("Begin to append CV bricks, number = " + to_string(bricks.size()));
  fcv.seekp(0, ios_base::end);
  for (auto &it : bricks) {
    appendCVbrick(it.first, it.second);
  }

  // insert other CV atoms into file
  theInfo("Begin to insert CV atoms number = " + to_string(atoms.size()));
  for (auto &it : atoms) {
    insertCVatom(it.first, it.second);
  }

  // sync with file
  fcv.flush();
};

void CVCacheHandle::appendCVbrick(const Kstr &ks, const CVatom &ci) {
  // append a new CVbrick into file
  CVbrick cvb(ks);
  cvb.vec[0] = ci;
  fcv.write(reinterpret_cast<const char *>(&cvb), sizeof(cvb));
};

void CVCacheHandle::insertCVatom(const CVKinfo &cbi, const CVatom &ci) {
  // write a CV atom into file
  size_t offset = sizeof(CVbrick) * cbi.ndx + sizeof(CVatom) * cbi.len;
  fcv.seekp(offset, ios_base::beg);
  fcv.write(reinterpret_cast<const char *>(&ci), sizeof(CVatom));
};

/********************************************************************************
 * @brief Section for write down the cv array
 ********************************************************************************/
void CVCacheHandle::writecva() {

  // file names data
  string dbcache = filename + ".data.cache";
  string dfile = filename + ".data";
  string kfile = filename + ".kstr";
  string gfile = filename + ".name";

  // write down the name list of cv
  ofstream fgn(gfile.c_str());
  sort(name.begin(), name.end());
  for (auto &it : name) {
    fgn << it << endl;
  }
  fgn.close();

  // write down CVatom and kstr index
  fcv.open(dbcache.c_str(), ios::in);
  ofstream fks(kfile.c_str(), ios::binary);
  ofstream fca(dfile.c_str(), ios::binary);
  for (auto &vcki : kstr) {
    size_t nItem = writeblock(vcki, fca);
    auto aks = make_pair(vcki.first, nItem);
    fks.write(reinterpret_cast<char *>(&aks), sizeof(aks));
  }
  fks.close();
  fca.close();
  fcv.close();
}

size_t CVCacheHandle::writeblock(const pair<Kstr, vector<CVKinfo>> &vcki,
                                 ostream &os) {
  size_t nItem(0);
  for (auto &cki : vcki.second) {
    nItem += cki.len;
    size_t offset = sizeof(CVbrick) * cki.ndx;
    size_t size = sizeof(CVatom) * cki.len;
    char buff[size];
    fcv.seekg(offset, ios_base::beg);
    fcv.read(buff, size);
    os.write(buff, size);
  }
  // write the kstr and seperate label for the K block
  KstrBlockSep ksep(vcki.first);
  os.write(reinterpret_cast<char *>(&ksep), sizeof(KstrBlockSep));
  return nItem;
};
