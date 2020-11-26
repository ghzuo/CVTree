/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-04-19 11:37:42
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-19 10:53:57
 */

#include "distmatrix.h"

Mdist::Mdist() : ng(0){};

void Mdist::init(const vector<string> &nmlist, bool chk) {
  ng = nmlist.size();
  name = nmlist;
  dist.resize(msize(), NAN);
  if (chk)
    cleanName();
};

// read the tranditional infile for the distance matrix
bool Mdist::readmtx(const string &file) {
  bool done = false;
  if (!fileExists(file)) {
#pragma omp critical
    { cerr << "Cannot find file: " << file << endl; }
    return done;
  }

#ifdef _HDF5
#pragma omp critical
  { done = readmtxh5(file); }
  if (done)
    return done;
#endif // _HDF5

#ifdef _NETCDF
#pragma omp critical
  { done = readmtxnc(file); }
  if (done)
    return done;
#endif // _NETCDF

  return readmtxtxt(file);
};

// select the distance matrix file type
void Mdist::writemtx(const string &file) {

#ifdef _NETCDF
  if (getsuffix(file) != "txt" && getsuffix(file) != "h5") {
#pragma omp critical
    { writemtxnc(file); }
    return;
  }
#endif // _NETCDF

#ifdef _HDF5
  if (getsuffix(file) != "txt") {
#pragma omp critical
    { writemtxh5(file); }
    return;
  }
#endif // _hDF5

  writemtxtxt(file);
  return;
};

// read the tranditional infile for the distance matrix
bool Mdist::readmtxtxt(const string &file) {
  ifstream dd(file.c_str());
  if (!dd.is_open()) {
    cerr << "Error opening file: " << file << endl;
    exit(1);
  }

  string line;
  getline(dd, line);
  resize(stoi(line));

  int idist(0);
  for (size_t i = 0; i < ng; ++i) {
    string line;
    if (getline(dd, line) && !line.empty()) {
      istringstream iss(line);
      iss >> name[i];
      for (int j = 0; j < i; ++j)
        iss >> dist[idist++];
      idist++;
    }
  }

  return true;
};

// write the tranditional infile for the distance matrix
void Mdist::writemtxtxt(const string &file) {
  ofstream dd(file.c_str());
  if (!dd.is_open()) {
    cerr << "Error opening file: " << file << endl;
    exit(1);
  }

  dd << ng << endl;
  for (size_t i = 0; i < ng; ++i) {
    dd << name[i] << " ";
    for (size_t j = 0; j < ng; ++j)
      dd << left << fixed << setprecision(10) << setw(13) << getdist(i, j);
    dd << endl;
  }
}

#ifdef _NETCDF
// read the netcdf file for the distance matrix
bool Mdist::readmtxnc(const string &file) {
  NcFile mtxFile(file.c_str(), NcFile::ReadOnly);
  if (!mtxFile.is_valid())
    return false;

  // get the dim information
  long lenWord = mtxFile.get_dim("word")->size();
  long ngenome = mtxFile.get_dim("genome")->size();
  long nmtx = mtxFile.get_dim("matrix")->size();

  // resize the matrix
  resize(ngenome);

  // read the space name
  NcVar *spname = mtxFile.get_var("spname");
  long counts[2] = {1, lenWord};
  char aname[lenWord];
  for (int i = 0; i < ngenome; ++i) {
    spname->set_cur(i, 0);
    spname->get(aname, counts);
    name[i] = aname;
  }

  // read the distance matrix
  float *fdist = new float[nmtx];
  NcVar *distance = mtxFile.get_var("distance");
  distance->get(fdist, nmtx);

  // set the distance
  if (nmtx == ngenome * (ngenome - 1) / 2) {
    for (int i = 0; i < nmtx; ++i)
      dist[i] = fdist[i];
  } else {
    // when the size is ng*(ng+1)/2
    size_t ndx(0);
    for (int i = 1; i < ng; ++i) {
      int m = (i + 1) * i / 2;
      for (int j = 0; j < i; ++j) {
        dist[ndx++] = fdist[m + j];
      }
    }
  }
  delete[] fdist;

  mtxFile.close();
  return true;
};

// write netcdf format distance matrix
void Mdist::writemtxnc(const string &file) {
  NcFile mtxFile(file.c_str(), NcFile::Replace);
  if (!mtxFile.is_valid()) {
    cout << "Couldn't open file!\n";
    cerr << NC_ERR << endl;
  }

  // output the space name
  size_t wlength(0);
  for (const auto &str : name)
    if (str.size() > wlength)
      wlength = str.size();
  ++wlength;

  NcDim *lenWord = mtxFile.add_dim("word", wlength);
  NcDim *ngenome = mtxFile.add_dim("genome", name.size());
  NcVar *spname = mtxFile.add_var("spname", ncChar, ngenome, lenWord);
  for (long i = 0; i < name.size(); ++i) {
    long nameSize = name[i].size() + 1;
    long counts[2] = {1, nameSize};
    spname->set_cur(i, 0);
    spname->put(name[i].c_str(), counts);
  }

  // output the distance matrix
  NcDim *mtxdist = mtxFile.add_dim("matrix", dist.size());
  NcVar *distance = mtxFile.add_var("distance", ncFloat, mtxdist);
  float *vdm = new float[dist.size()];
  for (int i = 0; i < dist.size(); ++i)
    vdm[i] = dist[i];
  distance->put(vdm, dist.size());
  delete[] vdm;

  mtxFile.close();
}
#endif // NETCDF

#ifdef _HDF5
// read hdf5 format distance matrix
bool Mdist::readmtxh5(const string &fname) {
  // check file
  if (!H5::H5File::isHdf5(fname.c_str()))
    return false;

  H5::H5File file(fname, H5F_ACC_RDONLY);
  // get the name list
  H5::DataSet nmset = file.openDataSet("Name");
  H5::DataSpace nmspace = nmset.getSpace();
  hsize_t nmdims[1];
  nmspace.getSimpleExtentDims(nmdims, NULL);
  H5::StrType nmtype(H5::PredType::C_S1, H5T_VARIABLE);
  H5::DataSpace mnspace(1, nmdims);
  vector<const char *> cstrs(nmdims[0]);
  nmset.read(cstrs.data(), nmtype, mnspace, nmspace);

  // set matrix size and the name
  resize(cstrs.size());
  for (size_t i = 0; i < cstrs.size(); ++i)
    name[i] = cstrs[i];

  // for distance
  H5::DataSet dset = file.openDataSet("Distance");
  H5::DataSpace dspace = dset.getSpace();
  hsize_t ddims[1];
  dspace.getSimpleExtentDims(ddims, NULL);
  H5::DataSpace mdspace(1, ddims);

  dset.read(dist.data(), H5::PredType::NATIVE_DOUBLE, mdspace, dspace);

  return true;
};

// write hdf5 format distance matrix
void Mdist::writemtxh5(const string &fname) {
  // open a h5 file to write
  H5::H5File file(fname, H5F_ACC_TRUNC);

  // for the name list
  vector<const char *> cstrs;
  for (auto &str : name)
    cstrs.push_back(str.c_str());
  hsize_t nmdim[1]{cstrs.size()};
  H5::DataSpace nmspace(1, nmdim);
  H5::StrType nmtype(H5::PredType::C_S1, H5T_VARIABLE);
  H5::DataSet nmset = file.createDataSet("Name", nmtype, nmspace);
  nmset.write(cstrs.data(), nmtype);

  // for the distance
  hsize_t ddim[1]{dist.size()};
  H5::DataSpace dspace(1, ddim);
  H5::FloatType dtype(H5::PredType::NATIVE_FLOAT);
  H5::DataSet dset = file.createDataSet("Distance", dtype, dspace);
  dset.write(dist.data(), H5::PredType::NATIVE_DOUBLE);

  file.close();
};
#endif // _HDF5

// readjust the distance matrix and their name by the index list
void Mdist::reduce(const vector<size_t> &ndx) {

  // put the name and distance in the tmp vector
  vector<double> tmpDist(dist);
  vector<string> tmpName(name);

  // resize the matrix
  resize(ndx.size());
  for (int i = 0; i < ng; ++i) {
    if (ndx[i] >= tmpName.size()) {
      cerr << "The index " << ndx[i] << " is out the size of matrix "
           << tmpName.size() << endl;
      exit(3);
    }

    name[i] = tmpName[ndx[i]];
    for (int j = 0; j <= i; ++j) {
      size_t itmp = ndx[i] < ndx[j] ? ndx[i] + ndx[j] * (ndx[j] - 1) / 2
                                    : ndx[j] + ndx[i] * (ndx[i] - 1) / 2;
      _setdist(j, i, tmpDist[itmp]);
    }
  }
}

// readjust the distance matrix and their name by the name list
void Mdist::reduce(const vector<string> &nlist) {

  // get the index map of the name
  unordered_map<string, size_t> mapMI;
  for (size_t i = 0; i < name.size(); ++i)
    mapMI[name[i]] = i;

  // resize the matrix
  vector<size_t> ndx;
  for (const auto &str : nlist)
    if (mapMI.find(str) == mapMI.end()) {
      cerr << "\nCannot find the genome " << str << "in the matrix file\n"
           << endl;
      exit(3);
    } else
      ndx.emplace_back(mapMI[str]);
  reduce(ndx);
};

// ... extend the matrix
void Mdist::extend(const vector<string> &nlist, const vector<double> &dd) {

  int orgSize = ng;
  int orgDistSize = dist.size();
  int newSize = ng + nlist.size();
  int newDistsize = newSize * (newSize + 1) / 2;

  if (dd.size() != newDistsize - orgDistSize) {
    cerr << "Error: the number of the append vector is unmatched!!" << endl;
    exit(4);
  }

  // resize the matrix
  resize(newSize);

  // add the name
  for (int i = 0; i < nlist.size(); ++i)
    name[orgSize + i] = nlist[i];

  // add the distance
  for (int i = 0; i < dd.size(); ++i)
    dist[orgDistSize + i] = dd[i];
};

// ... extend matrix by name list
void Mdist::extend(const vector<string> &nlist) {
  int orgSize = ng;
  int newSize = ng + nlist.size();

  resize(newSize);
  for (int i = 0; i < nlist.size(); ++i)
    name[orgSize + i] = nlist[i];
};

// ... translate genome to the index of the matrix
bool Mdist::name2ndx(const vector<string> &nm, vector<size_t> &ndx) const {

  // get the index map of the name
  unordered_map<string, size_t> mapMI;
  for (size_t i = 0; i < name.size(); ++i)
    mapMI[name[i]] = i;

  // resize the matrix
  ndx.resize(name.size());
  for (int i = 0; i < nm.size(); ++i) {
    if (mapMI.find(nm[i]) == mapMI.end()) {
      cerr << "\nCannot find the genome " << nm[i] << "in the matrix file\n"
           << endl;
      exit(3);
    } else {
      ndx[i] = mapMI[nm[i]];
    }
  }
  return true;
};

bool Mdist::ndx2name(const vector<size_t> &ndx, vector<string> &nm) const {
  nm.resize(ndx.size());
  for (int i = 0; i < ndx.size(); ++i)
    nm[i] = getname(i);
  return true;
};

// for check NAN distance
int Mdist::chkNAN(const vector<string> &gname,
                  vector<pair<size_t, size_t>> &nandist) const {

  vector<size_t> ndx;
  name2ndx(gname, ndx);
  return chkNAN(ndx, nandist);
};

int Mdist::chkNAN(const vector<size_t> &ndx,
                  vector<pair<size_t, size_t>> &nandist) const {

  for (size_t i = 0; i < ndx.size(); ++i) {
    if (ndx[i] >= ng) {
      cerr << "The index " << ndx[i] << " is out the size of matrix " << ng
           << endl;
      exit(3);
    }

    for (size_t j = 0; j < i; ++j) {
      if (std::isnan(getdist(ndx[i], ndx[j])))
        nandist.emplace_back(ndx[i], ndx[j]);
    }
  }

  return nandist.size();
};

int Mdist::chkAllNAN(vector<pair<size_t, size_t>> &nandist) const {
  for (size_t i = 0; i < ng; ++i) {
    for (size_t j = 0; j < i; ++j) {
      if (std::isnan(_getdist(j, i)))
        nandist.emplace_back(i, j);
    }
  }
  return nandist.size();
};

//.... the number of NAN of ith line
int Mdist::nNAN(size_t i) const {
  int n(0);
  size_t j = 0;
  for (; j < i; ++j) {
    if (std::isnan(_getdist(j, i)))
      ++n;
  }

  for (++j; j < ng; ++j) {
    if (std::isnan(_getdist(i, j)))
      ++n;
  }
  return n;
};

//.... number of NAN in the matrix
int Mdist::nNAN() const {
  int n(0);
  for (size_t i = 0; i < ng; ++i) {
    for (size_t j = i + 1; j < ng; ++j) {
      if (std::isnan(_getdist(i, j)))
        ++n;
    }
  }
  return n;
};

bool Mdist::isNAN(size_t i, size_t j) const {
  return std::isnan(getdist(i, j));
};

bool Mdist::isNAN(size_t m) const { return std::isnan(dist[m]); };

//.... has NAN in the ith line
bool Mdist::hasNAN(size_t i) const {
  size_t j = 0;
  for (; j < i; ++j) {
    if (std::isnan(_getdist(j, i)))
      return true;
  }

  for (++j; j < ng; ++j) {
    if (std::isnan(_getdist(i, j)))
      return true;
  }
  return false;
};

//.... has NAN in the matrix
bool Mdist::hasNAN() const {
  size_t m = msize();
  for (size_t i = 0; i < msize(); ++i) {
    if (std::isnan(dist[i]))
      return true;
  }
  return false;
};

// .. the infomation of matrix
string Mdist::info() const {

  string str("The dimension of the distance matrix is: ");
  str += to_string(ng);

  int n = nNAN();
  size_t m = msize();
  if (n == 0) {
    str += "\nAll distances are filled";
  } else if (n < m) {
    str += "\n" + to_string(n) + "/" + to_string(m) + " distances are blank";
  } else {
    str += "\nAll " + to_string(n) + " distances are blank";
  }

  return str;
}

// .. reset the size of the matrix
void Mdist::resize(size_t n) {
  ng = n;
  name.resize(n);
  dist.resize(n * (n - 1) / 2);
};

// option on sigle item
double Mdist::getdist(size_t i, size_t j) const {
  if (i < j) {
    if (j >= ng) {
      cerr << "Error: the index out of matrix" << endl;
      exit(4);
    }
    return _getdist(i, j);
  } else if (i > j) {
    if (i >= ng) {
      cerr << "Error: the index out of matrix" << endl;
      exit(4);
    }
    return _getdist(j, i);
  } else {
    if (j >= ng) {
      cerr << "Error: the index out of matrix" << endl;
      exit(4);
    }
    return 0.0;
  }
};

double Mdist::_getdist(size_t i, size_t j) const {
  return dist[i + (j - 1) * j / 2];
};

pair<size_t, size_t> Mdist::getIndex(size_t ndx) const {
  if (ndx >= dist.size()) {
    cerr << "Error: the index out of distance of matrix!" << endl;
    exit(4);
  }

  size_t i = 1;
  while (ndx > (i + 2) * (i - 1) / 2)
    i++;
  size_t j = ndx - (i - 1) * i / 2;

  return make_pair(j, i);
};

void Mdist::setdist(size_t i, size_t j, double d) {
  if (i < j) {
    if (j >= ng) {
      cerr << "Error: the index out of matrix" << endl;
      exit(4);
    }
    return _setdist(i, j, d);
  } else {
    if (i >= ng) {
      cerr << "Error: the index out of matrix" << endl;
      exit(4);
    }
    return _setdist(j, i, d);
  }
};

void Mdist::_setdist(size_t i, size_t j, double d) {
  dist[i + j * (j - 1) / 2] = d;
};

void Mdist::setdist(size_t i, double d) { dist[i] = d; };

string Mdist::getname(size_t i) const { return name[i]; };

void Mdist::setname(size_t i, const string &str) { name[i] = str; };

vector<string> Mdist::getNameList() const { return name; };

// format the genome name
void Mdist::cleanName() {
  for (int i = 0; i < ng; ++i) {
    name[i] = name[i].substr(name[i].find_last_of("/>") + 1);
  }
};

// assign distance by other distance matrix
void Mdist::assign(const Mdist &dm, vector<size_t> &hit) {

  // get the name-index map of the reference DM
  unordered_map<string, size_t> mapMI;
  for (size_t i = 0; i < dm.size(); ++i)
    mapMI[dm.getname(i)] = i;

  // find the match items
  vector<pair<size_t, size_t>> mlist;
  for (size_t i = 0; i < ng; ++i) {
    auto iter = mapMI.find(name[i]);
    if (iter != mapMI.end()) {
      pair<size_t, size_t> item(i, iter->second);
      mlist.emplace_back(item);
      hit.emplace_back(i);
    }
  }

  // renew the match distance
  for (auto ita = mlist.begin(); ita != mlist.end(); ++ita) {
    for (auto itb = ita + 1; itb != mlist.end(); ++itb) {
      setdist(ita->first, itb->first, dm.getdist(ita->second, itb->second));
    }
  }
};

// assign distance by other distance matrix
void Mdist::assign(const Mdist &dm) {

  // get the name-index map of the reference DM
  unordered_map<string, size_t> mapMI;
  for (size_t i = 0; i < dm.size(); ++i)
    mapMI[dm.getname(i)] = i;

  // find the match items
  vector<pair<size_t, size_t>> mlist;
  for (size_t i = 0; i < ng; ++i) {
    if (hasNAN(i)) {
      auto iter = mapMI.find(name[i]);
      if (iter != mapMI.end())
        mlist.emplace_back(make_pair(i, iter->second));
    }
  }

  // renew the match distance
  for (auto ita = mlist.begin(); ita != mlist.end(); ++ita) {
    for (auto itb = ita + 1; itb != mlist.end(); ++itb) {
      setdist(ita->first, itb->first, dm.getdist(ita->second, itb->second));
    }
  }
};

// assign distance by alist of distance matrix
void Mdist::assign(const string &refdm) {
  // assign distance by reference DM
  vector<string> dmlist;
  separateWord(dmlist, refdm);
  for (auto &fnm : dmlist) {
    Mdist xdm;
    bool readxdm;
    readxdm = xdm.readmtx(fnm);
    if (readxdm) {
      xdm.cleanName();
      assign(xdm);
      if (!hasNAN())
        return;
    }
  }
};

// global options
size_t Mdist::size() const { return ng; };

size_t Mdist::msize() const { return ng * (ng - 1) / 2; };

size_t Mdist::capacity() const {
  return dist.capacity() * sizeof(double) + sizeof(ng);
};
