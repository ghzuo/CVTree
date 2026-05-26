/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-17 09:31:30
 */

#include "rectangleMatrix.h"

// construct header by reading file
RectangleHeader::RectangleHeader(const string &fname) {
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL)
    throw runtime_error("Cannot open file for reading: " + gzfile);

  // read header
  read(fp);
  gzclose(fp);
};

void RectangleHeader::read(gzFile &fp) {
  // get the genome name
  gzline(fp, rowName);
  gzline(fp, colName);
  // get the matrix size
  gzread(fp, (char *)&(nrow), sizeof(nrow));
  gzread(fp, (char *)&(ncol), sizeof(ncol));
  gzread(fp, (char *)&(nsize), sizeof(nsize));
};

void RectangleHeader::write(gzFile &fp) const {
  // write the genome information
  string str = rowName + "\n" + colName + "\n";
  gzputs(fp, str.c_str());

  // write the size of CVArray
  gzwrite(fp, &(nrow), sizeof(nrow));
  gzwrite(fp, &(ncol), sizeof(ncol));
  gzwrite(fp, &(nsize), sizeof(nsize));
};

ostream &operator<<(ostream &os, const RectangleHeader &hd) {
  os << hd.rowName << ": " << hd.nrow << "\n" << hd.colName << ": " << hd.ncol;
  return os;
};

// set row name and col name
void Mrectangle::resetByHeader(const RectangleHeader &hd, float d0) {
  header = hd;
  vector<float> tmp(hd.nrow * hd.ncol, d0);
  data.swap(tmp);
}

// option on sigle item
float Mrectangle::get(size_t i, size_t j) const {
  if (i >= header.nrow || j >= header.ncol)
    throw out_of_range("Index " + outIndex(i, j) + " out of " +
                       outIndex(header.nrow, header.ncol) +
                       " in Mrectangle::get()");
  return _get(i, j);
};

float Mrectangle::_get(size_t i, size_t j) const { return data[index(i, j)]; };

void Mrectangle::set(size_t i, size_t j, float val) {
  if (i >= header.nrow || j >= header.ncol)
    throw out_of_range("Index " + outIndex(i, j) + " out of " +
                       outIndex(header.nrow, header.ncol) +
                       " in Mrectangle::set()");
  _set(i, j, val);
};

void Mrectangle::_set(size_t i, size_t j, float val) { data[index(i, j)] = val; };

void Mrectangle::add(size_t i, size_t j, float val) {
  if (i >= header.nrow || j >= header.ncol)
    throw out_of_range("Index " + outIndex(i, j) + " out of " +
                       outIndex(header.nrow, header.ncol) +
                       " in Mrectangle::add()");
  _add(i, j, val);
};

void Mrectangle::_add(size_t i, size_t j, float val) {
  data[index(i, j)] += val;
};

size_t Mrectangle::index(size_t i, size_t j) const {
  return header.ncol * i + j;
};

string Mrectangle::outIndex(size_t i, size_t j) const {
  return "(" + to_string(i) + ", " + to_string(j) + ")";
};

pair<size_t, size_t> Mrectangle::index(size_t ndx) const {
  pair<size_t, size_t> tmp;
  tmp.first = ndx / header.ncol;
  tmp.second = ndx % header.ncol;
  return tmp;
};

// .. the infomation of matrix
string Mrectangle::info() const {
  return "The dimension of the distance matrix is: " +
         std::to_string(header.nrow) + "x" + std::to_string(header.ncol);
}

void Mrectangle::write(const string &fname, float mindist) {
  // open and test file
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "wb")) == NULL)
    throw runtime_error("Cannot open file for reading: " + gzfile);

  // write data
  if (mindist < 0.0) {
    // for full matrix
    header.write(fp);
    gzwrite(fp, data.data(), data.size() * sizeof(data[0]));
  } else {
    // for threshold matrix
    vector<pair<size_t, float>> vec;
    for (size_t i = 0; i < data.size(); ++i)
      if (data[i] >= mindist)
        vec.emplace_back(i, data[i]);
    header.nsize = (long) vec.size();
    header.write(fp);
    gzwrite(fp, vec.data(), vec.size() * sizeof(vec[0]));
  }

  // close file
  gzclose(fp);
};

void Mrectangle::read(const string &fname) {
  try {
    // open file to read
    gzFile fp;
    string gzfile = addsuffix(fname, ".gz");
    if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL)
      throw runtime_error("Cannot open file for reading: " + gzfile);

    // read the header
    header.read(fp);

    // read data
    size_t dsize = header.nrow * header.ncol;
    data.resize(dsize);
    if (header.nsize < 0) {
      gzread(fp, (char *)data.data(), dsize * sizeof(data[0]));
    } else {
      vector<pair<size_t, float>> vec(header.nsize);
      gzread(fp, (char *)vec.data(), header.nsize * sizeof(vec[0]));
      for (auto &v : vec)
        data[v.first] = v.second;
    } 

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

ostream &operator<<(ostream &os, const Mrectangle &sm) {
  os << sm.header << "\n========= the similarity =====" << endl;
  for (size_t i = 0; i < sm.header.nrow; ++i) {
    for (size_t j = 0; j < sm.header.ncol; ++j)
      os << sm.get(i, j) << " ";
    os << endl;
  }
  return os;
}
