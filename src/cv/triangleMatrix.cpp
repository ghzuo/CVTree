/*
 * Copyright (c) 2026
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2026-05-15 22:54:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-16 13:26:19
 */


#include "triangleMatrix.h"

TriangleHeader::TriangleHeader(const string &fname) {
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL)
    throw runtime_error("Cannot open file for reading: " + gzfile);

  read(fp);
  gzclose(fp);
};

void TriangleHeader::read(gzFile &fp) {
  gzline(fp, name);
  gzread(fp, (char *)&(nsize), sizeof(nsize));
};

void TriangleHeader::write(gzFile &fp) const {
  string str = name + "\n";
  gzputs(fp, str.c_str());
  gzwrite(fp, &(nsize), sizeof(nsize));
};

ostream &operator<<(ostream &os, const TriangleHeader &hd) {
  os << hd.name << ": " << hd.nsize;
  return os;
};

void Mtriangle::reset(long n, float d0) {
  header.nsize = n;
  vector<float> tmp(n * (n + 1) / 2, d0);
  data.swap(tmp);
}

// Signed integer versions with negative value check
float Mtriangle::get(long i, long j) const {
  if (i < 0 || j < 0)
    throw out_of_range("Negative index (" + to_string(i) + ", " + to_string(j) +
                       ") in Mtriangle::get()");
  return get((size_t)i, (size_t)j);
};

void Mtriangle::set(long i, long j, float val) {
  if (i < 0 || j < 0)
    throw out_of_range("Negative index (" + to_string(i) + ", " + to_string(j) +
                       ") in Mtriangle::set()");
  set((size_t)i, (size_t)j, val);
};

void Mtriangle::add(long i, long j, float val) {
  if (i < 0 || j < 0)
    throw out_of_range("Negative index (" + to_string(i) + ", " + to_string(j) +
                       ") in Mtriangle::add()");
  add((size_t)i, (size_t)j, val);
};

// Unsigned integer versions with boundary check
float Mtriangle::get(size_t i, size_t j) const {
  return data[index(i, j)];
};

void Mtriangle::set(size_t i, size_t j, float val) {
  data[index(i, j)] = val;
};

void Mtriangle::add(size_t i, size_t j, float val) {
  data[index(i, j)] += val;
};

// Internal versions without boundary check
float Mtriangle::_get(size_t i, size_t j) const {
  return data[_index(i, j)];
};

void Mtriangle::_set(size_t i, size_t j, float val) {
  data[_index(i, j)] = val;
};

void Mtriangle::_add(size_t i, size_t j, float val) {
  data[_index(i, j)] += val;
};

// Index with size comparison and boundary check
size_t Mtriangle::index(size_t i, size_t j) const {
  if (i < j)
    std::swap(i, j);
  if (i >= header.nsize)
    throw out_of_range("Index (" + to_string(i) + ", " + to_string(j) +
                       ") out of " + to_string(header.nsize) + " in Mtriangle::index()");
  return _index(i, j);
};

// Internal index without any check, assumes i >= j
size_t Mtriangle::_index(size_t i, size_t j) const {
  return i * (i + 1) / 2 + j;
};

string Mtriangle::outIndex(size_t i, size_t j) const {
  return "(" + to_string(i) + ", " + to_string(j) + ")";
};

string Mtriangle::info() const {
  return "The dimension of the triangle matrix is: " +
         std::to_string(header.nsize) + "x" + std::to_string(header.nsize);
}

void Mtriangle::write(const string &fname, float mindist) {
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "wb")) == NULL)
    throw runtime_error("Cannot open file for writing: " + gzfile);

  if (mindist < 0.0) {
    header.write(fp);
    gzwrite(fp, data.data(), data.size() * sizeof(data[0]));
  } else {
    vector<pair<size_t, float>> vec;
    for (size_t i = 0; i < data.size(); ++i)
      if (data[i] >= mindist)
        vec.emplace_back(i, data[i]);
    header.write(fp);
    gzwrite(fp, vec.data(), vec.size() * sizeof(vec[0]));
  }

  gzclose(fp);
};

void Mtriangle::read(const string &fname) {
  try {
    gzFile fp;
    string gzfile = addsuffix(fname, ".gz");
    if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL)
      throw runtime_error("Cannot open file for reading: " + gzfile);

    header.read(fp);

    size_t dsize = header.nsize * (header.nsize + 1) / 2;
    data.resize(dsize);
    gzread(fp, (char *)data.data(), dsize * sizeof(data[0]));

    gzclose(fp);
  } catch (std::exception &e) {
    cerr << "Error reading file: " << fname << "\n" << e.what() << endl;
    exit(1);
  } catch (...) {
    cerr << "Error reading file: " << fname << "\n" << endl;
    exit(1);
  }
};

ostream &operator<<(ostream &os, const Mtriangle &sm) {
  os << sm.header << "\n========= the triangle matrix =====" << endl;
  for (size_t i = 0; i < sm.header.nsize; ++i) {
    for (size_t j = 0; j < sm.header.nsize; ++j)
      os << sm.get(i, j) << " ";
    os << endl;
  }
  return os;
}