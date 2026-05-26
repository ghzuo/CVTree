/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-05 Tuesday 19:41:23
 */

#include "kstring.h"

long Kstr::kmax() {
  float logbs = log2(Letter::nbase);
  float sz = sizeof(Kstr) * 8;
  return sz / logbs;
};

Kstr::Kstr() : ks(0){};

Kstr::Kstr(unsigned long s) { ks = s; };

Kstr::Kstr(const vector<Letter> &str) : ks(0) {
  for (auto &c : str)
    ks = ks * Letter::nbase + c.e;
};

string Kstr::decode() const {
  if(Letter::nbase == 0)
    return "code:" + to_string(ks);

  string str;
  for (unsigned long s = ks; s != 0; s /= Letter::nbase)
    str += Letter::decode(s % Letter::nbase);
  reverse(str.begin(), str.end());
  return str;
}

size_t Kstr::length() const {
  size_t n(0);
  unsigned long unit(1);
  while (unit <= ks) {
    ++n;
    unit *= Letter::nbase;
  }
  return n;
};

void Kstr::append(Letter c) {
  ks *= Letter::nbase;
  ks += c.e;
};

void Kstr::addhead(Letter c) {
  unsigned long unit(1);
  while (unit <= ks)
    unit *= Letter::nbase;
  ks += (unit * c.e);
};

void Kstr::choptail() { ks /= Letter::nbase; };

void Kstr::behead() {
  unsigned long unit(1);
  while (unit <= ks)
    unit *= Letter::nbase;
  unit /= Letter::nbase;
  ks %= unit;
};

void Kstr::forward(Letter c) {
  behead();
  append(c);
}

void Kstr::backward(Letter c) {
  choptail();
  addhead(c);
}

bool Kstr::operator<(const Kstr &r) const { return ks < r.ks; };

bool Kstr::operator>(const Kstr &r) const { return ks > r.ks; };

bool Kstr::operator==(const Kstr &r) const { return ks == r.ks; };

bool Kstr::operator!=(const Kstr &r) const { return ks != r.ks; };

bool Kstr::operator()(const Kstr &a, const Kstr &b) const { return a < b; };

ostream &operator<<(ostream &os, const Kstr &ks) {
  os << ks.decode();
  return os;
};


// operation for kstr vector
// used for the missing K string
string readvk(const string &filename, vector<Kstr> &vk) {
  gzFile fp;
  if ((fp = gzopen(filename.c_str(), "rb")) == NULL) {
    cerr << "CV file not found: \"" << filename << '"' << endl;
    exit(1);
  }

  size_t size;
  gzread(fp, (char *)&size, sizeof(size_t));

  vk.resize(size);
  gzread(fp, (char *)vk.data(), size * sizeof(Kstr));

  string charSet;
  gzline(fp, charSet);
  gzclose(fp);

  return charSet;
}

void writevk(const string &file, const vector<Kstr> &vk) {

  unsigned long size = vk.size();

  gzFile fp;
  if ((fp = gzopen(file.c_str(), "wb")) == NULL) {
    cerr << "Error happen on write cvfile: " << file << endl;
    exit(1);
  }

  gzwrite(fp, &size, sizeof(unsigned long));
  gzwrite(fp, vk.data(), size * sizeof(Kstr));

  gzputs(fp, Letter::decMap.c_str());
  gzclose(fp);
}


