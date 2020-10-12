/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-04-19 11:37:42
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 22:02:45
 */

#include "readgenome.h"

GeneType::GeneType(const string &str) { init(str); };

void GeneType::init(const string &str) {
  for (int i = 0; i < 128; ++i)
    mc[i] = 'A';
  if (str.compare("faa") == 0)
    aainit();
  else if (str.compare("ffn") == 0)
    nainit();
  else if (str.compare("fna") == 0)
    nainit();
  else {
    cerr << "unknown genome file type!\n" << endl;
    exit(1);
  }
};

void GeneType::aainit() {
  string aa = "ACDEFGHIKLMNPQRSTVWY";
  for (auto &c : aa) {
    mc[c] = c;
    letters.emplace_back(c);
  }

  // for the unfomrated aa
  mc['B'] = 'D';
  mc['U'] = 'C';
  mc['X'] = 'G';
  mc['Z'] = 'E';
};

void GeneType::nainit() {
  string na = "ACGT";
  for (auto &c : na) {
    mc[c] = c;
    letters.emplace_back(c);
  }

  // for the unfomrated na
  mc['R'] = 'A'; // AG
  mc['Y'] = 'T'; // CT
  mc['M'] = 'A'; // AC
  mc['K'] = 'T'; // TG
  mc['S'] = 'C'; // CG
  mc['W'] = 'A'; // AT

  mc['B'] = 'T'; // TCG
  mc['D'] = 'A'; // ATG
  mc['H'] = 'A'; // ATC
  mc['V'] = 'A'; // ACG

  mc['N'] = 'A'; // ACGT
  mc['X'] = 'A'; // ACGT
};

size_t GeneType::readgene(string &file, Genome &genome) const {
  ifstream infile(file.c_str());
  if (!infile) {
    cerr << "Cannot found the input file " << file << endl;
    exit(4);
  }
  // cout << " Read file: " << file << endl;

  for (string line; getline(infile, line);) {
    line = trim(line);
    if (line.empty() || line[0] == ';') {

    } else if (line[0] == '>') {
      genome.emplace_back();
    } else if(!genome.empty()){
      genome.back().append(line);
    }
  }
  infile.close();

  if(genome.size() <= 0){
    cerr << "The genome of " << file << " is empty!" << endl;
    exit(5);
  }

  size_t len(0);
  for (auto &gene : genome) {
    if (gene.size() == 0) {
      cerr << "Some empty gene in your genome file: " << file << endl;
    } else {
      checkgene(gene);
      len += gene.size();
    }
  }
  return len;
}

void GeneType::checkgene(string &str) const {
  if (*(str.rbegin()) == '*' || *(str.rbegin()) == '-')
    str.pop_back();
  str = toUpper(str);
  for (auto &c : str)
    c = mc[c];
};
