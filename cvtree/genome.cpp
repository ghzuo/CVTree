/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: Thu May 09 2024
 */

#include "genome.h"

char Letter::encMap[128];
string Letter::decMap;
long Letter::nbase;
vector<Letter> Letter::nonVoid;

void Letter::init(const string &str, const string &cgstr) {
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

  // use the upper case letters
  setLowercaseUpper();

  // for simplify
  if (!cgstr.empty())
    simplify(cgstr);

  // update the nbase and nonVoid
  renewAuxiliary();
};

void Letter::aainit() {
  // set decode and encMap
  decMap = "_ACDEFGHIKLMNPQRSTVWY";
  decMap4encMap();

  // for the unfomrated aa
  encMap['B'] = encMap['D'];
  encMap['U'] = encMap['C'];
  encMap['X'] = encMap['G'];
  encMap['Z'] = encMap['E'];
};

void Letter::nainit() {
  string decMap = "_ACGT";
  decMap4encMap();

  // for RNA Sequence
  encMap['U'] = encMap['T'];

  // for the unfomrated na
  encMap['R'] = encMap['A']; // AG
  encMap['Y'] = encMap['T']; // CT
  encMap['M'] = encMap['A']; // AC
  encMap['K'] = encMap['T']; // TG
  encMap['S'] = encMap['C']; // CG
  encMap['W'] = encMap['A']; // AT

  encMap['B'] = encMap['T']; // TCG
  encMap['D'] = encMap['A']; // ATG
  encMap['H'] = encMap['A']; // ATC
  encMap['V'] = encMap['A']; // ACG

  encMap['N'] = encMap['A']; // ACGT
  encMap['X'] = encMap['A']; // ACGT
};

void Letter::decMap4encMap() {
  // set default to 1, e.g. the first nonVoid
  for (auto &c : encMap)
    c = 1;

  // get the code of letter in decMap
  char ndx = 0;
  for (auto &c : decMap) {
    encMap[c] = ndx++;
  }
}

void Letter::simplify(const string &str) {
  // get the group list
  vector<string> cglist;
  separateWord(cglist, toUpper(str));

  //.... get the mapping of code: cgmap and the new cmap
  // get letter list in input info
  map<char, char> cgmap{make_pair(0, 0)};
  string newcmap = "_";
  char ndx = 1;
  for (auto &cg : cglist) {
    newcmap += cg.front();
    for (auto c : cg) {
      cgmap[encMap[c]] = ndx;
    }
    ndx++;
  }
  // add the lost letters in the old cmap
  for (auto c : decMap) {
    if (cgmap.find(encMap[c]) == cgmap.end()) {
      newcmap += c;
      cgmap[encMap[c]] = ndx++;
    }
  }

  // update encMap
  for (auto &c : encMap)
    c = cgmap[c];

  // update decMap
  decMap = newcmap;
};

void Letter::setLowercaseUpper(){
    // use the upper case letters
  for (long i = 97; i < 123; ++i)
    encMap[i] = encMap[i - 32];
}

void Letter::renewAuxiliary() {
  // set the nbase
  nbase = decMap.size();
  nonVoid.clear();
  for (char c : decMap) {
    if (c != '_')
      nonVoid.emplace_back(c);
  }
}

// initial letter by letters directly, only for decode method
void Letter::initByStr(const string &str) {
  decMap = toUpper(str);
  decMap4encMap();
  setLowercaseUpper();
  renewAuxiliary();
};

// decode the number
char Letter::decode(size_t n){
  return decMap.at(n);
};

//.. the non-static member method
Letter::Letter(char c) : e(encMap[c]){};

ostream &operator<<(ostream &os, const Letter c) {
  os << Letter::decMap[c.e];
  return os;
};

/********************************************************************************
 * @brief Struct Gene
 *
 ********************************************************************************/
Gene::Gene(const string &str) : head(str){};

size_t Gene::size() const { return code.size(); };

Letter Gene::operator[](size_t i) const { return code[i]; };

vector<Letter> Gene::substr(size_t pos, size_t length) const {
  vector<Letter> seg;
  seg.reserve(length);
  for (auto i = 0; i < length; ++i)
    seg.emplace_back(code[i]);
  return move(seg);
};

void Gene::translate() {
  code.reserve(seq.size());
  for (char c : seq)
    code.emplace_back(Letter(c));
}

ostream &operator<<(ostream &os, const Gene &g) {
  os << g.head << "\n";
  for (auto &c : g.code)
    os << c;
  return os;
};

/********************************************************************************
 * @brief Genome
 *
 ********************************************************************************/
size_t readFasta(const string &file, Genome &genome) {
  ifstream infile(file.c_str());
  if (!infile) {
    cerr << "Cannot found the input file " << file << endl;
    exit(4);
  }

  for (string line; getline(infile, line);) {
    line = trim(line);
    if (line.empty() || line[0] == ';') {

    } else if (line[0] == '>') {
      genome.emplace_back(Gene(line));
    } else if (!genome.empty()) {
      genome.back().seq.append(line);
    }
  }
  infile.close();

  // check read data
  if (genome.empty()) {
    cerr << "The format of " << file << " is not in Fasta!" << endl;
    exit(5);
  }

  size_t len(0);
  for (auto &gene : genome) {
    if (gene.seq.size() == 0) {
      cerr << "Some empty gene in your genome file: " << file << endl;
    } else {
      if (gene.seq.back() == '*' || gene.seq.back() == '-')
        gene.seq.pop_back();
      gene.translate();
      len += gene.size();
    }
  }

  return len;
}

void writeFasta(const string &file, const Genome &genome) {
  ofstream ofile(file.c_str());
  if (!ofile) {
    cerr << "Cannot open file to write: " << file << endl;
    exit(4);
  }

  for (auto &g : genome)
    ofile << g << endl;
};

Genome sampleGenome(const Genome &org, SampleMeth *smeth) {
  vector<long> ndxlist = (*smeth)(org.size());
  Genome gs;
  for (auto ndx : ndxlist) {
    gs.emplace_back(org[ndx]);
  }
  return gs;
};

/********************************************************************************
 * @brief for testing this module
 *
 ********************************************************************************/
#ifdef TESTING
int main(int argc, char *argv[]) {
  Letter::init("fna", "AG,CT");
  string file(argv[1]);

  Genome gn;
  readFasta(file, gn);
  for (auto &g : gn) {
    cout << g << "\n" << g.seq << endl;
  }
}
#endif // TESTING