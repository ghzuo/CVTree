/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-11 09:46:36
 */

#include "kstring.h"
#include "stringOpt.h"

void usage(string &program) {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << "  -i <cvfile>    input file name\n"
       << " [ -g ]          the type of genome file, only for old format cvfile\n"
       << " [ -n ]          output the number code, default: the letters\n"
       << " [ -h ]          Display this information\n"
       << endl;
  exit(1);
}

int main(int argc, char *argv[]) {

  string gtype;
  string cgstr;
  string infile;
  string program = argv[0];
  bool kstr = true;

  char ch;
  while ((ch = getopt(argc, argv, "i:g:n")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'g':
      gtype = optarg;
      break;
    case 'n':
      kstr = false;
      break;
    case 'h':
      usage(program);
    case '?':
      usage(program);
    }
  }

  // read cv
  CVvec cv;
  pair<double, string> info = readcv(infile, cv);

  // initial the char set
  if(!gtype.empty()){
    Letter::init(gtype, cgstr);
  }else if(!info.second.empty()){
    Letter::initByStr(info.second);
  }else{
    cerr << "Please input the type of genome: faa/fna/ffn" << endl;
    exit(2);
  }


  // output data
  cout << "The inner of CV: " << info.first << endl;
  cout << "The size  of CV: " << cv.size() << endl;
  cout << "The char set of CV: " << Letter::decMap << endl;

  for (const auto &cd : cv)
    cout << cd.first << "\t" << cd.second << endl;
}
