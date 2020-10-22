/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai, China.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-10-22 20:45:30
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-10-22 20:45:30
 */

#include "distmatrix.h"
#include "stringOpt.h"
#include <array>
#include <map>
#include <regex>
using namespace std;

void usage(string &program) {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -i input ]    Input matrix file, format determined by suffix\n"
       << " [ -o output ]   output matrix file, format determined by suffix\n"
       << " [ -h ]          Display this information\n"
       << endl;
  exit(1);
}

int main(int argc, char *argv[]) {

  // get the name of file
  string infile; 
  string outfile = "mdist.txt";
  string program = argv[0];

  char ch;
  while ((ch = getopt(argc, argv, "i:o:h")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'h':
      usage(program);
    case '?':
      usage(program);
    }
  }

  Mdist dm;
  dm.readmtx(infile);
  dm.writemtx(outfile);
}

