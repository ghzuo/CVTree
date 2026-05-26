/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-13 Saturday 20:55:11
 */

#include "distmatrix.h"
#include "kit.h"
#include <argparse/argparse.hpp>
#include <array>
#include <map>
#include <regex>
using namespace std;

struct Args {
  string infile;
  string outfile;
  string namelist;
  bool check;
  bool input;

  Args(int argc, char **argv)
      : outfile("mdist.txt"), check(false), input(false) {

    argparse::ArgumentParser parser("mconv", "1.0",
                                    argparse::default_arguments::help);

    parser.add_argument("-i", "--input")
        .help("Input matrix file, format determined by suffix")
        .required()
        .nargs(1)
        .store_into(infile);

    parser.add_argument("-o", "--output")
        .help("Output matrix file, format determined by suffix")
        .default_value(string("mdist.txt"))
        .nargs(1)
        .store_into(outfile);

    parser.add_argument("-n", "--namelist")
        .help("Output name list of dist matrix")
        .nargs(1)
        .store_into(namelist);

    parser.add_argument("-r", "--rename")
        .help("Set the name of matrix by the name list file")
        .nargs(1)
        .action([&](const auto &val) {
          namelist = val;
          input = true;
        });

    parser.add_argument("-C", "--check")
        .help("Check matrix information")
        .flag()
        .store_into(check);

    parser.add_description("Convert matrix format");

    try {
      parser.parse_args(argc, argv);
    } catch (const std::exception &e) {
      std::cerr << e.what() << std::endl;
      std::cout << parser;
      exit(1);
    }
  }
};

int main(int argc, char *argv[]) {

  Args args(argc, argv);

  Mdist dm;
  dm.readmtx(args.infile);
  theInfo("Read Distance Matrix from " + args.infile);
  if (args.check) {
    theInfo(dm.info());
  } else {
    theInfo("The dimension of distance matrix is: " + to_string(dm.size()));
  }

  if (!args.namelist.empty()) {
    // set the name of matrix by name list file
    vector<string> nmlist;
    if (args.input) {
      readlist(args.namelist, nmlist);
      if (nmlist.size() != dm.size()) {
        cerr << "The length of name list is not match size of matrix" << endl;
        exit(3);
      }
      for (size_t i = 0; i < nmlist.size(); ++i) {
        dm.setname(i, nmlist[i]);
      }
      dm.writemtx(args.outfile);
    } else {
      // read the name list from matrix
      ofstream fnm(args.namelist.c_str());
      nmlist = dm.getNameList();
      fnm << strjoin(nmlist.begin(), nmlist.end(), "\n") << endl;
      fnm.close();
    }
  } else {
    // convert the format of distance matrix
    dm.writemtx(args.outfile);
  }
}
