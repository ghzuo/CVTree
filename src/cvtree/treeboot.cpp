/*
 * Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-04-26 09:00:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-13 Saturday 22:34:44
 */

#include "treeboot.h"

int main(int argc, char *argv[]) {

  // get the input arguments
  Args myargs(argc, argv);

  /// read the first tree
  MarkNode mTree;
  mTree.innwk(myargs.maintree);

  /// boot the trees
  mTree.bootTree(myargs.boottree);

  // output result
  mTree.outnwk(myargs.outfile);
}

Args::Args(int argc, char **argv) : maintree(""), outfile(""), unroot(false) {

  string listfile = "treelist.txt";

  argparse::ArgumentParser parser("treeboot", "1.0",
                                  argparse::default_arguments::help);

  parser.add_argument("-i", "--main-tree")
      .help("The main tree file for bootstrap, \ndefault: the first at the "
            "treelist file")
      .nargs(1)
      .store_into(maintree);

  parser.add_argument("-o", "--output")
      .help("Output file for the bootstrap tree")
      .nargs(1)
      .store_into(outfile);

  parser.add_argument("-l", "--tree-list")
      .help("The tree list file")
      .default_value(listfile)
      .nargs(1)
      .store_into(listfile);

  parser.add_argument("-u", "--unroot")
      .help("Whether trees are unroot")
      .flag()
      .store_into(unroot);

  parser.add_description("Bootstrap tree from a list of trees");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  // If listfile was not provided via -l, read the default list
  readlist(listfile, boottree);

  if (maintree.empty()) {
    maintree = boottree.front();
    boottree.erase(boottree.begin());
  }

  if (outfile.empty())
    outfile = addnamelabel(maintree, "-bootstrap");
}
