/*
 * Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-04-26 09:00:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-12 16:07:09
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

Args::Args(int argc, char **argv)
    : maintree(""), outfile(""),unroot(false) {

  program = argv[0];
  string listfile = "treelist.txt";

  char ch;
  while ((ch = getopt(argc, argv, "i:o:l:uh")) != -1) {
    switch (ch) {
    case 'i':
      maintree = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'l':
      listfile = optarg;
      break;
    case 'u':
      unroot = true;
      break;
    case 'h':
      usage();
    case '?':
      usage();
    }
  }

  readlist(listfile, boottree);
  if(maintree.empty()){
    maintree = boottree.front();
    boottree.erase(boottree.begin());
  }

  if(outfile.empty())
    outfile = addnamelabel(maintree, "-bootstrap");
}

void Args::usage() {
  cerr
      << "\nProgram Usage: \n\n"
      << program << "\n"
      << " [ -i <maintree.nwk> ]  the main tree file for bootstrap.\n"
      << "                        default: the first at the treelist file\n"
      << " [ -l treelist.txt ]    the tree list file to, default: treelist.txt\n"
      << " [ -u ]                 whether trees are unroot, default: no\n"
      << " [ -h ]                 display this information\n"
      << endl;
  exit(1);
}
