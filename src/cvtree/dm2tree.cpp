/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-14 Sunday 14:37:42
 */

#include "dm2tree.h"

int main(int argc, char *argv[]) {
  // set a timer
  Timer mytimer;

  // get the input arguments
  Args myargs(argc, argv);

  // init the distance matrix and name list
  Mdist dm;
  dm.readmtx(myargs.distfile);
  theInfo("Read the distance matrix with " + to_string(dm.size()) + " items");

  // made the star tree by listfile
  // if no, all items in distance matrix are used
  if (!myargs.splist.empty()){
    dm.reduce(myargs.splist);
    theInfo(to_string(dm.size()) + " items are selected for tree");
  }

  // do the NJ algorithm and return the NJ tree
  Node* aTree = Node::initial();
  myargs.meth->execute(aTree, dm);

  // output the Tree
  ofstream nwk(myargs.outfile.c_str());
  (*aTree).outnwk(nwk);
  nwk.close();
}

Args::Args(int argc, char **argv) : distfile("infile"), outfile("Tree.nwk") {

  string wkdir("");
  string listfile;
  string methStr("NJ");
  string sclStr("");

  argparse::ArgumentParser parser("dm2tree", "1.0", 
                                  argparse::default_arguments::help);
  
  parser.add_argument("-i", "--input")
      .help("Selection index list of the distance matrix")
      .nargs(1)
      .store_into(listfile);
      
  parser.add_argument("-d", "--distance-matrix")
      .help("Input distance matrix")
      .default_value(distfile)
      .nargs(1)
      .store_into(distfile);
      
  parser.add_argument("-o", "--output")
      .help("Output newick tree")
      .default_value(outfile)
      .nargs(1)
      .store_into(outfile);
      
  parser.add_argument("-D", "--directory")
      .help("Working directory")
      .nargs(1)
      .store_into(wkdir);
      
  parser.add_argument("-m", "--method")
      .help("Select the tree method, NJ or UPGMA")
      .default_value(methStr)
      .nargs(1)
      .store_into(methStr);
    
  parser.add_argument("-T", "--dist-tranform")
      .help("Select the distance scaling method: log, sqrt, or cbrt")
      .choices("log", "sqrt", "cbrt")
      .nargs(1)
      .store_into(sclStr);
      
  parser.add_argument("-q", "--quiet")
      .help("Run command in quiet mode")
      .flag()
      .action([&](const auto &) { theInfo.quiet = true; });
      
  parser.add_description("Convert distance matrix to tree");
  
  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  if (wkdir != "") {
    addsuffix(wkdir, '/');
    distfile = wkdir + distfile;
    outfile = wkdir + outfile;
  }

  // read list file and name map
  if (!listfile.empty())
    readlist(listfile, splist, 1);

  meth = TreeMeth::create(methStr);
}


