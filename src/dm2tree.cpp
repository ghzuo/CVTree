/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-11-14 11:46:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 21:57:25
 */

#include "dm2tree.h"
Info theInfo;

int main(int argc, char *argv[]) {
  // set a timer
  Timer mytimer;

  // get the input arguments
  Args myargs(argc, argv);

  // init the distance matrix and name list
  Mdist dm;
  dm.readmtx(myargs.distfile, myargs.netcdf);

  // made the star tree by listfile
  // if no, all items in distance matrix are used
  if (!myargs.splist.empty())
    dm.reduce(myargs.splist);

  // do the NJ algorithm and return the NJ tree
  Node *aTree = neighborJoint(dm);

  // output the Tree
  ofstream nwk(myargs.outfile.c_str());
  (*aTree).outnwk(nwk);
  nwk.close();
}

Args::Args(int argc, char **argv)
    : distfile("infile"), outfile("Tree.nwk"), netcdf(false) {

  program = argv[0];
  string wkdir("");
  string listfile;

  char ch;
  while ((ch = getopt(argc, argv, "i:d:o:D:Cqh")) != -1) {
    switch (ch) {
    case 'i':
      listfile = optarg;
      break;
    case 'd':
      distfile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'D':
      wkdir = optarg;
      break;
    case 'C':
      netcdf = true;
      break;
    case 'q':
      theInfo.info = false;
      break;
    case 'h':
      usage();
    case '?':
      usage();
    }
  }

  if (wkdir != "") {
    addsuffix(wkdir, '/');
    distfile = wkdir + distfile;
    outfile = wkdir + outfile;
  }

  if (!listfile.empty())
    readlist(listfile, splist);
}

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -d infile ]     Input distance matrix, defaut: dist.matrix\n"
       << " [ -o Tree.nwk ]   Output newick tree, defaut: Tree.nwk\n"
       << " [ -i <list> ]     selection index list of the distance matrix,\n"
       << "                   if no defined, whole distance matrix are used\n"
       << " [ -C ]            Use the netcdf input format, default false\n"
       << " [ -q ]            Run command in queit mode\n"
       << " [ -h ]            Disply this information\n"
       << endl;
  exit(1);
}
