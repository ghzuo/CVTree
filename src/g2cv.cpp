/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-04-19 11:37:42
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-19 23:18:09
 */

#include "g2cv.h"
Info theInfo;

int main(int argc, char *argv[]) {

  // set the argures
  Args myargs(argc, argv);

#pragma omp parallel for
  for (int i = 0; i < myargs.flist.size(); ++i) {
    myargs.meth->execute(myargs.flist[i], myargs.klist);
  }
};

Args::Args(int argc, char **argv) {

  program = argv[0];
  string listfile("list");
  string listkval("5 6 7");
  string onefasta("");
  string gtype("faa");
  string gdir("");
  string cvdir("");
  string methStr("Hao");

  char ch;
  while ((ch = getopt(argc, argv, "G:i:k:V:g:m:f:qh")) != -1) {
    switch (ch) {
    case 'G':
      gdir = optarg;
      addsuffix(gdir, '/');
      break;
    case 'i':
      listfile = optarg;
      break;
    case 'V':
      cvdir = optarg;
      addsuffix(cvdir, '/');
      break;
    case 'g':
      gtype = optarg;
      break;
    case 'k':
      listkval = optarg;
      break;
    case 'm':
      methStr = optarg;
      break;
    case 'f':
      onefasta = optarg;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'h':
      usage();
    case '?':
      usage();
    }
  }

  // check the genome type
  if (gtype != "faa" && gtype != "ffn" && gtype != "fna") {
    cerr << "Only faa/ffn/fna are supported!\n" << endl;
    exit(1);
  }

  // set the method
  if (methStr == "Hao") {
    meth = new HaoMethod;
  } else if (methStr == "Count") {
    meth = new Counting;
  } else {
    cerr << "Unknow Method: " << methStr << endl;
    exit(3);
  }
  meth->init(cvdir, gtype);

  // get the kvalue
  vector<string> wd;
  separateWord(wd, listkval);
  for (auto &str : wd)
    klist.emplace_back(stoul(str));
  sort(klist.begin(), klist.end());
  uniqueWithOrder(klist);
  meth->checkK(klist);

  // get the input file name
  if (onefasta.empty()) {
    readlist(listfile, flist);
    uniqueWithOrder(flist);
  } else {
    flist.emplace_back(onefasta);
  }

  for (auto &gname : flist) {
    if (getsuffix(gname) == gtype)
      gname = delsuffix(gname);
  }

  if (!gdir.empty()) {
    for (auto &gname : flist) {
      gname = gdir + gname;
    }
  }
};

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -G <dir> ]      input genome file directory\n"
       << " [ -V <dir> ]      output cv directory\n"
       << " [ -i list ]       input species list, defaut: list\n"
       << " [ -f <Fasta> ]    get cv for only one fasta \n"
       << " [ -k '5 6 7' ]    values of k, defaut: K = 5 6 7\n"
       << " [ -g faa ]        the type of genome file, defaut: faa\n"
       << " [ -m Hao/Count ]  the method for cvtree, defaut: Hao\n"
       << " [ -q ]            Run command in queit mode\n"
       << " [ -h ]            disply this information\n"
       << endl;

  exit(1);
}
