/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-04-26 10:40:55
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-16 19:50:27
 */

#include "dist.h"
Info theInfo;

int main(int argc, char *argv[]) {
  // get the input arguments
  Args myargs(argc, argv);

  // init the distance matrix by species list
  Mdist dm;
  dm.init(myargs.glist);

  // assign the dm by reference DMs
  assignDM(myargs.refdm, myargs.netcdf, dm);
  theInfo.output(dm.info());

  if (dm.hasNAN()) {
    theInfo.output("Start distance calculate");

    do {
      // set the cvfile name
      vector<string> cvfile(myargs.glist);
      for (auto &str : cvfile)
        str += myargs.suffix + ".gz";

      // check the memory and divided steps
      IterStep theStep(dm, cvfile);
      theStep.checkSize(myargs.maxM, dm);
      theStep.execute(dm, myargs.method);
    } while (dm.hasNAN());
    theInfo.output("End the distance calculate");
  }

  // output the distance matrix
  dm.writemtx(myargs.outfile, myargs.netcdf);
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv)
    : outfile(""), netcdf(false), suffix(".faa.cv6") {

  program = argv[0];
  memorySize = getMemorySize() * 0.8;
  string listfile("list");
  string methStr("Hao");
  string cvdir("");

  char ch;
  while ((ch = getopt(argc, argv, "i:V:s:o:m:M:r:Cqh")) != -1) {
    switch (ch) {
    case 'i':
      listfile = optarg;
      break;
    case 'V':
      cvdir = optarg;
      addsuffix(cvdir, '/');
      break;
    case 'r':
      refdm = optarg;
      break;
    case 'M':
      str2number(optarg, memorySize);
      memorySize *= 1073741824;
      break;
    case 'm':
      methStr = optarg;
      break;
    case 's':
      suffix = optarg;
      if (*(suffix.begin()) != '.')
        suffix = '.' + suffix;
      break;
    case 'o':
      outfile = optarg;
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

  // read the genome list
  readlist(listfile, glist);
  uniqueWithOrder(glist);
  if (!cvdir.empty()) {
    for (auto &f : glist)
      f = cvdir + f;
  }

  // set the method
  method = Method::create(methStr);

  // set the outfile name
  if (outfile.empty()) {
    outfile = methStr + suffix;
  }

  //... Get The limit of memory size for cv
  float maxNameLen = 2048.0;
  float giga = 1073741824.0;
  size_t ng = glist.size();
  float bs = maxNameLen * ng + ng * (ng + 1) * sizeof(double) / 2 + giga;
  maxM = memorySize - bs;
}

void Args::usage() {
  cerr << "\nThe total physical memory of this computer is " << memorySize
       << " Byte\n"
       << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -o <dm> ]      Output distance matrix, defaut: "
          "<Method><Suffix>\n"
       << " [ -V <cvdir> ]   Super directory of extend cv files\n"
       << " [ -i list ]      Genome list for distance matrix, defaut: list\n"
       << " [ -s <Suffix> ]  Suffix of the cvfile, default: .faa.cv6\n"
       << " [ -r <matrix> ]  Reference distance matrixs, splite with ','\n"
       << " [ -M <N> ]       Runing memory size as G roughly, \n"
       << "                  default 80% of physical memory\n"
       << " [ -m Hao/InterList/InterSet] Method for cvtree, defaut: Hao\n"
       << " [ -C ]           Force use the netcdf compress distance matrix\n"
       << " [ -q ]           Run command in queit mode\n"
       << " [ -h ]           Disply this information\n"
       << endl;
  exit(1);
}