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

#include "cv2dm.h"
Info theInfo;

int main(int argc, char *argv[]) {
  // get the input arguments
  Args myargs(argc, argv);

  // init the distance matrix by species list
  Mdist dm;
  dm.init(myargs.glist);

  // assign the dm by reference DMs
  dm.assign(myargs.refdm);
  theInfo(dm.info());

  if (dm.hasNAN()) {
    theInfo("Start distance calculate");

    // set the cvfile name
    vector<string> cvfile(myargs.glist);
    for (auto &str : cvfile) {
      str += myargs.suffix + ".gz";
    }

    // do the calculation of distance
    myargs.meth->execute(cvfile, dm);

    theInfo("End the distance calculate");
  }

  // output the distance matrix
  dm.writemtx(myargs.outfile);
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv) : outfile(""), suffix(".faa.cv6") {

  program = argv[0];
  memorySize = getMemorySize() * 0.8;
  string listfile("list");
  string methStr("Cosine");
  string cvdir("");

  char ch;
  while ((ch = getopt(argc, argv, "i:V:s:o:m:M:r:qh")) != -1) {
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
    case 'q':
      theInfo.quiet = true;
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
  meth = DistMeth::create(methStr);

  // set the outfile name
  if (outfile.empty()) {
#ifdef _NETCDF
    outfile = methStr + suffix + ".nc";
#elif _HDF5
    outfile = methStr + suffix + ".h5";
#else
    outfile = methStr + suffix + ".txt";
#endif
  }

  //... Get The limit of memory size for cv
  meth->setMaxMem(memorySize, glist.size(), 1);
}

void Args::usage() {
  cerr << "\nThe total physical memory of this computer is " << memorySize
       << " Byte\n"
       << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -o <dm> ]      Output distance matrix, default: "
          "<Method><Suffix>.h5\n"
       << " [ -V <cvdir> ]   Super directory of extend cv files\n"
       << " [ -i list ]      Genome list for distance matrix, default: list\n"
       << " [ -s <Suffix> ]  Suffix of the cvfile, default: .faa.cv6\n"
       << " [ -r <matrix> ]  Reference distance matrices, split with ','\n"
       << " [ -M <N> ]       Running memory size as G roughly, \n"
       << "                  default 80% of physical memory\n"
       << " [ -m Cosine ]    Method for distance:\n"
       << "                  Cosine/Euclidean based on vector;\n"
       << "                  InterList/Min2Max based on count of kmers;\n"
       << "                  InterSet/Jaccard/Dice based set of kmers.\n"
       << "                  default: Cosine\n"
       << " [ -q ]           Run command in quiet mode\n"
       << " [ -h ]           Display this information\n"
       << endl;
  exit(1);
}