/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-13 Saturday 21:53:13
 */

#include "cv2dm.h"

int main(int argc, char *argv[]) {
  // get the input arguments
  Args myargs(argc, argv);
  theInfo("Get all options");

  // init the distance matrix by species list
  Mdist dm;
  dm.init(myargs.glist);
  theInfo("Initial the distance matrix");

  // assign the dm by reference DMs
  dm.assign(myargs.refdm);
  theInfo(dm.info());

  if (dm.hasNAN()) {
    theInfo("Start distance calculate");

    // set the cvfile name
    vector<string> cvfile(myargs.flist);
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

  memorySize = getMemorySize() * 0.8;
  string listfile("list");
  string methStr("Cosine");
  string cvdir("");
#ifdef _HDF5
  string dmftype = ".h5";
#else
  string dmftype = ".gz";
#endif

  // Define the available options and parameters
  argparse::ArgumentParser parser("cv2dm", "1.0",
                                  argparse::default_arguments::help);

  parser.add_argument("-i", "--input")
      .help("Input list file")
      .default_value(listfile)
      .nargs(1)
      .store_into(listfile);

  parser.add_argument("-V", "--cv-dir")
      .help("CV directory")
      .nargs(1)
      .action([&](const auto &val) { cvdir = addsuffix(val, '/'); });

  parser.add_argument("-r", "--reference-dm")
      .help("Reference distance matrix")
      .nargs(1)
      .store_into(refdm);

  parser.add_argument("-M", "--max-memory")
      .help("Memory size in GB, default is 0.8 of total memory")
      .nargs(1)
      .action([&](const auto &val) { memorySize = 1073741824 * stod(val); });

  parser.add_argument("-m", "--method")
      .help("Distance calculation method")
      .choices("Cosine", "Euclidean", "InterList", "Min2Max", "InterSet",
               "Jaccard", "Dice")
      .default_value(methStr)
      .nargs(1)
      .store_into(methStr);

  parser.add_argument("-s", "--suffix")
      .help("File suffix")
      .default_value(suffix)
      .nargs(1)
      .action([&](const auto &val) {
        suffix = val;
        if (*(suffix.begin()) != '.')
          suffix = '.' + suffix;
      });

  parser.add_argument("-o", "--output")
      .help("Output distance matrix file")
      .default_value(string("<Method>.<Suffix>") + dmftype)
      .nargs(1)
      .store_into(outfile);

  parser.add_argument("-q", "--quiet")
      .help("Run in quiet mode")
      .nargs(0)
      .action([](const auto &) { theInfo.quiet = true; });

  parser.add_description("Convert CV to distance matrix");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  // read the genome list
  map<string, string> nameMap;
  readNameMap(listfile, flist, nameMap);
  uniqueWithOrder(flist);

  // get the glist by flist and nameMap
  string suffstr = suffix.substr(1);
  for (auto &fname : flist) {
    auto iter = nameMap.find(fname);

    // delete the suffix of file
    if (getsuffix(fname) == suffstr)
      fname = delsuffix(fname);

    // set the genome name
    if (iter != nameMap.end()) {
      glist.emplace_back(iter->second);
    } else {
      glist.emplace_back(fname);
    }
  }

  // add the super folder
  if (!cvdir.empty()) {
    for (auto &f : flist)
      f = cvdir + f;
  }

  // set the method
  meth = DistMeth::create(methStr);

  // set the outfile name
  if (outfile.empty()) 
    outfile = methStr + suffix + dmftype;

  //... Get The limit of memory size for cv
  meth->setMaxMem(memorySize, flist.size(), 1);
}
