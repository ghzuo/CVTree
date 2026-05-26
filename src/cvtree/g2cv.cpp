/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-13 Saturday 20:46:19
 */

#include "g2cv.h"

int main(int argc, char *argv[]) {

  // set the argures
  Args myargs(argc, argv);

#pragma omp parallel for
  for (long i = 0; i < myargs.flist.size(); ++i) {
    myargs.meth->execute(myargs.flist[i], myargs.klist);
  }
};

Args::Args(int argc, char **argv) {

  string listfile("list");
  string listkval("5 6 7");
  string onefasta("");
  string gtype("faa");
  string gdir("");
  string cvdir("");
  string methStr("Hao");
  string cgstr("");

  // Define the available options and parameters
  argparse::ArgumentParser parser("g2cv", "1.0",
                                  argparse::default_arguments::help);

  parser.add_argument("-G", "--genome-dir")
      .help("Directory containing genome files")
      .default_value(gdir)
      .nargs(1)
      .action([&](const auto &val) { gdir = addsuffix(val, '/'); });

  parser.add_argument("-i", "--input")
      .help("Input list file")
      .default_value(listfile)
      .nargs(1)
      .store_into(listfile);

  parser.add_argument("-V", "--cv-dir")
      .help("Directory for CV files, default:  \nfor normal: <same to fasta "
            "file>\nfor bootstrap: ./resample/")
      .nargs(1)
      .action([&](const auto &val) { cvdir = addsuffix(val, '/'); });

  parser.add_argument("-g", "--genome-type")
      .help("Genome file type (faa/ffn/fna)")
      .choices("faa", "ffn", "fna")
      .default_value(gtype)
      .nargs(1)
      .store_into(gtype);

  parser.add_argument("-C", "--residue-group")
      .help("Residue groups string")
      .nargs(1)
      .store_into(cgstr);

  parser.add_argument("-k", "--k-values")
      .help("K values")
      .default_value(listkval)
      .nargs(1)
      .store_into(listkval);

  parser.add_argument("-m", "--method")
      .help("Method to build CV")
      .choices("Hao", "Count")
      .default_value(methStr)
      .nargs(1)
      .store_into(methStr);

  parser.add_argument("-f", "--file")
      .help("Single fasta file to CV")
      .nargs(1)
      .store_into(onefasta);

  parser.add_argument("-q", "--quiet")
      .help("Run in quiet mode")
      .nargs(0)
      .action([&](const auto &) { theInfo.quiet = true; });

  parser.add_description("Generate CV from genome files");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  // set the method
  meth = CVmeth::create(methStr, cvdir, gtype, cgstr);

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
    readlist(listfile, flist, 1);
    uniqueWithOrder(flist);
  } else {
    flist.emplace_back(onefasta);
  }

  // for (auto &gname : flist) {
  //   if (getsuffix(gname) == gtype)
  //     gname = delsuffix(gname);
  // }

  if (!gdir.empty()) {
    for (auto &gname : flist) {
      gname = gdir + gname;
    }
  }
}
