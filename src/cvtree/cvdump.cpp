/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-13 Saturday 21:08:39
 */

#include "cv.h"
#include "stringOpt.h"
#include <argparse/argparse.hpp>

int main(int argc, char *argv[]) {

  string gtype;
  string infile;
  bool kstr = true;

  // Define the available options and parameters
  argparse::ArgumentParser parser("cvdump", "1.0",
                                  argparse::default_arguments::help);
  
  parser.add_argument("-i", "--input")
      .help("Input CV file")
      .required()
      .nargs(1)
      .store_into(infile);
      
  parser.add_argument("-g", "--genome-type")
      .help("Type of genome file, only for old format CV file")
      .default_value(gtype)
      .nargs(1)
      .store_into(gtype);
      
  parser.add_argument("-n", "--numeric")
      .help("Output the number code instead of letters")
      .nargs(0)
      .action([&](const auto &) { kstr = false; });
      
  parser.add_description("Dump CV file contents");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  // read cv
  CVvec cv;
  string cgstr = cv.read(infile);

  // initial the char set
  if(!cgstr.empty()){
    Letter::initByStr(cgstr);
  }else if(!gtype.empty()){
    Letter::init(gtype);
  }else{
    cerr << "Please input the type of genome: faa/fna/ffn" << endl;
    exit(2);
  }


  // output data
  cout << "The inner of CV: " << cv.norm << endl;
  cout << "The size  of CV: " << cv.size() << endl;
  cout << "The char set of CV: " << Letter::decMap << endl;

  for (const auto &cd : cv.vec)
    cout << cd.first << "\t" << cd.second << endl;
}
