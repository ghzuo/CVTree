/*
 * Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-04-26 09:00:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-13 Saturday 22:20:45
 */

#include "treedist.h"

int main(int argc, char *argv[]) {

  // get the input arguments
  Args myargs(argc, argv);

  /// read the first tree
  MarkNode* aTree = MarkNode::initial();
  aTree->innwk(myargs.afile);
  vector<Node *> allLeafs;
  aTree->getLeafs(allLeafs);
  map<string, SetSym> mgi;
  size_t index(0);
  size_t nTotal = allLeafs.size();
  for (auto& leaf : allLeafs) {
    MarkNode* mleaf = dynamic_cast<MarkNode*>(leaf);
    mleaf->setOneContent(nTotal, index++);
    mgi[mleaf->name] = mleaf->content;
  }
  if(myargs.unroot){
    MarkNode* mleaf = dynamic_cast<MarkNode*>(allLeafs.front());
    mleaf->flipContent();
    mgi[mleaf->name] = mleaf->content;
  }
  aTree->setBranchContents();

  // read other tree
  MarkNode* bTree = MarkNode::initial();
  bTree->innwk(myargs.bfile);
  if(bTree->setAllContents(mgi)){
    set<SetSym> aset;
    aTree->getBranchContents(aset);
    set<SetSym> bset;
    bTree->getBranchContents(bset);
    
    for(auto& syb : aset)
      bset.erase(syb);
    long nTotal = aset.size();
    long nDiff  = bset.size();
    cout << nDiff * 2 << " " << float(nDiff)/nTotal << endl;
  }else{
    cerr << myargs.bfile << endl;
  }
}

Args::Args(int argc, char **argv)
    : afile("atree.nwk"), bfile("btree.nwk"), unroot(false) {

  argparse::ArgumentParser parser("treedist", "1.0", 
                                  argparse::default_arguments::help);
  
  parser.add_argument("-a", "--tree-a")
      .help("The tree file to compare")
      .default_value(afile)
      .nargs(1)
      .store_into(afile);
      
  parser.add_argument("-b", "--tree-b")
      .help("The other tree file to compare")
      .default_value(bfile)
      .nargs(1)
      .store_into(bfile);
      
  parser.add_argument("-u", "--unroot")
      .help("Whether this two tree are unroot")
      .flag()
      .store_into(unroot);
      
  parser.add_description("Calculate distance between two trees");
  
  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }
}


