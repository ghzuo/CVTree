/*
 * Copyright (c) 2024
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-05-11 14:51:30
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-12 14:55:37
 */

#include "sampling.h"

void sampleGenome(SampleMeth *smeth, const string &gfile,
                  const vector<string> &slist) {

  // genarate the resample fasta
  Genome gn;
  readFasta(gfile, gn, false);

  for (auto &sdir : slist) {
    string sfile = sdir + "fasta/" + getFileName(gfile);
    if (!fileExists(sfile)) {
      vector<long> ndxlist = (*smeth)(gn.size());
      Genome gs;
      for (auto ndx : ndxlist) {
        gs.emplace_back(gn[ndx]);
      }
      writeFasta(sfile, gs);
    }
  }
};
