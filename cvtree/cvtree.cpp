/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-12 17:54:46
 */

#include "cvtree.h"

int main(int argc, char *argv[]) {
  // get the input arguments
  Args myargs(argc, argv);

  // get the main trees
  maintree(myargs);

  // do the bootstrap
  if (myargs.smeth != NULL)
    doSampleTest(myargs);
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv) : treeName(""), dmName(""), smeth(NULL) {

  program = argv[0];
  memorySize = getMemorySize() * 0.8;

  string pdir("cvtree/");
  string listfile("list");
  string methStr("Hao");
  string gtype("faa");
  string cgstr;
  string gdir("");
  string cvdir(pdir + "cv/");
  string listkval("5 6 7");
  bool refself(false);
  string sdir;
  long nSample(10);

  char ch;
  while ((ch = getopt(argc, argv, "i:G:V:P:k:d:t:m:M:r:g:C:S:s:j:bRqh")) !=
         -1) {
    switch (ch) {
    case 'i':
      listfile = optarg;
      break;
    case 'V':
      cvdir = optarg;
      if (!cvdir.empty())
        addsuffix(cvdir, '/');
      break;
    case 'G':
      gdir = optarg;
      addsuffix(gdir, '/');
      break;
    case 'P':
      pdir = optarg;
      addsuffix(pdir, '/');
      break;
    case 'r':
      refdm = optarg;
      break;
    case 'R':
      refself = true;
      break;
    case 'M':
      str2number(optarg, memorySize);
      memorySize *= 1073741824;
      break;
    case 'm':
      methStr = optarg;
      break;
    case 'k':
      listkval = optarg;
      break;
    case 'd':
      dmName = optarg;
      break;
    case 't':
      treeName = optarg;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'g':
      gtype = optarg;
      break;
    case 'C':
      cgstr = optarg;
      break;
    case 'S':
      sdir = optarg;
      break;
    case 's':
      nSample = str2int(optarg);
      break;
    case 'b':
      smeth = SampleMeth::create("Bootstrap");
      break;
    case 'j':
      smeth = SampleMeth::create("Jackknife", str2float(optarg));
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

  // init genome type
  Letter::init(gtype, cgstr);

  // set the method
  if (methStr == "Hao" || methStr == "CVTree") {
    cmeth = CVmeth::create("Hao", cvdir, gtype);
    dmeth = DistMeth::create("Cosine");
    tmeth = TreeMeth::create("NJ");
  } else if (methStr == "InterSet") {
    cmeth = CVmeth::create("Count", cvdir, gtype);
    dmeth = DistMeth::create("InterSet");
    tmeth = TreeMeth::create("NJ");
  } else if (methStr == "InterList") {
    cmeth = CVmeth::create("Count", cvdir, gtype);
    dmeth = DistMeth::create("InterList");
    tmeth = TreeMeth::create("NJ");
  } else {
    vector<string> mlist;
    separateWord(mlist, methStr, ":");
    if (mlist.size() > 1) {
      cmeth = CVmeth::create(mlist[0], cvdir, gtype);
      dmeth = DistMeth::create(mlist[1]);
      if (mlist.size() > 2) {
        tmeth = TreeMeth::create(mlist[2]);
      } else {
        tmeth = TreeMeth::create("NJ");
      }
    } else {
      cerr << "Unknow Method: " << methStr << endl;
      exit(3);
    }
  }

  // get the kvalue and check
  vector<string> wd;
  separateWord(wd, listkval);
  for (auto &str : wd)
    klist.emplace_back(stoul(str));

  sort(klist.begin(), klist.end());
  uniqueWithOrder(klist);
  cmeth->checkK(klist);

  // get the input file name and genome name
  map<string, string> nameMap;
  readNameMap(listfile, flist, nameMap);
  uniqueWithOrder(flist);
  if (flist.size() < 3) {
    theInfo(
        "The number of species should be more than three.\nThere are only " +
        to_string(flist.size()) + " unique items in the list");
    exit(3);
  }

  // get the glist by flist and nameMap
  for (auto &fname : flist) {
    auto iter = nameMap.find(fname);

    // update gfile name
    if (fname.find_last_not_of(gtype) == string::npos &&
        fname.find_last_not_of("fasta") == string::npos)
      fname += "." + gtype;

    // set the genome name
    if (iter != nameMap.end()) {
      glist.emplace_back(iter->second);
    } else {
      glist.emplace_back(fname);
    }
  }

  // add the super folder
  if (!gdir.empty()) {
    for (auto &fname : flist) {
      fname = gdir + fname;
    }
  }

  // set the output tree name format
  if (treeName.empty()) {
    treeName = pdir + "tree/" + methStr + cmeth->cvsuff + "$.nwk";
  }

  // set the output dm name format
  if (dmName.empty()) {
#ifdef _HDF5
    dmName = pdir + "dm/" + methStr + cmeth->cvsuff + "$.h5";
#else
    dmName = pdir + "dm/" + methStr + cmeth->cvsuff + "$.gz";
#endif
  }

  // refer self distance matrix
  if (refself) {
    if (refdm.empty()) {
      refdm = dmName;
    } else {
      refdm += "," + dmName;
    }
  }

  //... Get The limit of memory size
  dmeth->setMaxMem(memorySize, flist.size(), klist.size());

  //... for resampling
  if (smeth != NULL) {
    // get the super folder
    if (sdir.empty())
      sdir = smeth->wkdir();
    addsuffix(sdir, '/');

    // get the sampling fold list
    for (long i = 0; i < nSample; ++i) {
      string tdir = sdir + int2lenStr(i, 4) + "/";
      blist.emplace_back(tdir);
    }
  }
}

void Args::usage() {
  cerr << "\nThe total physical memory of this computer is " << memorySize
       << " Byte\n"
       << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -d <dm> ]      Output distance matrix name, default: "
          "<Method>.<gtype>.<Suffix><K>\n"
       << " [ -t <nwk> ]     Output newick file name, default: "
          "<Method>.<gtype>.<Suffix><K>.nwk\n"
       << " [ -G <gdir> ]    Super directory of Input genome file, default: "
          "<current directory> \n"
       << " [ -g faa ]       Type of genome file [faa/ffn/fna], "
          "default: faa\n"
       << " [ -C <None> ]    Grouped letters, separated by ',', default: None\n"
       << " [ -V <cvdir> ]   Super directory of cv files\n"
       << " [ -i list ]      Genome list for calculating, default: list\n"
       << " [ -k '5 6 7' ]   Values of k, default: K = 5 6 7\n"
       << " [ -r <matrix> ]  Reference distance matrices, split with ','\n"
       << " [ -R ]           Refer the output distance matrix\n"
       << " [ -M <N> ]       Running memory size as G roughly,\n"
       << "                  default 80% of physical memory\n"
       << " [ -m Hao ]       Select CVTree Method from Hao/InterList/InterSet, "
          "default: Hao\n"
       << " [ -S resample/ ]  cache folder for resample, default: resample/\n"
       << " [ -s <n> ]        resample times, default: 10\n"
       << " [ -b ]            do bootstrap resampling\n"
       << " [ -j 0.8 ]        do jackknife resampling\n"
       << " [ -q ]           Run command in quiet mode\n"
       << " [ -h ]           Display this information\n"
       << endl;
  exit(1);
}

/********************************************************************************
 * @brief functions for obtaining main tree
 *
 * @param myargs
 ********************************************************************************/

void maintree(const Args &myargs) {
  // init distance matrixes
  vector<pair<size_t, Mdist>> dms(myargs.klist.size());
  initDMs(myargs.klist, myargs.glist, myargs.refdm, dms);

  // get the cvs for the NAN distances
  getCVs(myargs.cmeth, myargs.flist, dms);

  for (auto &kdm : dms) {
    // do the dm and tree
    getDM(myargs.dmeth, myargs.cmeth, myargs.flist, myargs.dmName, kdm);
    getTree(myargs.tmeth, myargs.treeName, kdm);
  }
}

/********************************************************************************
 * @brief Functions for resample
 *
 * @param myargs
 ********************************************************************************/

void doSampleTest(const Args &myargs) {
  theInfo("\n============ Start " + myargs.smeth->name + " Section ==========");
  // get all resample fasta files
  for (auto &gfile : myargs.flist) {
    sampleGenome(myargs.smeth, gfile, myargs.blist);
  }
  theInfo("Resample All genomes");

  // do resample jobs
  for (auto &sdir : myargs.blist) {
    theInfo("Start resampling job in " + sdir);
    // renew the parameters for the resample job
    vector<string> flist;
    for (auto nm : myargs.flist)
      flist.emplace_back(sdir + "fasta/" + getFileName(nm));
    myargs.cmeth->setCVdir(sdir + "cv/");
    string sdmName = sdir + "dm/" + getFileName(myargs.dmName);
    string streeName = sdir + "tree/" + getFileName(myargs.treeName);
    theInfo("Prepare all paremeters", 1);

    // init distance matrixes
    vector<pair<size_t, Mdist>> dms(myargs.klist.size());
    initDMs(myargs.klist, myargs.glist, sdmName, dms);

    // get the cvs for the NAN distances
    getCVs(myargs.cmeth, flist, dms);

    for (auto &kdm : dms) {
      // do the dm and tree
      getDM(myargs.dmeth, myargs.cmeth, flist, sdmName, kdm);
      getTree(myargs.tmeth, streeName, kdm);
    }
    theInfo("***", -1);
  }

  // boot the tree
  for(auto& k : myargs.klist){
    string maintree = nameWithK(myargs.treeName, k);
    bootTree(maintree, myargs.smeth->name, myargs.blist);
  }
  theInfo("Add bootstrap value on the main tree");

  theInfo("============ End Resampleing Section ==========");
}

/********************************************************************************
 * @brief the step functions
 *
 * @param myargs
 ********************************************************************************/
void initDMs(const vector<size_t> &klist, const vector<string> &glist,
             const string &refdm, vector<pair<size_t, Mdist>> &dms) {
#pragma omp parallel for ordered
  for (auto i = 0; i < klist.size(); ++i) {
    dms[i].first = klist[i];
    dms[i].second.init(glist);
    if (!refdm.empty()) {
      string str = nameWithK(refdm, dms[i].first);
      dms[i].second.assign(str);
    }
#pragma omp ordered
    theInfo(dms[i].second.info() + " for K=" + to_string(dms[i].first));
  }
}

void getCVs(CVmeth *cmeth, const vector<string> &flist,
            const vector<pair<size_t, Mdist>> &dms) {
  theInfo("Start check and calculate CVs for " + to_string(flist.size()) +
          " Genomes");
#pragma omp parallel for
  for (size_t i = 0; i < flist.size(); ++i) {
    // get the missing K list
    vector<size_t> klist;
    for (auto &kdm : dms) {
      if (kdm.second.hasNAN(i))
        klist.emplace_back(kdm.first);
    }

    // get the CVs
    cmeth->execute(flist[i], klist);
  }
  theInfo("CV Section: All CVs are obtained");
}

void getDM(DistMeth *dmeth, CVmeth *cmeth, const vector<string> &flist,
           const string &dnameTemple, pair<size_t, Mdist> &kdm) {
  // Calculate the NAN distance
  if (kdm.second.hasNAN()) {
    theInfo("Start calculate distance for K=" + to_string(kdm.first));
    // set the cvfile name
    vector<string> cvfile(flist);
    for (auto &str : cvfile) {
      str = cmeth->getCVname(str, kdm.first);
    }

    // do the calculation of distance
    dmeth->execute(cvfile, kdm.second);
    theInfo("End the calculate distance for K=" + to_string(kdm.first));
  }

  string dname = nameWithK(dnameTemple, kdm.first);
  if (!dname.empty())
    kdm.second.writemtx(dname);
}

void getTree(TreeMeth *tmeth, const string &tnameTemple,
             pair<size_t, Mdist> &kdm) {
  // do the NJ algorithm and return the NJ tree
  // initial the tree
  Node *aTree = Node::initial();
  theInfo("Start infer tree for K=" + to_string(kdm.first));
  tmeth->tree(aTree, kdm.second);
  theInfo("Get phylogenetic tree for K=" + to_string(kdm.first));

  // output the Tree
  string tname = nameWithK(tnameTemple, kdm.first);
  (*aTree).outnwk(tname);
  aTree->clear();
}

void bootTree(const string &maintree, const string &label,
              const vector<string> &blist) {
  /// read the first tree
  MarkNode mTree;
  mTree.innwk(maintree);

  /// boot the trees
  vector<string> btrflist;
  for (auto &bdir : blist)
    btrflist.emplace_back(bdir + "tree/" + getFileName(maintree));
  mTree.bootTree(btrflist);

  // output result
  mTree.outnwk(addnamelabel(maintree, label));
  mTree.clear();
};
