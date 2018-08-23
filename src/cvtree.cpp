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

#include "cvtree.h"

// for output infomation
Info theInfo;

int main(int argc, char *argv[]) {
  // get the input arguments
  Args myargs(argc, argv);

  // init domains
  vector<pair<size_t, Mdist>> dms(myargs.klist.size());
#pragma omp parallel for ordered
  for (size_t i = 0; i < myargs.klist.size(); ++i) {
    dms[i].first = myargs.klist[i];
    dms[i].second.init(myargs.glist);
    if (!myargs.refdm.empty()) {
      string str = nameWithK(myargs.refdm, dms[i].first);
      dms[i].second.assign(str, myargs.netcdf);
    }
#pragma omp ordered
    theInfo(dms[i].second.info() + " for K=" + to_string(dms[i].first));
  }

  // get the cvs for the NAN distances
  theInfo("Start check and calculate CVs for " +
                 to_string(myargs.glist.size()) + " Genomes");
#pragma omp parallel for
  for (size_t i = 0; i < myargs.glist.size(); ++i) {
    // get the missing K list
    vector<size_t> klist;
    for (size_t j = 0; j < myargs.klist.size(); ++j) {
      if (dms[j].second.hasNAN(i))
        klist.emplace_back(myargs.klist[j]);
    }

    // get the CVs
    myargs.cmeth->execute(myargs.glist[i], klist);
  }
  theInfo("All CVs are obtained");

  // Calculate the NAN distance
  for (size_t j = 0; j < myargs.klist.size(); ++j) {
    Mdist &dm = dms[j].second;
    size_t k = dms[j].first;

    if (dm.hasNAN()) {

      theInfo("Start calculate for K=" + to_string(k));

      // set the cvfile name
      vector<string> cvfile(myargs.glist);
      for (auto &str : cvfile) {
        str = myargs.cmeth->getCVname(str, k);
      }

      // do the calculation of distance
      myargs.dmeth->execute(cvfile, dm);

      theInfo("End the calculate distance for K=" + to_string(k));
    }

    // output the distance matrix
    string fname = nameWithK(myargs.dmName, k);
    mkpath(fname);
    dm.writemtx(fname, myargs.netcdf);
  }

// check the genome number bigger than 3
  if(myargs.glist.size() < 3){
      theInfo("There are only " + to_string(myargs.glist.size()) + " genomes, no tree will output");
      exit(2);
  }

// get the nwk tree
#pragma omp parallel for ordered
  for (size_t i = 0; i < dms.size(); ++i) {
    // do the NJ algorithm and return the NJ tree
    Node *aTree = neighborJoint(dms[i].second);

#pragma omp ordered
    theInfo("Get the Neighbor Joint tree for K=" +
                   to_string(dms[i].first));

    // output the Tree
    // output the distance matrix
    string fname = nameWithK(myargs.treeName, dms[i].first);
    mkpath(fname);
    ofstream nwk(fname.c_str());
    (*aTree).outnwk(nwk);
    nwk.close();
  }
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv) : treeName(""), dmName(""), netcdf(false) {

  program = argv[0];
  memorySize = getMemorySize() * 0.8;

  string listfile("list");
  string methStr("Hao");
  string gtype("faa");
  string gdir("");
  string cvdir("cv/");
  string listkval("5 6 7");

  char ch;
  while ((ch = getopt(argc, argv, "i:G:V:k:d:t:m:M:r:g:Cqh")) != -1) {
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
    case 'k':
      listkval = optarg;
      break;
    case 'd':
      dmName = optarg;
      break;
    case 't':
      treeName = optarg;
      break;
    case 'C':
      netcdf = true;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'g':
      gtype = optarg;
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
  if (methStr == "Hao" || methStr == "CVTree") {
    cmeth = new HaoMethod;
    dmeth = new Cosine;
  } else if (methStr == "Li" || methStr == "InterSet") {
    cmeth = new Counting;
    dmeth = new InterSet;
  } else if (methStr == "Zuo" || methStr == "InterList") {
    cmeth = new Counting;
    dmeth = new InterList;
  } else {
    cerr << "Unknow Method: " << methStr << endl;
    exit(3);
  }
  cmeth->init(cvdir, gtype);

  // get the kvalue and check
  vector<string> wd;
  separateWord(wd, listkval);
  for (auto &str : wd)
    klist.emplace_back(stoul(str));

  sort(klist.begin(), klist.end());
  uniqueWithOrder(klist);
  cmeth->checkK(klist);

  // get the input file name
  readlist(listfile, glist);
  uniqueWithOrder(glist);

  for (auto &gname : glist) {
    if (getsuffix(gname) == gtype)
      gname = delsuffix(gname);
  }

  if (!gdir.empty()) {
    for (auto &gname : glist) {
      gname = gdir + gname;
    }
  }

  // set the output tree name format
  if (treeName.empty()) {
    treeName = "tree/" + methStr + cmeth->cvsuff + "$.nwk";
  }

  // set the output dm name format
  if (dmName.empty()) {
    dmName = "dm/" + methStr + cmeth->cvsuff + "$.nc";
  }

  //... Get The limit of memory size
  dmeth->setMaxMem(memorySize, glist.size(), klist.size());

}

void Args::usage() {
  cerr << "\nThe total physical memory of this computer is " << memorySize
       << " Byte\n"
       << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -d <dm> ]         Output distance matrix format, defaut: "
          "<Method><Suffix><K>\n"
       << " [ -t <nwk> ]         Output distance matrix format, defaut: "
          "<Method><Suffix><K>.nwk\n"
       << " [ -G <dir> ]        Input genome file directory\n"
       << " [ -g faa ]          the type of genome file, defaut: faa\n"
       << " [ -V <cvdir> ]      Super directory of cv files\n"
       << " [ -i list ]         Genome list for distance matrix, defaut: "
          "list\n"
       << " [ -k '5 6 7' ]      values of k, defaut: K = 5 6 7\n"
       << " [ -r <matrix> ]     Reference distance matrixs, splite with ','\n"
       << " [ -M <N> ]          Runing memory size as G roughly,\n"
       << "                     default 80% of physical memory\n"
       << " [ -m Hao/InterList/InterSet] Method for cvtree, defaut: Hao\n"
       << " [ -C ]           Force use the netcdf compress distance matrix\n"
       << " [ -q ]           Run command in queit mode\n"
       << " [ -h ]           Disply this information\n"
       << endl;
  exit(1);
}

string nameWithK(const string &str, size_t k) {

  string kstr = to_string(k);
  if (str.empty())
    return kstr;

  if (str.find("$", 0) == std::string::npos)
    return str + "$";

  string sstr = str;
  size_t npos = 0;
  while ((npos = sstr.find("$", npos)) != std::string::npos) {
    sstr.replace(npos, 1, kstr);
    npos += kstr.length();
  }
  return sstr;
}

void mkpath(const string &nm) {
  size_t npos = 0;
  while ((npos = nm.find("/", npos)) != std::string::npos) {
    string dir = nm.substr(0, npos);
    mkdir(dir.c_str(), 0755);
    npos++;
  }
}
