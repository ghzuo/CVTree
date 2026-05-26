/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-13 Saturday 21:01:47
 */

#include "distmatrix.h"
#include "stringOpt.h"
#include <argparse/argparse.hpp>
#include <array>
#include <map>
#include <regex>
using namespace std;

struct Args {
  vector<string> dmfile;
  string spfile;
  string taxfile;
  bool html;
  bool list;

  Args(int argc, char **argv) : dmfile{"infile"}, html(false), list(false) {

    argparse::ArgumentParser parser("getdist", "1.0",
                                    argparse::default_arguments::help);

    parser.add_argument("-d", "--distance")
        .help("Distance matrix files separated by :")
        .default_value(string("infile"))
        .nargs(1)
        .action([this](const auto &val) { separateWord(dmfile, val, ":"); });

    parser.add_argument("-i", "--input")
        .help("The select genomes")
        .nargs(1)
        .store_into(spfile);

    parser.add_argument("-T", "--taxonomy")
        .help("Replace the file name by taxon name")
        .nargs(1)
        .store_into(taxfile);

    parser.add_argument("-H", "--html")
        .help("Output html table file")
        .flag()
        .store_into(html);

    parser.add_argument("-l", "--list")
        .help("Get the list of species of distance matrix")
        .flag()
        .store_into(list);

    parser.add_description("Get distances from distance matrix");

    try {
      parser.parse_args(argc, argv);
    } catch (const std::exception &e) {
      std::cerr << e.what() << std::endl;
      std::cout << parser;
      exit(1);
    }
  }
};

void getIndex(const string &file, const vector<string> &nmlist,
              vector<size_t> &ndxlist) {

  Mdist dm;
  dm.readmtx(file);
  map<string, int> names;

  for (long i = 0; i < dm.size(); ++i) {
    string str = dm.getname(i);
    size_t pos = str.find_last_of('>');
    if (pos != string::npos) {
      string sstr = str.substr(pos + 1);
      names[sstr] = i;
    } else {
      names[str] = i;
    }
  }

  for (auto &n : nmlist) {
    if (names.find(n) == names.end()) {
      cerr << "Cannot find genome named " << n << " in the distance matrix"
           << endl;
      exit(2);
    } else {
      ndxlist.emplace_back(names[n]);
    }
  }
}

int main(int argc, char *argv[]) {

  Args args(argc, argv);

  // for output the list
  if (args.list) {
    Mdist dm;
    dm.readmtx(args.dmfile[0]);
    for (size_t i = 0; i < dm.size(); ++i) {
      string str = dm.getname(i);
      size_t pos = str.find_last_of('>');
      if (pos != string::npos)
        str = str.substr(pos + 1);
      cout << str << endl;
    }
    return 0;
  }

  // get the specias list
  vector<string> splist;
  if (!args.spfile.empty()) {
    readlist(args.spfile, splist);
  } else if (stdin) {
    string str;
    while (cin >> str)
      splist.emplace_back(str);
  }

  // read the taxonomy map
  map<string, string> taxmap;
  if (!args.taxfile.empty()) {
    ifstream tax(args.taxfile.c_str());
    if (!tax) {
      cerr << "\nCannot found the input file " << args.taxfile << endl;
      exit(1);
    }

    regex matchReg("([^; ]+)[; ]+([^; ]+)");
    for (string line; getline(tax, line);) {
      if (!line.empty()) {
        smatch matchs;
        regex_search(line, matchs, matchReg);
        if (!matchs.empty()) {
          string theName = matchs[1].str();
          string theTax = matchs[2].str();
          auto pos = theTax.find("<T>");
          if (pos != string::npos)
            taxmap[theName] = theTax.substr(pos + 3);
        } else {
          cerr << "Error parse in " << line << endl;
        }
      }
    }
    tax.close();
  }

  if (splist.size() > 1) {
    // get the index
    vector<size_t> ndxlist;
    getIndex(args.dmfile[0], splist, ndxlist);

    // get the distance
    vector<Mdist> dmlist;
    for (auto &fname : args.dmfile) {
      Mdist dm;
      dm.readmtx(fname);
      dm.reduce(ndxlist);
      dmlist.emplace_back(dm);
    }

    // replace the file name by taxonomy name
    for (auto &gname : splist) {
      if (taxmap.find(gname) != taxmap.end())
        gname = taxmap[gname];
    }

    if (args.html) {
      for (long i = 0; i < splist.size(); ++i) {
        for (long j = i + 1; j < splist.size(); ++j) {
          cout << "<tr>"
               << "<td>" << splist[i] << "</td>"
               << "<td>" << splist[j] << "</td>";
          for (auto &dm : dmlist)
            cout << "<td>" << dm.getdist(i, j) << "</td>";
          cout << "</tr>" << endl;
        }
      }
    } else {
      for (long i = 0; i < splist.size(); ++i) {
        for (long j = i + 1; j < splist.size(); ++j) {
          cout << splist[i] << "\t" << splist[j] << "\t";
          for (auto &dm : dmlist)
            cout << dm.getdist(i, j) << "\t";
          cout << endl;
        }
      }
    }
  }
}
