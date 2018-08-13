/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-11-27 11:32:40
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2018-07-26 22:00:51
 */

#include <map>
#include <array>
#include "stringOpt.h"
#include "distmatrix.h"
using namespace std;

void usage(string &program)
{
    cerr << "\nProgram Usage: \n\n"
         << program
         << " [ -d infile ]       distance matrix files list separate by :, default: infile\n"
         << " [ -i selfile ]      the select genomes, defaut: none\n"
         << " [ -H ]              whether output html table, default: false\n"
         << " [ -l ]              get the list of species of distance matrix\n"
         << " [ -h ]              disply this information\n"
         << endl;
    exit(1);
}

void getIndex(const string &file, const vector<string> &nmlist, vector<size_t> &ndxlist)
{

    Mdist dm;
    dm.readmtx(file, true);
    map<string, int> names;

    for (int i = 0; i < dm.size(); ++i)
    {
        string str = dm.getname(i);
        size_t pos = str.find_last_of('>');
        if (pos != string::npos)
        {
            string sstr = str.substr(pos + 1);
            names[sstr] = i;
        }
        else
        {
            names[str] = i;
        }
    }

    for (auto &n : nmlist)
    {
        if (names.find(n) == names.end())
        {
            cerr << "Cannot find genome named " << n << " in the distance matrix" << endl;
            exit(2);
        }
        else
        {
            ndxlist.emplace_back(names[n]);
        }
    }
}

int main(int argc, char *argv[])
{

    // get the name of file
    vector<string> dmfile{"infile"};
    string spfile;
    string strlist;
    bool html(false);
    bool list(false);
    string program = argv[0];

    char ch;
    while ((ch = getopt(argc, argv, "d:i:Hlh")) != -1)
    {
        switch (ch)
        {
        case 'i':
            spfile = optarg;
            break;
        case 'd':
            separateWord(dmfile, optarg, ":");
            break;
        case 'H':
            html = true;
            break;
        case 'l':
            list = true;
            break;
        case 'h':
            usage(program);
        case '?':
            usage(program);
        }
    }

    // for output the list
    if (list)
    {
        Mdist dm;
        dm.readmtx(dmfile[0], true);
        for (size_t i = 0; i < dm.size(); ++i)
        {
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
    if (!spfile.empty())
    {
        readlist(spfile, splist);
    }
    else if (stdin)
    {
        string str;
        while (cin >> str)
            splist.emplace_back(str);
    }

    if (splist.size() > 1)
    {
        // get the index
        vector<size_t> ndxlist;
        getIndex(dmfile[0], splist, ndxlist);

        // get the distance
        vector<Mdist> dmlist;
        for (auto &fname : dmfile)
        {
            Mdist dm;
            dm.readmtx(fname, true);
            dm.reduce(ndxlist);
            dmlist.emplace_back(dm);
        }

        if (html)
        {
            for (int i = 0; i < splist.size(); ++i)
            {
                for (int j = i + 1; j < splist.size(); ++j)
                {
                    cout << "<tr>"
                         << "<td>" << splist[i] << "</td>"
                         << "<td>" << splist[j] << "</td>";
                    for (auto &dm : dmlist)
                        cout << "<td>" << dm.getdist(i, j) << "</td>";
                    cout << "</tr>" << endl;
                }
            }
        }
        else
        {
            for (int i = 0; i < splist.size(); ++i)
            {
                for (int j = i + 1; j < splist.size(); ++j)
                {
                    cout << splist[i] << "\t"
                         << splist[j] << "\t";
                    for (auto &dm : dmlist)
                        cout << dm.getdist(i, j) << "\t";
                    cout << endl;
                }
            }
        }
    }
}
