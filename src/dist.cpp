#include "dist.h"

int main(int argc, char* argv[]){
    // set a timer
    boost::timer mytimer;

    // get the input arguments
    Args myargs(argc,argv);

    // init the distance matrix and species list
    Mdist dm;
    
    // read the old matrix if exist
    if(!myargs.mtxfile.empty()){
        dm.readmtx(myargs.mtxfile, myargs.netcdf);

	if(!myargs.taxmap.empty()){
            for(int i=0; i<dm.size(); ++i){
                string str = dm.getname(i);
		size_t pos=str.find_last_of('>');
                if(pos!=string::npos){
		    string sstr = str.substr(pos+1);
		    if(myargs.taxmap.find(sstr) == myargs.taxmap.end())
			myargs.taxmap[sstr] = str;
                    dm.setname(i, sstr);
		}
            }
	}

	// adjust the distance matrix by orglist or orgIndex
	if(!myargs.orgIndex.empty())
            dm.reduce(myargs.orgIndex);
	else if(!myargs.orglist.empty())
            dm.reduce(myargs.orglist);
    }

    cerr << "*** Complete Distance Matrix of Exist Genomes, Time Elapsed: " 
	 << mytimer.elapsed() << "s" << endl;

    // extend the distance matrix
    if(!myargs.extlist.empty()){
	// resize the distance matrix and renew genome name list
	int ngorg = dm.size();
	int ngext = myargs.extlist.size();
        dm.extend(myargs.extlist);
	
	vector<IterStep> steplist;
	checkCVsize(myargs, ngorg, steplist);

	int index(0);
	for(const auto &iter : steplist){

	    // output the runing blocks
	    cerr << "*** Block " << ++index << ": from " 
		 << iter.first  - ngorg + 1 << " to " 
		 << iter.second - ngorg     << ", Total cv size of this block is " 
		 << iter.size/1073741824    << "G" << endl;

	    // read the extend cv
	    vector<CVvec> cvins(iter.second);
	    for(int i=iter.first; i<iter.second; ++i)
		for(int j=0; j<myargs.suflist.size(); ++j)
		    readcv(myargs.extdir+dm.getname(i)+myargs.suflist[j], cvins[i-iter.first]);

	    // for the distances between the extend genomes
	    vector<pair<int,int> > dolist;
	    for(int i=0; i<iter.second-iter.first; ++i){
		for(int j=0; j<i; ++j){
		    dolist.emplace_back(i,j);
		}
	    }

            #pragma omp parallel for schedule (dynamic)
	    for(int i=0; i<dolist.size(); ++i)
		dm.setdist(dolist[i].first+iter.first, 
		           dolist[i].second+iter.first, 
		           dist(cvins[dolist[i].first],cvins[dolist[i].second]));

	    // for the distance betwee genomes of next steps
            #pragma omp parallel for schedule (dynamic)
	    for(int i=iter.second; i<dm.size(); ++i){
		//read an orginal cv
		CVvec cv;
		for(int j=0; j<myargs.suflist.size(); ++j)
		    readcv(myargs.extdir+dm.getname(i)+myargs.suflist[j], cv);
	    
		// obtain the distances between the genome and the extend genomes
		for(int j=iter.first; j<iter.second; ++j)
		    dm.setdist(i, j, dist(cvins[j-iter.first], cv));
	    }
	
	    // for the distance between orginal genomes and extend genomes
            #pragma omp parallel for schedule (dynamic)
	    for(int i=0; i<ngorg; ++i){
		//read an orginal cv
		CVvec cv;
		for(int j=0; j<myargs.suflist.size(); ++j)
		    readcv(myargs.orgdir+dm.getname(i)+myargs.suflist[j], cv);
	    
		// obtain the distances between the genome and the extend genomes
		for(int j=iter.first; j<iter.second; ++j)
		    dm.setdist(i, j, dist(cvins[j-iter.first], cv));
	    }

	    cerr << "*** Complete Block " << index << ", Time Elapsed: " 
		 << mytimer.elapsed() << "s" << endl;
	}
    }

    // output the matrix
    if(myargs.outtax && !myargs.taxmap.empty()){
        for(int i=0; i<dm.size(); ++i){
            string orgname = dm.getname(i);
            dm.setname(i, myargs.taxmap[orgname]);
        }
    }
    
    dm.writemtxnc(myargs.outfile, myargs.netcdf);

    cerr << "*** Complete Program, Time Elapsed: " << mytimer.elapsed() << "s" << endl;
}

Args::Args(int argc, char** argv):outfile("dist.matrix"),orgdir("cv/"),extdir("cv/"),
				  outtax(true),netcdf(false){
    
    program = argv[0];
    memorySize = getMemorySize() * 0.8;
    string extfile, orgfile, ndxfile, taxfile;

    char ch;
    while ((ch = getopt(argc, argv, "i:I:n:e:E:s:o:m:t:M:TCh")) != -1){
      switch (ch){
      case 'i': 
	  extfile = optarg; break;
      case 'I':
	  extdir = optarg; 
	  addsuffix(extdir, '/'); 
	  break;
      case 'n': 
	  ndxfile = optarg; break;
      case 'e': 
	  orgfile = optarg; break;
      case 'E':
	  orgdir = optarg; 
	  addsuffix(orgdir, '/');
	  break;
      case 'm':
	  mtxfile = optarg; break;
      case 't':
	  taxfile = optarg; break;
      case 'T':
	  outtax = false; break;
      case 'M':
	  str2number(optarg, memorySize);
	  memorySize *= 1073741824;
	  break;
      case 's':
	  separateWord(suflist, optarg);
	  for(auto &str : suflist)
	      if(*(str.begin()) != '.')
		  str = '.' + str;
	  break;
      case 'o':
	  outfile = optarg; break;
      case 'C':
	  netcdf = true; break;
      case 'h':
	  usage();
      case '?':
	  usage();
      }
    }

    if(suflist.empty())
	suflist.emplace_back(".cv6.gz");

    // read the extend list
    if(!extfile.empty())
	readlist(extfile, extlist);

    // read the index or name list of origin genomes
    if(!ndxfile.empty())
    	readlist(ndxfile, orgIndex);
    else if(!orgfile.empty())
    	readlist(orgfile, orglist);

    // read the taxonomy map
    if(!taxfile.empty()){
	ifstream tax(taxfile.c_str());
	if(!tax){
	    cerr << "\nCannot found the input file " << taxfile << endl;
	    exit(1);
	}

        regex matchReg("([^; ]+)[; ]+([^; ]+)");
	for(string line; getline(tax, line); ){
	    if(!line.empty()){
                smatch matchs;
                regex_search(line, matchs, matchReg);
		if(!matchs.empty()){
		    taxmap[matchs[1].str()] = matchs[2].str() + "<T>" + matchs[1].str();
		}else{
		    cerr << "Error parse in " << line << endl;
		}
	    }
	}
	tax.close();

	/// check the list
	if(!extlist.empty()){
	    vector<string> vstr;
	    for(const auto &str : extlist)
		if(taxmap.find(str) == taxmap.end())
		    vstr.emplace_back(str);

	    if(!vstr.empty()){
		cerr << "Cannot find the taxonomy of: \n";
		for(const auto &str : vstr)
		    cerr << str << "\n";
		cerr << endl;
		exit(9);		
	    }
	}
    }
}


void Args::usage(){
    cerr << "\nThe total physical memory of this computer is " << memorySize << " Byte\n"
	 << "\nProgram Usage (VERSION: " << GHZ_VERSION << "): \n\n" 
	 << program  <<"\n" 
	 <<" [ -o dist.matrix ]   Output distance matrix, defaut: dist.matrix\n"
	 <<" [ -I extdir ]        Directory of extend cv files, defaut: cv\n"
	 <<" [ -i infile ]        Extend cv file list, defaut: no extend cv used\n"
	 <<" [ -s cv6.gz ]        Suffix of cv file, defaut: cv6.gz\n"
	 <<" [ -E orgdir  ]       Directory of the orginal cv files\n"
	 <<" [ -m indist.matrix ] Input distance matrix, default: no input matrix used\n"
	 <<" [ -e orglist ]       Name of selected genomes, which are in input distance matrix\n"
	 <<" [ -n ndxlist ]       Index of selected genomes, which are input distance matrix\n"
	 <<"                      The index file used first!\n"
	 <<" [ -t taxfile ]       input taxonomy information\n"
	 <<" [ -T ]               Do not output taxonomy information\n"
	 <<" [ -M <N> ]           Runing memory size as G roughly, default 80% physical memory\n"
	 <<" [ -C ]               Force use the netcdf compress distance matrix\n"
	 <<" [ -h ]               disply this information\n"
	 << endl;
    exit(1);
}

void checkCVsize(const Args& myargs, size_t ibeg, vector<IterStep>& slist){

    //... GettThe memory size for cv
    float maxNameLen = 2048.0;
    float giga = 1073741824.0;
    size_t ng = ibeg + myargs.extlist.size();
    float bs = maxNameLen*ng + maxNameLen*myargs.taxmap.size() + ng*(ng+1)*sizeof(double) + giga;
    float maxM = myargs.memorySize - float(bs);
    float totalsize(0.0);
    
    IterStep astep(ibeg);
    for(int i=0; i<myargs.extlist.size(); ++i){
	// get size size of cv file
	for(int j=0; j<myargs.suflist.size(); ++j){
	    float oneCVSize = cvsize(myargs.extdir+myargs.extlist[i]+myargs.suflist[j]);
	    astep.size += oneCVSize;
	    totalsize  += oneCVSize;
	}

	// check size
	if(astep.size > maxM){
	    astep.second = ibeg + i + 1;
	    slist.emplace_back(astep);
	    astep.first = astep.second;
	    astep.size = 0;
	}
    }

    if(astep.size != 0){
        astep.second = ibeg + myargs.extlist.size();
	slist.emplace_back(astep);
    }

    // output the size information
    cerr << "*** The memory of the program limited to: " 
	 << myargs.memorySize / giga 
	 << "G\n*** The total size of all "
	 << myargs.extlist.size() << " CVs is: " 
	 << totalsize / giga
	 << "G\n";

    if(slist.size() > 1)
	cerr << "*** The  new CVs will be divided into "
	     << slist.size() << " blocks for memory limit" << endl;
    else
	cerr << "*** All new CVs will complete in a block" << endl;
}
