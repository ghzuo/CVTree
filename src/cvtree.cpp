#include "cvtree.h"

int main(int argc, char* argv[]){
    // set a timer
    Timer mytimer;

    Args myargs(argc, argv);

    //get gene type to read and read gene file
    GeneType mygene(myargs.gtype);

    //get the genome letters map to check the sequcne
    int kmax = Kstr::init(mygene.letters);
    myargs.kcheck(kmax);

    #pragma omp parallel for schedule (dynamic)
    for(int i=0; i<myargs.flist.size(); ++i){
        string fname = myargs.flist[i].first;
	string infile = myargs.indir + fname;
	Genome genome;
	mygene.readgene(infile, genome);

	//get the cv of all K of the genome
	map<int, CVmap> mvc;
	map<int, double> nstr;
	for(auto &k : myargs.slist){
	    CVmap cv;
	    nstr[k] = count(genome, k, cv);
	    mvc[k] = cv;
	}

	for(auto &k : myargs.klist){
	    CVmap cv;
	    if(myargs.method == 1){
                double factor = nstr[k]*nstr[k-2]/(nstr[k-1]*nstr[k-1]);
	    	subtract(mvc[k], mvc[k-1], mvc[k-2], factor, cv);
	    }else
		cv = mvc[k];

	    char kbuf[5];
	    sprintf(kbuf, "%d", k);
	    string outfile = myargs.outdir + fname + myargs.suffix + kbuf + ".gz";
	    writecv(cv, outfile);
	}
    }

    cerr << "*** Time Elapsed: " << mytimer.elapsed() << "s" << endl;
};


Args::Args(int argc, char** argv):gtype("faa"),method(1),indir(""),outdir(""){

    program = argv[0];
    string listfile("list");
    string listkval("3 4 5 6 7");
    string onefasta("");

    char ch;
    while((ch = getopt(argc, argv, "I:i:k:O:g:S:s:f:h")) != -1){
      switch(ch){
      case 'I':
	  indir = optarg; addsuffix(indir, '/'); break;
      case 'i': 
	  listfile = optarg; break;
      case 'O':
	  outdir = optarg; addsuffix(outdir, '/'); break;
      case 'g':
	  gtype = optarg; break;
      case 'k':
	  listkval = optarg; break;
      case 'S':
	  method = str2int(optarg); break;
      case 's':
          suffix = optarg;
          if(*(suffix.begin()) != '.')
              suffix = '.' + suffix;
	  break;
      case 'f':
          onefasta = optarg; break;
      case 'h':
	  usage();
      case '?':
	  usage();
      }
    }

    //check the genome type
    if(gtype != "faa" && gtype != "ffn" && gtype != "fna"){
	cerr << "Only faa/ffn/fna are supported!\n" << endl;
	exit(1);
    }

    //get the input file name
    if(onefasta.empty()){
        ifstream list(listfile.c_str());
        if(!list){
            cerr << "Cannot found the input file " << listfile << endl;
            exit(1);
        }

        set<string> tmpSet;
        for(string line; getline(list, line); ){
            string fname = trim(line);
            if(!fname.empty()){
                if(getsuffix(fname) != gtype)
                    fname += "." + gtype;
                string infile = indir + fname;
                if(tmpSet.insert(fname).second)
                    flist.emplace_back(fname,getFileSize(infile));
            }
        }
        list.close();
        sort(flist.begin(),flist.end(),bySecond);
    }else{
        flist.emplace_back(onefasta, 0);
    }    

    // get the kvalue
    vector<string> wd;
    separateWord(wd, listkval);
    for(auto &str : wd){
	size_t kval;
	str2number(str, kval);
	klist.insert(kval);
    }

    //get the statistic list
    if(method == 1){
	for(auto &k : klist){
	    slist.insert(k);
	    slist.insert(k-1);
	    slist.insert(k-2);
	}
    }else{
	slist = klist;
    }

    // determine the cv suffix
    if(suffix.empty()){
        if(method == 0)
            suffix = ".ncv";
        else
            suffix = ".cv";
    }
}

void Args::kcheck(int kmax){
    // check k value range
    if(*(slist.begin()) < 1 || *(--slist.end())> kmax){
	cerr << "The range of k value is ["
	     << ((method == 1) ? 3 : 1)
	     << "," << kmax << "]" << endl;
       exit(3);
    }
};


void Args::usage(){
    cerr << "\nProgram Usage: \n\n" 
	 << program  << "\n"
	 <<" [ -I <dir> ]        input genome file directory\n"
	 <<" [ -i list ]         input species list, defaut: list\n"
         <<" [ -f <Fasta> ]      get cv for only one fasta \n"
	 <<" [ -k '3 4 5 6 7' ]  values of k, defaut: N = 3 4 5 6 7\n"
	 <<" [ -g faa ]          the type of genome file, defaut: faa\n"
	 <<" [ -O <dir> ]        output cv directory\n"
	 <<" [ -S 0/1]           whethe do the subtract, defaut: 1\n"
         <<" [ -s .ncv/.cv ]     Suffix of cv file, \n"
         <<"                     defaut: .ncv for method 0, .cv for method 1\n"
	 <<" [ -h ]              disply this information\n"
	 << endl;

    exit(1);
}

size_t count(const Genome& genome, size_t k, CVmap& cv){
    size_t n(0);
    for(const auto &gene : genome){
	//the number of kstring of the gene
        n += (gene.size()-k+1);

	//get the first k string
	Kstr ks(gene.substr(0,k));
	CVmap::iterator iter=cv.find(ks);
	if(iter == cv.end()){
	    cv[ks] = 1.0;
	}else{
	    ++(iter->second);
	}

	//get the next kstring
	for(int i=k; i<gene.size(); ++i){
	    ks.behead();
	    ks.append(gene[i]);
	    CVmap::iterator iter=cv.find(ks);
	    if(iter == cv.end()){
		cv[ks] = 1.0;
	    }else{
		++(iter->second);
	    }
	}
    }
    return n;
};


void subtract(const CVmap& mck, const CVmap& mckM1, const CVmap& mckM2, 
	      double factor, CVmap& cv){

    CVmap::const_iterator iter;
    for(const auto &cd : mckM2){
	Kstr ksM2 = cd.first;
	double nksM2 = cd.second;
	for(auto &c : ksM2.charSet){
	    Kstr ksM1A = ksM2;
	    ksM1A.append(c);
	    iter = mckM1.find(ksM1A);
	    if(iter != mckM1.end()){
		double nksM1A = iter->second;

		for(auto &d : ksM2.charSet){
		    Kstr ksM1B = ksM2;
		    ksM1B.addhead(d);
		    iter = mckM1.find(ksM1B);
		    if(iter != mckM1.end()){
			double nksM1B = iter->second;
			double nks0 = factor * nksM1B*nksM1A/nksM2;

			Kstr ks = ksM1B;
			ks.append(c);
			iter = mck.find(ks);
			double nks = 0;
			if(iter != mck.end())
			    nks = iter->second;
			cv[ks] = (nks - nks0)/nks0;
		    }
		}
	    }
	}
    }
    
};

bool bySecond(const pair<string,long>& a, const pair<string,long>& b){
    return a.second > b.second;
};

