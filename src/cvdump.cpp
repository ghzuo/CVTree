#include "global.h"
#include "stringOpt.h"
#include "readgenome.h"
#include "kstring.h"

void usage(string& program){
    cerr << "\nProgram Usage: \n\n" 
	 << program  
	 <<"  -i <cvfile>        input file name\n"
     <<" [ -g faa|ffn ]      the type of genome file, defaut: faa\n"
	 <<" [ -n ]              output the number code, default: the letters\n"
	 <<" [ -h ]              disply this information\n"
	 << endl;
    exit(1);
}

int main(int argc, char* argv[]){

    string gtype = "faa";
    string infile;
    string program = argv[0];
    bool kstr = true;

    char ch;    
    while ((ch = getopt(argc, argv, "i:g:nh")) != -1){
      switch (ch){
      case 'i': 
	  infile = optarg; break;
      case 'g':
	  gtype = optarg; break;
      case 'n':
	  kstr = false; break;
      case 'h':
	  usage(program);
      case '?':
	  usage(program);
      }
    }

    //get gene type to read and read gene file
    GeneType mygene(gtype);
    Kstr::init(mygene.letters);
    
    CVvec cv;
    cout << "The inner of CV: " << readcv(infile, cv) <<endl;
    cout << "The size  of CV: " << cv.size() << endl;
    
    for(const auto& cd : cv)
        cout << cd.first <<"\t" <<cd.second << endl;
}
