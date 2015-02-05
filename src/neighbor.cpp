#include "neighbor.h"

int main(int argc, char* argv[]){
    // set a timer
    boost::timer mytimer;

    // get the input arguments
    Args myargs(argc,argv);

    // init the distance matrix and name list
    Mdist dm;
    if(myargs.netcdf)
	dm.readmtxnc(myargs.distfile);
    else
	dm.readmtx(myargs.distfile);
    
    // made the star tree by listfile
    // if no, all items in distance matrix are used
    if(!myargs.splist.empty())
        dm.reduce(myargs.splist);
        
    // do the NJ algorithm and return the NJ tree
    Node *aTree = neighborJoint(dm);
    //delete the outgrouup
    (*aTree).children.pop_back();
    
    // output the Tree
    ofstream nwk(myargs.outfile.c_str());
    (*aTree).outnwk(nwk);
    nwk.close();

    cerr << "*** Time Elapsed: " << mytimer.elapsed() << "s" << endl;
}

Args::Args(int argc, char** argv):distfile("infile"),outfile("FullTaxonomy.nwk"),netcdf(false){

    program = argv[0];
    string wkdir("");
    string listfile;

    char ch;
    while ((ch = getopt(argc, argv, "i:d:o:D:Ch")) != -1){
      switch (ch){
      case 'i': 
	  listfile = optarg; break;
      case 'd':
	  distfile = optarg; break;
      case 'o':
	  outfile = optarg; break;
      case 'D':
	  wkdir = optarg; break;
      case 'C':
	  netcdf = true; break;
      case 'h':
	  usage();
      case '?':
	  usage();
      }
    }

    if(wkdir != ""){
	addsuffix(wkdir, '/');
	distfile = wkdir + distfile;
	outfile = wkdir + outfile;
    }

    if(!listfile.empty())
        readlist(listfile, splist);
}

void Args::usage(){
    cerr << "\nProgram Usage (VERSION: " << GHZ_VERSION << "): \n\n" 
	 << program  <<"\n" 
	 <<" [ -d infile ]        input distance matrix, defaut: dist.matrix\n"
	 <<" [ -o NJtree.nwk ]    output newick tree, defaut: NJtree.nwk\n"
	 <<" [ -i <listfile> ]    the selection index list of the distance matrix\n"
	 <<"                      if no defined, whole distance matrix are used\n"
	 <<" [ -C ]               use the netcdf input format, default false\n"
	 <<" [ -h ]               disply this information\n"
	 << endl;
    exit(1);
}

// select the leafs to of the output tree
void selectLeafs(const Mdist& dm, const string& listfile, vector<Node*>& nodes){
    if(listfile.empty()){
	nodes.resize(dm.size());
	for(size_t i=0; i<dm.size(); ++i){
	    nodes[i] = new Node(i);
	}
    }else{
	ifstream inf(listfile.c_str());
	if(!inf.is_open()){
	    cerr << "Error opening file: " << listfile << endl;
	    exit(3);
	}
	
	while(inf.good()){
	    Node* np = new Node;
	    inf >> (*np).id;
	    nodes.emplace_back(np);
	}
    }
}


