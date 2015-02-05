#include "distmatrix.h"

Mdist::Mdist():ng(0){};

// read the tranditional infile for the distance matrix
void Mdist::readmtx(const string& file){
    if(getsuffix(file) == "nc")
	readmtxnc(file);
    else
	readmtxtxt(file);
};

// read the tranditional infile for the distance matrix
void Mdist::readmtxtxt(const string& file){
    ifstream dd(file.c_str());
    if(!dd.is_open()){
	cerr << "Error opening file: " << file << endl;
	exit(1);
    }

    string line;
    getline(dd, line);
    resize(stoi(line));
    
    int idist(0);
    for(size_t i=0; i<ng; ++i){
        string line;
        if(getline(dd, line) && !line.empty()){
            istringstream iss(line);
            iss >> name[i];
            for(int j=0; j<i; ++j)
                iss >> dist[idist++];
            idist++;
        }
    }
};

// read the netcdf file for the distance matrix
void Mdist::readmtxnc(const string& file){
    NcFile mtxFile(file.c_str(), NcFile::ReadOnly);
    if (!mtxFile.is_valid()){
	cout << "Couldn't open file!\n";
	cerr <<  NC_ERR << endl;
    }

    //get the dim information
    long lenWord = mtxFile.get_dim("word")->size();
    long ngenome = mtxFile.get_dim("genome")->size();
    long nmtx = mtxFile.get_dim("matrix")->size();

    //resize the matrix
    resize(ngenome);
    
    // read the space name
    NcVar *spname=mtxFile.get_var("spname");
    long counts[2] = {1, lenWord};
    char aname[lenWord];
    for(int i=0; i<ngenome; ++i){
	spname->set_cur(i,0);
	spname->get(aname, counts);
        name[i] = aname;
    }

    // read the distance matrix
    float *fdist = new float[nmtx];
    NcVar *distance=mtxFile.get_var("distance");
    distance->get(fdist, nmtx);
    for(int i=0; i<nmtx; ++i)
        dist[i] = fdist[i];
    delete [] fdist;

    mtxFile.close();
};


// select the distance matrix file type
void Mdist::writemtx(const string& file){
    if(getsuffix(file) == "nc")
	writemtxnc(file);
    else
	writemtxtxt(file);
};


// write the tranditional infile for the distance matrix
void Mdist::writemtxtxt(const string& file){
    ofstream dd(file.c_str());
    if(!dd.is_open()){
	cerr << "Error opening file: " << file << endl;
	exit(1);
    }
	
    dd << ng << endl;
    for(size_t i=0; i<ng; ++i){
	dd << name[i] << " ";
	for(size_t j=0; j<ng; ++j)
	    dd << left << fixed << setprecision(10) 
	       << setw(13) << getdist(i, j);
	dd << endl;
    }
}

// write the tranditional infile for the distance matrix
void Mdist::writemtxnc(const string& file){
    NcFile mtxFile(file.c_str(), NcFile::Replace);
    if (!mtxFile.is_valid()){
	cout << "Couldn't open file!\n";
	cerr <<  NC_ERR << endl;
    }

    //output the space name
    size_t wlength(0);
    for(const auto &str : name)
	if(str.size()>wlength)
	    wlength = str.size();
    ++wlength;

    NcDim* lenWord = mtxFile.add_dim("word", wlength);
    NcDim* ngenome = mtxFile.add_dim("genome", name.size());
    NcVar *spname = mtxFile.add_var("spname", ncChar, ngenome, lenWord);
    for(long i=0; i<name.size(); ++i){
	long nameSize = name[i].size()+1;
    	long counts[2] = {1, nameSize};
    	spname->set_cur(i,0);
    	spname->put(name[i].c_str(),counts);
    }

    // output the distance matrix
    NcDim *mtxdist=mtxFile.add_dim("matrix", dist.size());
    NcVar *distance=mtxFile.add_var("distance", ncFloat, mtxdist);
    float *vdm = new float[dist.size()];
    for(int i=0; i<dist.size(); ++i)
        vdm[i] = dist[i];
    distance->put(vdm, dist.size());
    delete [] vdm;

    mtxFile.close();
}

// readjust the distance matrix and their name by the index list
void Mdist::reduce(const vector<size_t>& ndx){

    // put the name and distance in the tmp vector
    vector<double> tmpDist(dist);
    vector<string> tmpName(name);

    // resize the matrix
    resize(ndx.size());
    for(int i=0; i<ng; ++i){
        if(ndx[i] >= tmpName.size()){
            cerr << "The index " << ndx[i]
                 << " is out the size of matrix "
                 << tmpName.size() << endl;
            exit(3);
        }
            
        name[i] = tmpName[ndx[i]];
	for(int j=0; j<=i; ++j){
            size_t itmp = ndx[i]<ndx[j] ? ndx[i]+ndx[j]*(ndx[j]+1)/2 : ndx[j]+ndx[i]*(ndx[i]+1)/2;
            _setdist(j,i,tmpDist[itmp]);
        }
    }
}

// readjust the distance matrix and their name by the name list
void Mdist::reduce(const vector<string>& nlist){
    
    // get the index map of the name
    unordered_map<string, size_t> mapMI;
    for(size_t i=0; i<name.size(); ++i)
	mapMI[name[i]] = i;

    // resize the matrix
    vector<size_t> ndx;
    for(const auto &str : nlist)
	if(mapMI.find(str) == mapMI.end()){
	    cerr << "\nCannot find the genome " << str << "in the matrix file\n" << endl;
	    exit(3);
	}else
	    ndx.emplace_back(mapMI[str]);
    reduce(ndx);
};

// ... extend the matrix
void Mdist::extend(const vector<string>& nlist, const vector<double>& dd){

    int orgSize     = ng;
    int orgDistSize = dist.size();
    int newSize     = ng + nlist.size();
    int newDistsize = newSize * (newSize + 1) / 2;
        
    if( dd.size() != newDistsize - orgDistSize){
	cerr << "Error: the number of the append vector is unmatched!!" << endl;
	exit(4);
    }

    // resize the matrix
    resize(newSize);

    // add the name
    for(int i=0; i<nlist.size(); ++i)
        name[orgSize+i] = nlist[i];

    // add the distance
    for(int i=0; i<dd.size(); ++i)
        dist[orgDistSize+i] = dd[i];
};

// ... extend matrix by name list
void Mdist::extend(const vector<string>& nlist){
    int orgSize = ng;
    int newSize = ng + nlist.size();
    
    resize(newSize);
    for(int i=0; i<nlist.size(); ++i)
        name[orgSize+i] = nlist[i];
};

// .. reset the size of the matrix
void Mdist::resize(size_t n){
    ng = n;
    name.resize(n);
    dist.resize(n*(n+1)/2);
};

// option on sigle item 
double Mdist::getdist(size_t i, size_t j) const{
    if(i<j)
	return _getdist(i,j);
    return _getdist(j,i);
};

double Mdist::_getdist(size_t i, size_t j) const{
    if(j>=ng){
	cerr << "Error: the index out of matrix" << endl;
	exit(4);
    }

    return dist[i+(j+1)*j/2];
};

void Mdist::setdist(size_t i, size_t j, double d){
    if(i<j)
	return _setdist(i,j, d);
    return _setdist(j,i, d);
};

void Mdist::_setdist(size_t i, size_t j, double d){
    if(j>=ng){
	cerr << "Error: the index out of matrix" << endl;
	exit(4);
    }

    dist[i+j*(j+1)/2] = d;
};

string Mdist::getname(size_t i) const {
    return name[i];
};

void Mdist::setname(size_t i, const string& str){
    name[i] = str;
}

// global options
size_t Mdist::size() const{
    return ng;
};

size_t Mdist::msize() const{
    return ng*(ng+1)/2;
};

size_t Mdist::capacity() const{
    return dist.capacity()*sizeof(double) + sizeof(ng);
};

