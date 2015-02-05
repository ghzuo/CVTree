#include "kstring.h"

size_t Kstr::nbase;
char Kstr::cmap[128];
vector<char> Kstr::charSet;

int Kstr::init(const vector<char>& letters){
    nbase = letters.size()+1;
    charSet = letters;
    for(int i=1; i<nbase; ++i){
	cmap[letters[i-1]] = i;
	cmap[i] = letters[i-1];
    }

    float logbs = log2(nbase);
    float sz = sizeof(unsigned long)*8;
    return sz/logbs;
};

Kstr::Kstr():ks(0){};

Kstr::Kstr(unsigned long s){
    ks = s;
};

Kstr::Kstr(const string& str):ks(0){
    for(auto &c : str)
	ks = ks*nbase + cmap[c];
};

string Kstr::decode() const{
    string str;
    for(unsigned long s=ks; s!=0; s /= nbase)
	str += char(cmap[s%nbase]);
    reverse(str.begin(), str.end());
    return str;
}

size_t Kstr::length() const{
    size_t n(0);
    unsigned long unit(1);
    while(unit <= ks){
	++n;
	unit *= nbase;
    }
    return n;
};

void Kstr::append(char c){
    ks *= nbase;
    ks += cmap[c];
};

void Kstr::addhead(char c){
    unsigned long unit(1);
    while(unit <= ks)
	unit *= nbase;
    ks += (unit*cmap[c]);
};

void Kstr::choptail(){
    ks /= nbase;
};

void Kstr::behead(){
    unsigned long unit(1);
    while(unit <= ks)
	unit *= nbase;
    unit /= nbase;
    ks%=unit;
};

void Kstr::forward(char c){
    behead();
    append(c);    
}

void Kstr::backward(char c){
    choptail();
    addhead(c);    
}

bool Kstr::operator<(const Kstr& r) const{
    return ks < r.ks;
};

bool Kstr::operator>(const Kstr& r) const{
    return ks > r.ks;
};

bool Kstr::operator==(const Kstr& r) const{
    return ks == r.ks;
};

bool Kstr::operator!=(const Kstr& r) const{
    return ks != r.ks;
};

bool Kstr::operator()(const Kstr& a, const string& b) const{
    return a<b;    
};

ostream& operator<<(ostream& os, const Kstr &ks){
    os << ks.decode();
    return os;
};

/// output cv in binary file and gzip it
void writecv(const CVmap& cv, const string& file){

    CVvec kslist(cv.size());
    int index(-1);
    for(const CVdim& cdim : cv)
        kslist[++index] = cdim;
    sort(kslist.begin(), kslist.end(),
         [](const CVdim& a, const CVdim& b){return a.first < b.first;});
    writecv(kslist, file);
};

// output cv vector file and gzip it
void writecv(const CVvec& cv, const string& file){
    double inner(0);
    unsigned long size = cv.size();

    for(const CVdim& c : cv)
	inner += c.second * c.second;

    gzFile fp;
    if((fp=gzopen(file.c_str(),"wb"))==NULL){
	cerr << "Error happen on write cvfile: " << file << endl;
	exit(1);
    }

    gzwrite(fp, &inner, sizeof(double));
    gzwrite(fp, &size,  sizeof(unsigned long));

    for(const CVdim& s : cv)
    	gzwrite(fp, &s, sizeof(CVdim));
    gzclose(fp);
    
}

//Read the cv from binary gz file
//read the cv from operation on the cv vector
//used for the distance calculate
double readcv(const string& filename, CVvec& cv){
    gzFile fp;
    if((fp=gzopen(filename.c_str(),"rb"))==NULL){
	cerr<<"CV file not found: \""<< filename << '"' << endl;
	exit(1);
    }

    double inner;
    gzread(fp,(char*)&inner, sizeof(double));
    double m = sqrt(inner);
    
    mlong size;
    gzread(fp,(char*)&size, sizeof(mlong));

    CVdim cdim;
    while(gzread(fp,(char*)&cdim, sizeof(CVdim))>0){
	cdim.second /= m;
	cv.emplace_back(cdim);
    }
    gzclose(fp);

    return inner;
}

size_t cvsize(const string& filename){
    gzFile fp;
    if((fp=gzopen(filename.c_str(),"rb"))==NULL){
	cerr<<"CV file not found: \""<< filename << '"' << endl;
	exit(1);
    }

    double inner;
    gzread(fp,(char*)&inner, sizeof(double));
    size_t size;
    gzread(fp,(char*)&size, sizeof(size_t));
    gzclose(fp);
    
    return size*(sizeof(long)+sizeof(double));
};


double module(const CVvec& cv){
    CVvec::const_iterator iter = cv.begin();
    CVvec::const_iterator iterEnd = cv.end();
    double mm(0);
    for(;iter!=iterEnd;++iter)
	mm += ((*iter).second)*((*iter).second);
    return sqrt(mm);
}

void normalize(CVvec& cv){
    double m = module(cv);
    CVvec::iterator iter = cv.begin();
    CVvec::iterator iterEnd = cv.end();
    for(;iter!=iterEnd;++iter)
	(*iter).second /= m;
}

double dist(const CVvec& cv1, const CVvec& cv2){
    double d(1);
    CVvec::const_iterator iter1 = cv1.begin();
    CVvec::const_iterator iter1End = cv1.end();
    CVvec::const_iterator iter2 = cv2.begin();
    CVvec::const_iterator iter2End = cv2.end();

    for(;;){
	if((*iter1).first < (*iter2).first){
	    if(++iter1 == iter1End) break;
	}else if((*iter1).first > (*iter2).first){
	    if(++iter2 == iter2End) break;
	}else{
	    d -= (*iter1).second * (*iter2).second;
	    if(++iter1 == iter1End) break;
	    if(++iter2 == iter2End) break;
	}
    }
    return 0.5*d;
}

// operation for kstr vector
// used for the missing K string
void readvk(const string& filename, vector<Kstr>& vk){
    gzFile fp;
    if((fp=gzopen(filename.c_str(),"rb"))==NULL){
	cerr<<"CV file not found: \""<< filename << '"' << endl;
	exit(1);
    }

    double inner;
    gzread(fp,(char*)&inner, sizeof(double));
    
    mlong size;
    gzread(fp,(char*)&size, sizeof(mlong));

    CVdim cdim;
    while(gzread(fp,(char*)&cdim, sizeof(CVdim))>0)
	vk.emplace_back(cdim.first);
    gzclose(fp);
}

void writevk(const string& file, const vector<Kstr>& vk){
    
    unsigned long size = vk.size();
    double inner(1.0);

    gzFile fp;
    if((fp=gzopen(file.c_str(),"wb"))==NULL){
	cerr << "Error happen on write cvfile: " << file << endl;
	exit(1);
    }

    gzwrite(fp, &inner, sizeof(double));
    gzwrite(fp, &size,  sizeof(unsigned long));

    for(const auto &s : vk){
	CVdim cv(s, 1.0);
    	gzwrite(fp, &cv, sizeof(CVdim));
    }
    gzclose(fp);
}


