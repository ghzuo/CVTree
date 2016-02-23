#include "readgenome.h"

GeneType::GeneType(const string& str){
    for(int i=0; i<128; ++i)
	mc[i] = 'A';

    if(str.compare("faa") == 0)
	aainit();
    else if(str.compare("ffn") == 0)
	nainit();
    else if(str.compare("fna") == 0)
	nainit();
    else{
	cerr << "unknown genome file type!\n" << endl;
	exit(1);
    }
};

void GeneType::aainit(){
    string aa="ACDEFGHIKLMNPQRSTVWY";
    for(auto &c : aa){
	mc[c] = c;
	letters.emplace_back(c);
    }

    // for the unfomrated aa
    mc['B']='D';
    mc['U']='C';
    mc['X']='G';
    mc['Z']='E';
};

void GeneType::nainit(){
    string na="ACGT";
    for(auto &c : na){
	mc[c] = c;
	letters.emplace_back(c);
    }

    // for the unfomrated na
    mc['R']='A';//AG
    mc['Y']='T';//CT
    mc['M']='A';//AC
    mc['K']='T';//TG
    mc['S']='C';//CG
    mc['W']='A';//AT

    mc['B']='T';//TCG
    mc['D']='A';//ATG
    mc['H']='A';//ATC
    mc['V']='A';//ACG

    mc['N']='A';//ACGT
    mc['X']='A';//ACGT
};

size_t GeneType::readgene(string& file, Genome& genome) const{
    ifstream infile(file.c_str());
    if(! infile){
	cerr << "Cannot found the input file " << file << endl;
	exit(4);
    }
    //cout << " Read file: " << file << endl;

    for(string line; getline(infile, line); ){
        line = trim(line);
        if(line.empty()){
            
        }else if(line[0] == '>'){
            genome.emplace_back();
        }else{
            genome.back().append(line);
        }
    }
    infile.close();

    size_t len(0);
    for(auto &gene : genome){
	if (gene.size() == 0){
	    cerr << "Some empty gene in your genome file: " << file <<  endl;
	} else {
	    checkgene(gene);
	    len += gene.size();
	}
    }
    return len;
}

void GeneType::checkgene(string& str) const{
    if(*(str.rbegin()) == '*' ||*(str.rbegin()) == '-')
        str.pop_back();
    str = toUpper(str);
    for(auto &c : str)
	c = mc[c];
};

 
