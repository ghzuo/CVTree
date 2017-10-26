# CVTree
 
CVTree stands for `Composition Vector Tree` which is the implementation
of an alignment-free algorithm to generate a dissimilarity matrix from
comparatively large collection of DNA or Amino Acid sequences,
preferably whole-genome data, for phylogenetic studies.

Please visit our webserve of CVTree, where you can use the cvtree tools 
more convenient.  The CVTree3 Web Server have two identical but independent 
installations at cvtree.big.ac.cn (Beijing Institute of Genomics, Beijing) 
and tlife.fudan.edu.cn/cvtree3 (Fudan University, Shanghai).

Currently due to management issue only the one at 
http://cvtree.big.ac.cn can be used.

#### Main Programs
* cv:  Get the composition vector based on the fasta file of the genome.
* tree:  Get the phylogeny tree based on the composition vectors and
  neighbor-joint method.
* runCVTree.pl: run the two programs in one time.

## Installation

#### Preparation
* cmake >= 2.6
* g++ >= 4.8 or other compiler supporting C++11 standard
* Library: libz, netcdf, netcdf_cpp

#### Compiling
1. unzip the package file and change into it
2. mkdir build and change into it
3. cmake .. or add some options you wanted
4. make
5. make install

## Reference
* Ji Qi, Bin Wang, Bailin Hao (2004) Whole proteome prokaryote phylogeny
  without sequence alignment: a K-string composition approach, J Mol
  Evol, 58: 1–11
* Guanghong Zuo, Bailin Hao (2015) CVTree3 web server for
  whole-genome-based and alignment-free prokaryotic phylogeny and
  taxonomy, Genomics Proteomics & Bioinformatics, 13: 321-331

## License

This software is free for non-commercial use. For commercial use,
a software agreement is required.
