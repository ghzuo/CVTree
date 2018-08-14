# CVTree
 
CVTree stands for `Composition Vector Tree` which is the implementation
of an alignment-free algorithm to generate a dissimilarity matrix from
comparatively large collection of DNA or Amino Acid sequences,
preferably whole-genome data, for phylogenetic studies.

Please visit our webserve of CVTree, where you can use the cvtree tools 
more convenient.  The CVTree3 Web Server have two identical but independent 
installations at http://cvtree.big.ac.cn (Beijing Institute of Genomics, Beijing) 
and http://tlife.fudan.edu.cn/cvtree (Fudan University, Shanghai).

#### Main Programs
* cvtree: the main program, it get the phylogeny tree based from the fasta file of genomes.
* g2dv:  Get the composition vector based on the fasta file of the genome.
* cv2dm: Get the distance matrix based on the compostion vector
* dm2tree: Get the phylogeny tree from the distance matrix by neighbor-joint method.
* getdist: Show select distances from the distamce matrix
* cvdump: Show composition vector
* diffmtx: compare two distance matrixes


## Installation

#### Preparation
* cmake >= 2.6
* g++ >= 4.8 or other compiler supporting C++11 standard
* compiler with support openmp for parallel
* Library: libz, netcdf, netcdf_cpp

#### Compiling
1. unzip the package file and change into it
2. mkdir build and change into it
3. cmake .. or add some options you wanted
4. make
5. make manual
6. make install

## Run Programs with Example
If this is the first time you use CVTree package, please go to the
"example" folder. Edit "list" to include the genome names, and run 
the cvtree command to get the phylogeny tree by: `../build/cvtree -G faa`

## Todo
* The docuemt of the software.

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
