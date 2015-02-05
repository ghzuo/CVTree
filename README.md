# CVTree
A whole-genome and alignment-free prokaryotic phylogeny tool

### cv: get the composition vector

Program Usage:

cv [ -I faa ]          input genome file directory, defaut: faa
   [ -i list ]         input species list, defaut: list
   [ -k '3 4 5 6 7' ]  values of k, defaut: N = 3 4 5 6 7
   [ -g faa ]          the type of genome file, defaut: faa
   [ -O cv ]           output cv directory, defaut: cv
   [ -S 0/1]           whethe do the subtract, defaut: 1
   [ -h ]              disply this information

### tree: get the phylogeny tree
tree [ -o dist.matrix ]   Output distance matrix, defaut: dist.matrix
     [ -I extdir ]        Directory of extend cv files, defaut: cv
     [ -i infile ]        Extend cv file list, defaut: no extend cv used
     [ -s cv6.gz ]        Suffix of cv file, defaut: cv6.gz
     [ -E orgdir  ]       Directory of the orginal cv files
     [ -m indist.matrix ] Input distance matrix, default: no input matrix used
     [ -e orglist ]       Name of selected genomes, which are in input distance matrix
     [ -n ndxlist ]       Index of selected genomes, which are input distance matrix
                          The index file used first!
     [ -t taxfile ]       input taxonomy information
     [ -T ]               Do not output taxonomy information
     [ -M <N> ]           Runing memory size as G roughly, default 80% physical memory
     [ -C ]               Force use the netcdf compress distance matrix
     [ -h ]               disply this information

### runCVTree.pl
run the two program in one time!
