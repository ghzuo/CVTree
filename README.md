# CVTree

CVTree stands for `Composition Vector Tree` which is the implementation
of an alignment-free algorithm to generate a dissimilarity matrix from
comparatively large collection of DNA or Amino Acid sequences,
preferably whole-genome data, for phylogenetic studies.

Please visit our webserve of CVTree, where you can use the cvtree tools
more convenient. The CVTree3 Web Server have two identical but independent
installations at http://cvtree.big.ac.cn (Beijing Institute of Genomics, Beijing)
and http://tlife.fudan.edu.cn/cvtree (Fudan University, Shanghai).

## Main Programs

- cvtree: the main program, it get the phylogeny tree based from the fasta file of genomes.
- g2dv: Get the composition vector based on the fasta file of the genome.
- cv2dm: Get the distance matrix based on the compostion vector
- dm2tree: Get the phylogeny tree from the distance matrix by neighbor-joint method.
- getdist: Show select distances from the distamce matrix
- cvdump: Show composition vector
- diffmtx: compare two distance matrixes

## Installation

### Compile with CMake

#### Preparation

- cmake >= 3.0
- g++ >= 4.8 or other compiler supporting C++11 standard
- require ligrary: libz
- compiler with support openmp for parallel (_option_)
- Library (_option_): netcdf, netcdf_cpp
- Library (_option_): libhdf5 for c++

#### Compiling

1. unzip the package file and change into it
2. mkdir build and change into it
3. cmake .. or add some options you wanted
4. make
5. make manual (_option_)
6. make install

### Run Programms in Docker

Docker allows users run programs on both Windows and Linux/MacOS.
You can download docker free and reference https://docs.docker.com/install/
to install it. After install docker, basic usages for CVTree are:

1. Build/download docker image: `docker build -t="xxxxx/cvtree" .`
   or `docker pull ghzuo/cvtree`. In this step, a image with cvtree 
   programs will obtained. A small docker image, named "ghzuo/cvtree-lite"
   and have the size only 100M, can also obtained from docker hub. It 
   contains only the binary cvtree tools.
2. Start container from image:
   `docker run --rm -it -v $PWD/example:/root/data xxxxx/cpplearn`
   In this step, you will enter the cvtree container, and the "example" folder
   of this project will be find in the "data" folder. Enter the data folder,
   and run `cvtree -G faa`. You will get the result for eight genomes in the "list"
   file. You can change the path "$PWD/example" to your own data directory.
3. Exit and stop container: `exit` in docker terminal.
4. More usage for docker can reference https://docs.docker.com/.

## Run Programs with Example

If this is the first time you use CVTree package, please go to the
"example" folder. Edit "list" to include the genome names, and run
the cvtree command to get the phylogeny tree by:

    ../build/cvtree -G faa

More detail of the command usage can be obtaion by `-h` option.

## Todo

- More detail docuemt of the software.

## Reference

- Ji Qi, Bin Wang, Bailin Hao (2004) Whole proteome prokaryote phylogeny
  without sequence alignment: a K-string composition approach, J Mol
  Evol, 58: 1–11
- Guanghong Zuo, Bailin Hao (2015) CVTree3 web server for
  whole-genome-based and alignment-free prokaryotic phylogeny and
  taxonomy, Genomics Proteomics & Bioinformatics, 13: 321-331

## License

This software is free for non-commercial use. For commercial use,
a software agreement is required.
