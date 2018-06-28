#!/bin/bash
# Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2018-06-28 11:50:14
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2018-06-28 11:50:14

######################################
#####   The parameters for CVTree
######################################
## path of cv command
CVcomm="../build/bin/cv"
## path of tree command
TREEcomm="../build/bin/tree"
## list of K values
klist="3 4 5 6 7"
## folder include the fasta files
gdir="faa/"
## the type of genome: faa or ffn
gtype="faa"
## file include all name of genomes(one genome per line)
gfile="list";
## folder to save cache files (cv and distance matrix)
cache="cache/"
## output folder (the nwk files)
nwkdir="nwk/"
## type of cv file
cvmth=1
## set the number of cores
export OMP_NUM_THREADS=4

###############################################
#####   Run programs
#####   Please don't change the code below 
#####   if you don't sure what you do
################################################

## run the cv command
cvdir=$cache"cv/"
mkdir -p $cvdir
$CVcomm -i $gfile -k "$klist" -I $gdir -O $cvdir -S $cvmth -g $gtype

## run the tree command
dmdir=$cache"dm/"
mkdir -p $dmdir
mkdir -p $nwkdir

for k in $klist; do
    suffix=$gtype".cv"$k".gz"
    dmfile=$dmdir"K"$k"."$gtype".nc"
    outtre=$nwkdir"K"$k"."$gtype".nwk"
    $TREEcomm -C -o $outtre -I $cvdir -i $gfile -d $dmfile
done

