#!/usr/bin/perl -w
use Tie::File;
use strict;
use warnings;

# A perl script to run the CVTree stand alone
# Author: G.H.Zuo
(my $COMMDIR = $0) =~ s/[^\/]+$//;

### parameter for cv
my @klist  = qw/3 4 5 6 7/;
my $gdir   = "faa/";
my $gtype  = "faa";
my $gfile  = "list";
my $cvdir  = "cv/";
my $tredir = "result/";
my $cvmth  = 1;

### options for excution
my $MaxNP = 4;
my $quit  = 1;

### the commands
my $CVTREE   = "OMP_NUM_THREADS=$MaxNP " . $COMMDIR . "cv ";
$CVTREE .= "-S " . $cvmth . " ";
$CVTREE .= "-g " . $gtype . " ";

my $BUIDTREE = "OMP_NUM_THREADS=$MaxNP " . $COMMDIR . "tree ";

### Get the workdir and the datadir

### Obtain the cv files
mkdir $cvdir unless -d $cvdir;
my $kstr = '"' . (join " ", @klist) . '"';
my $CVTRCOMM = "$CVTREE -i $gfile -k $kstr -I $gdir -O $cvdir";
sysrun($CVTRCOMM);

### Obtain the tree
mkdir $tredir unless -d $tredir;
foreach my $k(@klist){
    my $suffix = '.' . $gtype . ".cv" . $k . ".gz";
    my $outtre = $tredir . "K" . $k . "." . $gtype . ".nwk";

    my $TREECOMM = $BUIDTREE . " -o $outtre -I $cvdir -i $gfile";

    ## do the update
    sysrun($TREECOMM);
}

###### subroutines to run command #########
sub sysrun{
    my $command = shift(@_);
    print ">>>>>>>>>>>>\$ ", $command, "\n\n" unless $quit;
    system($command) == 0 || die;
    return 0;
}
#### End of sysrun



