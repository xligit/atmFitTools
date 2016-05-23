#! /usr/bin/perl
#

use strict;
use warnings;
use Cwd;
use File::Basename;

if (@ARGV != 3) {
  print "Usage: ./genweights_BANFF.pl [input ROOT ntuple file path&name, wildcards incl.] [ banff covariance matrix ] [output file path]\n";
  exit;
}

system("hostname");
system("date");

my $wd = Cwd::getcwd();
print "In: ", $wd, "\n";

my $BINDIR="/home/xiaoyue/atmFitTools/genWeights";

my $EXEC="$BINDIR/genWeightsFromSK_BANFF14c.exe";

my $INFPATH=$ARGV[0];
my $BANFFFILE=$ARGV[1];
my $OUTPATH=$ARGV[2];

my $srcname = basename($INFPATH);
my $srcdir = dirname($INFPATH);

chdir($srcdir);
$srcdir = Cwd::getcwd();
my @list = glob $srcname;
my $nfiles=@list;
chdir($wd);

print "Source dir.: ", $srcdir, "\n";
print "Source file: ", $srcname, "\n";
print " Processing ", $nfiles, " input files...", "\n";

for (my $ifile=0; $ifile<$nfiles; $ifile++) {
#for (my $ifile=0; $ifile<1; $ifile++) {  
  
  my $fnamnoext=$list[$ifile];
  $fnamnoext =~ s|.root||;# remove extension
  system("$EXEC -s $srcdir/$list[$ifile] -o ${OUTPATH}/${fnamnoext}_wgt.root -i ${BANFFFILE} >& ${OUTPATH}/${fnamnoext}_wgt.log");
}

print "Done!\n";
system("date");
