#! /usr/bin/perl
#

use strict;
use warnings;
use Cwd;

if (@ARGV != 5) {
  print "Usage: ./genweights_throws.pl [input ROOT ntuple directory] [+1 for nu-mode, -1 for antinu-mode] [1 for appearance sample, 0 otherwise] [job name] [t2k or sk]\n";
  exit;
}

my $wd = Cwd::getcwd();
print $wd, "\n";

my $scriptdir;
my $banffwgtdir;

if ($ARGV[4] eq "t2k" && $ARGV[1] eq "-1") {
    system("mkdir -p /disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/banffwgt_scripts/");
    system("mkdir -p /disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/banffwgt/");
    $scriptdir="/disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/banffwgt_scripts/";
    $banffwgtdir="/disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/banffwgt/";
}

elsif ($ARGV[4] eq "t2k" && $ARGV[1] eq "+1") { 
    system("mkdir -p /disk2/usr5/xiaoyue/t2kmc14a/nu_mode/banffwgt_scripts/");
    system("mkdir -p /disk2/usr5/xiaoyue/t2kmc14a/nu_mode/banffwgt/");
    $scriptdir="/disk2/usr5/xiaoyue/t2kmc14a/nu_mode/banffwgt_scripts/";
    $banffwgtdir="/disk2/usr5/xiaoyue/t2kmc14a/nu_mode/banffwgt/";
}

elsif ($ARGV[4] eq "sk") {
    system("mkdir -p /disk2/usr5/xiaoyue/skmc/test_fc/banffwgt_scripts/");
    system("mkdir -p /disk2/usr5/xiaoyue/skmc/test_fc/banffwgt/");
    $scriptdir="/disk2/usr5/xiaoyue/skmc/test_fc/banffwgt_scripts/";
    $banffwgtdir="/disk2/usr5/xiaoyue/skmc/test_fc/banffwgt/";
}

else {
    print "Invalid input!\n";
    exit;
}

print $scriptdir, "\n", $banffwgtdir;

# location of T2KReWeight environment setup script
my $setup_sh="/home/xiaoyue/atmFitTools/genWeights/my_setup.sh";

my $parfile="/disk2/sklb/t2kmc/14a/postfit_data_1p1h_biascorrection_20160310/postfit_data_1p1h_biascorrection_20160310.root";

my $srcdir=$ARGV[0]; # directory where the SK ROOT ntuples are

my $iHorn=$ARGV[1]; # +1 for nu-mode, -1 for antinu-mode
my $iApp=$ARGV[2]; # Set to 1 for nue/nueb appearance sample

my $jobname=$ARGV[3];

my $nthrows=0; # number of systematic throws (first element in the weight array is based on the BANFF fit central value)

# Set any other genWeights options here (--use-prefit --drop-flux --drop-xsec)
#my $OtherOptions="";
my $OtherOptions;
if ($ARGV[4] eq "t2k") {
    $OtherOptions="-dslist /disk2/sklb/t2kmc/14a/postfit_data_1p1h_biascorrection_20160310/disable_list.txt";
} else {
    $OtherOptions="-dslist /disk2/sklb/t2kmc/14a/postfit_data_1p1h_biascorrection_20160310/disable_list.txt --drop-flux";
}

my $nMaxThreads=4; # maximum number of threads

#my $jobname;

# random number seed used for systematic throws. make sure to use the same value for all MC files!
#my $randsd = int(rand(2147483646));
my $randsd = 41855552;

chdir(${srcdir});
my @list = <*.root>;
my $nfls=@list;

print $nfls, "\n";

chdir($scriptdir);

if ($nfls > 0) {
  
  my $nseqs=0;
  my $nthreads=0;
  for (my $ifile=0; $ifile<$nfls; $ifile++) {
    if ($nthreads == 0) {#new job
      $nthreads = 1;
#      $jobname="${ifile}";
      open (SHS,">gen_${jobname}.sh");
      print SHS "#!/bin/bash

#scl enable devtoolset-2 bash

source ${setup_sh}

cd ${banffwgtdir}

date

\(
";
    }
    
    $nseqs++;
    
    $list[$ifile] =~ s/.root//;
    
    print SHS "/home/xiaoyue/atmFitTools/genWeights/genWeights_2015c.exe -i $srcdir/$list[$ifile].root -horn $iHorn -app $iApp -p $parfile -o $list[$ifile]_wgt.root -t $nthrows -r $randsd $OtherOptions
/disk/sklb/software/others/WeightOut/WeightOut $list[$ifile]
";
    
    if ($nseqs == 1000) {
      print SHS "\)&

\(
";
      $nseqs=0;
      $nthreads++;
    }
    
    if ($nthreads > $nMaxThreads || $ifile==$nfls-1) {
      print SHS "echo Hello!
\)&
      
################################################
wait
date
";
      close SHS;
      
      chmod 0755,"gen_${jobname}.sh";
      $nthreads=0;
    }
  }
  
}
else {
  print "Cannot find ROOT ntuple! \n";
}

