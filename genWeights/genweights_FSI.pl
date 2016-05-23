#! /usr/bin/perl
#

use strict;
use warnings;
use Cwd;
use File::Basename;

if (@ARGV < 3) {
  print "Usage: ./genweights_FSI.pl [input ROOT ntuple file path&name, wildcards incl.] [ Reweight option: 0:FSI&SI, 1:FSI only, 2:SI only ] [t2k or sk] [+1 or -1 for horn current]\n";
  exit;
}

system("hostname");
system("date");

my $wd = Cwd::getcwd();
print "In: ", $wd, "\n";

my $BINDIR="/home/xiaoyue/atmFitTools/genWeights";

my $EXEC="$BINDIR/genWeightsFromSK_NEUTFSI.exe";

my $INFPATH=$ARGV[0];

my $strOpt="";
my $outDir="";

my $OptFSI=$ARGV[1];
my $det=$ARGV[2];
my $horn="";

if (@ARGV > 3) {
    $horn=$ARGV[3];
}

if ($OptFSI==0) {
  $strOpt="FSISI";
  if ($horn eq "+1" && $det eq "t2k") {
      $outDir="/disk2/usr5/xiaoyue/t2kmc14a/nu_mode/fsisi/";
  } elsif ($horn eq "-1" && $det eq "t2k") {
      $outDir="/disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/fsisi/";
  } elsif ($det eq "sk") {
      $outDir="/disk2/usr5/xiaoyue/skmc/test_fc/fsisi/";
  }
}
elsif ($OptFSI==1) {
  $strOpt="FSI";
  if ($horn eq "+1" && $det eq "t2k") {
      $outDir="/disk2/usr5/xiaoyue/t2kmc14a/nu_mode/fsi/";
  } elsif ($horn eq "-1" && $det eq "t2k") {
      $outDir="/disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/fsi/";
  } elsif ($det eq "sk") {
      $outDir="/disk2/usr5/xiaoyue/skmc/test_fc/fsi/"; 
  }
}
elsif ($OptFSI==2) {
  $strOpt="SI";
  if ($horn eq "+1" && $det eq "t2k") {
      $outDir="/disk2/usr5/xiaoyue/t2kmc14a/nu_mode/si/";
  } elsif ($horn eq "-1" && $det eq "t2k") {
      $outDir="/disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/si/";
  } elsif ($det eq "sk") {
      $outDir="/disk2/usr5/xiaoyue/skmc/test_fc/si/"; 
  }
}
else {
  print "Invalid option : ", $OptFSI, "\n";
  exit(-1);
}

print "Reweight mode:", $strOpt, "\n";

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

#system("ln -s $BINDIR/pi_fsi_tbl.txt .");

for (my $ifile=0; $ifile<$nfiles; $ifile++) {
#for (my $ifile=0; $ifile<1; $ifile++) {
  
  my $fnamnoext=$list[$ifile];
  $fnamnoext =~ s|.root||;# remove extension

  system("$EXEC -s $srcdir/$list[$ifile] -o ${outDir}${fnamnoext}_${strOpt}_wgt.root -r $OptFSI >& ${outDir}${fnamnoext}_${strOpt}.log1");
  system("/disk/sklb/software/others/WeightOut/WeightOut ${outDir}${fnamnoext}_${strOpt} >& ${outDir}${fnamnoext}_${strOpt}.log2");
  system("cat ${outDir}${fnamnoext}_${strOpt}.log1 ${outDir}${fnamnoext}_${strOpt}.log2 > ${outDir}${fnamnoext}_${strOpt}.log");
  system("rm ${outDir}${fnamnoext}_${strOpt}.log1 ${outDir}${fnamnoext}_${strOpt}.log2");
}

print "Done!\n";
system("date");
