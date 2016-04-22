#! /usr/bin/perl
#

use strict;
use warnings;
use Cwd;
use File::Basename;

if (@ARGV != 2) {
  print "Usage: ./genweights_throws.pl [input ROOT ntuple file path&name, wildcards incl.] [ Reweight option: 0:FSI&SI, 1:FSI only, 2:SI only ] \n";
  exit;
}

system("hostname");
system("date");

my $wd = Cwd::getcwd();
print "In: ", $wd, "\n";

my $BINDIR="/disk2/usr5/shimpeit/FSI_RW/code";

my $EXEC="$BINDIR/genWeightsFromSK_NEUTFSI.exe";

my $INFPATH=$ARGV[0];

my $strOpt="";
my $OptFSI=$ARGV[1];
if ($OptFSI==0) {
  $strOpt="FSISI";
}
elsif ($OptFSI==1) {
  $strOpt="FSI";
}
elsif ($OptFSI==2) {
  $strOpt="SI";
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

system("ln -s $BINDIR/pi_fsi_tbl.txt .");

for (my $ifile=0; $ifile<$nfiles; $ifile++) {
  
  my $fnamnoext=$list[$ifile];
  $fnamnoext =~ s|.root||;# remove extension
  
  system("$EXEC -s $srcdir/$list[$ifile] -o ${fnamnoext}_${strOpt}_wgt.root -r $OptFSI >& ${fnamnoext}_${strOpt}.log1");
  system("/disk/sklb/software/others/WeightOut/WeightOut ${fnamnoext}_${strOpt} >& ${fnamnoext}_${strOpt}.log2");
  system("cat ${fnamnoext}_${strOpt}.log1 ${fnamnoext}_${strOpt}.log2 > ${fnamnoext}_${strOpt}.log");
  system("rm ${fnamnoext}_${strOpt}.log1 ${fnamnoext}_${strOpt}.log2");
  
}

print "Done!\n";
system("date");
