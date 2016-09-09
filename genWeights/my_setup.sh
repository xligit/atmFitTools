#!/bin/bash -f
#
# Environment setup script for T2KReWeight (called by example_config.sh)
#

export LANG C
export LC_ALL C

gcc -v |& grep 4.8 > /dev/null
if [ $? != 0 ];
then
   echo " Please run 'scl enable devtoolset-2 bash' at first. "
   exit 1
fi

export CXX=g++
export CPP="gcc -E"
export CC=gcc
export FC=gfortran
export F77=gfortran

export ROOTSYS=/home/skofl/sklib_gcc4.8.2/root_v5.28.00h

export GLOBALANATOOLS=/disk/sklb/software/gcc_482/GlobalAnalysisTools

# The following are the standard libs and env variables needed when 
# compiling T2KReWeight with oaAnalysis enabled
export T2KREWEIGHT=${GLOBALANATOOLS}/T2KReWeight_v1r23

export PATH=$T2KREWEIGHT/bin:$PATH:$T2KREWEIGHT/app:$ROOTSYS/bin:$PATH;
export LD_LIBRARY_PATH=$T2KREWEIGHT/lib:$ROOTSYS/lib:$LD_LIBRARY_PATH;


# ND280 oaAnalysis libs (See README for generating this)
# Only needed with ./configure --enable-oaanalysis
export OAANALYSISLIBS=$T2KREWEIGHT/libReadoaAnalysis
export LD_LIBRARY_PATH=$OAANALYSISLIBS:$LD_LIBRARY_PATH;


# GENIE and dependencies
# Only needed with ./configure --enable-genie
#export GENIE=/cygdrive/c/T2K/lib/genie/genie_2.6.2
#export PYTHIA6_LIB=/cygdrive/c/T2K/lib/pythia6/v6_412/lib
#export LIBXML_LIB=/usr/lib64
#export LOG4CPP_LIB=/usr/lib64
#export LHAPDF_LIB=/usr/local/lib
#export LD_LIBRARY_PATH=$LOG4CPP_LIB:$LIBXML_LIB:$LHAPDF_LIB:$PYTHIA6_LIB:$ROOTSYS/lib:$GENIE/lib:$LD_LIBRARY_PATH;
#export LIBXML_INC=/usr/include/libxml2
#export LOG4CPP_INC=/usr/include/log4cpp
#export LHAPDF_INC=/cygdrive/c/T2K/lib/lhapdf-5.8.5/include
#export LHAPATH=/usr/local/share/lhapdf
#export PATH=$GENIE/bin:$PATH;


# NEUT and CERNLIB
# Only needed with ./configure --enable-neut
export CERN=/home/skofl/sklib_gcc4.8.2/cern
export CERN_LEVEL=2005
export CERN_ROOT=${CERN}/${CERN_LEVEL}
export PATH=${CERN_ROOT}/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CERN_ROOT/lib
export CVSCOSRC=$CERN/$CERN_LEVEL/src

export NEUT_ROOT=/disk/sklb/software/gcc_482/neut_5.3.3_v1r21
export PATH=$NEUT_ROOT/src/neutsmpl/bin:$PATH
export LD_LIBRARY_PATH=$NEUT_ROOT/src/neutclass:$NEUT_ROOT/src/reweight:$LD_LIBRARY_PATH;


# JNuBeam
# Only needed with ./configure --enable-jnubeam
export JNUBEAM=${GLOBALANATOOLS}/JReWeight_v1r11
export LD_LIBRARY_PATH=${JNUBEAM}:$LD_LIBRARY_PATH;
export JREWEIGHT_INPUTS=${JNUBEAM}/inputs



# NIWG
# Only needed with ./configure --enable-niwg
export NIWG=${GLOBALANATOOLS}/NIWGReWeight_v1r19
export LD_LIBRARY_PATH=${NIWG}:$LD_LIBRARY_PATH;
export NIWGREWEIGHT_INPUTS=${NIWG}/inputs



# GEANT
# Only needed with ./configure --enable-geant
export GEANTRW=${GLOBALANATOOLS}/GEANTReWeight_v1r1
export LD_LIBRARY_PATH=${GEANTRW}:$LD_LIBRARY_PATH;
export GEANTREWEIGHT_INPUTS=${GEANTRW}/inputs

