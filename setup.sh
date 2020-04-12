#!/bin/bash

# setupATLAS
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
#source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh
lsetup git
lsetup asetup

#asetup 21.2,AnalysisBase,latest,here
asetup 21.2.39,AnalysisBase,here
rm -f CMakeLists.txt

#lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"
