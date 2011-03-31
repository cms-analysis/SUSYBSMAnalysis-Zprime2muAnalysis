#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V00-05-00 MuonAnalysis/Examples
cvs co -r V01-13-00 MuonAnalysis/MuonAssociators
cvs co -r V00-03-00 -d UserCode/Examples UserCode/AEverett/Examples
cvs co -r V00-05-00 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer

# Temporary hack until Sam provides .xml BuildFiles
scram b -c

scramv1 b -j 4
popd
