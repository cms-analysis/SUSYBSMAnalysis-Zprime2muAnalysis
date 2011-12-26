#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V03-03-12 RecoLuminosity/LumiDB
cvs co -r V00-05-00 MuonAnalysis/Examples
cvs co -r V01-13-00 MuonAnalysis/MuonAssociators
cvs co -r V00-03-00 -d UserCode/Examples UserCode/AEverett/Examples
cvs co -r V00-05-00 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer

# Patch to turn off a rogue cout.
addpkg RecoVertex/VertexTools
cvs update -r 1.3 RecoVertex/VertexTools/src/InvariantMassFromVertex.cc

# Temporary hack until Sam provides .xml BuildFile.
scram b -c

scramv1 b -j 8
popd
