#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V03-03-12 RecoLuminosity/LumiDB
cvs co -r V01-13-00 MuonAnalysis/MuonAssociators
cvs co -r V00-07-00 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer

# Patch to turn off a rogue cout.
addpkg RecoVertex/VertexTools
cvs update -r 1.3 RecoVertex/VertexTools/src/InvariantMassFromVertex.cc

scramv1 b -j 8
popd
