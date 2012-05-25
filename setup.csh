#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V03-05-00 RecoLuminosity/LumiDB
cvs co -r V02-00-01 MuonAnalysis/MuonAssociators
cvs co -r V00-08-07 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer

scramv1 b -j 8
popd
