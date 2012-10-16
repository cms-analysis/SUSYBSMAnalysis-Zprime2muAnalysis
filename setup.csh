#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V04-01-09 RecoLuminosity/LumiDB
cvs co -r V00-08-11 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer

scramv1 b -j 8
popd
