#!/bin/tcsh

pushd $CMSSW_BASE/src

#cvs co -r V03-05-00 RecoLuminosity/LumiDB # Might still be needed?
cvs co -r V00-08-11 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer

scramv1 b -j 8
popd
