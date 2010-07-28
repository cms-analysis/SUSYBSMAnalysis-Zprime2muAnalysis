#!/bin/tcsh

pushd $CMSSW_BASE/src
cvs co -r V00-03-02 MuonAnalysis/Examples
scramv1 b -j 4
popd

