#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V00-03-04 MuonAnalysis/Examples
cvs co -r V01-09-01 MuonAnalysis/MuonAssociators
#cvs co -r V02-04-00 RecoMuon/TrackingTools
cvs co -r V02-06-01 SimMuon/MCTruth
cvs co -r V01-08-17 SimTracker/TrackAssociation
cvs co -r V00-01-03-01 SimTracker/TrackAssociatorESProducer
cvs co -r V00-04-00 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer
cvs co -r V00-00-02 -d UserCode/Examples UserCode/AEverett/Examples

scramv1 b -j 4
popd
