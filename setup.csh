#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V00-03-04 MuonAnalysis/Examples
cvs co -r V01-09-01 MuonAnalysis/MuonAssociators
cvs co -r V00-05-00 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer
cvs co -r V00-00-02 -d UserCode/Examples UserCode/AEverett/Examples
cvs co -r lumi2010-Oct12 RecoLuminosity/LumiDB
cvs co -r V01-06-09 FWCore/PythonUtilities

echo $CMSSW_VERSION | grep -v CMSSW_3_8_ > /dev/null
if (! $?) then
  cvs co -r V02-06-01 SimMuon/MCTruth
  cvs co -r V01-08-17 SimTracker/TrackAssociation
  cvs co -r V00-01-03-01 SimTracker/TrackAssociatorESProducer
endif

# fix InputTag.h location
#find . -type f \! -name setup.csh | xargs grep -l "interface/InputTag.h" | xargs -I {} sed -i {} -e "s#FWCore/ParameterSet/interface/InputTag.h#FWCore/Utilities/interface/InputTag.h#g"
sed -i SimMuon/MCTruth/test/testAssociatorRecoMuon.cc -e "s#FWCore/ParameterSet/interface/InputTag.h#FWCore/Utilities/interface/InputTag.h#g"

scramv1 b -j 4
popd
