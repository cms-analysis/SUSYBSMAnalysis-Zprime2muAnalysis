#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V00-03-04 MuonAnalysis/Examples
cvs co -r V01-09-01 MuonAnalysis/MuonAssociators
cvs co -r V00-05-00 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer
cvs co -r V00-00-02 -d UserCode/Examples UserCode/AEverett/Examples
cvs co -r V02-01-01 RecoLuminosity/LumiDB
cvs co -r V01-07-01 FWCore/PythonUtilities
cvs co -D "11/19/2010" HeavyFlavorAnalysis/Onia2MuMu/certification/7TeV/Collisions10/Reprocessing/Cert_132440-144114_7TeV_Sep17ReReco_Collisions10_JSON_BPAG.txt  
cvs co -D "11/19/2010" HeavyFlavorAnalysis/Onia2MuMu/certification/7TeV/Collisions10/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_BPAG_v2.txt

# fix InputTag.h location
#find . -type f \! -name setup.csh | xargs grep -l "interface/InputTag.h" | xargs -I {} sed -i {} -e "s#FWCore/ParameterSet/interface/InputTag.h#FWCore/Utilities/interface/InputTag.h#g"
sed -i SimMuon/MCTruth/test/testAssociatorRecoMuon.cc -e "s#FWCore/ParameterSet/interface/InputTag.h#FWCore/Utilities/interface/InputTag.h#g"

scramv1 b -j 4
popd
