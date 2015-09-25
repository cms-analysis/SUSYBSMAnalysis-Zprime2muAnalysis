#!/bin/tcsh

pushd $CMSSW_BASE/src

git clone https://github.com/cms-sw/RecoLuminosity-LumiDB.git RecoLuminosity/LumiDB
cd $CMSSW_BASE/src/RecoLuminosity/LumiDB
git checkout V04-02-10
cd ../..

git clone -b 72XRelease git@github.com:Sam-Harper/usercode.git SHarper
#addpkg DataFormats/MuonReco
#cvs update -r V09-06-00 DataFormats/MuonReco/interface/MuonCocktails.h
#cvs update -r V09-06-00 DataFormats/MuonReco/src/MuonCocktails.cc

##- Needed for MET filters and MET
##- Taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETRecipe53X
##- Works for CMSSW_5_3_11 to CMSSW_5_3_11_patch3
#addpkg DataFormats/PatCandidates V06-05-06-12
#addpkg PhysicsTools/PatAlgos V08-09-62
#addpkg PhysicsTools/PatUtils V03-09-28
#addpkg RecoMET/METAnalyzers V00-00-08
#addpkg DataFormats/METReco V03-03-11-01
#addpkg JetMETCorrections/Type1MET V04-06-09-02

scramv1 b -j 8
popd
