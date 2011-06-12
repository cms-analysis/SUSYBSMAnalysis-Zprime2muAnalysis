#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V00-05-00 MuonAnalysis/Examples
cvs co -r V01-13-00 MuonAnalysis/MuonAssociators
cvs co -r V00-03-00 -d UserCode/Examples UserCode/AEverett/Examples
cvs co -r V00-05-00 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer
cvs co -r V06-00-06 PhysicsTools/RooStatsCms

# Patch to turn off a rogue cout.
addpkg RecoVertex/VertexTools
cvs update -r 1.3 RecoVertex/VertexTools/src/InvariantMassFromVertex.cc

# Temporary hack until Sam provides .xml BuildFile.
scram b -c

# Overwrite the stock rescalc.cc with our own. (Could just make a bin/
# dir of our own but this is good enough for now.)
ln -sf $CMSSW_BASE/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/ResonanceCalculator/rescalc.cc $CMSSW_BASE/src/PhysicsTools/RooStatsCms/bin/rescalc.cc

scramv1 b -j 4
popd
