#ifndef TEVMUHELPER_H
#define TEVMUHELPER_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

class TeVMuHelper {
public:
  TeVMuHelper() : ptCut(20), isoCut(10) {}

  enum CutResult {PASS, PT, ISO};

  unsigned leptonIsCut(const reco::Candidate& lepton) const {
    unsigned retval = 0;

    // pT cut.
    if (lepton.pt() <= ptCut)
      retval |= PT;

    // Sum pT in cone of dR=0.3 cut
    // Try to cast to reco::Muon.
    const reco::Muon* muon = toConcretePtr<reco::Muon>(lepton);
    if (muon != 0 && muon->isIsolationValid() &&
	muon->getIsolationR03().sumPt > isoCut)
      retval |= ISO;

    return retval;
  }

private:
  double ptCut;
  double isoCut;
};

#endif
