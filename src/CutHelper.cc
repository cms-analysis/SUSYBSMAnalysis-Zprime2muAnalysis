#include "DataFormats/MuonReco/interface/Muon.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/CutHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

bool leptonIsCut(const reco::Candidate& lepton) {
  if (lepton.pt() < 20)
    return true;

  const reco::Muon* muon = toConcretePtr<reco::Muon>(lepton);
  if (muon && muon->isolationR03().sumPt > 10)
    return true;

  return false;
}
