// This module is deprecated.

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

bool leptonIsCut(const reco::Candidate& lepton) {
  if (lepton.pt() < 20)
    return true;

  const reco::Muon* muon = toConcretePtr<reco::Muon>(lepton);
  if (muon && muon->isolationR03().sumPt > 10)
    return true;

  return false;
}

using namespace std;
using namespace edm;
using namespace reco;

class DileptonPicker : public EDProducer {
public:
  explicit DileptonPicker(const ParameterSet&);
  
private:
  virtual void produce(Event&, const EventSetup&);
  
  InputTag src;
  unsigned maxDileptons;
};

DileptonPicker::DileptonPicker(const ParameterSet& cfg)
  : src(cfg.getParameter<InputTag>("src")),
    maxDileptons(cfg.getParameter<unsigned>("maxDileptons"))
{
  produces<CompositeCandidateCollection>();
}

void DileptonPicker::produce(Event& event, const EventSetup& eSetup) {
  static const bool debug = false;

  Handle<CompositeCandidateCollection> hdileptons;
  event.getByLabel(src, hdileptons);

  // Make the output collection.
  auto_ptr<CompositeCandidateCollection> dileptons(new CompositeCandidateCollection);

  ostringstream out;

  if (!hdileptons.failedToGet()) {
    if (hdileptons->size() > 0) {
      if (debug) out << hdileptons->size() << " dileptons before selection and overlap removal; ";

      for (CompositeCandidateCollection::const_iterator dil = hdileptons->begin(); dil != hdileptons->end(); ++dil)
	if (!leptonIsCut(*dil->daughter(0)->masterClone()) && !leptonIsCut(*dil->daughter(1)->masterClone()))
	  dileptons->push_back(*dil);

      if (debug) out << dileptons->size() << " dileptons before overlap removal; ";
      
      // Remove dileptons of lower invariant mass that are comprised of a
      // lepton that has been used by a higher invariant mass one.
      if (dileptons->size() > 1)
	removeDileptonOverlap(*dileptons, debug);

      if (dileptons->size() > maxDileptons)
	dileptons->erase(dileptons->begin() + maxDileptons, dileptons->end());

      if (debug) out << dileptons->size() << " dileptons kept.";
    }
    else {
      if (debug) out << "No dileptons in source collection.";
    }
    
    if (debug) LogInfo("DileptonPicker") << out.str();
  }
  else
    LogWarning("DileptonPicker")
      << "no collection " << src << " in event; producing empty collection";

  event.put(dileptons);
}

DEFINE_FWK_MODULE(DileptonPicker);
