#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TeVMuHelper.h"

using namespace std;
using namespace edm;
using namespace reco;

class DileptonPicker : public EDProducer {
public:
  explicit DileptonPicker(const ParameterSet&);
  ~DileptonPicker() {}
  
private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() {}
  
  InputTag src;
  bool doingGen;
  unsigned maxDileptons;
  vector<int> pdgIds;
};

DileptonPicker::DileptonPicker(const ParameterSet& cfg)
  : src(cfg.getParameter<InputTag>("src")),
    doingGen(cfg.getParameter<bool>("doingGen")),
    maxDileptons(cfg.getParameter<unsigned>("maxDileptons")),
    pdgIds(cfg.getParameter<vector<int> >("pdgIds"))
{
  produces<CompositeCandidateCollection>();
}

void DileptonPicker::produce(Event& event, const EventSetup& eSetup) {
  static const bool debug = false;

  Handle<CompositeCandidateCollection> hdileptons;
  event.getByLabel(src, hdileptons);

  // Make the output collection.
  auto_ptr<CompositeCandidateCollection>
    dileptons(new CompositeCandidateCollection);

  ostringstream out;

  if (!hdileptons.failedToGet()) {
    if (hdileptons->size() > 0) {
      if (debug) out << hdileptons->size() << " dileptons before selection and overlap removal; ";

      if (doingGen)
	genDileptonsOnly(*hdileptons, *dileptons, debug);
      else {
	// Default cuts are pT < 20 GeV and isolation sumPt < 10.
	const unsigned cuts = TeVMuHelper::PT | TeVMuHelper::ISO;
	cutDileptons(event, *hdileptons, *dileptons, cuts,
		     pdgIds[0], pdgIds[1], debug);

	if (debug) out << dileptons->size() << " dileptons before overlap removal; ";

	// Remove dileptons of lower invariant mass that are comprised of a
	// lepton that has been used by a higher invariant mass one.
	if (dileptons->size() > 1)
	  removeDileptonOverlap(*dileptons, debug);
      }

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
