#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

using namespace std;
using namespace edm;
using namespace reco;

// It is possible that the modules which produce the HLT collections
// of leptons may not be run in the HLT path (if L1 did not fire the
// L2 module will not run, same for L2 and L3), in which case no
// collection is put in the event. If that is the case, put an empty
// one so we don't crash later down the analysis path, or copy the
// existing one if it exists.  (We are already doing something special
// for L3 muons in L3MuonSanitizer, which includes this
// functionality.)

template <typename CollectionType>
class HLTLeptonSanitizer : public EDProducer {
public:
  explicit HLTLeptonSanitizer(const ParameterSet&);
  ~HLTLeptonSanitizer() {}
  
private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() {}
  
  InputTag src;
};


template <typename CollectionType>
HLTLeptonSanitizer<CollectionType>::HLTLeptonSanitizer(const ParameterSet& cfg) {
  src = cfg.getParameter<InputTag>("src");
  produces<CollectionType>();
}

template <typename CollectionType>
void HLTLeptonSanitizer<CollectionType>::produce(Event& event,
						 const EventSetup& eSetup) {
  // Try to get the HLT lepton collection from the event.
  Handle<CollectionType> leptons;
  event.getByLabel(src, leptons);

  // Make the output collection.
  auto_ptr<CollectionType> cands(new CollectionType);

  if (!leptons.failedToGet())
    // Copy the existing one.
    *cands = *leptons;

  // Fail silently to prevent spamming messages. The user will know when
  // something is amiss from empty plots.

  //else
  //  edm::LogWarning("HLTLeptonSanitizer")
  //    << "no collection " << src << " in event; producing empty collection";

  event.put(cands);
}

typedef HLTLeptonSanitizer<RecoChargedCandidateCollection> L2MuonSanitizer;
typedef HLTLeptonSanitizer<RecoEcalCandidateCollection> L2ElectronSanitizer;
typedef HLTLeptonSanitizer<ElectronCollection> L3ElectronSanitizer;

DEFINE_FWK_MODULE(L2MuonSanitizer);
DEFINE_FWK_MODULE(L2ElectronSanitizer);
DEFINE_FWK_MODULE(L3ElectronSanitizer);
