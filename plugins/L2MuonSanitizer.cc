#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

using namespace std;
using namespace edm;
using namespace reco;

// It is possible that the L2MuonCandidateProducer may not be run in
// the HLT path, in which case no hltL2MuonCandidates collection is
// put in the event. If that is the case, put an empty one so we don't
// crash later down the analysis path, else copy the existing
// one. (The same is true for the L3 muons, but we are already doing
// something special for them in L3MuonSanitizer.)

class L2MuonSanitizer : public EDProducer {
public:
  explicit L2MuonSanitizer(const ParameterSet&);
  ~L2MuonSanitizer() {}

private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() {}

  InputTag src;
};

L2MuonSanitizer::L2MuonSanitizer(const ParameterSet& cfg) {
  src = cfg.getParameter<InputTag>("src");
  produces<RecoChargedCandidateCollection>();
}

void L2MuonSanitizer::produce(Event& event,
				  const EventSetup& eSetup) {
  // Try to get the L2 muon collection from the event.
  Handle<RecoChargedCandidateCollection> muons;
  bool ok = true;
  try {
    event.getByLabel(src, muons);
  } catch (...) {
    ok = false;
  }

  // Make the output collection.
  auto_ptr<RecoChargedCandidateCollection>
    cands(new RecoChargedCandidateCollection);

  if (ok && !muons.failedToGet())
    // Copy the existing one.
    *cands = *muons;
  else
    edm::LogWarning("L2MuonSanitizer")
      << "no collection " << src << " in event; producing empty collection";

  event.put(cands);
}

DEFINE_FWK_MODULE(L2MuonSanitizer);
