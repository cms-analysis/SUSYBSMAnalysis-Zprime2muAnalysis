#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

using namespace std;
using namespace edm;

class HLTLeptonsFromTriggerEvent : public EDProducer {
public:
  explicit HLTLeptonsFromTriggerEvent(const ParameterSet&);
  ~HLTLeptonsFromTriggerEvent() {}
  
private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() {}

  bool isWanted(const InputTag& tag) const;
  vector<InputTag> src;
};

HLTLeptonsFromTriggerEvent::HLTLeptonsFromTriggerEvent(const ParameterSet& cfg) {
  src = cfg.getParameter<vector<InputTag> >("src");
  produces<vector<reco::RecoChargedCandidate> >();
}

bool HLTLeptonsFromTriggerEvent::isWanted(const InputTag& tag) const {
  for (vector<InputTag>::const_iterator stag = src.begin(); stag != src.end(); ++stag)
    if (*stag == tag) return true;
  return false;
}

void HLTLeptonsFromTriggerEvent::produce(Event& event,
					 const EventSetup& eSetup) {
  Handle<trigger::TriggerEvent> trigEvent;
  event.getByLabel("hltTriggerSummaryAOD", trigEvent);

  auto_ptr<vector<reco::RecoChargedCandidate> > cands(new vector<reco::RecoChargedCandidate>);

  if (trigEvent.isValid()) {
    vector<pair<int, int> > inds;
    int ind_prev = 0;
    for (trigger::size_type iC = 0; iC < trigEvent->sizeCollections(); ++iC) {
      if (isWanted(trigEvent->collectionTag(iC)))
	inds.push_back(make_pair(ind_prev, trigEvent->collectionKey(iC) - 1));
      ind_prev = trigEvent->collectionKey(iC);
    }

    for (unsigned i = 0; i < inds.size(); ++i) {
      int ind_first = inds[i].first, ind_last = inds[i].second;
      if (ind_first >= 0 && ind_last >= ind_first) {
	const trigger::TriggerObjectCollection& TOC(trigEvent->getObjects());
	for (trigger::size_type iO = ind_first; iO <= ind_last; ++iO) {
	  const trigger::TriggerObject& TO(TOC.at(iO));

	  // TriggerObjects store the charge in their PDG ids, but
	  // they don't bother to have a charge() method or to fill
	  // the particle's charge below properly.  Also, for L2
	  // electrons (which do not have a charge), the id member of
	  // the TriggerObject is set to 0. Hardcode the id and the
	  // charge to that of an electron.
	  int id = TO.id();
	  reco::Particle::Charge q;
	  if (id == 0) {
	    id = 11;
	    q = -1;
	  }
	  else
	    q = -TO.id()/abs(TO.id());

	  reco::Particle::LorentzVector p4(TO.px(), TO.py(), TO.pz(), TO.energy());
	  reco::RecoChargedCandidate cand(q, p4);
	  cand.setPdgId(id);
	  cands->push_back(cand);
	}
      }
    }  
  }
  // Fail silently to prevent spamming messages. The user will know when
  // something is amiss from empty plots.
  //else
  //  edm::LogWarning("HLTLeptonsFromTriggerEvent")
  //    << "no collection " << src << " in event; producing empty collection";

  event.put(cands);
}

DEFINE_FWK_MODULE(HLTLeptonsFromTriggerEvent);
