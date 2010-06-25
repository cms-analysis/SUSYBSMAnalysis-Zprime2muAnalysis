#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class Zprime2muCompositeCandidatePicker : public edm::EDProducer {
public:
  explicit Zprime2muCompositeCandidatePicker(const edm::ParameterSet&);
  
private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  struct reverse_mass_sort {
    bool operator()(const pat::CompositeCandidate& lhs, const pat::CompositeCandidate& rhs) {
      return lhs.mass() > rhs.mass();
    }
  };

  void remove_overlap(pat::CompositeCandidateCollection&) const;

  edm::InputTag src;
  StringCutObjectSelector<pat::CompositeCandidate> selector;
  unsigned max_candidates;
};

Zprime2muCompositeCandidatePicker::Zprime2muCompositeCandidatePicker(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
    selector(cfg.getParameter<std::string>("cut")),
    max_candidates(cfg.getParameter<unsigned>("max_candidates"))
{
  produces<pat::CompositeCandidateCollection>();
}

void Zprime2muCompositeCandidatePicker::remove_overlap(pat::CompositeCandidateCollection& cands) const {
  // Don't bother doing anything if there's just one candidate.
  if (cands.size() < 2) return; 

  // Sort candidates so we keep the ones with larger invariant
  // mass. Could make configurable to choose other sorting.
  sort(cands.begin(), cands.end(), reverse_mass_sort());
                                                 
  pat::CompositeCandidateCollection::iterator p, q;
  for (p = cands.begin(); p != cands.end() - 1; ) {
    for (q = p + 1; q != cands.end(); ++q) {         
      // Check to see if any of the leptons in p is in q also. If so,
      // remove q (e.g. the one with lower invariant mass since we
      // have sorted the vector already), reset pointers and restart.

      // To do this we need the unique ids of the daughters, i.e. the
      // refs into the original lepton collections.
      typedef std::vector<reco::CandidateBaseRef> refs;
      refs prefs, qrefs;
      for (size_t i = 0; i < p->numberOfDaughters(); ++i)
	prefs.push_back(p->daughter(i)->masterClone());
      for (size_t i = 0; i < q->numberOfDaughters(); ++i)
	qrefs.push_back(q->daughter(i)->masterClone());

      // Compare every pair of (pref, qref) to check for any lepton
      // being shared.
      bool any_shared = false;
      refs::const_iterator pr = prefs.begin(), pre = prefs.end(),
	qr = qrefs.begin(), qre = qrefs.end();
      for ( ; pr != pre && !any_shared; ++pr)
	for ( ; qr != qre && !any_shared; ++qr)
	  if (pr == qr)
	    any_shared = true;

      if (any_shared) {
	cands.erase(q);
	p = cands.begin();
      }
      else
	++p;
    }
  }
}

void Zprime2muCompositeCandidatePicker::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<pat::CompositeCandidateCollection> cands;
  event.getByLabel(src, cands);

  std::auto_ptr<pat::CompositeCandidateCollection> new_cands(new pat::CompositeCandidateCollection);

  // Copy all the candidates that pass the specified cuts into the new
  // output vector.
  for (pat::CompositeCandidateCollection::const_iterator c = cands->begin(), ce = cands->end(); c != ce; ++c)
    if (selector(*c))
      new_cands->push_back(*c);
  
  // Remove cands of lower invariant mass that are comprised of a
  // lepton that has been used by a higher invariant mass one.
  remove_overlap(*new_cands);
  
  // Only return the maximum number of candidates specified.
  if (new_cands->size() > max_candidates)
    new_cands->erase(new_cands->begin() + max_candidates, new_cands->end());
  
  event.put(new_cands);
}

DEFINE_FWK_MODULE(Zprime2muCompositeCandidatePicker);
