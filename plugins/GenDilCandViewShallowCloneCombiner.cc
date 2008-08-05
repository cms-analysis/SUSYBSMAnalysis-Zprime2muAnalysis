#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "PhysicsTools/CandAlgos/interface/CandCombiner.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "PhysicsTools/UtilAlgos/interface/ParameterAdapter.h"
 
struct GenDilPairSelector {
  template<typename T1, typename T2>
  bool operator()(const T1 &t1, const T2 &t2) const {
    // Make sure the leptons come from the same resonance.
    const reco::Candidate* mom1 = t1.mother();
    const reco::Candidate* mom2 = t2.mother();

    // Go up the brem chain.
    while (mom1 && mom1->pdgId() == t1.pdgId()) mom1 = mom1->mother();
    while (mom2 && mom2->pdgId() == t2.pdgId()) mom2 = mom2->mother();

    if (mom1 == 0 || mom2 == 0) return false;

    return mom1 == mom2;
  }
};

namespace reco {
  namespace modules {
    template<>
    struct ParameterAdapter<GenDilPairSelector> {
      static GenDilPairSelector make(const edm::ParameterSet& cfg) {
	return GenDilPairSelector();
      }
    };

    // CandCombiner is in a funny state -- some interface changes
    // (specificially the dropped template parameters below) went into
    // 2_0_9, but not into any 2_1_X version. For now, comment out
    // until things settle.
    typedef CandCombiner<
      //      reco::CandidateView,
      StringCutObjectSelector<reco::Candidate>,
      //      reco::CompositeCandidateCollection,
      GenDilPairSelector,
      combiner::helpers::ShallowClone
      > GenDilCandViewShallowCloneCombiner;
  
    DEFINE_FWK_MODULE(GenDilCandViewShallowCloneCombiner);
  } 
}
