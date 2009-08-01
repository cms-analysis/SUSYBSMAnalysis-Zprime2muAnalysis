#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/CandAlgos/interface/CandCombiner.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

struct GenDilPairSelector {
  template<typename T1, typename T2>
  bool operator()(const T1 &t1, const T2 &t2) const {
    // Make sure the leptons come from the same resonance.
    const reco::Candidate* mom1 = t1.mother();
    const reco::Candidate* mom2 = t2.mother();

    // Go up the brem chain.
    while (mom1 && mom1->pdgId() == t1.pdgId()) mom1 = mom1->mother();
    while (mom2 && mom2->pdgId() == t2.pdgId()) mom2 = mom2->mother();

    if (mom1 == 0 || mom2 == 0 || mom1 != mom2) return false;

    int pdgId = mom1->pdgId();
    return pdgId == 32 || pdgId == 23 || pdgId == 39 || pdgId == 5000039;
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

    typedef CandCombiner<
      StringCutObjectSelector<reco::Candidate>,
      GenDilPairSelector,
      combiner::helpers::ShallowClone
      > GenDilCandViewShallowCloneCombiner;
  
    DEFINE_FWK_MODULE(GenDilCandViewShallowCloneCombiner);
  } 
}
