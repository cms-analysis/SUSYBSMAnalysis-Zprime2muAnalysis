#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/CandAlgos/interface/CandMatcher.h"
#include "PhysicsTools/HepMCCandAlgos/interface/MCTruthPairSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ParameterAdapter.h"
#include "PhysicsTools/Utilities/interface/AnyPairSelector.h"

// Used to filter gen/sim dilepton creation.

struct PtMinPairSelector {
   template<typename T1, typename T2>
   bool operator()(const T1 &t1, const T2 &t2) const {
     static const double ptMin = 1;
     return t1.pt() > ptMin && t2.pt() > ptMin;
   }
};

namespace reco {
  namespace modules {
    template<>
    struct ParameterAdapter<PtMinPairSelector> {
      static PtMinPairSelector make(const edm::ParameterSet& cfg) {
        return PtMinPairSelector();
      }
    };
  }
}

typedef reco::modules::CandMatcher<
  PtMinPairSelector,
  reco::CandidateView
  > PtMinDeltaRViewMatcher;

DEFINE_FWK_MODULE(PtMinDeltaRViewMatcher);

typedef reco::modules::CandMatcher<
  helpers::MCTruthPairSelector<reco::Candidate>,
  reco::CandidateView
  > MCTruthDeltaRViewMatcher;

DEFINE_FWK_MODULE(MCTruthDeltaRViewMatcher);
