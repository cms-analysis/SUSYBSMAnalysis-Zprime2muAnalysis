#include "CommonTools/CandAlgos/interface/CandCombiner.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/PatUtils/interface/StringParserTools.h"

struct LeptonPairSelector {
  bool passes(const reco::Candidate& c) const {
    const std::type_info& type = typeid(c);
    if (type == typeid(pat::Muon)) {
      const pat::Muon& mu = static_cast<const pat::Muon&>(c);
      if (!mu.hasUserInt("cutFor"))
	throw cms::Exception("BadCandidateForLeptonSelector") << "lepton does not have cutFor userInt!\n";
      return mu.userInt("cutFor") == 0;
    }
    else if (type == typeid(pat::Electron)) {
      const pat::Electron& el = static_cast<const pat::Electron&>(c);
      if (!el.hasUserInt("cutFor"))
        throw cms::Exception("BadCandidateForLeptonSelector") << "lepton does not have cutFor userInt!\n";
      return el.userInt("cutFor") == 0;
    }
    else
      throw cms::Exception("BadCandidateForLeptonSelector") << "typeid name is " << type.name() << "\n";
  }

  bool operator()(const reco::Candidate& c1, const reco::Candidate& c2) const {
    return passes(c1) && passes(c2);
  }
};

namespace reco {
  namespace modules {
    template<>
    struct ParameterAdapter<PATStringCutObjectSelector> {
      static PATStringCutObjectSelector make(const edm::ParameterSet& cfg) {
	return PATStringCutObjectSelector(cfg.getParameter<std::string>("cut"));
      }
    };

    template <>
    struct ParameterAdapter<LeptonPairSelector> {
      static LeptonPairSelector make(const edm::ParameterSet& cfg) {
	return LeptonPairSelector(); // can imagine passing in cut masks from the cfg to be used in the cutFor part above
      }
    };

    typedef CandCombiner<
      //StringCutObjectSelector<pat::CompositeCandidate>,
      PATStringCutObjectSelector,

      //AnyPairSelector,
      LeptonPairSelector,

      combiner::helpers::ShallowClone,
      pat::CompositeCandidateCollection
      > PATCandViewShallowCloneCombiner;
  
    DEFINE_FWK_MODULE(PATCandViewShallowCloneCombiner);
  } 
}
