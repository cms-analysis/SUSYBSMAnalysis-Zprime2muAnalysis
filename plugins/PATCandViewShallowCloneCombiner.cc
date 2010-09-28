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

struct LooseTightPairSelector {
  StringCutObjectSelector<reco::Candidate, true> loose;
  StringCutObjectSelector<reco::Candidate, true> tight;
  
  LooseTightPairSelector(const std::string& loose_cut, const std::string& tight_cut) : loose(loose_cut), tight(tight_cut) {}

  bool electron_ok(const reco::Candidate& cel) const {
    const pat::Electron& el = static_cast<const pat::Electron&>(cel);
    if (!el.hasUserInt("cutFor"))
      throw cms::Exception("BadCandidateForLeptonSelector") << "electron does not have cutFor userInt!\n";
    return el.userInt("cutFor") == 0;
  }

  bool operator()(const reco::Candidate& c1, const reco::Candidate& c2) const {
    // If one is an electron, fall back to requiring the muon pass
    // both the loose and tight cuts.
    const bool e1 = typeid(c1) == typeid(pat::Electron);
    const bool e2 = typeid(c2) == typeid(pat::Electron);
    if (e1 && e2)
      return electron_ok(c1) && electron_ok(c2);
    else if (e1)
      return electron_ok(c1) && loose(c2) && tight(c2);
    else if (e2)
      return electron_ok(c2) && loose(c1) && tight(c1);

    return loose(c1) && loose(c2) && (tight(c1) || tight(c2));
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

    template<>
    struct ParameterAdapter<LooseTightPairSelector> {
      static LooseTightPairSelector make(const edm::ParameterSet& cfg) {
	return LooseTightPairSelector(cfg.getParameter<std::string>("loose_cut"), cfg.getParameter<std::string>("tight_cut"));
      }
    };

    typedef CandCombiner<
      PATStringCutObjectSelector,
      LeptonPairSelector,
      combiner::helpers::ShallowClone,
      pat::CompositeCandidateCollection
      > PATCandViewShallowCloneCombiner;
  
    DEFINE_FWK_MODULE(PATCandViewShallowCloneCombiner);

    typedef CandCombiner<
      PATStringCutObjectSelector,
      AnyPairSelector,
      combiner::helpers::ShallowClone,
      pat::CompositeCandidateCollection
      > PATRawCandViewShallowCloneCombiner;
  
    DEFINE_FWK_MODULE(PATRawCandViewShallowCloneCombiner);

    typedef CandCombiner<
      PATStringCutObjectSelector,
      LooseTightPairSelector,
      combiner::helpers::ShallowClone,
      pat::CompositeCandidateCollection
      > LooseTightCandViewShallowCloneCombiner;

    DEFINE_FWK_MODULE(LooseTightCandViewShallowCloneCombiner);
  } 
}
