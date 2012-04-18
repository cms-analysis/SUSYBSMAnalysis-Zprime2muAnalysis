#include "CommonTools/CandAlgos/interface/CandCombiner.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/PatUtils/interface/StringParserTools.h"

struct Zprime2muPairSelector {
  StringCutObjectSelector<reco::Candidate, true> loose;
  StringCutObjectSelector<reco::Candidate, true> tight;
  const unsigned electron_cut_mask;
  const std::string module_label;

  Zprime2muPairSelector(const std::string& loose_cut,
			const std::string& tight_cut,
			const unsigned el_mask,
			const std::string& label)
    : loose(loose_cut),
      tight(tight_cut),
      electron_cut_mask(el_mask),
      module_label(label)
  {}

  bool electron_ok(const reco::Candidate& cel) const {
    const pat::Electron& el = static_cast<const pat::Electron&>(cel);
    if (!el.hasUserInt("cutFor"))
      throw cms::Exception("BadCandidateForLeptonSelector") << "electron does not have cutFor userInt!\n";
    return (el.userInt("cutFor") & electron_cut_mask) == 0;
  }

  bool operator()(const reco::Candidate& c1, const reco::Candidate& c2) const {
    // Candidates c1 and c2 can either both be muons, both electrons,
    // or one muon and one electron. For electrons, the cut codes from
    // HEEPSelector are checked, respecting the electron_cut_mask to
    // turn off cuts, so an electron passes if cutFor is 0. Muons must
    // pass the loose cuts, while at least one muon must pass the
    // tight selection.

    const bool e1 = typeid(c1) == typeid(pat::Electron);
    const bool e2 = typeid(c2) == typeid(pat::Electron);
    if (e1 && e2)
      return electron_ok(c1) && electron_ok(c2);
    else if (e1)
      return electron_ok(c1) && loose(c2) && tight(c2);
    else if (e2)
      return electron_ok(c2) && loose(c1) && tight(c1);

    //printf("in mumu Zprime2muPairSelector %s with c1 pt %f c2 pt %f  loose1 %i loose2 %i tight1 %i tight2 %i\n", module_label.c_str(), c1.pt(), c2.pt(), loose(c1), loose(c2), tight(c1), tight(c2));
    return loose(c1) && loose(c2) && (tight(c1) || tight(c2));
  }
};

namespace reco {
  namespace modules {
    template<>
    struct ParameterAdapter<StringCutObjectSelector<reco::Candidate, true> > {
      static StringCutObjectSelector<reco::Candidate, true> make(const edm::ParameterSet& cfg) {
	return StringCutObjectSelector<reco::Candidate, true>(cfg.getParameter<std::string>("cut"));
      }
    };

    template<>
    struct ParameterAdapter<Zprime2muPairSelector> {
      static Zprime2muPairSelector make(const edm::ParameterSet& cfg) {
	return Zprime2muPairSelector(cfg.getParameter<std::string>("loose_cut"),
				     cfg.getParameter<std::string>("tight_cut"),
				     cfg.existsAs<unsigned>("electron_cut_mask") ? cfg.getParameter<unsigned>("electron_cut_mask") : 0xFFFFFFFF,
				     cfg.getParameter<std::string>("@module_label"));
      }
    };

    typedef CandCombiner<
      StringCutObjectSelector<reco::Candidate, true>,
      Zprime2muPairSelector,
      combiner::helpers::ShallowClone,
      pat::CompositeCandidateCollection
      > Zprime2muCombiner;

    DEFINE_FWK_MODULE(Zprime2muCombiner);
  } 
}
