#include <algorithm>

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/MuonCocktails.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"

namespace patmuon {
  reco::TrackRef pmcTrack(const pat::Muon& mu) {
    return tevOptimized(mu.innerTrack(), mu.globalTrack(), mu.tpfmsMuon(), mu.pickyMuon());
  }

  reco::TrackRef tmrTrack(const pat::Muon& mu) {
    return TMR(mu.innerTrack(), mu.tpfmsMuon());
  }

  reco::TrackRef sigmaSwitchTrack(const pat::Muon& mu) {
    return sigmaSwitch(mu.globalTrack(), mu.innerTrack());
  }

  TrackType trackNameToType(std::string name) {
    std::transform(name.begin(), name.end(), name.begin(), tolower);
    static const std::string names[nTrackTypes] = { "global", "inner", "outer", "tpfms", "picky", "pmc", "tmr", "sigmaswitch" };
    for (size_t i = 0; i < nTrackTypes; ++i)
      if (names[i] == name)
	return TrackType(i);
    throw cms::Exception("BadNameForTrackType") << "invalid name: " << name << "\n";
  }

  reco::TrackRef trackByType(const pat::Muon& mu, TrackType t) {
    switch (t) {
    case TkGlobal: return mu.globalTrack();
    case TkInner: return mu.innerTrack();
    case TkOuter: return mu.outerTrack();
    case TkTPFMS: return mu.tpfmsMuon();
    case TkPicky: return mu.pickyMuon();
    case TkPMC: return pmcTrack(mu);
    case TkTMR: return tmrTrack(mu);
    case TkSigmaSwitch: return sigmaSwitchTrack(mu);
    case nTrackTypes: default: throw cms::Exception("nTrackTypesNotATrackType");
    }
  }

  reco::TrackRef trackByName(const pat::Muon& mu, const std::string& name) {
    return trackByType(mu, trackNameToType(name));
  }

  reco::TrackRef getPickedTrack(const pat::Muon& mu) {
    return trackByType(mu, TrackType(mu.userInt("trackUsedForMomentum")));
  }
}
