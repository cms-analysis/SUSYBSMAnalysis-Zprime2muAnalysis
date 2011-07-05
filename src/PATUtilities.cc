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
    if (name == std::string("tkonly"))
      return TkInner;
    for (size_t i = 0; i < nTrackTypes; ++i)
      if (track_names[i] == name)
	return TrackType(i);
    return nTrackTypes;
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
    case nTrackTypes: default: return reco::TrackRef();
    }
  }

  reco::TrackRef trackByName(const pat::Muon& mu, const std::string& name) {
    return trackByType(mu, trackNameToType(name));
  }

  reco::TrackRef userDataTrack(const pat::Muon& mu, const std::string& name) {
    if (mu.hasUserData(name))
      return reco::TrackRef(mu.userData<reco::TrackCollection>(name), 0);
    else
      return reco::TrackRef();
  }
  
  TrackType getPickedTrackType(const pat::Muon& mu) {
    if (!mu.hasUserInt("trackUsedForMomentum"))
      throw cms::Exception("Zp2muPATUtilities") << "getPickedTrackType called on muon without trackUsedForMomentum userInt!\n";
    return TrackType(mu.userInt("trackUsedForMomentum"));
  }
  
  reco::TrackRef getPickedTrack(const pat::Muon& mu) {
    return trackByType(mu, getPickedTrackType(mu));
  }
  
  bool wasCocktailUsed(const TrackType type) {
    return type == TkPMC || type == TkTMR || type == TkSigmaSwitch;
  }

  bool wasCocktailUsed(const pat::Muon& mu) {
    return wasCocktailUsed(getPickedTrackType(mu));
  }

  TrackType resolveCocktail(const pat::Muon& mu) {
    TrackType type = getPickedTrackType(mu);
    if (!wasCocktailUsed(type)) // could throw but leave up to caller to check
      return type;
    reco::TrackRef picked = getPickedTrack(mu);
    if (picked == mu.globalTrack())
      return TkGlobal;
    else if (picked == mu.innerTrack())
      return TkInner;
    else if (picked == mu.tpfmsMuon())
      return TkTPFMS;
    else if (picked == mu.pickyMuon())
      return TkPicky;
    else
      return nTrackTypes;
  }
}
