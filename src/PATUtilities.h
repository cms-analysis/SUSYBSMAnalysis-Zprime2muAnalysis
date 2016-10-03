#ifndef Zp2mu_PATUtilities_h
#define Zp2mu_PATUtilities_h

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

namespace patmuon {
  enum TrackType { TkGlobal, TkInner, TkOuter, TkTPFMS, TkPicky, TkDYT, TkTuneP, TkTMR, TkTunePNew, TkStartup, nTrackTypes };
  const std::string track_names[nTrackTypes+1] = { "global", "inner", "outer", "tpfms", "picky","dyt", "tunep", "tmr", "tunepnew", "startup", "invaild" };
//startup: tunep in barrel and tracker-only in endcaps
    
  TrackType trackNameToType(std::string name);
  reco::TrackRef trackByType(const pat::Muon& mu, TrackType t);
  reco::TrackRef trackByName(const pat::Muon& mu, const std::string& name);
  reco::TrackRef userDataTrack(const pat::Muon& mu, const std::string& name);
  TrackType getPickedTrackType(const pat::Muon& mu);
  reco::TrackRef getPickedTrack(const pat::Muon& mu);
  bool wasCocktailUsed(const TrackType type);
  bool wasCocktailUsed(const pat::Muon& mu);
  TrackType whichTrack(const pat::Muon& mu, const reco::TrackRef& tk);
  TrackType resolveCocktail(const pat::Muon& mu);
}

#endif
