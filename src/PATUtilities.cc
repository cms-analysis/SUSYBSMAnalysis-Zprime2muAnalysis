#include <algorithm>

#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"

namespace patmuon {
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
    bool hasTeVMuons = true;
    if (mu.hasUserInt("hasTeVMuons")) hasTeVMuons = mu.userInt("hasTeVMuons");
    switch (t) {
    case TkGlobal: return mu.globalTrack();
    case TkInner:{
//        std::cout<<"case TkInner: muon eta "<<mu.eta()<<" muon pt "<<mu.pt()<<std::endl;
        return mu.innerTrack();
    }
    case TkOuter: return mu.outerTrack();
    case TkTPFMS:{
	if (hasTeVMuons) return mu.tpfmsMuon();
	else return mu.tunePMuonBestTrack();
   }
   case TkPicky:{
	if (hasTeVMuons) return mu.pickyMuon();
	else mu.tunePMuonBestTrack();
    }
    case TkDYT:{
	   if (hasTeVMuons) return mu.dytTrack();
	   else mu.tunePMuonBestTrack();
    }
    case TkTuneP: {
	if (hasTeVMuons) muon::tevOptimized(mu, 200, 4, 6, -1).first;
	else return mu.tunePMuonBestTrack();
    }
    case TkTMR: return muon::TMR(mu.innerTrack(), mu.tpfmsMuon()).first;
    case TkTunePNew: {
//        std::cout<<"case TkTunePNew: muon eta "<<mu.eta()<<" muon pt "<<mu.pt()<<std::endl;
         if (hasTeVMuons) return muon::tevOptimized(mu, 200, 17, 40, 0.25).first;
         else return mu.tunePMuonBestTrack();
        
    }
    case TkStartup:{
        if(fabs(mu.eta())<1.2){
//            std::cout<<"case TkStartup: muon eta "<<mu.eta()<<" muon pt "<<mu.pt()<<std::endl;
//            std::cout<<" tunep "<<std::endl;
                if (hasTeVMuons) return muon::tevOptimized(mu, 200, 17, 40, 0.25).first;
         	else return mu.tunePMuonBestTrack();

	 }
        else {
//            std::cout<<"case TkStartup: muon eta "<<mu.eta()<<" muon pt "<<mu.pt()<<std::endl;
//            std::cout<<" tracker "<<std::endl;
            return mu.innerTrack();
        }
    }
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
    return type == TkTuneP || type == TkTMR || type == TkTunePNew;
  }

  bool wasCocktailUsed(const pat::Muon& mu) {
    return wasCocktailUsed(getPickedTrackType(mu));
  }

  TrackType whichTrack(const pat::Muon& mu, const reco::TrackRef& tk) {
    if (tk == mu.globalTrack())
      return TkGlobal;
    else if (tk == mu.innerTrack())
      return TkInner;
    else if (tk == mu.outerTrack())
      return TkOuter;
    else if (tk == mu.tpfmsMuon())
      return TkTPFMS;
    else if (tk == mu.dytTrack())
      return TkDYT;
    else if (tk == mu.pickyMuon())
      return TkPicky;
    else
      return nTrackTypes;
  }

  TrackType resolveCocktail(const pat::Muon& mu) {
    TrackType type = getPickedTrackType(mu);
    if (!wasCocktailUsed(type)) // could throw but leave up to caller to check
      return type;
    return whichTrack(mu, getPickedTrack(mu));
  }
}
