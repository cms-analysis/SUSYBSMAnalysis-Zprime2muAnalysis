#ifndef Zp2mu_MuonCocktails_h
#define Zp2mu_MuonCocktails_h

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

reco::TrackRef TMR(const reco::TrackRef& trackerTrack,
		   const reco::TrackRef& fmsTrack,
		   const double cut=4.);

reco::TrackRef TMR(const reco::Muon& muon,
		   const reco::TrackToTrackMap& fmsMap,
		   const double cut=4.);

reco::TrackRef sigmaSwitch(const reco::TrackRef& combinedTrack,
			   const reco::TrackRef& trackerTrack,
			   const double nSigma=2.,
			   const double ptThreshold=200.);

reco::TrackRef sigmaSwitch(const reco::Muon& muon,
			   const double nSigma=2.,
			   const double ptThreshold=200.);

reco::TrackRef tevOptimized(const reco::TrackRef& trackerTrack,
			    const reco::TrackRef& gmrTrack,
			    const reco::TrackRef& fmsTrack,
			    const reco::TrackRef& pmrTrack);

reco::TrackRef tevOptimized(const reco::TrackRef& combinedTrack,
			    const reco::TrackRef& trackerTrack,
			    const reco::TrackToTrackMap tevMap1,
			    const reco::TrackToTrackMap tevMap2,
			    const reco::TrackToTrackMap tevMap3);

reco::TrackRef tevOptimized(const reco::Muon& muon,
			    const reco::TrackToTrackMap tevMap1,
			    const reco::TrackToTrackMap tevMap2,
			    const reco::TrackToTrackMap tevMap3);

#endif
