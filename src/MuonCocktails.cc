#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/MuonCocktails.h"

reco::TrackRef TMR(const reco::TrackRef& trackerTrack,
		   const reco::TrackRef& fmsTrack,
		   const double cut) {
  double probTK  = 0;
  double probFMS = 0;
  
  if (trackerTrack.isNonnull() && trackerTrack->numberOfValidHits())
    probTK = muon::trackProbability(trackerTrack);
  if (fmsTrack.isNonnull() && fmsTrack->numberOfValidHits())
    probFMS = muon::trackProbability(fmsTrack);
  
  bool TKok  = probTK > 0;
  bool FMSok = probFMS > 0;

  if (TKok && FMSok) {
    if (probFMS - probTK > cut)
      return trackerTrack;
    else
      return fmsTrack;
  }
  else if (FMSok)
    return fmsTrack;
  else if (TKok)
    return trackerTrack;
  else
    return reco::TrackRef();
}

reco::TrackRef TMR(const reco::Muon& muon,
		   const reco::TrackToTrackMap& fmsMap,
		   const double cut) {
  reco::TrackToTrackMap::const_iterator fmsTrack = fmsMap.find(muon.globalTrack());
  return TMR(muon.innerTrack(), fmsTrack != fmsMap.end() ? fmsTrack->val : reco::TrackRef(), cut);
}

reco::TrackRef sigmaSwitch(const reco::TrackRef& combinedTrack,
			   const reco::TrackRef& trackerTrack,
			   const double nSigma,
			   const double ptThreshold) {
  if (combinedTrack->pt() < ptThreshold || trackerTrack->pt() < ptThreshold)
    return trackerTrack;

  double delta = fabs(trackerTrack->qoverp() - combinedTrack->qoverp());
  double threshold = nSigma * trackerTrack->qoverpError();

  return delta > threshold ? trackerTrack : combinedTrack;
}

reco::TrackRef sigmaSwitch(const reco::Muon& muon,
			   const double nSigma,
			   const double ptThreshold) {
  return sigmaSwitch(muon.globalTrack(), muon.innerTrack(), nSigma, ptThreshold);
}

reco::TrackRef tevOptimized(const reco::TrackRef& trackerTrack,
			    const reco::TrackRef& gmrTrack,
			    const reco::TrackRef& fmsTrack,
			    const reco::TrackRef& pmrTrack) {
  const reco::TrackRef refit[4] = { 
    trackerTrack, 
    gmrTrack, 
    fmsTrack, 
    pmrTrack 
  }; 

  double prob[4] = {0.}; 
  int chosen = 3; 
 
  for (unsigned int i = 0; i < 4; ++i) 
    if (refit[i].isNonnull() && refit[i]->numberOfValidHits()) 
      prob[i] = muon::trackProbability(refit[i]); 
 
  if (prob[3] == 0.) { 
    if (prob[2] > 0.) chosen=2; 
    else if (prob[1] > 0.) chosen=1; 
    else if (prob[0] > 0.) chosen=0; 
  } 
 
  if (prob[0] > 0. && prob[3] > 0. && prob[3] - prob[0] > 30.) chosen = 0; 
  if (prob[2] > 0. && prob[chosen] - prob[2] > 0.) chosen = 2;
 
  return refit[chosen]; 
} 

reco::TrackRef tevOptimized(const reco::TrackRef& combinedTrack,
			    const reco::TrackRef& trackerTrack,
			    const reco::TrackToTrackMap tevMap1,
			    const reco::TrackToTrackMap tevMap2,
			    const reco::TrackToTrackMap tevMap3) {
  reco::TrackToTrackMap::const_iterator gmrTrack = tevMap1.find(combinedTrack);
  reco::TrackToTrackMap::const_iterator fmsTrack = tevMap2.find(combinedTrack);
  reco::TrackToTrackMap::const_iterator pmrTrack = tevMap3.find(combinedTrack);
  
  reco::TrackRef gmr = gmrTrack != tevMap1.end() ? gmrTrack->val : reco::TrackRef();
  reco::TrackRef fms = fmsTrack != tevMap2.end() ? fmsTrack->val : reco::TrackRef();
  reco::TrackRef pmr = pmrTrack != tevMap3.end() ? pmrTrack->val : reco::TrackRef();
  
  return tevOptimized(trackerTrack, gmr, fms, pmr);
}

reco::TrackRef tevOptimized(const reco::Muon& muon,
			    const reco::TrackToTrackMap tevMap1,
			    const reco::TrackToTrackMap tevMap2,
			    const reco::TrackToTrackMap tevMap3) {
  return muon::tevOptimized(muon.globalTrack(), muon.innerTrack(), tevMap1, tevMap2, tevMap3);
}
