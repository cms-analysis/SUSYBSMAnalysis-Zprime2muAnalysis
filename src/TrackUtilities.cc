#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"

const reco::Track* getMainTrack(const reco::CandidateBaseRef& cand) {
  const reco::RecoCandidate* rcand =
    toConcretePtr<reco::RecoCandidate>(cand);
  if (rcand != 0)
    return rcand->bestTrack();
  else
    return 0;
}

double pError(const reco::Track* tk) {
  // dumb identity:
  // p = q / qoverp
  // dp_dqoverp = - q/qoverp^2

  const double qoverp = tk->qoverp();
  const double q = tk->charge();
  
  double sig_p = tk->qoverpError() * q/qoverp/qoverp;
  if (sig_p < 0) sig_p = -sig_p;
  return sig_p;
}

double ptError(const reco::Track* tk) {
  return tk->ptError();
}

double invPtError(const reco::Track* tk) {
  return invError(tk->pt(), tk->ptError());
}
  
double invPError(const reco::Track* tk) {
  return invError(tk->p(), pError(tk));
}

double ptError(const reco::CandidateBaseRef& cand) {
  return getMainTrack(cand)->ptError();
}

double pError(const reco::CandidateBaseRef& cand) {
  return pError(getMainTrack(cand));
}

double invPtError(const reco::CandidateBaseRef& cand) {
  return invPtError(getMainTrack(cand));
}

double invPError(const reco::CandidateBaseRef& cand) {
  return invPError(getMainTrack(cand));
}

bool matchTracks(const reco::CandidateBaseRef& cand1,
		 const reco::CandidateBaseRef& cand2) {
  return fabs(deltaPhi(cand1->phi(), cand2->phi())) < .5 &&
    fabs(cand1->eta() - cand2->eta()) < .5;
}
