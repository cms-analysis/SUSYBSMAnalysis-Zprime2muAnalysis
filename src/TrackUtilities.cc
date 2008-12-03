#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"

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

bool matchTracks(const reco::CandidateBaseRef& cand1,
		 const reco::CandidateBaseRef& cand2) {
  return fabs(deltaPhi(cand1->phi(), cand2->phi())) < .5 &&
    fabs(cand1->eta() - cand2->eta()) < .5;
}
