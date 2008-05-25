#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneratorUtilities.h"

void SetP4M(reco::Particle::LorentzVector& v,
	    double pt, double phi, double p, double theta, double m) {
  v.SetCoordinates(pt*cos(phi), pt*sin(phi), p*cos(theta), sqrt(p*p + m*m));
}

const reco::Candidate* mother(const reco::CandidateBaseRef& cand) {
  int pId = cand->pdgId();
  const reco::Candidate* mom = cand->mother();
  while (mom != 0 && mom->pdgId() == pId)
    mom = mom->mother();
  return mom;
}
  
int motherId(const reco::CandidateBaseRef& cand) {
  const reco::Candidate* mom = mother(cand);
  if (mom != 0) return mom->pdgId();
  else return 0;
  /*    throw cms::Exception("motherId")
        << "+++ cannot find mother for particle! pdgId: " << pId
        << " status: " << cand->status() << " +++\n";*/
}

int grandmotherId(const reco::CandidateBaseRef& cand) {
  const reco::Candidate *mom = mother(cand);
  if (mom == 0) return 0;
  const reco::Candidate *gm = mom->mother();
  if (gm != 0) return gm->pdgId();
  else return 0;
  /* throw cms::Exception("grandmotherId")
     << "+++ cannot find grandmother for particle! pdgId: " << p->pdgId()
     << " status: " << cand->status() 
     << " mother pdgId: " << mId << " status: " << mom->status() << " +++\n";*/
}
 
const reco::Candidate* sameMother(const reco::CandidateBaseRef& cand1,
				  const reco::CandidateBaseRef& cand2) {
  const reco::Candidate* m1 = mother(cand1);
  const reco::Candidate* m2 = mother(cand2);

  bool sameMother = m1 != 0 && m2 != 0 && m1 == m2;
  //m1->pdgId() == m2->pdgId() &&
  //m1->status() == m2->status() && m1->p4() == m2->p4();

  if (sameMother) return m1;
  else return 0;
}
