#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"

using namespace std;

void removeDileptonOverlap(reco::CompositeCandidateCollection& dileptons,
			   const bool debug) {
  if (dileptons.size() < 2) return;

  // Sort dileptons so we keep the ones with larger invariant mass.
  sort(dileptons.begin(), dileptons.end(), reverse_mass_sort());

  ostringstream out;
  if (debug)
    out << "removeDileptonOverlap: starting size " << dileptons.size() << endl;

  reco::CompositeCandidateCollection::iterator pdi, qdi;
  for (pdi = dileptons.begin(); pdi != dileptons.end() - 1;) {
    for (qdi = pdi+1; qdi != dileptons.end(); qdi++) {
      // Get the unique ids of the dilepton daughters, i.e. the
      // indices into the original lepton collections.
      int p0 = dileptonDaughter(*pdi, 0).key();
      int p1 = dileptonDaughter(*pdi, 1).key();
      int q0 = dileptonDaughter(*qdi, 0).key();
      int q1 = dileptonDaughter(*qdi, 1).key();

      if (debug) out << "  this pair: " << p0 << " " << p1 << " "
		     << q0 << " " << q1 << " " << endl;
      if (p0 == q0 || p0 == q1 || p1 == q1 || p1 == q0) {
	// If either lepton is shared, remove the second dilepton
	// (i.e. the one with lower invariant mass since we have
	// sorted the vector already), reset pointers and restart.
	dileptons.erase(qdi);
	pdi = dileptons.begin();
	break;
      }
      else {
	pdi++;
      }
    }
  }

  if (debug) {
    out << " done. end size: " << dileptons.size();
    edm::LogInfo("removeDileptonOverlap") << out.str();
  }
}

int numDaughtersInAcc(const reco::CompositeCandidate& dil,
		      const double etaCut) {
  int count = 0;
  for (unsigned ilep = 0; ilep < dil.numberOfDaughters(); ilep++)
    if (fabs(dil.daughter(ilep)->eta()) < etaCut)
      count++;
  return count;
}

const reco::CandidateBaseRef
dileptonDaughter(const reco::CompositeCandidate& dil,
		 const unsigned i) {
  if (i < 0 || i >= dil.numberOfDaughters())
    return reco::CandidateBaseRef(); // an invalid reference
  return dil.daughter(i)->masterClone();
}

const reco::CandidateBaseRef
dileptonDaughterByCharge(const reco::CompositeCandidate& dil,
			 const int charge) {
  for (unsigned ilep = 0; ilep < dil.numberOfDaughters(); ilep++)
    if (dil.daughter(ilep)->charge() == charge)
      return dil.daughter(ilep)->masterClone();
  return reco::CandidateBaseRef();
}

reco::Particle::LorentzVector resonanceP4(const pat::CompositeCandidate& cand) {
  // Start with the p4 of the candidate, which is the sum of the
  // daughter p4s by construction.
  reco::Particle::LorentzVector p4 = cand.p4();

  std::vector<int> used_already;

  for (pat::CompositeCandidate::const_iterator dau = cand.begin(), daue = cand.end(); dau != daue; ++dau) {
    // Only muons can have photons matched to them. Also, remember
    // that the CompositeCandidates were formed using shallow clones,
    // so need to get the daughters' masterClones to be able to get
    // the original pat::Muon*.
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>(&*dau->masterClone());
    if (mu == 0)
      continue;

    if (mu->hasUserInt("photon_index")) {
      int ndx = mu->userInt("photon_index");

      // Don't double count, in the case where two daughter muons
      // ended up matching to the same photon.
      if (std::find(used_already.begin(), used_already.end(), ndx) == used_already.end()) {
	used_already.push_back(ndx);
	const reco::Particle::LorentzVector* pho_p4 = mu->userData<reco::Particle::LorentzVector>("photon_p4");
	if (pho_p4 == 0)
	  throw cms::Exception("BadlyFormedPhotonUserData") << "candidate daughter had userInt for photon_index but no p4!\n";
	p4 += *pho_p4;
      }
    }
  }

  return p4;
}
