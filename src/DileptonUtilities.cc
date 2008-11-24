#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneratorUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TeVMuHelper.h"

using namespace std;

// Sorting functor for removeDileptonOverlap().
struct reverse_mass_sort {
  bool operator()(const reco::Candidate& lhs, const reco::Candidate& rhs) {
    return lhs.mass() > rhs.mass();
  }
  bool operator()(const reco::CompositeCandidate& lhs, const reco::CompositeCandidate& rhs) {
    return lhs.mass() > rhs.mass();
  }
};

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

void genDileptonsOnly(const reco::CompositeCandidateCollection& dils,
		      reco::CompositeCandidateCollection& newDils,
		      const bool debug) {
  ostringstream out;
  if (debug)
    out << "cutDileptons: starting size " << dils.size() << endl;

  reco::CompositeCandidateCollection::const_iterator dil;
  for (dil = dils.begin(); dil != dils.end(); dil++) {
    // Make sure we pick up the actual dilepton from the resonance
    // -- i.e. the one comprised of leptons with the same mother
    // that has a PDG id of one of the resonances we care about.
    const reco::Candidate* mom  = sameMother(dileptonDaughter(*dil, 0),
					     dileptonDaughter(*dil, 1));

    if (debug) {
      out << "  filtering at generator level: mom pointer: " << mom;
      if (mom) out << " mom's id: " << mom->pdgId();
      out << endl;
    }

    if (mom && HardInteraction::IsResonance(mom->pdgId())) {
      if (debug) out << "  dilepton not cut!" << endl;
      newDils.push_back(*dil);
    }
  }

  if (debug) {
    out << " done. end size: " << newDils.size();
    edm::LogInfo("cutDileptons") << out.str();
  }
}

void cutDileptons(const edm::Event& event,
		  const reco::CompositeCandidateCollection& dils,
		  reco::CompositeCandidateCollection& newDils,
		  unsigned cuts, int pdgId1, int pdgId2,
		  const bool debug) {
  ostringstream out;
  if (debug)
    out << "cutDileptons: starting size " << dils.size() << endl;

  TeVMuHelper tmh;
  tmh.initEvent(event);

  reco::CompositeCandidateCollection::const_iterator dil;
  for (dil = dils.begin(); dil != dils.end(); dil++) {
    // If requested (by both PDG ids being nonzero), cut out the
    // dileptons that are not made up of the requested
    // leptons. (Useful for separating mu+mu+ from mu-mu- in the
    // output of CandCombiner.)
    if (pdgId1 != 0 && pdgId2 != 0) {
      int p1 = dil->daughter(0)->pdgId();
      int p2 = dil->daughter(1)->pdgId();
      if (debug) out << "  filtering by PDG ids: requested: " << pdgId1 << " "
		     << pdgId2 << " got: " << p1 << " " << p2 << endl;
      if (!((p1 == pdgId1 || p1 == pdgId2) && (p2 == pdgId1 || p2 == pdgId2)))
	continue;
    }

    // If the lepton is cut with a code that is present in the cuts
    // bitmask, cut this dilepton.
    if (!tmh.dileptonIsCut(*dil, cuts)) {
      if (debug) out << "  dilepton not cut!" << endl;
      newDils.push_back(*dil);
    }
  }

  if (debug) {
    out << " done. end size: " << newDils.size();
    edm::LogInfo("cutDileptons") << out.str();
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
