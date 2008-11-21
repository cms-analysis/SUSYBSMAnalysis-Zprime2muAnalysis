#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

using namespace std;
using namespace reco;

HardInteraction::HardInteraction(int lepFlavor, //vector<int> resIds,
				 bool allowFakeRes)
  : leptonFlavor(lepFlavor), //resonanceIds(resIds),
    allowFakeResonance(allowFakeRes)
{
  Clear();
  // Don't fill the structure here, since we want to be able to throw
  // on error.
}

HardInteraction::~HardInteraction() {
  // If we faked the resonance by adding the l+ and l-, we own its
  // pointer, so delete it. Caution: we own none of the other
  // pointers, since they come directly from the genParticles
  // collection!
  if (resonanceIsFake)
    delete resonance;
}

/*
void SetResonanceIds(const std::vector<int>& resIds) {
  resonanceIds = resIds;
}
*/

bool HardInteraction::IsResonance(int pdgId) {
  /*
  vector<int>::const_iterator p = resonanceIds.begin();
  for (; p != resonanceIds.end(); p++)
    if (*p == pdgId)
      return true;
  return false;
  */

  // Z', Z0/gamma*, G, or G*
  return pdgId == 32 || pdgId == 23 || pdgId == 39 || pdgId == 5000039;
}

void HardInteraction::Clear() {
  quark = resonance = lepPlus = lepMinus = lepPlusNoIB = lepMinusNoIB = 0;
  bremPhotons.clear();
  resonanceIsFake = false;
}

void HardInteraction::Fill(const edm::Event& event) {
  edm::Handle<GenParticleCollection> genParticles;
  event.getByLabel("genParticles", genParticles);
  Fill(*genParticles);
}

void HardInteraction::Fill(const GenParticleCollection& genParticles) {
  // Reset everything before filling.
  Clear();

  // Look in the doc lines for the hard-interaction resonance and
  // leptons.
  GenParticleCollection::const_iterator genp = genParticles.begin();
  for (; genp != genParticles.end(); genp++) {
    if (genp->status() == 3) {
      int pdgId = genp->pdgId();

      if (IsResonance(pdgId)) {
	// We found the resonance (Z0/Z'/etc.). Make sure we didn't
	// find a second one.
	if (resonance != 0)
	  throw cms::Exception("HardInteraction")
	    << "Found second resonance (pdgId: " << pdgId << ") in event!\n";

	resonance = &*genp;
      }
      else if (pdgId == leptonFlavor) {
	// We found the l-. Make sure we didn't find a second one.
	if (lepMinus != 0)
	  throw cms::Exception("HardInteraction")
	    << "Found second l- in event!\n";

	lepMinusNoIB = &*genp;
      }
      else if (pdgId == -leptonFlavor) {
	// We found the l+. Make sure we didn't find a second one.
	if (lepPlus != 0)
	  throw cms::Exception("HardInteraction")
	    << "Found second l+ in event!\n";

	lepPlusNoIB = &*genp;
      }
    }
  }

  // We should always find the l+ and l-.
  if (lepMinusNoIB == 0 || lepPlusNoIB == 0)
    throw cms::Exception("HardInteraction")
      << "Couldn't find at least one of the pre-brem l-l+! l- = "
      << lepMinusNoIB << " l+ = " << lepPlusNoIB << endl;

  // We get the quark/antiquark from either the resonance or directly
  // from one of the leptons (see below). Start by assuming it is from
  // the resonance (we'll check for the null pointer later if this is
  // not true).
  const Candidate* mothersAreQuarks = resonance;

  // Some COMPHEP samples do not put the resonance back into the HepMC
  // record. So, if we didn't find a resonance, and are allowed to,
  // build one (we really only care about its four vector and charge
  // at this point). Use the doc-line l+l- we found to do so.
  if (resonance == 0) {
    if (allowFakeResonance) {
      Particle::Charge q = lepPlusNoIB->charge() + lepMinusNoIB->charge();
      Particle::LorentzVector p4 = lepPlusNoIB->p4() + lepMinusNoIB->p4();
      // Don't know the creation vertex of the resonance. Set to (0,0,0)
      // by default.
      Particle::Point vtx;
      // We also don't know what the pdgId of the resonance was. Pass
      // one in?
      int pdgId = 99;      
      int status = 3;
    
      // We have ownership of this pointer, and we will delete it in our
      // destructor.
      resonance = new GenParticle(q, p4, vtx, pdgId, status, true);
      resonanceIsFake = true;

      // In these resonanceless COMPHEP samples, the quark/antiquark are
      // direct mothers of the leptons. So now, assume we can get them
      // from the l-.
      mothersAreQuarks = lepMinusNoIB;
    }
    else
      throw cms::Exception("HardInteraction")
	<< "Did not find the resonance in the event, and forbidden"
	<< " from faking it!\n";
  }

  // Did we successfully identify a Candidate that has the
  // quark/antiquark as its mothers? Try to get them.
  if (mothersAreQuarks != 0 && mothersAreQuarks->numberOfMothers() == 2) {
    for (unsigned m = 0; m < 2; m++) {
      const Candidate* mom = mothersAreQuarks->mother(m);
      int momId = mom->pdgId();

      // For now, don't count gluons (id = 21), but we should
      // implement them later here and in all the above (for gg->G*).
      if (abs(momId) > 6 || momId == 0)
	throw cms::Exception("HardInteraction")
	  << "Mother " << m << " is not a quark! its pdgId = "
	  << momId << endl;
      
      if (momId > 0)
	quark = mom;
      else
	antiquark = mom;
    }
  }
  
  // If we didn't find the quark or antiquark, bomb.
  if (quark == 0 || antiquark == 0)
    throw cms::Exception("HardInteraction")
      << "Couldn't find at least one of the quark/antiquark! quark = "
      << quark << " antiquark = " << antiquark << endl;

  // Now, pick up the leptons after brem. They are daughter leptons of
  // the doc-line leptons and have status = 1, i.e. they are
  // final-state. Also grab the brem photons.
  for (Candidate::const_iterator dau = lepMinusNoIB->begin();
       dau != lepMinusNoIB->end(); dau++) {
    if (dau->status() == 1) {
      int pdgId = dau->pdgId();
      if (pdgId == leptonFlavor)
	lepMinus = &*dau;
      else if (pdgId == 22)
	bremPhotons.push_back(&*dau);
    }
  }

  for (Candidate::const_iterator dau = lepPlusNoIB->begin();
       dau != lepPlusNoIB->end(); dau++) {
    if (dau->status() == 1) {
      int pdgId = dau->pdgId();
      if (pdgId == -leptonFlavor)
	lepPlus = &*dau;
      else if (pdgId == 22)
	bremPhotons.push_back(&*dau);
    }
  }

  // If we didn't find the final-state l- or l+, bomb.
  if (lepMinus == 0 || lepPlus == 0)
    throw cms::Exception("HardInteraction")
      << "Couldn't find at least one of the final-state l-l+! l- = "
      << lepMinus << " l+ = " << lepPlus << endl;
}
