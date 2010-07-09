#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

HardInteraction::HardInteraction(const edm::ParameterSet cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
    doingElectrons(cfg.getParameter<bool>("doingElectrons")),
    leptonFlavor(doingElectrons ? 11 : 13),
    leptonMass(doingElectrons ? 0.000511 : 0.10566),
    allowFakeResonance(cfg.getParameter<bool>("allowFakeResonance")),
    resonanceIds(cfg.getParameter<std::vector<int> >("resonanceIds"))
{
  Clear();
}

HardInteraction::HardInteraction(bool doingElec, bool allowFakeRes)
  : src(edm::InputTag("genParticles")),
    doingElectrons(doingElec),
    leptonFlavor(doingElectrons ? 11 : 13),
    leptonMass(doingElectrons ? 0.000511 : 0.10566),
    allowFakeResonance(allowFakeRes)
{
  Clear();
}

HardInteraction::~HardInteraction() {
  // If we faked the resonance by adding the l+ and l-, we own its
  // pointer, so delete it. Caution: we own none of the other
  // pointers, since they come directly from the genParticles
  // collection!
  if (resonanceIsFake)
    delete resonance;
}

bool HardInteraction::IsResonance(int id) const {
  if (resonanceIds.size() > 0) {
    for (std::vector<int>::const_iterator p = resonanceIds.begin(), pe = resonanceIds.end(); p != resonanceIds.end(); p++)
      if (*p == id)
	return true;
    return false;
  }
  else 
    // Default: Z', Z0/gamma*, G, or G*.
    return id == 32 || id == 23 || id == 39 || id == 5000039;
}

void HardInteraction::Clear() {
  quark = resonance = lepPlus = lepMinus = lepPlusNoIB = lepMinusNoIB = 0;
  bremPhotons.clear();
  resonanceIsFake = false;
}

bool HardInteraction::IsValid() const {
  return quark != 0 && resonance != 0 && lepPlus != 0 && lepMinus != 0 && lepPlusNoIB != 0 && lepMinusNoIB != 0;
}

void HardInteraction::Fill(const edm::Event& event) {
  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel(src, genParticles);
  Fill(*genParticles);
}

void HardInteraction::Fill(const reco::GenParticleCollection& genParticles) {
  // Reset everything before filling.
  Clear();

  // Look in the doc lines for the hard-interaction resonance and
  // leptons.
  reco::GenParticleCollection::const_iterator genp = genParticles.begin();
  for (; genp != genParticles.end(); genp++) {
    if (genp->status() == 3) {
      int pdgId = genp->pdgId();

      if (IsResonance(pdgId)) {
	// We found the resonance (Z0/Z'/etc.). Make sure we didn't
	// find a second one.
	if (resonance != 0)
	  edm::LogWarning("HardInteraction")
	    << "Found second resonance (pdgId: " << pdgId << ") in event!";
	else
	  resonance = &*genp;
      }
      else if (pdgId == leptonFlavor) {
	// We found the l-. Make sure we didn't find a second one.
	if (lepMinus != 0)
	  edm::LogWarning("HardInteraction") << "Found second l- in event!";
	else
	  lepMinusNoIB = &*genp;
      }
      else if (pdgId == -leptonFlavor) {
	// We found the l+. Make sure we didn't find a second one.
	if (lepPlus != 0)
	  edm::LogWarning("HardInteraction") << "Found second l+ in event!";
	else
	  lepPlusNoIB = &*genp;
      }
    }
  }

  // We should always find the l+ and l-.
  if (lepMinusNoIB == 0 || lepPlusNoIB == 0)
    edm::LogWarning("HardInteraction")
      << "Couldn't find at least one of the pre-brem l-l+! l- = "
      << lepMinusNoIB << " l+ = " << lepPlusNoIB;

  // We get the quark/antiquark from either the resonance or directly
  // from one of the leptons (see below). Start by assuming it is from
  // the resonance (we'll check for the null pointer later if this is
  // not true).
  const reco::Candidate* mothersAreQuarks = resonance;

  // Some COMPHEP samples do not put the resonance back into the HepMC
  // record. So, if we didn't find a resonance, and are allowed to,
  // build one (we really only care about its four vector and charge
  // at this point). Use the doc-line l+l- we found to do so.
  if (resonance == 0) {
    if (allowFakeResonance && lepPlusNoIB && lepMinusNoIB) {
      reco::Particle::Charge q = lepPlusNoIB->charge() + lepMinusNoIB->charge();
      reco::Particle::LorentzVector p4 = lepPlusNoIB->p4() + lepMinusNoIB->p4();
      // Don't know the creation vertex of the resonance. Set to (0,0,0)
      // by default.
      reco::Particle::Point vtx;
      // We also don't know what the pdgId of the resonance was. Pass
      // one in?
      int pdgId = 99;      
      int status = 3;
    
      // We have ownership of this pointer, and we will delete it in our
      // destructor.
      resonance = new reco::GenParticle(q, p4, vtx, pdgId, status, true);
      resonanceIsFake = true;

      // In these resonanceless COMPHEP samples, the quark/antiquark are
      // direct mothers of the leptons. So now, assume we can get them
      // from the l-.
      mothersAreQuarks = lepMinusNoIB;
    }
    else
      edm::LogWarning("HardInteraction")
	<< "Did not find the resonance in the event, and forbidden"
	<< " from faking it!";
  }

  // Did we successfully identify a Candidate that has the
  // quark/antiquark as its mothers? Try to get them.
  if (mothersAreQuarks != 0 && mothersAreQuarks->numberOfMothers() == 2) {
    for (unsigned m = 0; m < 2; m++) {
      const reco::Candidate* mom = mothersAreQuarks->mother(m);
      int momId = mom->pdgId();

      // For now, don't count gluons (id = 21), but we should
      // implement them later here and in all the above (for gg->G*).
      if (abs(momId) > 6 || momId == 0)
	edm::LogWarning("HardInteraction")
	  << "Mother " << m << " is not a quark! its pdgId = " << momId;
      
      if (momId > 0)
	quark = mom;
      else
	antiquark = mom;
    }
  }
  
  // If we didn't find the quark or antiquark, bomb.
  if (quark == 0 || antiquark == 0)
    edm::LogWarning("HardInteraction")
      << "Couldn't find at least one of the quark/antiquark! quark = "
      << quark << " antiquark = " << antiquark;

  // Now, pick up the leptons after brem. They are daughter leptons of
  // the doc-line leptons and have status = 1, i.e. they are
  // final-state. Also grab the brem photons.
  if (lepMinusNoIB) {
    for (reco::Candidate::const_iterator dau = lepMinusNoIB->begin();
	 dau != lepMinusNoIB->end(); dau++) {
      if (dau->status() == 1) {
	int pdgId = dau->pdgId();
	if (pdgId == leptonFlavor)
	  lepMinus = &*dau;
	else if (pdgId == 22)
	  bremPhotons.push_back(&*dau);
      }
    }
  }

  if (lepPlusNoIB) {
    for (reco::Candidate::const_iterator dau = lepPlusNoIB->begin();
	 dau != lepPlusNoIB->end(); dau++) {
      if (dau->status() == 1) {
	int pdgId = dau->pdgId();
	if (pdgId == -leptonFlavor)
	  lepPlus = &*dau;
	else if (pdgId == 22)
	  bremPhotons.push_back(&*dau);
      }
    }
  }

  // If we didn't find the final-state l- or l+, bomb.
  if (lepMinus == 0 || lepPlus == 0)
    edm::LogWarning("HardInteraction")
      << "Couldn't find at least one of the final-state l-l+! l- = "
      << lepMinus << " l+ = " << lepPlus;
}
