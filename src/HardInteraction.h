#ifndef HARDINTERACTION_H
#define HARDINTERACTION_H

#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/Framework/interface/Event.h"

// A class for extracting and storing the generator-level particles
// that were in the hard interaction (+ the final-state leptons after
// bremsstrahlung) from the MC record, assuming the event was of the
// form q qbar -> resonance -> l+l-, and that there was only one hard
// interaction per event. We store pointers to all these particles.

struct HardInteraction {
  HardInteraction(int lepFlavor, // std::vector<int> resIds,
		  bool allowFakeRes=false);
  ~HardInteraction();

  // Set the list of qualifying resonance ids.
  //static void SetResonanceIds(const std::vector<int>& resIds);

  // Determine whether the passed in pdgId is one of the resonances of
  // interest (currently one of Z0 (inc. DY), Z', or G*).
  static bool IsResonance(int pdgId);

  // Clear out the structure: reset pointers to null, empty
  // bremPhotons, reset flags.
  void Clear();

  // Store pointers to all the particles from the genParticles collection.
  void Fill(const reco::GenParticleCollection& genParticles);
  
  // The same, but get the genParticles from the event.
  void Fill(const edm::Event& event);

  // Pointers to the particles in the hard interaction (as well as the
  // final-state leptons after brem.)
  const reco::Candidate* quark;
  const reco::Candidate* antiquark;
  const reco::Candidate* resonance;
  const reco::Candidate* lepPlus;
  const reco::Candidate* lepMinus;
  const reco::Candidate* lepPlusNoIB;
  const reco::Candidate* lepMinusNoIB;
  std::vector<const reco::Candidate*> bremPhotons;

  // The lepton flavor to look for (e.g. 13 for muons, 11 for electrons).
  int leptonFlavor;

  // The PDG ids of the resonances of interest.
  // JMTBAD Hard-coded for now in isResonance().
  //static std::vector<int> resonanceIds;

  // Whether to allow the building of the resonance from the doc-line
  // leptons if the resonance is not found in the event record (this
  // is the case for at least some COMPHEP samples).
  bool allowFakeResonance;

  // Flag declaring whether we built the resonance ourselves, and
  // therefore own and should delete its pointer at destruction.
  bool resonanceIsFake;
};

#endif
