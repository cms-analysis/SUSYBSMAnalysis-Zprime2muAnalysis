#ifndef Zprime2muHelper_h
#define Zprime2muHelper_h

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

namespace pat {
  class CompositeCandidate;
}

class Zprime2muHelper {
 public:
  Zprime2muHelper(const edm::ParameterSet&);

  void initEvent(const edm::Event&, const edm::EventSetup&);
  reco::Particle::LorentzVector resonanceP4(const pat::CompositeCandidate&) const;

  // Keep this many highest-invariant-mass dileptons.
  const unsigned maxDileptons;

  // whether we are looking at electrons instead of muons;
  const bool doingElectrons;

  // whether to look at generator-level information;
  const bool useGen;

  // whether to look at GEANT tracks;
  const bool useSim;

  // whether to look at trigger quantities;
  const bool useTrigger;

  // whether to expect quantities from raw data;
  const bool useRaw;

  // whether to look at reconstructed quantities;
  const bool useReco;

  // basic quantities for the chosen lepton (muon or electron);
  const unsigned int leptonFlavor; // PDG ID
  const double       leptonMass;   // in GeV/c^2

  ////////////////////////////////////////////////////////////////////
  // Event data
  ////////////////////////////////////////////////////////////////////

  // Helper object to extract the important generator-level objects:
  // the resonance and its decay products (including brem. photons).
  HardInteraction hardInteraction;
};

#endif
