#ifndef GenEventTopology_h
#define GenEventTopology_h

#include "DataFormats/Candidate/interface/Particle.h"
#include "FWCore/Framework/interface/Event.h"

enum TopologyStatus {
  topo_ttbar   = 1<<0,
  topo_Wp      = 1<<1,
  topo_Wm      = 1<<2,
  topo_Z       = 1<<3,
  topo_Z2      = 1<<4,
  topo_Zee     = 1<<5,
  topo_Zmm     = 1<<6,
  topo_Ztt     = 1<<7
};

unsigned GenEventTopology(const edm::Event& event, int pdgId[2], reco::Particle::LorentzVector p4[2]);

#endif // GenEventTopology_h
