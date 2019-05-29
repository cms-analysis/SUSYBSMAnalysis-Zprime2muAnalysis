#ifndef LRWEIGHTPRODUCER_H
#define LRWEIGHTPRODUCER_H

#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include <complex>

namespace edm {
  class ParameterSet;
}

struct LRWeightProducer {
  LRWeightProducer(const edm::ParameterSet cfg);
  ~LRWeightProducer();

  double calculateWeight(const edm::Event& event, HardInteraction*& hardInteraction, double alpEM);


  // Whether doing electrons.
  const bool doingElectrons;
  const int lambda;
  const bool doingLR;
  const bool calculate;
  // The lepton flavor to look for (e.g. 13 for muons, 11 for electrons).
  const int leptonFlavor;

  // The lepton mass, for convenience later.
  const double leptonMass;


const double ef[20] = { 0., -1./3., 2./3., -1./3., 2./3., -1./3., 
  2./3., -1./3., 2./3., 0., 0., -1., 0., -1., 0., -1., 0., -1., 0., 0.};
const double af[20] = { 0., -1., 1., -1., 1., -1., 1., -1., 1., 
0., 0., -1., 1., -1., 1., -1., 1., -1., 1., 0.};
};

#endif
