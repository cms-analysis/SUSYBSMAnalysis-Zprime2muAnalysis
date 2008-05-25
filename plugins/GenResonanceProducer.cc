#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

using namespace std;
using namespace edm;
using namespace reco;

class GenResonanceProducer : public EDProducer {
public:
  explicit GenResonanceProducer(const ParameterSet&);
  ~GenResonanceProducer() {}

private:
  virtual void produce(Event&, const EventSetup&);

  int leptonFlavor;
};

GenResonanceProducer::GenResonanceProducer(const ParameterSet& cfg) 
  : leptonFlavor(cfg.getParameter<int>("leptonFlavor"))
{  
  produces<CompositeCandidateCollection>();
}

void GenResonanceProducer::produce(Event& event,
				   const EventSetup& eSetup) {
  // Set up the output collection.
  auto_ptr<CompositeCandidateCollection>
    resonances(new CompositeCandidateCollection);

  // Use a HardInteraction object to extract the particles we want.
  HardInteraction hi(leptonFlavor, true);
  hi.Fill(event);

  // Make a CompositeCandidate out of the resonance object, setting as
  // its daughters the final-state leptons and brem photons.
  // JMTBAD Do we need shallow clones?
  CompositeCandidate newRes(*hi.resonance);
  newRes.addDaughter(*hi.lepMinus);
  newRes.addDaughter(*hi.lepPlus);
  for (unsigned i = 0; i < hi.bremPhotons.size(); i++)
    newRes.addDaughter(*hi.bremPhotons[i]);
  
  // Store what we made in the event.
  resonances->push_back(newRes);
  event.put(resonances);
}

DEFINE_FWK_MODULE(GenResonanceProducer);
