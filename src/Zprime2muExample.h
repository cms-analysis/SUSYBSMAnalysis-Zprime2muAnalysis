#ifndef ZP2MUEXAMPLE_H
#define ZP2MUEXAMPLE_H

#include "TH1F.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

class Zprime2muExample : public Zprime2muAnalysis {
 public:
  explicit Zprime2muExample(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  TH1F* hDilMass;
};

#endif // ZP2MUEXAMPLE_H
