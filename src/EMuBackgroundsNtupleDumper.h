#ifndef EmuBackgroundsNtupleDumper_h
#define EmuBackgroundsNtupleDumper_h

#include "TTree.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

class EMuBackgroundsNtupleDumper : public edm::EDAnalyzer {
 public:
  explicit EMuBackgroundsNtupleDumper(const edm::ParameterSet&);
  virtual void analyze(const edm::Event& ev, const edm::EventSetup& setup);

 private:
  int selfProcId;
  TTree* tree;
};

#endif
