#ifndef ZP2MUTRIGCOMPARISON_H
#define ZP2MUTRIGCOMPARISON_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muRecLevelAnalysis.h"

class Zprime2muTrigComparison : public Zprime2muRecLevelAnalysis {
 public:
  explicit Zprime2muTrigComparison(const edm::ParameterSet& config)
    : Zprime2muRecLevelAnalysis(config) {}
  void analyze(const edm::Event& event, const edm::EventSetup& eSetup);

 private:
  // Function to translate the algorithm algo defined for trigger
  // level lvl, requiring nmu muons (e.g. single or dimuon trigger),
  // and return whether we judge the event to pass the trigger at this
  // level for this algorithm.
  bool TriggerTranslator(const std::string& algo, const unsigned int lvl, 
			 const unsigned int nmu) const;

  // Use TriggerTranslator to compare the "official" decisions stored
  // by store*Decision with what we calculate.
  void compareTrigDecision(const edm::Event& event) const;

};

#endif // ZP2MUTRIGCOMPARISON_H
