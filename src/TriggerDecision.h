#ifndef TRIGGERDECISION_H
#define TRIGGERDECISION_H

#include <string>
#include <vector>

#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"

class TriggerDecision {
 public:
  // Set up per-lifetime things, including which trigger paths to
  // check.
  void init(const edm::ParameterSet& config, const bool dbg=false);

  // Process the event: get the L1/HLT decisions. Return the result of
  // storeHLTDecision().
  bool initEvent(const edm::Event& event, const bool ignoreSubLevels);

  // Get Level-1 decisions for trigger paths we are interested in,
  // storing them in a bitmap.
  void storeL1Decision(const edm::Event& event);

  // Same idea, but for levels 2 and 3 of the HLT. Returns whether the
  // extracted L2 and L3 decisions agree with the overall HLT
  // decision. If ignoreSubLevels, set L2 and L3 decisions equal to
  // the overall HLT one.
  bool storeHLTDecision(const edm::Event& event, const bool ignoreSubLevels);

  // Return the trigger bits as we have packed them in
  // storeL1/HLTDecision().
  unsigned getWord(const int irec) const;
  
  // Check if we passed the trigger at this rec level.
  bool pass(const int irec) const;

  // Check if we passed the trigger overall.
  bool pass() const;

  // Allow getting the L1 particle map.
  const l1extra::L1ParticleMapCollection& getL1ParticleMap() const
    { return *l1MapColl; }

  // Allow reading which L1/HLT paths we use.
  const std::vector<l1extra::L1ParticleMap::L1TriggerType>& getL1Paths() const
    { return l1paths; }
  const std::vector<std::string>& getHLTModules(unsigned which) const
    { return hltModules[which]; }

 private:
  bool debug;

  // If we're doing electrons, modify what trigger paths are checked.
  bool doingElectrons;

  // If useTrigger is set false, ignore all the trigger information
  // and pass at every level.
  bool useTrigger;

  // The input tags for the L1 and HLT decision objects.
  edm::InputTag l1ParticleMap;
  edm::InputTag hltResults;

  // Which trigger paths to use.
  std::vector<l1extra::L1ParticleMap::L1TriggerType> l1paths;
  edm::Handle<l1extra::L1ParticleMapCollection> l1MapColl;
  std::vector<std::string> hltModules[2]; // in order: L2, L3
  std::vector<std::string> hltPaths;

  // The decision we calculated at each level for the current event.
  bool passTrig[TRIG_LEVELS];
  
  // The results of each trigger path, packed as a bitmap for the
  // current event.
  unsigned trigWord[TRIG_LEVELS];
};

#endif
