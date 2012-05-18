#ifndef TRIGGERDECISION_H
#define TRIGGERDECISION_H

#include <string>
#include <vector>
#include "TString.h" 

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
  void initEvent(const edm::Event& event);

  // Get Level-1 decisions for trigger paths we are interested in,
  // storing them in a bitmap.
  void storeL1Decision(const edm::Event& event);

  // Same idea, but for levels 2 and 3 of the HLT. Returns whether the
  // extracted L2 and L3 decisions agree with the overall HLT
  // decision. If ignoreSubLevels, set L2 and L3 decisions equal to
  // the overall HLT one.
  void storeHLTDecision(const edm::Event& event);

  // Return the trigger bits as we have packed them in
  // storeL1/HLTDecision().
  unsigned getWord(const int irec) const;
  
  // Check if we passed the trigger at this rec level.
  bool pass(const int irec) const;

  // Check if we passed the trigger overall.
  bool pass() const;

  // If a trigger level is queried, return its result; otherwise if an
  // offline level is passed in, return the overall trigger result.
  bool pass_all(const int irec) const;

  const std::vector<std::string>& l1_paths() const { return l1Paths; }
  const std::vector<std::string>& hlt_paths() const { return hltPaths; }

  bool l1_path_pass(size_t ipath)  const { assert(ipath < l1Paths.size());  return bool(trigWord[1] & (1 << ipath)); }
  bool hlt_path_pass(size_t ipath) const { assert(ipath < hltPaths.size()); return bool(trigWord[3] & (1 << ipath)); }

    bool pass_hlt_path(edm::Event const& event,std::string shltpath);
    void const dumpPaths(edm::Event const& event) const;

 private:
  bool debug;

  // If we're doing electrons, modify what trigger paths are checked.
  bool doingElectrons;

  // If useTrigger is set false, ignore all the trigger information
  // and pass at every level.
  bool useTrigger;

  // The input tags for the L1 and HLT decision objects.
  edm::InputTag l1GtObjectMap;
  edm::InputTag hltResults;

  // Which trigger paths to use.
  std::vector<std::string> l1Paths;
  std::vector<std::string> hltPaths;

  // The decision we calculated at each level for the current event.
  bool passTrig[TRIG_LEVELS];
  
  // The results of each trigger path, packed as a bitmap for the
  // current event.
  unsigned trigWord[TRIG_LEVELS];
};

#endif
