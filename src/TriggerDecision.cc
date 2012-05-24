#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerDecision.h"

using namespace std;

void TriggerDecision::init(const edm::ParameterSet& config, const bool dbg) {
  debug = dbg;

  doingElectrons = config.getParameter<bool>("doingElectrons");
  useTrigger     = config.getParameter<bool>("useTrigger");
  l1GtObjectMap  = config.getParameter<edm::InputTag>("l1GtObjectMap");
  hltResults     = config.getParameter<edm::InputTag>("hltResults");

  if (useTrigger) {
    l1Paths  = config.getParameter<vector<string> >("l1Paths");
    hltPaths = config.getParameter<vector<string> >("hltPaths");
  }
}

void TriggerDecision::initEvent(const edm::Event& event) {
  for (int i_rec = 0; i_rec < TRIG_LEVELS; i_rec++) {
    passTrig[i_rec] = true;
    trigWord[i_rec] = 0;
  }

  // If we're to ignore trigger info, leave passTrig as true.
  if (useTrigger) {
    storeL1Decision(event);
    storeHLTDecision(event);
  }
}

void TriggerDecision::storeL1Decision(const edm::Event& event) {
  // Get Level-1 decisions for trigger paths we are interested in.
  edm::Handle<L1GlobalTriggerObjectMapRecord> l1Map;
  event.getByLabel(l1GtObjectMap, l1Map);

  static int ifois = 0;
  if (!l1Map.isValid()) {
    if (ifois == 0) {    
      edm::LogWarning("storeL1Decision")
	<< "  ++ This warning will be printed only once ++\n"
	<< "L1GlobalTriggerObjectMapRecord with label "
	<< l1GtObjectMap.encode() << " not found!";
      ifois = 1;
    }
    return;
  }

  ostringstream out;
  if (debug) out << "storeL1Decision for run # " << event.id().run()
		 << " event # " << event.id().event() << ":\n";

  // Loop over chosen paths, check trigger decisions, and pack them
  // into trigbbits.
  unsigned int trigbits = 0;
  for (int ipath = 0; ipath < int(l1Paths.size()); ipath++) {
    const L1GlobalTriggerObjectMap* map = l1Map->getObjectMap(l1Paths[ipath]);
    bool fired = map->algoGtlResult();
    if (fired) trigbits = trigbits | (1 << ipath);
    if (debug) out << "  " << l1Paths[ipath] << " (index "
		   << map->algoBitNumber() << "): decision " << fired << endl;
  }

  if (debug) {
    out << " L1 official trigbits: " << trigbits;
    cout << out.str();
  }

  trigWord[lL1] = trigbits;
  passTrig[lL1] = trigbits != 0;
}

void TriggerDecision::storeHLTDecision(const edm::Event& event) {
  // Get the HLT TriggerResults object, from which the official HLT
  // path decisions can be extracted.
  edm::Handle<edm::TriggerResults> hltRes;
  event.getByLabel(hltResults, hltRes);

  // Get the map between path names and indices.
  const edm::TriggerNames& hltTrigNames = event.triggerNames(*hltRes);

  ostringstream out;

  // Dump the full list of trigger paths if needed.
  static int ifois = 0;
  if (debug && ifois == 0) {
    cout << "Full list of available trigger paths:" << endl;
    for (unsigned int i=0; i<hltTrigNames.size(); i++)
      cout << "trigger path " << i << ": "
           << hltTrigNames.triggerName(i) << endl;
    ifois = 1;
  }

  if (debug) out << "storeHLTDecision for run # " << event.id().run()
		 << " event # " << event.id().event() << ":\n";

  // Get the official HLT results, and pack them.
  unsigned int hlt_trigbits = 0;
  int paths_defined = hltRes->size();
  for (unsigned int i = 0; i < hltPaths.size(); i++) {
    int ndx = hltTrigNames.triggerIndex(hltPaths[i]);
    if (ndx >= paths_defined) {
      edm::LogWarning("storeHLTDecision")
	<< "+++ HLT path " << hltPaths[i]
	<< " is not available in TriggerResults object; skipping it... +++";
    }
    else {
      bool fired = hltRes->accept(ndx);
      if (debug)
	out << "  " << hltPaths[i] << " (index " << ndx << "): " 
	    << " decision " << fired << endl;
      if (fired) hlt_trigbits |= 1 << i;
    }
  }

  // In the new HLT data model, filter decisions and objects firing
  // the trigger are stored in a TriggerEvent object. However, the L2
  // filter module results for muons are not saved (saveTag = False in
  // HLTMuonL2PreFilter in the HLT table) separately, e.g. we do not
  // have both SingleMuNoIsoL2PreFiltered and
  // SingleMuonNoIsoL3PreFiltered saved. (The L2 decisions and
  // TriggerObjects are in HLTDEBUG.) So, skip extracting L2 and L3
  // decisions separately for now.

  trigWord[lL2] = trigWord[lL3] = hlt_trigbits;
  passTrig[lL2] = passTrig[lL3] = hlt_trigbits != 0;
  for (int l = lL2; l <= lL3; l++) 
    out << "  trigWord[l" << l << "]: " << trigWord[l] << endl;

  if (debug) cout << out.str();
}

unsigned TriggerDecision::getWord(const int irec) const {
  if (irec < 0 || irec > lL3)
    throw cms::Exception("trigWord")
      << "L" << irec << " trigger is unknown!\n";

  return trigWord[irec];
}

bool TriggerDecision::pass(const int irec) const {
  // If we're ignoring trigger info, everything passes.
  if (!useTrigger) return true;

  if (irec < 0 || irec > lL3)
    throw cms::Exception("passTrigger")
      << "L" << irec << " trigger is unknown!\n";

  return passTrig[irec];
}

bool TriggerDecision::pass() const {
  // Again, if we're ignoring trigger info, everything passes.
  if (!useTrigger) return true;

  unsigned int decision = 1;

  for (int itrig = lL1; itrig <= lL3; itrig++)
    decision *= trigWord[itrig];
  
  return decision != 0;
}

bool TriggerDecision::pass_all(const int irec) const {
  if (irec <= lL3) return pass(irec);
  else return pass();
}
