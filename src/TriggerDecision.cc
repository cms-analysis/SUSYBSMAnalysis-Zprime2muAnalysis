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

bool TriggerDecision::pass_hlt_path(edm::Event const& event,std::string shltpath){

    bool retval = false;
    edm::Handle<edm::TriggerResults> hltRes;
    event.getByLabel(hltResults, hltRes);
    const edm::TriggerNames& hltTrigNames = event.triggerNames(*hltRes);

  // Get the map between path names and indices.
    int paths_defined = hltRes->size();
    int ndx = hltTrigNames.triggerIndex(shltpath);
//    std::cout<<hltPaths[i]<<"\t"<<ndx<<std::endl;
    if (ndx >= paths_defined) return retval; 
    else retval = hltRes->accept(ndx);

    return retval;
} 
/*
bool TriggerDecision::pass_hlt_path_like(edm::Event const& event,std::string shltpath){

    bool retval = false;
    edm::Handle<edm::TriggerResults> hltRes;
    event.getByLabel(hltResults, hltRes);

  // Get the map between path names and indices.
    const edm::TriggerNames& hltTrigNames = event.triggerNames(*hltRes);
    int paths_defined = hltRes->size();
    int ndx = hltTrigNames.triggerIndex(shltpath);
//    std::cout<<hltPaths[i]<<"\t"<<ndx<<std::endl;
    if (ndx >= paths_defined) return retval; 
    else retval = hltRes->accept(ndx);

    return retval;
} 
*/
void TriggerDecision::storeL1Decision(const edm::Event& event) {
  // Get Level-1 decisions for trigger paths we are interested in.
  edm::Handle<L1GlobalTriggerObjectMapRecord> l1Map;
  event.getByLabel(l1GtObjectMap, l1Map);
  if (!l1Map.isValid()) {
    edm::LogWarning("storeL1Decision")
      << "L1GlobalTriggerObjectMapRecord with label "
      << l1GtObjectMap.encode() << " not found!";
    return;
  }

  ostringstream out;
  if (debug) out << "storeL1Decision:\n";

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
    edm::LogVerbatim("storeL1Decision") << out.str();
  }

  trigWord[lL1] = trigbits;
  passTrig[lL1] = trigbits != 0;
}
//
// store the HLT decisions
//
void TriggerDecision::storeHLTDecision(const edm::Event& event) {
  // Get the HLT TriggerResults object, from which the official HLT
  // path decisions can be extracted.
  edm::Handle<edm::TriggerResults> hltRes;
  event.getByLabel(hltResults, hltRes);

  // Get the map between path names and indices.
  const edm::TriggerNames& hltTrigNames = event.triggerNames(*hltRes);

  ostringstream out;
  if (debug) out << "storeHLTDecision:\n";

  // Get the official HLT results, and pack them.
  unsigned int hlt_trigbits = 0;
  int paths_defined = hltRes->size();
  for (unsigned int i = 0; i < hltPaths.size(); i++) {
    int ndx = hltTrigNames.triggerIndex(hltPaths[i]);
//    std::cout<<hltPaths[i]<<"\t"<<ndx<<std::endl;
    if (ndx >= paths_defined) {
      edm::LogWarning("storeHLTDecision")
	<< "+++ HLT path " << hltPaths[i]
	<< " is not available in TriggerResults object; skipping it... +++";
    }
    else {
      bool fired = hltRes->accept(ndx);
      if (debug||true)
	out << "  " << hltPaths[i] << " (index " << ndx << "): " 
//	std::cout << "  " << hltPaths[i] << " (index " << ndx << "): " 
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

  if (debug) edm::LogVerbatim("storeHLTDecision") << out.str();
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
//
// dump available paths
//
void const TriggerDecision::dumpPaths(edm::Event const& event) const {
        std::cout<<"---- dumping the trigger names"<<std::endl;

    edm::Handle<edm::TriggerResults> hltRes;
    event.getByLabel(hltResults, hltRes);
    const edm::TriggerNames& hltTrigNames = event.triggerNames(*hltRes);

    std::vector<std::string>  hlNames_;  // name of each HLT algorithm
    hlNames_=hltTrigNames.triggerNames();
    std::cout<<"\t num Triggers: "<<hlNames_.size()<<std::endl;
/*
  unsigned int hlt_trigbits = 0;
  int paths_defined = hltRes->size();
  for (unsigned int i = 0; i < hltPaths.size(); i++) {
    int ndx = hltTrigNames.triggerIndex(hltPaths[i]);
    std::cout<<hltPaths[i]<<"\t"<<ndx<<std::endl;
*/
/*
    edm::Handle<edm::TriggerResults> triggerEventHandle_;
    event.getByLabel(hltResults, triggerEventHandle_);
    edm::Handle<edm::TriggerResults> triggerResultsHandle_;
  // Get the map between path names and indices.
//    if (!triggerEventHandle_.isValid()) {
        const edm::TriggerNames & triggerNames = event.triggerNames(*triggerResultsHandle_);
        std::vector<std::string>  hlNames_;  // name of each HLT algorithm
        hlNames_=triggerNames.triggerNames();
        std::cout<<"\t num Triggers: "<<hlNames_.size()<<std::endl;
        for (unsigned int i=0; i!=hlNames_.size(); ++i) {
            TString _s = hlNames_[i];
            if (!_s.Contains("Mu")) continue;
            std::cout<<i<<"\t"<<_s.Data()<<std::endl;
        }
        std::cout<<"---- done looping on triggers"<<std::endl;

        return;
    } else {
        std::cout<<"The triggerEvent is broken, dude."<<std::endl;
    }
*/
}
