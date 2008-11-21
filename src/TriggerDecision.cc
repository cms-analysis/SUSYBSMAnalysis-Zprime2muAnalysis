#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "FWCore/Framework/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerDecision.h"

using namespace std;

void TriggerDecision::init(const edm::ParameterSet& config,
			   const bool dbg) {
  debug = dbg;

  doingElectrons     = config.getParameter<bool>("doingElectrons");
  useTrigger         = config.getParameter<bool>("useTrigger");
  l1GtObjectMap      = config.getParameter<edm::InputTag>("l1GtObjectMap");
  hltResults         = config.getParameter<edm::InputTag>("hltResults");

  if (useTrigger) {
    if (!doingElectrons) {
      // Level-1 paths we want to use for the trigger decision.
      l1Paths.push_back("L1_SingleMu7");
      l1Paths.push_back("L1_DoubleMu3");

      // Level-2 paths (actually, the names of the modules ran)
      hltModules[0].push_back("SingleMuNoIsoL2PreFiltered");
      hltModules[0].push_back("DiMuonNoIsoL2PreFiltered");
      // Level-3 paths (module names)
      hltModules[1].push_back("SingleMuNoIsoL3PreFiltered");
      hltModules[1].push_back("DiMuonNoIsoL3PreFiltered");

      // HLT paths (the logical ANDs of L2 and L3 single/dimuon paths
      // above) in 2E30 menu.
      hltPaths.push_back("HLT_Mu15");      // former HLT1MuonNonIso
      hltPaths.push_back("HLT_DoubleMu3"); // former HLT2MuonNonIso
    }
    else {
      l1Paths.push_back("L1_SingleEG15");

      // For now, just look at the overall HLT decision for electrons.
      // JMTBAD trigger names have changed; which are the "right" ones for electrons?
      hltPaths.push_back("HLT_EM80");
      hltPaths.push_back("HLT_EM200");
      //hltPaths.push_back("HLT1ElectronRelaxed");
    }
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
    LogTrace("storeL1Decision") << out.str();
  }

  trigWord[l1] = trigbits;
  passTrig[l1] = trigbits != 0;
}

void TriggerDecision::storeHLTDecision(const edm::Event& event) {
  // Get the HLT TriggerResults object, from which the official HLT
  // path decisions can be extracted.
  edm::Handle<edm::TriggerResults> hltRes;
  event.getByLabel(hltResults, hltRes);

  // Get the map between path names and indices.
  edm::TriggerNames hltTrigNames(*hltRes);

  ostringstream out;
  if (debug) out << "storeHLTDecision:\n";

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

  trigWord[l2] = trigWord[l3] = hlt_trigbits;
  passTrig[l2] = passTrig[l3] = hlt_trigbits != 0;
  for (int l = l2; l <= l3; l++) 
    out << "  trigWord[l" << l << "]: " << trigWord[l] << endl;

  edm::Handle<trigger::TriggerEvent> trigEvent;
  event.getByLabel("hltTriggerSummaryAOD", trigEvent);
  if (trigEvent.isValid()) {
    // Get indices of L2 and L3 muon objects.
    int l2ind_first = -99, l2ind_last = -99;
    int l3ind_first = -99, l3ind_last = -99;
    int ind_prev = 0;
    const trigger::size_type nC(trigEvent->sizeCollections());
    for (trigger::size_type iC = 0; iC < nC; ++iC) {
      std::string encodedCollectionTag(trigEvent->collectionTag(iC).encode());
      if (strstr(encodedCollectionTag.c_str(), "hltL2MuonCandidates")) {
	l2ind_first = ind_prev;
	l2ind_last  = trigEvent->collectionKey(iC)-1;
      }
      else if (strstr(encodedCollectionTag.c_str(), "hltL3MuonCandidates")) {
	l3ind_first = ind_prev;
	l3ind_last  = trigEvent->collectionKey(iC)-1;
      }
      ind_prev = trigEvent->collectionKey(iC);
    }

    // Access collections of L2 and L3 muons.  Shall we convert them to
    // reco::muons, similarly to what is done in L3MuonSanitizer.cc, and
    // save into allLeptons?
    const trigger::TriggerObjectCollection& TOC(trigEvent->getObjects());
    if (!doingElectrons) {
      if (l2ind_first >= 0 && l2ind_last >= l2ind_first) {
	out << "  L2 muon(s): #, id, pt, eta, phi\n";
	for (trigger::size_type iO = l2ind_first; iO <= l2ind_last; ++iO) {
	  const trigger::TriggerObject& TO(TOC[iO]);
	  out << "  " << iO << " " << TO.id() << " " << TO.pt()
	      << " " << TO.eta() << " " << TO.phi() << "\n";
	}
      }
      if (l3ind_first >= 0 && l3ind_last >= l3ind_first) {
	out << "  L3 muon(s): #, id, pt, eta, phi\n";
	for (trigger::size_type iO = l3ind_first; iO <= l3ind_last; ++iO) {
	  const trigger::TriggerObject& TO(TOC[iO]);
	  out << "  " << iO << " " << TO.id() << " " << TO.pt()
	      << " " << TO.eta() << " " << TO.phi() << "\n";
	}
      }
    }
  }

  if (debug) LogTrace("storeHLTDecision") << out.str();
}

unsigned TriggerDecision::getWord(const int irec) const {
  if (irec < 0 || irec > l3)
    throw cms::Exception("trigWord")
      << "L" << irec << " trigger is unknown!\n";

  return trigWord[irec];
}

bool TriggerDecision::pass(const int irec) const {
  // If we're ignoring trigger info, everything passes.
  if (!useTrigger) return true;

  if (irec < 0 || irec > l3)
    throw cms::Exception("passTrigger")
      << "L" << irec << " trigger is unknown!\n";

  return passTrig[irec];
}

bool TriggerDecision::pass() const {
  // Again, if we're ignoring trigger info, everything passes.
  if (!useTrigger) return true;

  unsigned int decision = 1;

  for (int itrig = l1; itrig <= l3; itrig++)
    decision *= trigWord[itrig];
  
  return decision != 0;
}

