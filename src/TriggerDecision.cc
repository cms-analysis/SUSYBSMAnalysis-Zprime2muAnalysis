#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"

#include "FWCore/Framework/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerDecision.h"

using namespace std;

void TriggerDecision::init(const edm::ParameterSet& config,
			   const bool dbg) {
  debug = dbg;

  doingElectrons = config.getParameter<bool>("doingElectrons");
  useTrigger     = config.getParameter<bool>("useTrigger");
  l1ParticleMap  = config.getParameter<edm::InputTag>("l1ParticleMap");
  hltResults     = config.getParameter<edm::InputTag>("hltResults");

  if (useTrigger) {
    if (!doingElectrons) {
      // Level-1 paths we want to use for the trigger decision.
      l1paths.push_back(l1extra::L1ParticleMap::kSingleMu7);
      l1paths.push_back(l1extra::L1ParticleMap::kDoubleMu3);
      
      // Level-2 paths (actually, the names of the modules ran)
      hltModules[0].push_back("SingleMuNoIsoL2PreFiltered");
      hltModules[0].push_back("DiMuonNoIsoL2PreFiltered");
      // Level-3 paths (module names)
      hltModules[1].push_back("SingleMuNoIsoL3PreFiltered");
      hltModules[1].push_back("DiMuonNoIsoL3PreFiltered");
      
      // HLT paths (the logical ANDs of L2 and L3 single/dimuon paths
      // above)
      hltPaths.push_back("HLT1MuonNonIso");
      hltPaths.push_back("HLT2MuonNonIso");
    }
    else {
      l1paths.push_back(l1extra::L1ParticleMap::kSingleEG15);
      
      // For now, just look at the overall HLT decision for electrons.
      hltPaths.push_back("HLT1EMHighEt");
      hltPaths.push_back("HLT1EMVeryHighEt");
      hltPaths.push_back("HLT1ElectronRelaxed");
    }
  }
}

bool TriggerDecision::initEvent(const edm::Event& event,
				const bool ignoreSubLevels) {
  for (int i_rec = 0; i_rec < TRIG_LEVELS; i_rec++) {
    passTrig[i_rec] = true;
    trigWord[i_rec] = 0;
  }

  if (useTrigger) {
    event.getByLabel(l1ParticleMap, l1MapColl);
    if (!l1MapColl.isValid()) {
      edm::LogWarning("storeL1Decision")
	<< "L1ParticleMapCollection with label [" << l1ParticleMap.encode()
	<< "] not found!" << endl;
      return false;
    }
    storeL1Decision(event);
    return storeHLTDecision(event, ignoreSubLevels);
  }

  return true;
}

void TriggerDecision::storeL1Decision(const edm::Event& event) {
  // Get Level-1 decisions for trigger paths we are interested in.
  ostringstream out;
  if (debug) out << "storeL1Decision:\n";

  // Loop over chosen paths, check trigger decisions, and save them into
  // "trigbits".
  unsigned int trigbits = 0;
  int nl1paths = l1paths.size();
  for (int ipath = 0; ipath < nl1paths; ipath++) {
    const l1extra::L1ParticleMap& thisMap = (*l1MapColl)[l1paths[ipath]];
    bool fired = thisMap.triggerDecision();
    if (fired) trigbits = trigbits | (1 << ipath);
    if (debug) out << "  " << thisMap.triggerName() << " (index "
		   << l1paths[ipath] << "): decision " << fired << endl;
  }

  if (debug) {
    out << " L1 official trigbits: " << trigbits;
    LogTrace("storeL1Decision") << out.str();
  }

  trigWord[l1] = trigbits;
  passTrig[l1] = (trigbits != 0);
}

bool TriggerDecision::storeHLTDecision(const edm::Event& event,
				       const bool ignoreSubLevels) {
  // Getting the result of the entire HLT path is done easily with the
  // TriggerResults object, however to get at the separate decisions
  // for L2 and L3 we have to do a little hacky magic.

  // Try to get the HLT TriggerResults object now, before
  // trying to get the HLTFilterObjectWithRefs below
  // so that if there is no HLT information in the file, 
  // getByLabel will go ahead and throw an exception.
  edm::Handle<edm::TriggerResults> hltRes;
  event.getByLabel(hltResults, hltRes);
  edm::TriggerNames hltTrigNames;
  hltTrigNames.init(*hltRes);

  ostringstream out;
  if (debug) out << "storeHLTDecision:\n";

  // Get the official HLT results, and pack them.
  unsigned int hlt_trigbits = 0;
  for (unsigned int i = 0; i < hltPaths.size(); i++) {
    int ndx = hltTrigNames.triggerIndex(hltPaths[i]);
    bool fired = hltRes->accept(ndx);
    if (debug)
      out << " HLT path #" << ndx << ": " << hltPaths[i]
	  << " decision = " << fired << endl;
    if (fired) hlt_trigbits |= 1 << i;
  }

  if (ignoreSubLevels) {
    trigWord[l2] = trigWord[l3] = hlt_trigbits;
    passTrig[l2] = passTrig[l3] = hlt_trigbits != 0;
    for (int l = l2; l <= l3; l++) 
      out << "  trigWord[l" << l << "]: " << trigWord[l] << endl;
    if (debug) LogTrace("storeHLTDecision") << out.str();
    return true;
  }
  
  // Extract L2 and L3 decisions by seeing if the corresponding
  // HLTFilterObjectWithRefs exists and seeing how many muons it holds.
  for (unsigned int lvl = 0; lvl < 2; lvl++) {
    unsigned int l = l2 + lvl;
    unsigned int trigbits = 0;
    for (unsigned int ipath = 0; ipath < hltModules[lvl].size(); ipath++) {
      const string& trigName = hltModules[lvl][ipath];
      edm::Handle<reco::HLTFilterObjectWithRefs> hltFilterObjs;

      bool fired = true;
      unsigned failAt = 0;
      try {
	event.getByLabel(trigName, hltFilterObjs);
      }
      catch (const cms::Exception& e) {
	fired = false;
	failAt = 1;
      }

      // There may be an HLTFilterObject in the event even if the
      // trigger did not accept; the real check is to make sure that
      // the minimum number of muons for the trigger was met.
      const unsigned minNMuons = ipath + 1;
      unsigned nmu = 0;
      if (fired) {
	nmu = hltFilterObjs->size();
	if (nmu < minNMuons) {
	  fired = false;
	  failAt = 2;
	}
      }
	
      if (debug)
	out << "  " << trigName << ": decision = " << fired
	    << " (#mu: " << nmu << "; failAt: " << failAt << ")\n";

      if (fired) {
	trigbits = trigbits | (1 << ipath);

	if (debug) {
	  out << "  " << trigName << " filter result muons:\n";
	  reco::HLTFilterObjectWithRefs::const_iterator muItr;
	  int imu = 0;
	  for (muItr = hltFilterObjs->begin(); muItr != hltFilterObjs->end();
	       muItr++) {
	    out << "    #" << imu++ << " q = " << muItr->charge()
		<< " p = (" << muItr->px() << ", " << muItr->py()
		<< ", " << muItr->pz() << ", " << muItr->energy() << ")\n"
		<< "     pt = " << muItr->pt() << " eta = " << muItr->eta()
		<< " phi = " << muItr->phi() << endl;
	  }
	}
      }
    }

    trigWord[l] = trigbits;
    passTrig[l] = (trigbits != 0);
    out << "  trigWord[l" << l << "]: " << trigWord[l] << endl;
  }
   
  bool retval = true;
  // Check that the official full HLT path decisions agree with
  // what we extracted for the official L3 decision.
  if (hlt_trigbits != trigWord[l3]) {
    edm::LogWarning("storeHLTDecision")
      << "+++ Warning: official HLT"
      << " decision disagrees with extracted L3 decision:"
      << " official HLT: " << hlt_trigbits
      << " extracted L3: " << trigWord[l3] << " +++";
    retval = false;
  }

  if (debug) LogTrace("storeHLTDecision") << out.str();

  return retval;
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

  // JMTBAD keep electron and muon trigger info separate for every
  // event, so we can look at opposite-flavor events. Possibly just
  // shift the bits over?

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

