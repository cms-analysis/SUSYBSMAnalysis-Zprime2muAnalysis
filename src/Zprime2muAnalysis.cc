#include "TStyle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

using namespace std;

////////////////////////////////////////////////////////////////////
// Configuration and driver methods
////////////////////////////////////////////////////////////////////

Zprime2muAnalysis::Zprime2muAnalysis(const edm::ParameterSet& config) 
  : verbosity(Verbosity(config.getUntrackedParameter<int>("verbosity", 0))),
    maxDileptons(config.getParameter<unsigned>("maxDileptons")),
    doingElectrons(config.getParameter<bool>("doingElectrons")),
    useGen(config.getParameter<bool>("useGen")),
    useSim(config.getParameter<bool>("useSim")),
    useReco(config.getParameter<bool>("useReco")),
    useOtherMuonRecos(config.getParameter<bool>("useOtherMuonRecos")),
    usingAODOnly(config.getParameter<bool>("usingAODOnly")),
    eventNum(-1),
    eventsDone(0),
    genDils(config.getParameter<edm::InputTag>("genDileptons")),
    hltDils(config.getParameter<edm::InputTag>("hltDileptons")),
    recDils(config.getParameter<edm::InputTag>("recDileptons")),
    bestDils(config.getParameter<edm::InputTag>("bestDileptons"))
{
  InitROOT();

  // Hidden option for me (JMT) to make tkdiffing outputs when
  // validating new versions easier.
  if (!config.getUntrackedParameter<bool>("dateHistograms", true))
    gStyle->SetOptDate(0);

  if (!doingElectrons) {
    leptonFlavor = 13;
    leptonMass = 0.10566;
  }
  else {
    leptonFlavor = 11;
    leptonMass = 0.000511;
  }

  trigDecision.init(config, verbosity >= VERBOSITY_SIMPLE);

  // For muons, only apply the cuts that were used for the 2007 AN for
  // now.  JMTBAD The result of performing cuts on electrons is the
  // "best" electrons collection produced by HEEPSelector; so there is
  // kind of an asymmetry between the muon code path and the electron
  // one.  Perhaps use HEEPHelper here instead of TeVMuHelper to make
  // a cut code?
  unsigned cutMask = doingElectrons ? 0 : (TeVMuHelper::PT | TeVMuHelper::ISO);
  tevMuHelper.setCutMask(cutMask);
}

void Zprime2muAnalysis::analyze(const edm::Event& event,
				const edm::EventSetup& eSetup) {
  // We could store the whole event/eventsetup object if needed, but
  // for now just store the event number.
  eventNum = event.id().event();
 
  // Keep track of how many events we run over total.
  eventsDone++;

  // Get the trigger decision from the event. If there was some funny
  // business, print out the trigger collections.  For now, don't
  // bother looking at sub-levels for electrons. Also, the AOD does
  // not store the L2 path result separately, so ignore the sub-levels
  // if we're running on AOD. JMTBAD For L2, we could take the result
  // of TriggerTranslator()...
  trigDecision.initEvent(event, doingElectrons || usingAODOnly);

  // Get the main dilepton collections: gen, HLT, default offline, and
  // "best" offline.
  if (useGen)
    event.getByLabel(genDils,  genDileptons);
  if (useReco) {
    event.getByLabel(hltDils,  hltDileptons);
    event.getByLabel(recDils,  recDileptons);
    event.getByLabel(bestDils, bestDileptons);
  }
}

DEFINE_FWK_MODULE(Zprime2muAnalysis);
