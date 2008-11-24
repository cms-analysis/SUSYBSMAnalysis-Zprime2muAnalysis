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
    useTrigger(config.getParameter<bool>("useTrigger")),
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
  tevMuHelper.setCutMask(doingElectrons ? 0 : TeVMuHelper::tevMuCuts);
}

void Zprime2muAnalysis::analyze(const edm::Event& event,
				const edm::EventSetup& eSetup) {
  // We could store the whole event/eventsetup object if needed, but
  // for now just store the event number.
  eventNum = event.id().event();
 
  // Keep track of how many events we run over total.
  eventsDone++;

  // Get the trigger decision from the event. For now, don't bother
  // looking at sub-levels since L2 decisions are no longer stored in
  // the event. JMTBAD For L2, we could take the result of
  // TriggerTranslator()...
  trigDecision.initEvent(event);

  // Store the pointer to the current event in the helper object.
  tevMuHelper.initEvent(event);

  // Get the main dilepton collections: gen, HLT, default offline, and
  // "best" offline.
  if (useGen)
    event.getByLabel(genDils,  genDileptons);
  if (useTrigger)
    event.getByLabel(hltDils,  hltDileptons);
  if (useReco) {
    event.getByLabel(recDils,  recDileptons);
    event.getByLabel(bestDils, bestDileptons);
  }
}

double getIso(const reco::CandidateBaseRef& cand) {
  const reco::Muon* mu = toConcretePtr<reco::Muon>(cand);
  return mu == 0 ? 0 : mu->isolationR03().sumPt;
}

void Zprime2muAnalysis::dumpLepton(ostream& output,
				   const reco::CandidateBaseRef& cand) const {
  output << " pdgId: " << cand->pdgId()
	 << " charge: " << cand->charge() << " pt: " << cand->pt()
	 << " eta: " << cand->eta() << " phi: " << cand->phi()
	 << " sumPt: " << getIso(cand) << endl;
}

void Zprime2muAnalysis::dumpDilepton(ostream& output,
				     const reco::CompositeCandidate& cand) const {
  output << "Dilepton: charge: " << cand.charge()
	 << " pt: " << cand.pt() << " eta: " << cand.eta()
	 << " phi: " << cand.phi() << " mass: " << cand.mass() << endl;

  int larger = dileptonDaughter(cand, 0)->p() > dileptonDaughter(cand, 1)->p() ? 0 : 1;
  int smaller = larger == 0 ? 1 : 0;
  const reco::CandidateBaseRef& cand1 = dileptonDaughter(cand, larger);
  const reco::CandidateBaseRef& cand2 = dileptonDaughter(cand, smaller);

  output << "Higher momentum daughter:\n";
  dumpLepton(output, cand1);
  output << "Lower momentum daughter:\n";
  dumpLepton(output, cand2);
}

DEFINE_FWK_MODULE(Zprime2muAnalysis);

