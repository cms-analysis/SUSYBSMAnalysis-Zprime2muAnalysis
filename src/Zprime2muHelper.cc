#include "TStyle.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muHelper.h"

Zprime2muHelper::Zprime2muHelper(const edm::ParameterSet& config) 
  : verbosity(Verbosity(config.getUntrackedParameter<int>("verbosity", 0))),
    maxDileptons(config.getParameter<unsigned>("maxDileptons")),
    doingElectrons(config.getParameter<bool>("doingElectrons")),
    useGen(config.getParameter<bool>("useGen")),
    useSim(config.getParameter<bool>("useSim")),
    useTrigger(config.getParameter<bool>("useTrigger")),
    useRaw(config.getParameter<bool>("useRaw")),
    useReco(config.getParameter<bool>("useReco")),
    leptonFlavor(doingElectrons ? 11 : 13),
    leptonMass(doingElectrons ? 0.000511 : 0.10566)
{
  InitROOT();

  // Hidden option for me (JMT) to make tkdiffing outputs when
  // validating new versions easier.
  if (!config.getUntrackedParameter<bool>("dateHistograms", true))
    gStyle->SetOptDate(0);

  if (useGen) hardInteraction.init(leptonFlavor, true);
  if (useTrigger) trigDecision.init(config, verbosity >= V_SIMPLE);
}

void Zprime2muHelper::initEvent(const edm::Event& event, const edm::EventSetup& eSetup) {
  // Get the hard interaction from the MC record.
  if (useGen) hardInteraction.Fill(event);

  // Get the trigger decision from the event. For now, don't bother
  // looking at sub-levels since L2 decisions are no longer stored in
  // the event. JMTBAD For L2, we could take the result of
  // TriggerTranslator()...
  if (useTrigger) trigDecision.initEvent(event);
}

reco::Particle::LorentzVector Zprime2muHelper::resonanceP4(const pat::CompositeCandidate& cand) const {
  // JMTBAD todo, using the stored photon four-vectors in each daughter muon.
  return cand.p4();
}
