#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muHelper.h"

Zprime2muHelper::Zprime2muHelper(const edm::ParameterSet& config) 
  : maxDileptons(config.getParameter<unsigned>("maxDileptons")),
    doingElectrons(config.getParameter<bool>("doingElectrons")),
    useGen(config.getParameter<bool>("useGen")),
    useSim(config.getParameter<bool>("useSim")),
    useTrigger(config.getParameter<bool>("useTrigger")),
    useRaw(config.getParameter<bool>("useRaw")),
    useReco(config.getParameter<bool>("useReco")),
    leptonFlavor(doingElectrons ? 11 : 13),
    leptonMass(doingElectrons ? 0.000511 : 0.10566)
{
  // Hidden option for me (JMT) to make tkdiffing outputs when
  // validating new versions easier.
  InitROOT(config.getUntrackedParameter<bool>("dateHistograms", true));

  if (useGen) hardInteraction.init(leptonFlavor, true);
  if (useTrigger) trigDecision.init(config);
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
  // Start with the p4 of the candidate, which is the sum of the
  // daughter p4s by construction.
  reco::Particle::LorentzVector p4 = cand.p4();

  std::vector<int> used_already;

  for (pat::CompositeCandidate::const_iterator dau = cand.begin(), daue = cand.end(); dau != daue; ++dau) {
    // Only muons can have photons matched to them. Also, remember
    // that the CompositeCandidates were formed using shallow clones,
    // so need to get the daughters' masterClones to be able to get
    // the original pat::Muon*.
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>(&*dau->masterClone());
    if (mu == 0)
      continue;

    if (mu->hasUserInt("photon_index")) {
      int ndx = mu->userInt("photon_index");

      // Don't double count, in the case where two daughter muons
      // ended up matching to the same photon.
      if (std::find(used_already.begin(), used_already.end(), ndx) == used_already.end()) {
	used_already.push_back(ndx);
	const reco::Particle::LorentzVector* pho_p4 = mu->userData<reco::Particle::LorentzVector>("photon_p4");
	if (pho_p4 == 0)
	  throw cms::Exception("BadlyFormedPhotonUserData") << "candidate daughter had userInt for photon_index but no p4!\n";
	p4 += *pho_p4;
      }
    }
  }

  return p4;
}
