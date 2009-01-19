#include "TStyle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"

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
    lBest(config.getParameter<int>("bestRecLevel")),
    eventNum(-1),
    eventsDone(0)
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
  cutHelper.setCutMask(config.getParameter<unsigned>("cutMask"));
  hardInteraction.init(leptonFlavor, true);
  recLevelHelper.init(config);
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
  cutHelper.initEvent(event);

  // Get the hard interaction from the MC record.
  if (useGen) hardInteraction.Fill(event);

  // Store per-rec-level information from the event (includes getting
  // all the match maps from the event).
  recLevelHelper.initEvent(event);

  // At each rec level (including the "best" leptons), store the
  // lepton and dilepton collections, including the "resonances"
  // (dileptons + closest photons added in).
  for (int irec = 0; irec < MAX_LEVELS; irec++) {
    recLevelHelper.getLeptons(event, irec, allLeptons[irec]);
    recLevelHelper.getDileptons(event, irec, allDileptons[irec]);
  }
  
  // Dump the event if appropriate.
  if (verbosity >= VERBOSITY_SIMPLE)
    dumpEvent();
}

bool Zprime2muAnalysis::skipRecLevel(const int level) const {
  return
    (level == lGN && !useGen && !useSim) ||
    (!useTrigger && level >= lL1 && level <= lL3) ||
    (!useReco && level >= lGR) ||
    (!useOtherMuonRecos && level >= lFS && level <= lPR);
}

double Zprime2muAnalysis::resonanceMass(const reco::CompositeCandidate& dil) const {
  // Gen level muons don't have photons matched to them; we can take
  // the mass of the resonance directly from the MC record.
  if (recLevelHelper.recLevel(dil) == lGN) {
    if (hardInteraction.resonance)
      return hardInteraction.resonance->mass();
    else
      return dil.mass();
  }

  const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
  const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);

  reco::Particle::LorentzVector p4 = dil.p4(); // == lep0->p4() + lep1->p4() by construction.

  // If the leptons are electrons, or there were no photons in the
  // cone dR < 0.1, these four-vectors will be 0.
  reco::Particle::LorentzVector ph0 = recLevelHelper.closestPhoton(lep0);
  reco::Particle::LorentzVector ph1 = recLevelHelper.closestPhoton(lep1);

  p4 += ph0;
  // Don't double-count.
  if (ph1 != ph0) p4 += ph1;

  return p4.mass();
}

void Zprime2muAnalysis::dumpEvent(const bool trigOnly) const {
  unsigned imu, idil;
  int irec;
  ostringstream out;

  out << "\n******************************** Event " << eventNum
      << " (" << eventsDone << ")\n";

  int imax = trigOnly ? lL3+1 : MAX_LEVELS;
  for (irec = trigOnly ? lL1 : lGN; irec < imax; irec++) {
    if (irec >= lGR)
      out << endl;
    if (irec >= lOP)
      out << levelName(irec) << " cocktail muons:\n";
    for (imu = 0; imu < allLeptons[irec].size(); imu++)
      dumpLepton(out, allLeptons[irec][imu]);
    for (idil = 0; idil < allDileptons[irec].size(); idil++)
      dumpDilepton(out, allDileptons[irec][idil]);
    if (irec == lGN)
      out << endl;
  }
  
  out << "\nDilepton masses at levels  ";
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    if (i_rec >= lL1 && i_rec <= lL3) continue;
    out << " " << levelName(i_rec);
    if (i_rec < MAX_LEVELS-1) out << "     ";
  }

  out << setw(6) << setprecision(5);
  for (unsigned int i_dil = 0; i_dil < maxDileptons; i_dil++) {
    out << "\n Dilepton # " << i_dil << "; dimu mass: ";
    for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
      if (i_rec >= lL1 && i_rec <= lL3) continue;
      if (allDileptons[i_rec].size() > i_dil)
	out << allDileptons[i_rec][i_dil].mass();
      else
	out << " ---  ";
      if (i_rec < MAX_LEVELS-1) out << ", ";
    }
    out << "\n";

    out << "                res mass: ";
    for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
      if (i_rec >= lL1 && i_rec <= lL3)
	continue;
      if (allDileptons[i_rec].size() > i_dil)
	out << resonanceMass(allDileptons[i_rec][i_dil]);
      else
	out << " ---  ";
      if (i_rec < MAX_LEVELS-1) out << ", ";
    }
  }

  out << "\n********************************\n\n";
  edm::LogVerbatim("dumpEvent") << out.str();
}

void Zprime2muAnalysis::dumpLepton(ostream& output,
					   reco::CandidateBaseRef cand) const {
  // Make sure we're looking at the master ref (to handle shallow
  // clones which are daughters of dileptons).
  if (cand->hasMasterClone())
    cand = cand->masterClone().castTo<reco::CandidateBaseRef>();

  const int level = recLevelHelper.recLevel(cand);

  output << levelName(recLevelHelper.originalRecLevel(cand))
	 << " #"        << setw(1) << recLevelHelper.id(cand)
	 << " (-> GN #" << recLevelHelper.genMatchId(cand) << ")"
	 << "  Q: "     << setw(2) << cand->charge();

  if (level == lGN) {
    output << " Origin: " << setw(4) << motherId(cand)
	   << " Phi: "    << setw(7) << setprecision(4) << cand->phi()
	   << " Eta: "    << setw(7) << setprecision(4) << cand->eta()
	   << " Pt: "     << setw(7) << setprecision(5) << cand->pt()
	   << " P: "      << setw(7) << setprecision(5) << cand->p()
	   << endl;
  }
  else if (level == lL1) {
    int quality = -1;
    if (!doingElectrons) {
      const l1extra::L1MuonParticle& l1mu
	= toConcrete<l1extra::L1MuonParticle>(cand);
      quality = l1mu.gmtMuonCand().quality();
    }
    output << " Quality: "<< setw(2) << quality
	   << " Phi: "    << setw(7) << setprecision(4) << cand->phi()
	   << " Eta: "    << setw(7) << setprecision(4) << cand->eta()
	   << " Pt: "     << setw(7) << setprecision(5) << cand->pt()
	   << " P: "      << setw(7) << setprecision(5) << cand->p()
	   << endl;
  }
  else if (level == lL2 || level == lL3) {
    const reco::RecoCandidate& lep = toConcrete<reco::RecoCandidate>(cand);
    const reco::Track* tk = lep.track().get();
    output << " Phi: "      << setw(8) << setprecision(4) << cand->phi()
	   << "      Eta: " << setw(8) << setprecision(4) << cand->eta();
    if (tk && !doingElectrons) {
      output << endl
	     << "   Muhits: "    << setw(3) << tk->hitPattern().numberOfValidMuonHits()
	     << "   Chi2/Ndof: " << setw(8) << setprecision(4) << tk->chi2()
	     << "/"              << setw(2) << tk->ndof() << endl;
      output << "   P:     " << setw(7) << setprecision(5) << cand->p()
	   << " +/- " << pError(tk)
	     << "   "
	     << "   Pt:    " << setw(7) << setprecision(5) << cand->pt()
	     << " +/- " << ptError(tk)
	     << endl;
    }
    else
      output << "   Pt:    " << setw(7) << setprecision(5) << cand->pt()
	     << "   P:     " << setw(7) << setprecision(5) << cand->p()
	     << endl;
  }
  else {
    const reco::Track* globalTrack = 0;
    const reco::Track* trackerTrack = 0;
    const reco::Track* standAloneTrack = 0;
    double sumptr03 = -999;

    const reco::RecoCandidate& lep = toConcrete<reco::RecoCandidate>(cand);
    if (doingElectrons) {
      globalTrack = lep.gsfTrack().get();
    }
    else {
      globalTrack = lep.combinedMuon().get();
      trackerTrack = lep.track().get();
      standAloneTrack = lep.standAloneMuon().get();

      const reco::Muon* mu = toConcretePtr<reco::Muon>(cand);
      if (mu && mu->isIsolationValid()) sumptr03 = mu->isolationR03().sumPt;
    }
    
    output << " Phi: "      << setw(8) << setprecision(4) << cand->phi()
	   << "      Eta: " << setw(8) << setprecision(4) << cand->eta()
	   << endl;
 
    if (globalTrack) {
      const reco::HitPattern& hp = globalTrack->hitPattern();
      int px = hp.numberOfValidPixelHits();
      int si = hp.numberOfValidTrackerHits() - px;
      int pxl = hp.numberOfLostPixelHits();
      int sil = hp.numberOfLostTrackerHits() - pxl;
      int nh = hp.numberOfValidHits();
      int nhl = hp.numberOfLostHits();
      output << "   Pixhits: " << setw(3) << px + pxl 
	     << " Silhits: "   << setw(3) << si + sil
	     << " Rechits: "   << setw(3) << nh + nhl
	     << endl;
    }

    output << "   P:          " << setw(7) << setprecision(5) << cand->p();
    if (globalTrack != 0) 
      output << " +/- "    << pError(globalTrack);
    output << endl;

    output << "   Pt:         " << setw(7) << cand->pt();
    if (globalTrack != 0)
      output << " +/- "    << ptError(globalTrack)
	     << "   Chi2/Ndof: "  << setw(11) << globalTrack->chi2()
	     << "/"        << setw(2) << globalTrack->ndof();
    output << endl;

    double tkpt, tkpterr, tkchi2;
    if (trackerTrack == 0 || doingElectrons)
      tkpt = tkpterr = tkchi2 = -999;
    else {
      tkpt = trackerTrack->pt();
      tkpterr = ptError(trackerTrack);
      tkchi2 = trackerTrack->chi2();
    }
    output << "   Tracker Pt: " << setw(7) << tkpt << " +/- " << tkpterr
	   << "   Tracker Chi2: " << setw(8) << tkchi2 << endl;

    double stapt, stapterr, stachi2;
    if (standAloneTrack == 0 || doingElectrons)
      stapt = stapterr = stachi2 = -999;
    else {
      stapt = standAloneTrack->pt();
      stapterr = ptError(standAloneTrack);
      stachi2 = standAloneTrack->chi2();
    }
    output << "   MuonFit Pt: " << setw(7) << stapt << " +/- " << stapterr
	   << "   MuonFit Chi2: " << setw(8) << stachi2 << endl;

    if (globalTrack != 0) {
      output << "   Vertex position: " << setw(11) << cand->vx() 
	     << " "                    << setw(11) << cand->vy() 
	     << " "                    << setw(11) << cand->vz() << endl;
      if (globalTrack->extra().isAvailable())
	output << "   Track. position: " << setw(11) << globalTrack->innerPosition().X()
	       << " "                    << setw(11) << globalTrack->innerPosition().Y()
	       << " "                    << setw(11) << globalTrack->innerPosition().Z()
	       << endl;
      output << "   Track. position: " << setw(11) << globalTrack->vx()
	     << " "                    << setw(11) << globalTrack->vy()
	     << " "                    << setw(11) << globalTrack->vz()
	     << endl;
    }

    output << "   Closest photon: " << recLevelHelper.closestPhoton(cand)
	   << "   Sum pT (dR<0.3): " << sumptr03
	   << "   Cut code  " << cutHelper.leptonIsCut(*cand) << endl;

    output << "   Loc. code: " << whereIsLepton(cand, doingElectrons)
	   << "   p4 (p, E): " << cand->p4() << endl;

    if (!doingElectrons && (globalTrack != 0 && standAloneTrack != 0 && trackerTrack != 0)) {
      output << "   Combined track: charge: " << setw(2) << globalTrack->charge()
	     << " p: " << globalTrack->momentum() << endl
	     << "   Standalone mu : charge: " << setw(2) << standAloneTrack->charge()
	     << " p: " << standAloneTrack ->momentum() << endl
	     << "   Tracker track : charge: " << setw(2) << trackerTrack->charge()
	     << " p: " << trackerTrack->momentum() << endl;
    }
  }
}

void Zprime2muAnalysis::dumpDilepton(ostream& output,
					     const reco::CompositeCandidate& cand,
					     bool dumpLeptons) const {
  output << "Dilepton: charge: " << cand.charge()
	 << " pt: " << cand.pt() << " eta: " << cand.eta()
	 << " phi: " << cand.phi() << " mass: " << cand.mass() << endl;

  int larger = dileptonDaughter(cand, 0)->p() > dileptonDaughter(cand, 1)->p() ? 0 : 1;
  const reco::CandidateBaseRef& cand1 = dileptonDaughter(cand, larger);
  const reco::CandidateBaseRef& cand2 = dileptonDaughter(cand, !larger);

  if (dumpLeptons) {
    output << "Higher momentum daughter:\n";
    dumpLepton(output, cand1);
    output << "Lower momentum daughter:\n";
    dumpLepton(output, cand2);
  }
  else
    output << "  higher-p daughter: " << recLevelHelper.id(cand1)
	   << " pdgId: " << cand1->pdgId()
	   << " charge: " << cand1->charge() << " pt: " << cand1->pt()
	   << " eta: " << cand1->eta() << " phi: " << cand1->phi() << endl
	   << "   lower-p daughter: " << recLevelHelper.id(cand2)
	   << " pdgId: " << cand2->pdgId()
	   << " charge: " << cand2->charge() << " pt: " << cand2->pt()
	   << " eta: " << cand2->eta() << " phi: " << cand2->phi() << endl;
}

DEFINE_FWK_MODULE(Zprime2muAnalysis);

