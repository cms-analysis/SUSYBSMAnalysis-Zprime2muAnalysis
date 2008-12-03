#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muRecLevelAnalysis.h"

using namespace std;

Zprime2muRecLevelAnalysis::Zprime2muRecLevelAnalysis(const edm::ParameterSet& config)
  : Zprime2muAnalysis(config)
{
  recLevelHelper.init(config);
}

void Zprime2muRecLevelAnalysis::analyze(const edm::Event& event,
					const edm::EventSetup& eSetup) {
  Zprime2muAnalysis::analyze(event, eSetup);

  // Store per-rec-level information from the event (includes getting
  // all the match maps from the event).
  recLevelHelper.initEvent(event);

  // At each rec level (including the "best" leptons), store the
  // lepton and dilepton collections, including the "resonances"
  // (dileptons + closest photons added in).
  for (int irec = 0; irec < MAX_LEVELS; irec++) {
    recLevelHelper.getLeptons(event, irec, allLeptons[irec]);
    recLevelHelper.getDileptons(event, irec, RecLevelHelper::DIL,
				allDileptons[irec]);
    recLevelHelper.getDileptons(event, irec, RecLevelHelper::RES,
				allResonances[irec]);
  }

  // Dump the event if appropriate.
  if (verbosity >= VERBOSITY_SIMPLE)
    dumpEvent();
}

bool Zprime2muRecLevelAnalysis::skipRecLevel(const int level) const {
  return
    (level == lgen && !useGen && !useSim) ||
    (!useTrigger && level >= l1 && level <= l3) ||
    (!useReco && level >= lgmr) ||
    (!useOtherMuonRecos && level > lgmr && level < lbest);
}

void Zprime2muRecLevelAnalysis::dumpEvent(const bool trigOnly) const {
  unsigned imu, idil;
  int irec;
  ostringstream out;

  out << "\n******************************** Event " << eventNum
      << " (" << eventsDone << ")\n";

  int imax = trigOnly ? l3 : lbest;
  for (irec = trigOnly ? l1 : lgen; irec <= imax; irec++) {
    if (irec >= l3)
      out << endl;
    if (irec == lbest)
      out << "Best off-line muons: \n";
    for (imu = 0; imu < allLeptons[irec].size(); imu++)
      dumpLepton(out, allLeptons[irec][imu]);
    for (idil = 0; idil < allDileptons[irec].size(); idil++)
      dumpDilepton(out, allDileptons[irec][idil]);
  }
  
  out << "\nDilepton masses at levels  ";
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    out << " " << i_rec;
    if (i_rec < MAX_LEVELS-1) out << "      ";
  }

  out << setw(6) << setprecision(5);
  for (unsigned int i_dil = 0; i_dil < maxDileptons; i_dil++) {
    out << "\n Dilepton # " << i_dil << "; dimu mass: ";
    for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
      if (allDileptons[i_rec].size() > i_dil)
	out << allDileptons[i_rec][i_dil].mass();
      else
	out << " ---  ";
      if (i_rec < MAX_LEVELS-1) out << ", ";
    }
    out << "\n";

    out << "                res mass: ";
    for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
      if (allResonances[i_rec].size() > i_dil)
	out << allResonances[i_rec][i_dil].mass();
      else
	out << " ---  ";
      if (i_rec < MAX_LEVELS-1) out << ", ";
    }
  }

  LogTrace("dumpEvent") << out.str();
}

void Zprime2muRecLevelAnalysis::dumpLepton(ostream& output,
					   reco::CandidateBaseRef cand) const {
  // Make sure we're looking at the master ref (to handle shallow
  // clones which are daughters of dileptons).
  if (cand->hasMasterClone())
    cand = cand->masterClone().castTo<reco::CandidateBaseRef>();

  const int level = recLevelHelper.recLevel(cand);

  output << levelName(recLevelHelper.originalRecLevel(cand), true)
	 << " #"          << setw(1) << recLevelHelper.id(cand)
	 << "  Q: "       << setw(2) << cand->charge();

  if (level == lgen) {
    //const reco::GenParticleCandidate& genmu
    //  = toConcrete<reco::GenParticleCandidate>(cand);
    // JMTBAD fake motherline with pointer value
    output << " Origin: " << setw(4) << motherId(cand)
	   << "/"         << setw(4) << mother(cand)
	   << " Phi: "    << setw(7) << setprecision(4) << cand->phi()
	   << " Eta: "    << setw(7) << setprecision(4) << cand->eta()
	   << " Pt: "     << setw(7) << setprecision(5) << cand->pt()
	   << " P: "      << setw(7) << setprecision(5) << cand->p()
	   << endl;
  }
  else if (level == l1) {
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
  else if (level == l2) {
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
	//	   << " Seed: "      << setw(3) << seedIndex(cand)
	     << endl;
    }
    output << "   Closest  :";
    for (int i = 0; i < MAX_LEVELS-1; i++) {
      output << setw(5) << recLevelHelper.closestLeptonId(cand, i)
	     << "(" << levelName(i, true) << ")";
    }
    output << endl;
    output << "   Same-seed:";
    for (int i = 0; i < MAX_LEVELS-1; i++) {
      output << setw(5) << recLevelHelper.sameSeedLeptonId(cand, i)
	     << "(" << levelName(i, true) << ")";
    }
    output << endl;
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
    output << "   Forward Pt: " << setw(7) << 0 //rhs.forwardPt()
	   << " +/- "    << 0 //rhs.errForwardPt()
	   << "   Back Chi2: "  << setw(11) << 0 //rhs.backChi2()
	   << endl;
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
    output << "   MuonFit Pt: " << setw(7) << stapt
	   << " +/- " << stapterr
	   << "   MuonFit Chi2: " << setw(8) << stachi2 << endl;
    output << "   Tracker Chi2 diff.: " << setw(7) << 0 // rhs.trackerChi2Diff()
	   << "   MuonFit Chi2 diff.: " << setw(7) << 0 // rhs.muonFitChi2Diff()
	   << endl;
    if (globalTrack != 0) {
      output << "   Vertex position: " << setw(11) << cand->vx() 
	     << " "                    << setw(11) << cand->vy() 
	     << " "                    << setw(11) << cand->vz() << endl;
      if (globalTrack->extra().isAvailable())
	output << "   Track. position: " << setw(11) << globalTrack->innerPosition().X()
	       << " "                    << setw(11) << globalTrack->innerPosition().Y()
	       << " "                    << setw(11) << globalTrack->innerPosition().Z()
	       << endl;
    }

    if (!doingElectrons)
      output << "   Closest photon: " << recLevelHelper.closestPhoton(cand)
	     << "   Sum pT (dR<0.3): " << sumptr03
	     << "   Cut code  " << tevMuHelper.leptonIsCut(*cand) << endl;

    output << "   Loc. code: " << whereIsLepton(cand, doingElectrons)
	   << "   p4 (p, E): " << cand->p4() << endl;

    if (doingElectrons && level > l3) {
      const reco::GsfElectron& el =
	toConcrete<reco::GsfElectron>(cand);
      const reco::SuperClusterRef& sc = el.superCluster();
      output << "   Et: " << el.et() << "   Eta^{sc}: " << sc->eta()
	     << "   El. classification: " << el.classification() << endl
	     << "   Delta Eta_in: " << el.deltaEtaSuperClusterTrackAtVtx()
	     << "   Delta Phi_in: " << el.deltaPhiSuperClusterTrackAtVtx() 
	     << "   HoE: " << el.hadronicOverEm() << endl;
    }

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

void Zprime2muRecLevelAnalysis::dumpDilepton(ostream& output,
					     const reco::CompositeCandidate& cand,
					     bool dumpLeptons) const {
  output << "Dilepton: charge: " << cand.charge()
	 << " pt: " << cand.pt() << " eta: " << cand.eta()
	 << " phi: " << cand.phi() << " mass: " << cand.mass() << endl;

  int larger = dileptonDaughter(cand, 0)->p() > dileptonDaughter(cand, 1)->p() ? 0 : 1;
  int smaller = larger == 0 ? 1 : 0;
  const reco::CandidateBaseRef& cand1 = dileptonDaughter(cand, larger);
  const reco::CandidateBaseRef& cand2 = dileptonDaughter(cand, smaller);

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

DEFINE_FWK_MODULE(Zprime2muRecLevelAnalysis);
