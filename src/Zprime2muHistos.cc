/**
  \class    Zprime2muHistos
  \brief    Plots basic lepton and dilepton quantities for each rec level.

  \author   Jordan Tucker, Slava Valuev
  \version  $Id$
*/

#include "TString.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muHistos.h"

using namespace std;

Zprime2muHistos::Zprime2muHistos(const edm::ParameterSet& config) : Zprime2muAnalysis(config) {
  // Get the parameters specific to the data sample on which we are running.
  string dataSet = config.getParameter<string>("dataSet");
  edm::ParameterSet dataSetConfig = config.getParameter<edm::ParameterSet>(dataSet);

  peakMass     = dataSetConfig.getParameter<double>("peakMass");
  lowerMassWin = dataSetConfig.getParameter<double>("lowerMassWin");
  upperMassWin = dataSetConfig.getParameter<double>("upperMassWin");
  binSize      = dataSetConfig.getParameter<int>   ("binSize");
  maxTrigMass  = dataSetConfig.getParameter<double>("maxTrigMass");

  bookBasicLeptonHistos();
  bookOfflineLeptonHistos();
  bookDileptonHistos();
}

const double pi = TMath::Pi();

TString makeName(const char* name, Long_t i) {
  return TString(name) + i;
}

void Zprime2muHistos::bookBasicLeptonHistos() {
  // Basic kinematics at all levels (gen, L1/HLT, offline reco)..
  for (int rec = 0; rec < MAX_LEVELS; ++rec) {
    TString level(levelName(rec, true));

    // Lepton multiplicity.
    NLeptons[rec] = fs->make<TH1F>(makeName("NLeptons", rec), "#/event, " + level, 10, 0, 10);

    // Lepton eta, y, phi.
    LeptonEta[rec] = fs->make<TH1F>(makeName("LeptonEta", rec), level + " #eta", 100,  -5,  5);
    LeptonRap[rec] = fs->make<TH1F>(makeName("LeptonRap", rec), level + " y",    100,  -5,  5);
    LeptonPhi[rec] = fs->make<TH1F>(makeName("LeptonPhi", rec), level + " #phi", 100, -pi, pi);
    
    // Lepton momenta: p, p_T, p_z.
    // Special binning for momenta for L1 leptons.
    bool isL1 = rec == l1;

    LeptonPt[rec] = fs->make<TH1F>(makeName("LeptonPt", rec), level + " p_{T}", 100,  0, isL1 ? 150 :   peakMass);
    LeptonPz[rec] = fs->make<TH1F>(makeName("LeptonPz", rec), level + " p_{z}", 100,  0, isL1 ? 800 : 2*peakMass);
    LeptonP[rec]  = fs->make<TH1F>(makeName("LeptonP",  rec), level + " p",     100,  0, isL1 ? 800 : 2*peakMass);

    // Lepton momenta versus pseudorapidity.
    LeptonPVsEta[rec]  = fs->make<TProfile>(makeName("LeptonPVsEta",  rec), level + " p vs. #eta",     100, -6, 6);
    LeptonPtVsEta[rec] = fs->make<TProfile>(makeName("LeptonPtVsEta", rec), level + " p_{T} vs. #eta", 100, -6, 6);
  }
}

void Zprime2muHistos::bookOfflineLeptonHistos() {
  // Offline reconstructed quantities.
  for (int rec = 0; rec < MAX_LEVELS; ++rec) {
    TString level(levelName(rec, true));

    // Delta R < 0.3 isolation variables.
    IsoSumPt[rec]   = fs->make<TH1F>(makeName("IsoSumPt",   rec), level + " Iso. (#Delta R < 0.3) #Sigma p_{T}", 30, 0, 30);
    IsoEcal[rec]    = fs->make<TH1F>(makeName("IsoEcal",    rec), level + " Iso. (#Delta R < 0.3) ECAL",         30, 0, 30);
    IsoHcal[rec]    = fs->make<TH1F>(makeName("IsoHcal",    rec), level + " Iso. (#Delta R < 0.3) HCAL",         30, 0, 30);
    IsoNTracks[rec] = fs->make<TH1F>(makeName("IsoNTracks", rec), level + " Iso. (#Delta R < 0.3) nTracks",      10, 0, 10);
    IsoNJets[rec]   = fs->make<TH1F>(makeName("IsoNJets",   rec), level + " Iso. (#Delta R < 0.3) nJets",        10, 0, 10);
    
    // Track hit counts.
    NPxHits[rec] = fs->make<TH1F>(makeName("NPxHits", rec), level + " # pixel hits",    8, 0,  8);
    NStHits[rec] = fs->make<TH1F>(makeName("NStHits", rec), level + " # strip hits",   26, 0, 26);
    NTkHits[rec] = fs->make<TH1F>(makeName("NTkHits", rec), level + " # tracker hits", 33, 0, 33);
    NMuHits[rec] = fs->make<TH1F>(makeName("NMuHits", rec), level + " # muon hits",    51, 0, 51);
    NHits[rec]   = fs->make<TH1F>(makeName("NHits",   rec), level + " # hits",         78, 0, 78);

    // Other track variables.
    Chi2dof[rec] = fs->make<TH1F>(makeName("Chi2dof", rec), level + " #chi^{2}/dof", 100, 0, 5);
    TrackD0[rec] = fs->make<TH1F>(makeName("TrackD0", rec), level + " |d_{0}|",      100, 0, 0.1);
    TrackDz[rec] = fs->make<TH1F>(makeName("TrackDz", rec), level + " |dz|",         100, 0, 5);
  }
}

void Zprime2muHistos::bookDileptonHistos() {
  // Dilepton quantities.
  for (int rec = lgen; rec < MAX_LEVELS; ++rec) {
    // No trigger-level dileptons.
    if (rec >= l1 && rec <= l3) continue;

    TString level(levelName(rec, true));

    // Dilepton multiplicity.
    NDileptons[rec] = fs->make<TH1F>(makeName("NDileptons", rec), "# dileptons/event, " + level, maxDileptons + 1, 0, maxDileptons + 1);

    // Dilepton eta, y, phi.
    DileptonEta[rec] = fs->make<TH1F>(makeName("DileptonEta", rec), level + " dil. #eta", 100,  -5,  5);
    DileptonRap[rec] = fs->make<TH1F>(makeName("DileptonRap", rec), level + " dil. y",    100,  -5,  5);
    DileptonPhi[rec] = fs->make<TH1F>(makeName("DileptonPhi", rec), level + " dil. #phi", 100, -pi, pi);
    
    // Dilepton momenta: p, p_T, p_z.
    DileptonPt[rec] = fs->make<TH1F>(makeName("DileptonPt", rec), level + " dil. p_{T}", 100,  0,        0.5*peakMass);
    DileptonPz[rec] = fs->make<TH1F>(makeName("DileptonPz", rec), level + " dil. p_{z}", 100,  0, 1000 + upperMassWin);
    DileptonP[rec]  = fs->make<TH1F>(makeName("DileptonP",  rec), level + " dil. p",     100,  0, 1000 + upperMassWin);

    // Dilepton momenta versus pseudorapidity.
    DileptonPVsEta[rec]  = fs->make<TProfile>(makeName("DileptonPVsEta",  rec), level + " dil. p vs. #eta",     100, -6, 6);
    DileptonPtVsEta[rec] = fs->make<TProfile>(makeName("DileptonPtVsEta", rec), level + " dil. p_{T} vs. #eta", 100, -6, 6);

    // Dilepton invariant mass.
    DileptonMass[rec]            = fs->make<TH1F>(makeName("DileptonMass",            rec), level + " dil. mass",       binSize, lowerMassWin, upperMassWin);
    DileptonWithPhotonsMass[rec] = fs->make<TH1F>(makeName("DileptonWithPhotonsMass", rec), level + " dil+#gamma mass", binSize, lowerMassWin, upperMassWin);

    // Opposite/like-sign counts.
    DileptonSigns[rec] = fs->make<TH1F>(makeName("DileptonSigns", rec), level + " dil. sign", 3, 0, 3);

    // Plots comparing the daughter lepton momenta.
    DileptonDeltaPt[rec]  = fs->make<TH1F>(makeName("DileptonDeltaPt",  rec), level + " dil. |p_{T}^{1}| - |p_{T}^{2}|",                100, 0, 100);
    DileptonDeltaP[rec]   = fs->make<TH1F>(makeName("DileptonDeltaP",   rec), level + " dil. |p^{1}| - |p^{2}|",                        100, 0, 100);
    DileptonPtErrors[rec] = fs->make<TH2F>(makeName("DileptonPtErrors", rec), level + " dil. #sigma_{p_{T}}^{1} v. #sigma_{p_{T}}^{2}", 100, 0, 100, 100, 0, 100);
  }
}

void Zprime2muHistos::fillBasicLeptonHistos(const reco::CandidateBaseRef& lep, const int rec) {
  LeptonEta[rec]->Fill(lep->eta());
  LeptonRap[rec]->Fill(lep->rapidity());
  LeptonPhi[rec]->Fill(lep->phi());

  LeptonPt[rec]->Fill(lep->pt());
  LeptonPz[rec]->Fill(fabs(lep->pz()));
  LeptonP[rec] ->Fill(lep->p());

  LeptonPtVsEta[rec]->Fill(lep->eta(), lep->pt());
  LeptonPVsEta[rec] ->Fill(lep->eta(), lep->p());
}

void Zprime2muHistos::fillOfflineMuonHistos(const reco::Muon* lep, const int rec) {
  const reco::MuonIsolation& iso = lep->isolationR03();
  IsoSumPt[rec]  ->Fill(iso.sumPt);
  IsoEcal[rec]   ->Fill(iso.emEt);
  IsoHcal[rec]   ->Fill(iso.hadEt + iso.hoEt);
  IsoNTracks[rec]->Fill(iso.nTracks);
  IsoNJets[rec]  ->Fill(iso.nJets);

  const reco::TrackRef& globalTrack = lep->globalTrack();
  if (globalTrack.isAvailable()) {
    Chi2dof[rec]->Fill(globalTrack->normalizedChi2());
    TrackD0[rec]->Fill(fabs(globalTrack->d0()));
    TrackDz[rec]->Fill(fabs(globalTrack->dz()));

    const reco::HitPattern& hp = globalTrack->hitPattern();
    NPxHits[rec]->Fill(hp.numberOfValidPixelHits());
    NStHits[rec]->Fill(hp.numberOfValidStripHits());
    NTkHits[rec]->Fill(hp.numberOfValidTrackerHits());
    NMuHits[rec]->Fill(hp.numberOfValidMuonHits());
    NHits[rec]  ->Fill(hp.numberOfValidHits());
  }
}

void Zprime2muHistos::fillOfflineElectronHistos(const reco::GsfElectron* lep, const int rec) {
  // Can add electron quantities here.
}

void Zprime2muHistos::fillDileptonSigns(const int rec) {
  int totalCharge = 0;
  for (reco::CandidateBaseRefVector::const_iterator lep = allLeptons[rec].begin(); lep != allLeptons[rec].end(); ++lep)
    totalCharge += (*lep)->charge();
  
  // If total charge is equal to number of leptons, then positive
  // dilepton could be found. Else if total charge is equal to
  // -(number of leptons), then negative dilepton could be found.
  // Otherwise opposite sign dilepton could be found.
  int nLeptons = int(allLeptons[rec].size());
  if (totalCharge == nLeptons)
    DileptonSigns[rec]->Fill(2);
  else if (totalCharge == -nLeptons)
    DileptonSigns[rec]->Fill(1);
  else
    DileptonSigns[rec]->Fill(0);
}

void Zprime2muHistos::fillBasicDileptonHistos(const reco::CompositeCandidate& dil, const int rec) {
  DileptonEta[rec]->Fill(dil.eta());
  DileptonRap[rec]->Fill(dil.rapidity());
  DileptonPhi[rec]->Fill(dil.phi());

  DileptonPt[rec]->Fill(dil.pt());
  DileptonPz[rec]->Fill(fabs(dil.pz()));
  DileptonP[rec] ->Fill(dil.p());

  DileptonPtVsEta[rec]->Fill(dil.eta(), dil.pt());
  DileptonPVsEta[rec] ->Fill(dil.eta(), dil.p());

  DileptonMass[rec]->Fill(dil.mass());
  DileptonWithPhotonsMass[rec]->Fill(resonanceMass(dil));
}

void Zprime2muHistos::fillDileptonDaughterHistos(const reco::CompositeCandidate& dil, const int rec) {
  const reco::CandidateBaseRef& lep_minus = dileptonDaughterByCharge(dil, -1);
  const reco::CandidateBaseRef& lep_plus  = dileptonDaughterByCharge(dil, +1);

  DileptonDeltaPt[rec]->Fill(fabs(lep_minus->pt()) - fabs(lep_plus->pt()));
  DileptonDeltaP[rec] ->Fill(fabs(lep_minus->p())  - fabs(lep_plus->p()));

  const reco::Muon* mu_minus = toConcretePtr<reco::Muon>(lep_minus);
  const reco::Muon* mu_plus  = toConcretePtr<reco::Muon>(lep_plus);
  if (mu_minus && mu_plus) {
    const reco::Track* tk_minus = mu_minus->globalTrack().get();
    const reco::Track* tk_plus  = mu_plus ->globalTrack().get();
    if (tk_minus && tk_plus) {
      DileptonPtErrors[rec]->Fill(ptError(tk_minus), ptError(tk_plus));
    }
  }
}

void Zprime2muHistos::fillLeptonHistos(const int rec) {
  // Leptons.
  NLeptons[rec]->Fill(allLeptons[rec].size());

  for (reco::CandidateBaseRefVector::const_iterator lep = allLeptons[rec].begin(); lep != allLeptons[rec].end(); ++lep) {
    fillBasicLeptonHistos(*lep, rec);
  
    const reco::Muon* muon = toConcretePtr<reco::Muon>(*lep);
    if (muon) fillOfflineMuonHistos(muon, rec);

    const reco::GsfElectron* electron = toConcretePtr<reco::GsfElectron>(*lep);
    if (electron) fillOfflineElectronHistos(electron, rec);
  }
}

void Zprime2muHistos::fillDileptonHistos(const int rec) {
  if (allLeptons[rec].size() >= 2)
    fillDileptonSigns(rec);

  NDileptons[rec]->Fill(allDileptons[rec].size());

  unsigned j = 0;
  for (reco::CompositeCandidateCollection::const_iterator dil = allDileptons[rec].begin(); dil != allDileptons[rec].end(); ++dil, ++j) {
    fillBasicDileptonHistos(*dil, rec);
    fillDileptonDaughterHistos(*dil, rec);
  }
}

void Zprime2muHistos::analyze(const edm::Event& event, const edm::EventSetup& eSetup) {
  // Delegate filling our lepton vectors to the parent class.
  Zprime2muAnalysis::analyze(event, eSetup);

  for (int rec = 0; rec < MAX_LEVELS; ++rec) {
    fillLeptonHistos(rec);
    if (rec == lgen || rec > l3)
      fillDileptonHistos(rec);
  }
}

DEFINE_FWK_MODULE(Zprime2muHistos);
