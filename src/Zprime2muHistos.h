#ifndef ZP2MUHISTOS_H
#define ZP2MUHISTOS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

class Zprime2muHistos : public Zprime2muAnalysis {
 public:
  explicit Zprime2muHistos(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void bookTriggerHistos();
  void bookBasicLeptonHistos();
  void bookOfflineLeptonHistos();
  void bookDileptonHistos();

  void fillTriggerHistos();
  void fillBasicLeptonHistos(const reco::CandidateBaseRef& lep, const int rec);
  void fillOfflineMuonHistos(const reco::Muon* lep, const int rec);
  void fillOfflineElectronHistos(const reco::GsfElectron* lep, const int rec);

  void fillDileptonSigns(const int rec);
  void fillBasicDileptonHistos(const reco::CompositeCandidate& dil, const int rec);
  void fillDileptonDaughterHistos(const reco::CompositeCandidate& dil, const int rec);

  void fillLeptonHistos(const int rec);
  void fillDileptonHistos(const int rec);

  // Parameters specified in the config file.
  double peakMass;
  double lowerMassWin;
  double upperMassWin;
  int    binSize;
  double maxTrigMass;

  // Trigger quantity histograms; level gen will not be booked and is
  // a placeholder.
  TH1F* TriggerBits[TRIG_LEVELS];
  TH1F* NLeptonsTriggered[TRIG_LEVELS];
  TH1F* NLeptonsFailed[TRIG_LEVELS];

  // Histograms of basic quantities for leptons, booked for all
  // levels.
  TH1F* NLeptons[MAX_LEVELS];
  TH1F* LeptonEta[MAX_LEVELS];
  TH1F* LeptonRap[MAX_LEVELS];
  TH1F* LeptonPhi[MAX_LEVELS];
  TH1F* LeptonPt[MAX_LEVELS];
  TH1F* LeptonPz[MAX_LEVELS];
  TH1F* LeptonP[MAX_LEVELS];
  TProfile* LeptonPVsEta[MAX_LEVELS];
  TProfile* LeptonPtVsEta[MAX_LEVELS];

  // For these offline quantities, histograms for levels gen-L3 will
  // not be booked.
  TH1F* IsoSumPt[MAX_LEVELS];
  TH1F* IsoEcal[MAX_LEVELS];   
  TH1F* IsoHcal[MAX_LEVELS];   
  TH1F* IsoNTracks[MAX_LEVELS];
  TH1F* IsoNJets[MAX_LEVELS];  
  TH1F* NPxHits[MAX_LEVELS];
  TH1F* NStHits[MAX_LEVELS];
  TH1F* NTkHits[MAX_LEVELS];
  TH1F* NMuHits[MAX_LEVELS];
  TH1F* NHits[MAX_LEVELS];  
  TH1F* Chi2dof[MAX_LEVELS];
  TH1F* TrackD0[MAX_LEVELS];
  TH1F* TrackDz[MAX_LEVELS];

  // For the dileptons, histograms for levels L1-L3 won't be booked
  // (but lGN will).
  TH1F* NDileptons[MAX_LEVELS];
  TH1F* DileptonEta[MAX_LEVELS];
  TH1F* DileptonRap[MAX_LEVELS];
  TH1F* DileptonPhi[MAX_LEVELS];
  TH1F* DileptonPt[MAX_LEVELS];
  TH1F* DileptonPz[MAX_LEVELS];
  TH1F* DileptonP[MAX_LEVELS];
  TProfile* DileptonPVsEta[MAX_LEVELS];
  TProfile* DileptonPtVsEta[MAX_LEVELS];
  TH1F* DileptonMass[MAX_LEVELS];
  TH1F* DileptonWithPhotonsMass[MAX_LEVELS];
  TH1F* DileptonSigns[MAX_LEVELS];
  TH1F* DileptonDeltaPt[MAX_LEVELS];
  TH1F* DileptonDeltaP[MAX_LEVELS];
  TH2F* DileptonPtErrors[MAX_LEVELS];
};

#endif // ZP2MUHISTOS_H
