#ifndef ZP2MURESOLUTION_H
#define ZP2MURESOLUTION_H

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TString.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

class Zprime2muResolution : public Zprime2muAnalysis {
 public:
  explicit Zprime2muResolution(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void bookGenLevelHistos();
  void bookEfficiencyHistos();
  void bookLeptonResolutionHistos();
  void bookChargeResolutionHistos();
  void bookDileptonResolutionHistos();

  int encodeLeptonOrigin(const int id) const;

  void fillGenLevelHistos();

  void fillLeptonEfficiencyHistos();
  void fillTriggerEfficiencyHistos();
  void fillDileptonEfficiencyHistos();

  void fillLeptonResolution(const reco::CandidateBaseRef& gen_lep, const reco::CandidateBaseRef& lep, const int rec);
  void fillLeptonExtraMomentumResolution(const reco::CandidateBaseRef& gen_lep, const reco::CandidateBaseRef& lep, const int rec);
  void fillChargeResolution(const reco::CandidateBaseRef& gen_lep, const reco::CandidateBaseRef& lep, const int rec);
  void fillDileptonMassResolution(const reco::CompositeCandidate& gen_dil, const reco::CompositeCandidate& dil, const int rec);

  void fillLeptonHistos(const reco::CandidateBaseRef& lep, const int rec);
  void fillLeptonHistos(const int rec);
  void fillDileptonHistos(const int rec);

  // Parameters specified in the config file.
  bool   leptonsFromDileptons;
  double peakMass;
  double lowerMassWin;
  double upperMassWin;
  int    binSize;
  double maxTrigMass;

  TH1F* LeptonOrigin[2];

  TH1F* GenMassAllEvents;
  TH1F* GenMassInAccept;

  TH1F* EffVsEta[MAX_LEVELS];
  TH1F* EffVsPhi[MAX_LEVELS];
  TH1F* EffVsPt[MAX_LEVELS];

  TH1F* TrigEffVsDilMass[TRIG_LEVELS][3];
  TH1F* DilRecEffVsMass[MAX_LEVELS][3];

  TH1F* LeptonEtaDiff[MAX_LEVELS];
  TH1F* LeptonPhiDiff[MAX_LEVELS];

  TH1F* LeptonPtDiff[MAX_LEVELS];

  TH1F* LeptonPtRes[MAX_LEVELS];
  TH1F* LeptonPRes[MAX_LEVELS];

  TH1F* LeptonInvPtRes[MAX_LEVELS];
  TH1F* LeptonInvPRes[MAX_LEVELS];

  TProfile* LeptonInvPtResVPtGen[MAX_LEVELS];
  TProfile* LeptonInvPResVPGen[MAX_LEVELS];
  
  TH1F* LeptonInvPtPull[MAX_LEVELS];
  TH1F* LeptonInvPPull[MAX_LEVELS];
  
  TH1F* LeptonInvPtResBarrel[MAX_LEVELS];
  TH1F* LeptonInvPResBarrel[MAX_LEVELS];
    
  TH1F* LeptonInvPtPullBarrel[MAX_LEVELS];
  TH1F* LeptonInvPPullBarrel[MAX_LEVELS];

  TH1F* LeptonInvPtResEndcap[MAX_LEVELS];
  TH1F* LeptonInvPResEndcap[MAX_LEVELS];
    
  TH1F* LeptonInvPtPullEndcap[MAX_LEVELS];
  TH1F* LeptonInvPPullEndcap[MAX_LEVELS];

  TH1F* ChargeDiff[MAX_LEVELS];

  TH1F* ChargeRightVInvPt[MAX_LEVELS];
  TH1F* ChargeWrongVInvPt[MAX_LEVELS];

  TH1F* DileptonMassRes[MAX_LEVELS];
  TH1F* DileptonResMassRes[MAX_LEVELS];
  TH1F* ResonanceMassRes[MAX_LEVELS];

  TProfile* DileptonMassResVMass[MAX_LEVELS];
  TProfile* DileptonResMassResVMass[MAX_LEVELS];
  TProfile* ResonanceMassResVMass[MAX_LEVELS];
};

#endif // ZP2MURESOLUTION_H
