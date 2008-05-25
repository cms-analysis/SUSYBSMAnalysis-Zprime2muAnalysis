#ifndef ZP2MURESOLUTION_H
#define ZP2MURESOLUTION_H

#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muRecLevelAnalysis.h"

class Zprime2muResolution : public Zprime2muRecLevelAnalysis {
 public:
  explicit Zprime2muResolution(const edm::ParameterSet&);
  ~Zprime2muResolution();
  void beginJob(const edm::EventSetup& eSetup);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

 private:
  void BookResHistos();
  void BookEffHistos();
  void BookPtResHistos();
  void BookDilResHistos();
  void BookChargeResHistos();

  void calcResolution(const bool debug = false);

  void fillEffHistos();
  void fillPtResHistos(const bool debug = false);
  void fillChargeResHistos(const bool debug = false);
  void fillMuonHistos(const int rec, const bool debug = false);
  void fillDilResHistos(const bool debug = false);
  void fillSignOfDilepton(const int rec, const bool debug = false);

  int getOrigin(const int motherId);

  void getHistosFromFile();

  void DrawResHistos();
  void DrawAcceptance();

  // roughly sorted by order booked, and in what booking fcn
  // booked in BookResHistos()
  TH1F *GenMassAllEvents, *GenMassInAccept;
  TH1F *Origin[2];
  TH1F *TrigResult[TRIG_LEVELS][3], *TrigMass[TRIG_LEVELS][3];
  TH1F *EventsInAccFailed, *L1TrigPassSingleMu, *L1TrigFailSingleMu;
  TH2F *L1TrigFailMu2VsMu1, *L1TrigPassMu2VsMu1;
  TH1F *L2MuonHits, *L3TrackerHits, *GMRMuonHits[3], *GMRChi2dof[3];
  TH1F *SumPtR03[MAX_LEVELS][2];
  TH1F *NMuons[TRIG_LEVELS][4], *NumDilVsRec, *SignOfDil[TRIG_LEVELS];
  TH1F *MuonEta[MAX_LEVELS], *MuonRap[MAX_LEVELS], *MuonPhi[MAX_LEVELS];
  TH1F *MuonPt[MAX_LEVELS], *MuonPz[MAX_LEVELS], *MuonP[MAX_LEVELS];
  TH2F *MuonPVsEta[MAX_LEVELS], *MuonPtVsEta[MAX_LEVELS];
  TH1F *AllDilMass[MAX_LEVELS], *DilMass[MAX_LEVELS];
  TH2F *DilMassVsEta[MAX_LEVELS], *DilMassVsY[MAX_LEVELS];
  TH2F *MuMVsMuP[MAX_LEVELS][3];
  TH1F *Phi[MAX_LEVELS][3], *Eta[MAX_LEVELS][3], *Rapidity[MAX_LEVELS][3];
  TProfile *PVsEta[MAX_LEVELS][3], *PtVsEta[MAX_LEVELS][3];
  TH1F *Pt[MAX_LEVELS][3], *Pz[MAX_LEVELS][3], *P[MAX_LEVELS][3];
  TH1F *ZonDilMass[TRIG_LEVELS], *ZofDilMass[TRIG_LEVELS];
  TH1F *Eta4muons, *Pt4muons;
  TH1F *EtaRes[MAX_LEVELS], *PhiRes[MAX_LEVELS], *PtDiff[MAX_LEVELS];
  TH1F *GenInvPtRes[3], *PRes[MAX_LEVELS], *GenInvPRes[3];
  TProfile *GenPhiResVsPhi[3], *PResVsP[MAX_LEVELS], *GenPResVsPt[3];
  TProfile *GenInvPResVsPt[3], *GenInvPtResVsPt[3];
  TH2F *GenEtaResScat[3], *GenPhiResScat[3], *GenPtResScat[3];
  TH1F *AllDilMassRes, *GenDilMassRes[3], *GenDilMassFrRes[3];
  TProfile *GenDilMassResScat[3], *GenDilMassFrResScat[3];
  TH1F *L1EtaRes[2], *L1PhiRes[2], *L1PtDiff[2];
  TH2F *L1EtaResScat[2], *L1PhiResScat[2], *L1PtResScat[2];
  TH1F *L2EtaRes, *L2PhiRes, *L2PtDiff;
  TH2F *L2EtaResScat, *L2PhiResScat, *L2PtResScat;
  // booked in BookEffHistos()
  TH1F *EffVsEta[MAX_LEVELS], *EffVsPhi[MAX_LEVELS], *EffVsPt[MAX_LEVELS];
  TH1F *RecMass[4][3];
  TH1F *RecMassByLoc[2][10];
  // booked in BookPtResHistos()
  TProfile *ResidualPt[6];
  TH1F *MuInvPRes[2][3], *MuInvPtRes[2][3], *MuPtDiff[2];
  TH1F *TotInvPtRes[4], *TotInvPtPull[4];
  TH1F *InvPtRes[4][2], *InvPtPull[4][2];
  // booked in BookDilResHistos()
  TH1F *DilMassComp[MAX_LEVELS][2];
  TH1F *DilMassRes[MAX_LEVELS+1][6], *DilPtRes[MAX_LEVELS];
  TH2F *MuPVsMuM[2];
  // booked in BookChargeResHistos()
  TH1F *QRes[MAX_LEVELS];
  TH1F *QResVsPt[3][2], *QResVsInvPt[3][2];
  TH1F *QResVsP[3][2],  *QResVsInvP[3][2];

  // config file parameters
  double peakMass;
  double lowerMassWin;
  double upperMassWin;
  int binSize;
  double maxTrigMass;

  std::string outputFile;
  TFile* histoFile;
  bool useHistosFromFile;
};

#endif // ZP2MURESOLUTION_H
