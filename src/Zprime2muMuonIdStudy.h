#ifndef ZPRIME2MUMUIDSTUDY_H
#define ZPRIME2MUMUIDSTUDY_H

#include "TH1F.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muRecLevelAnalysis.h"

class Zprime2muMuonIdStudy : public Zprime2muRecLevelAnalysis {
 public:
  explicit Zprime2muMuonIdStudy(const edm::ParameterSet&);
  virtual void analyze(const edm::Event& ev, const edm::EventSetup& setup);
  void endJob();

 private:
  void bookHistos();
  TH1F* makeHist(const char* baseName, int irec,
		 int nbins, double min, double max) const;
  void drawHistos() const;

  TH1F *h_OK, *h_eOK, *h_isoOK, *h_timeOK, *h_matchOK;

  TH1F* h_calE_em   [2];
  TH1F* h_calE_emS9 [2];
  TH1F* h_calE_had  [2];
  TH1F* h_calE_hadS9[2];
  TH1F* h_calE_ho   [2];
  TH1F* h_calE_hoS9 [2];

  TH1F* h_iso03_emEt   [2];
  TH1F* h_iso03_hadEt  [2];
  TH1F* h_iso03_hoEt   [2];
  TH1F* h_iso03_nJets  [2];
  TH1F* h_iso03_nTracks[2];
  TH1F* h_iso03_sumPt  [2];

  TH1F* h_iso05_emEt   [2];
  TH1F* h_iso05_hadEt  [2];
  TH1F* h_iso05_hoEt   [2];
  TH1F* h_iso05_nJets  [2];
  TH1F* h_iso05_nTracks[2];
  TH1F* h_iso05_sumPt  [2];

  TH1F* h_time_freeInverseBeta   [2];
  TH1F* h_time_freeInverseBetaErr[2];
  TH1F* h_time_inverseBeta       [2];
  TH1F* h_time_inverseBetaErr    [2];
  TH1F* h_time_nStations         [2];
  TH1F* h_time_timeAtIpInOut     [2];
  TH1F* h_time_timeAtIpInOutErr  [2];
  TH1F* h_time_timeAtIpOutIn     [2];
  TH1F* h_time_timeAtIpOutInErr  [2];
};

#endif
