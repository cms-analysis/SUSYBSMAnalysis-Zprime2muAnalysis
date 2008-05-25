#ifndef ZPRIME2MUBACKGROUNDS_H
#define ZPRIME2MUBACKGROUNDS_H

#include "TH1F.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muRecLevelAnalysis.h"

class Zprime2muMatchStudy : public Zprime2muRecLevelAnalysis {
 public:
  explicit Zprime2muMatchStudy(const edm::ParameterSet&);
  virtual void analyze(const edm::Event& ev, const edm::EventSetup& setup);
  void endJob();

 private:
  void bookHistos();
  TH1F* makeHist(const char* baseName, int type, int irec, int jrec,
		 const char* title, int nbins, double min, double max) const;
  void drawHistos() const;

  TH1F* hDeltaR[2][MAX_LEVELS][MAX_LEVELS-1];
  TH1F* hDeltaPt[2][MAX_LEVELS][MAX_LEVELS-1];
  TH1F* hDeltaQ[2][MAX_LEVELS][MAX_LEVELS-1];
};

#endif
