#ifndef ZP2MUASYM_H
#define ZP2MUASYM_H

#include <string>

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFitData.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TFMultiD.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

const int FIT_ARRAY_SIZE = 50000;// max size of arrays for unbinned fits. 

class Zprime2muAsymmetry : public Zprime2muAnalysis {
 public:
  explicit Zprime2muAsymmetry(const edm::ParameterSet&);

  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

 private:
  static const double nSigma;
  static const unsigned int nNormPoints;
  static const unsigned int nSmearPoints;

  void bookFitHistos();

  void makeFakeData(int nEvents, double A_FB, double b);
  void smearGenData(const AsymFitData& data);

  void fillFitData(const edm::Event& event);
  void dumpFitData();

  void countAsymmetry(const int which, double& Afb, double& eAfb);
  void asymmetryByEventWeighting(const int which, double& Afb, double& eAfb);
  void asymmetryByEventWeightingWithAngularInfo(const int which, double& Afb, double& eAfb, bool use_h=true, bool bin_x=true);

  TFMultiD* getNewFitFcn(int which);
  void fitAsymmetry();

  void drawFitHistos();

  TH1F* h_genCosNoCut;
  TH2F* AsymFitSmearHisto[6];
  TH1F* AsymFitHistoGen[6];
  TH1F* AsymFitHistoRec[6];
  TH1F* AsymFitSmearHistoDif[6];
  TH1F* AsymFitHistoGenSmeared[6];
  TH1F* AsymFitHistoGenByType[2][6];
  TH1F* AsymFitHistoRecByType[2][6];

  TH1F* h_cos_theta_true;
  TH1F* h_cos_theta_cs_acc;
  TH1F* h_cos_theta_cs_acc_fixed;
  TH1F* h_cos_theta_cs;
  TH1F* h_cos_theta_cs_fixed;
  TH1F* h_cos_theta_cs_rec;
  TH1F* h_b_mass;
  TH1F* h_f_mass;
  TH1F* h_b_smass;
  TH1F* h_f_smass;
  TH1F* h_gen_sig[2];
  TH2F* h2_rap_cos_d[2];
  TH2F* h2_rap_cos_d_uncut[2];
  TH2F* h2_rap_cos_d_rec;
  TH1F* mistagProbEvents[3];

  // Data arrays for unbinned fits.  Fixed size arrays for now.
  // Order of the arrays: generated events, reconstructed events,
  //                      smeared generated events
  int nfit_used[3];   // number of events
  double pt_dil_data[3][FIT_ARRAY_SIZE];
  double phi_dil_data[3][FIT_ARRAY_SIZE];
  double mass_dil_data[3][FIT_ARRAY_SIZE];
  double rap_dil_data[3][FIT_ARRAY_SIZE];
  double cos_theta_cs_data[3][FIT_ARRAY_SIZE];
  double phi_cs_data[3][FIT_ARRAY_SIZE];
  // weights to keep track of which data point is from gg/qqbar
  double qqbar_weights[3][FIT_ARRAY_SIZE];
  double gg_weights[3][FIT_ARRAY_SIZE];

  // arrays for generating fake data from distributions
  double fake_cos_true[FIT_ARRAY_SIZE];
  double fake_cos_cs[FIT_ARRAY_SIZE];
  double fake_rap[FIT_ARRAY_SIZE];
  int fake_mistag_true[FIT_ARRAY_SIZE];
  int fake_mistag_cs[FIT_ARRAY_SIZE];

  // config file parameters
  bool noFit;
  int numFits;
  std::string paramCacheFile;
  bool internalBremOn;
  bool fixbIn1DFit;
  bool useCosTrueInFit;
  bool artificialCosCS;
};

#endif // ZP2MUASYM_H
