#ifndef ZP2MUASYM_H
#define ZP2MUASYM_H

#include <string>

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "AnalysisExamples/Zprime2muAnalysis/src/AsymFitData.h"
#include "AnalysisExamples/Zprime2muAnalysis/src/TFMultiD.h"
#include "AnalysisExamples/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

const int FIT_ARRAY_SIZE = 50000;// max size of arrays for unbinned fits. 

// a class for passing around the needed generator-level momenta easily
struct GenMomenta {
  HepLorentzVector res;
  HepLorentzVector mum;
  HepLorentzVector mup;
  HepLorentzVector mum_noib;
  HepLorentzVector mup_noib;
  HepLorentzVector quark;
};

class Zprime2muAsymmetry : public Zprime2muAnalysis {
 public:
  explicit Zprime2muAsymmetry(const edm::ParameterSet&);
  ~Zprime2muAsymmetry();

  // public:
  void beginJob(const edm::EventSetup&);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

 private:
  static const double nSigma;
  static const unsigned int nNormPoints;
  static const unsigned int nSmearPoints;

  enum FITTYPE { ASYM, GRAV, GRAV_QQBAR, GRAV_GG, GRAV_THEORY };

  void bookFrameHistos();
  void bookFitHistos();

  void calcHistoAsymmetry(const TH1F* IdF, const TH1F* IdB, TH1F* IdA,
			  fstream& out);
  void calcFrameAsym();
  void fillFrameHistos();

  void fillFitData();
  void dumpFitData();

  // JMTBAD GenKineAna methods... to be moved
  void printGenParticle(const HepMC::GenParticle* p) const;
  bool getGenerated4Vectors(const edm::HepMCProduct& mcp,
			    unsigned int eventNum,
			    GenMomenta& genMom);
  bool computeFitQuantities(const edm::HepMCProduct& mcp, 
			    bool internalBremOn, int eventNum,
			    AsymFitData& data);
  void bookParamHistos();
  void fillParamHistos(bool internalBremOn, bool fakeData);
  void getAsymParams(bool internalBremOn);

  TFMultiD* getNewFitFcn(int fitType);
  void evalLikelihoods();
  void fitAsymmetry();

  void fitCosCS(TH1F* cos_hist, int params);
  void drawFrameHistos();
  void drawFitHistos();

  void deleteHistos();

  TH1F *h_genCosNoCut;
  TH2F *AsymFitSmearHisto[6];
  TH1F *AsymFitHistoGen[6], *AsymFitHistoRec[6], *AsymFitSmearHistoDif[6];
  // look at the angular distributions separately by type
  TH1F *AsymFitHistoGenByType[2][6], *AsymFitHistoRecByType[2][6];

  TH1F *cosGJ[NUM_REC_LEVELS][2], *cosCS[NUM_REC_LEVELS][2];
  TH1F *cosBoost[NUM_REC_LEVELS], *cosW[NUM_REC_LEVELS];
  TH1F *cosCSRes[NUM_REC_LEVELS-1];
  TH2F *rap_vs_cosCS[NUM_REC_LEVELS], *rap3_vs_rap0;
  TH1F *FMassGJ[NUM_REC_LEVELS][2], *FMassCS[NUM_REC_LEVELS][2];
  TH1F *BMassGJ[NUM_REC_LEVELS][2], *BMassCS[NUM_REC_LEVELS][2];
  TH1F *AMassGJ[NUM_REC_LEVELS][2], *AMassCS[NUM_REC_LEVELS][2];
  TH1F *FMassBoost[NUM_REC_LEVELS], *FMassW[NUM_REC_LEVELS];
  TH1F *BMassBoost[NUM_REC_LEVELS], *BMassW[NUM_REC_LEVELS];
  TH1F *AMassBoost[NUM_REC_LEVELS], *AMassW[NUM_REC_LEVELS];
  TH1F *FRapGJ[NUM_REC_LEVELS][2],  *FRapCS[NUM_REC_LEVELS][2];
  TH1F *BRapGJ[NUM_REC_LEVELS][2],  *BRapCS[NUM_REC_LEVELS][2];
  TH1F *ARapGJ[NUM_REC_LEVELS][2],  *ARapCS[NUM_REC_LEVELS][2];
  TH1F *FRapBoost[NUM_REC_LEVELS],  *FRapW[NUM_REC_LEVELS];
  TH1F *BRapBoost[NUM_REC_LEVELS],  *BRapW[NUM_REC_LEVELS];
  TH1F *ARapBoost[NUM_REC_LEVELS],  *ARapW[NUM_REC_LEVELS];
  TH1F *FPseudGJ[NUM_REC_LEVELS],   *FPseudCS[NUM_REC_LEVELS];
  TH1F *BPseudGJ[NUM_REC_LEVELS],   *BPseudCS[NUM_REC_LEVELS];
  TH1F *FPseudBoost[NUM_REC_LEVELS],*FPseudW[NUM_REC_LEVELS];
  TH1F *BPseudBoost[NUM_REC_LEVELS],*BPseudW[NUM_REC_LEVELS];
  TProfile *cosCS3_diffsq_vs_cosCS0;
  TH1F *FMBoostCut[NUM_REC_LEVELS][6];
  TH1F *BMBoostCut[NUM_REC_LEVELS][6];
  TH1F *AsymMBoostCut[NUM_REC_LEVELS][6];

  // histos used to get the mistag parameterizations
  TH1F *h_mass_dil[2];
  TH1F *h_rap_mistag, *h_rap_nomistag;
  TH1F *h_mistag[6][3];
  TH2F *h2_mistag[4];
  TH2F *h2_pTrap, *h2_pTrap_mistag;
  TH2F *h2_rap_cos_mistag, *h2_rap_cos_nomistag;
  TH2F *h2_rap_cos_p;
  TH1F *h_cos_mistag, *h_cos_cs, *h_cos_mistag_prob;
  TH1F *h_rap_dil[2], *h_pt_dil;
  TH1F *h_phi_cs;

  // Unbinned fits of cos(theta) vs Y.  Fixed size arrays for now.
  int nfit_ycos_used; // number of events
  int nfit_used[2];   // number of events
  double fit_data1[FIT_ARRAY_SIZE]; // cos(theta_CS)
  double fit_data2[FIT_ARRAY_SIZE]; // Y
  double pt_dil_data[2][FIT_ARRAY_SIZE];
  double phi_dil_data[2][FIT_ARRAY_SIZE];
  double mass_dil_data[2][FIT_ARRAY_SIZE];
  double rap_dil_data[2][FIT_ARRAY_SIZE];
  double cos_theta_cs_data[2][FIT_ARRAY_SIZE];
  double phi_cs_data[2][FIT_ARRAY_SIZE];
  // weights to keep track of which data point is from gg/qqbar
  double qqbar_weights[2][FIT_ARRAY_SIZE];
  double gg_weights[2][FIT_ARRAY_SIZE];

  bool ibOnParams;

  // config file parameters
  VERBOSITY verbosity;
  std::string outputFileBase;
  std::string genSample;
  double peakMass;
  bool onPeak;
  bool noFit;
  bool onlyEvalLLR;
  FITTYPE fitType;
  int numFits;
  int maxParamEvents;
};

#endif // ZP2MUASYM_H
