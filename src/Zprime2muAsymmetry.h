#ifndef ZP2MUASYM_H
#define ZP2MUASYM_H

#include <string>

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFitData.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TFMultiD.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

const int FIT_ARRAY_SIZE = 50000;// max size of arrays for unbinned fits. 

// a class for passing around the needed generator-level momenta easily
struct GenMomenta {
  math::XYZTLorentzVector res;
  math::XYZTLorentzVector mum;
  math::XYZTLorentzVector mup;
  math::XYZTLorentzVector mum_noib;
  math::XYZTLorentzVector mup_noib;
  math::XYZTLorentzVector quark;
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

  void makeFakeData(int nEvents, double A_FB, double b);
  void smearGenData(const AsymFitData& data);

  double calcAFBError(double f, double b);
  void calcAsymmetry(double f, double b, double& A_FB, double& e_A_FB);
  void calcAsymmetry(const TH1F* h_cos, double& A_FB, double& e_A_FB);
  void calcAsymmetry(const TH1F* IdF, const TH1F* IdB, TH1F* IdA,
		     fstream& out);

  void calcFrameAsym();
  void fillFrameHistos();

  void fillGenData(const reco::CandidateCollection& mcp);
  void fillFitData(const reco::CandidateCollection& mcp);
  void dumpFitData();

  // JMTBAD GenKineAna methods... to be moved
  bool getGenerated4Vectors(const reco::CandidateCollection& ,
			    unsigned int eventNum,
			    GenMomenta& genMom);
  bool computeFitQuantities(const reco::CandidateCollection&, 
			    int eventNum, AsymFitData& data);

  void bookParamHistos();
  void fillParamHistos(bool fakeData);
  void getAsymParams();

  TFMultiD* getNewFitFcn(int fitType);
  void evalLikelihoods();
  void fitAsymmetry();

  void fitCosCS(TH1F* cos_hist, int params);
  void drawFrameHistos();
  void drawFitHistos();

  void deleteHistos();

  std::string print(const reco::Candidate& par) const;

  // store the event number so we don't have to get it from the edm
  // object every time
  int eventNum;
  TH1F *h_genCosNoCut;
  TH2F *AsymFitSmearHisto[6];
  TH1F *AsymFitHistoGen[6], *AsymFitHistoRec[6], *AsymFitSmearHistoDif[6];
  TH1F *AsymFitHistoGenSmeared[6];
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

  TH1F *h_cos_theta_true, *h_cos_theta_cs_acc;
  TH1F *h_cos_theta_cs, *h_cos_theta_cs_fixed;
  TH1F *h_b_mass, *h_f_mass, *h_b_smass, *h_f_smass;
  TH1F *h_gen_sig[2];
  TH2F *h2_rap_cos_d[2];

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

  std::vector<double> angDist;

  // config file parameters
  edm::InputTag genMuons;
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
  bool useCachedParams;
  bool internalBremOn;
  bool fixbIn1DFit;
  bool useCosTrueInFit;
  bool artificialCosCS;
};

#endif // ZP2MUASYM_H
