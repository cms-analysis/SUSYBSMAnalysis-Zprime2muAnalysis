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

  void fillFitData(const edm::Event& event);
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

  TH1F *cosGJ[MAX_LEVELS][2], *cosCS[MAX_LEVELS][2];
  TH1F *cosBoost[MAX_LEVELS], *cosW[MAX_LEVELS];
  TH1F *cosCSRes[MAX_LEVELS-1];
  TH2F *rap_vs_cosCS[MAX_LEVELS], *rap3_vs_rap0;
  TH1F *FMassGJ[MAX_LEVELS][2], *FMassCS[MAX_LEVELS][2];
  TH1F *BMassGJ[MAX_LEVELS][2], *BMassCS[MAX_LEVELS][2];
  TH1F *AMassGJ[MAX_LEVELS][2], *AMassCS[MAX_LEVELS][2];
  TH1F *FMassBoost[MAX_LEVELS], *FMassW[MAX_LEVELS];
  TH1F *BMassBoost[MAX_LEVELS], *BMassW[MAX_LEVELS];
  TH1F *AMassBoost[MAX_LEVELS], *AMassW[MAX_LEVELS];
  TH1F *FRapGJ[MAX_LEVELS][2],  *FRapCS[MAX_LEVELS][2];
  TH1F *BRapGJ[MAX_LEVELS][2],  *BRapCS[MAX_LEVELS][2];
  TH1F *ARapGJ[MAX_LEVELS][2],  *ARapCS[MAX_LEVELS][2];
  TH1F *FRapBoost[MAX_LEVELS],  *FRapW[MAX_LEVELS];
  TH1F *BRapBoost[MAX_LEVELS],  *BRapW[MAX_LEVELS];
  TH1F *ARapBoost[MAX_LEVELS],  *ARapW[MAX_LEVELS];
  TH1F *FPseudGJ[MAX_LEVELS],   *FPseudCS[MAX_LEVELS];
  TH1F *BPseudGJ[MAX_LEVELS],   *BPseudCS[MAX_LEVELS];
  TH1F *FPseudBoost[MAX_LEVELS],*FPseudW[MAX_LEVELS];
  TH1F *BPseudBoost[MAX_LEVELS],*BPseudW[MAX_LEVELS];
  TProfile *cosCS3_diffsq_vs_cosCS0;
  TH1F *FMBoostCut[MAX_LEVELS][6];
  TH1F *BMBoostCut[MAX_LEVELS][6];
  TH1F *AsymMBoostCut[MAX_LEVELS][6];

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
  TH1F *h_cos_theta_cs_rec;
  TH1F *h_b_mass, *h_f_mass, *h_b_smass, *h_f_smass;
  TH1F *h_gen_sig[2];
  TH2F *h2_rap_cos_d[2];
  TH2F *h2_rap_cos_d_uncut[2];
  TH2F *h2_rap_cos_d_rec;
  TH1F *mistagProbEvents[2];

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
  std::vector<std::string> genSampleFiles;
  double peakMass;
  bool onPeak;
  bool noFit;
  bool onlyEvalLLR;
  FITTYPE fitType;
  int numFits;
  int maxParamEvents;
  bool useCachedParams;
  std::string paramCacheFile;
  bool calcParamsOnly;
  bool internalBremOn;
  bool fixbIn1DFit;
  bool useCosTrueInFit;
  bool artificialCosCS;
};

#endif // ZP2MUASYM_H
