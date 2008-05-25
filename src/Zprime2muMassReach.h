#ifndef ZP2MUMASSREACH_H
#define ZP2MUMASSREACH_H

#include <string>

#include "TH1F.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

// max number of files for mass fits
const int NF = 2;
// max size of arrays for unbinned fits
const int MASS_FIT_ARRAY_SIZE = 20000;

// Auxilliary structure containing event information: invariant mass and weight
struct evt {
  double mass, weight;

  evt& operator = (const evt& rhs) {
    if (this != &rhs) {
      mass    = rhs.mass;
      weight  = rhs.weight;
    }
    return *this;
  }

  bool operator > (const evt& rhs) const {return (mass > rhs.mass);}
  bool operator < (const evt& rhs) const {return (mass < rhs.mass);}
  friend ostream& operator << (ostream& output, const evt& rhs) {
    output << " mass = "   << std::setw(7) << rhs.mass
	   << " weight = " << std::setw(7) << rhs.weight;
    return output;
  }
};

class Zprime2muMassReach : public Zprime2muAnalysis {
 public:
  explicit Zprime2muMassReach(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

 private:
  void bookMassHistos();
  void dilMassPlots(bool debug=false);
  void fillMassArrays();
  void scaleMassHistos();
  void drawMassHistos();
  void bookMassFitHistos();
  void fitMass();
  double unbinnedMassFit(int nentries, double *mass_data,
			 double *weight, const std::string fittype,
			 const bool bckgonly=false);
  int unbinnedFitExec(const Char_t *funcname, const Option_t *option,
		      const std::vector<evt>& events, double& log_ML);
  void analyzeUnbinnedMassFits(const double sumw,
			       const std::vector<evt>& events,
			       const int npars, const double *pars,
			       const std::string fittype, const bool bckgonly,
			       const double mass_min, const double mass_max);
  void ks(std::vector<evt> data, TF1 *func, const double& mass_min,
	  const double& mass_max, double& d, double& prob);
  void drawUnbinnedMassFitsHistos();

  enum FILETYPE { SIGNAL, BKG1, BKG2 };

  int eventCount; // sequential event counter (not the official event number!)
  int fileNum;

  // Configuration parameters:
  std::string psFile;

  bool kDYEvents;
  bool kBackgroundFit;
  bool kGenuineEvents;
  bool kGMR;
  bool useL3Muons;
  bool kFixedMass;
  bool kFixedFWHM;
  bool kSmoothedSample;
  bool kRandomSeed;
  bool kExpPlots;
  bool kBinnedFit;
  unsigned int kNexp;
  double intLumi;

  bool fitGenMass;
  bool fitRecMass;

  std::string resModel;
  unsigned int resMassId;
  unsigned int nBins;
  std::vector<double> massWin;
  std::vector<double> lowerGenMass;
  std::vector<double> upperGenMass;
  std::vector<double> XSec;
  std::vector<double> KFactor;
  std::vector<double> effFilter;
  std::vector<unsigned int> nGenEvents;
 
  // Fixed size arrays for now.
  // number of events
  static int nfit_genmass_used[NF];
  static int nfit_recmass_used[NF];
  // mass
  static double fit_genmass[NF][MASS_FIT_ARRAY_SIZE];
  static double fit_recmass[NF][MASS_FIT_ARRAY_SIZE];
  // generated and reconstructed event number
  static int fit_genevent[NF][MASS_FIT_ARRAY_SIZE];
  static int fit_recevent[NF][MASS_FIT_ARRAY_SIZE];
  // sample ids
  static int fit_sid[NF];

  static const double fwhm_over_sigma;

  // Histograms filled by dilMassPlots
  // three of each (signal, main bkgnd, aux bkgnd)
  TH1F *GenDilMass[3], *HltDilMass[3], *OffDilMass[3];

  // Histograms for the mass fit
  TH1F *GenMassFitData, *RecMassFitData, *RecMassFitDataSmoothed;
  TH1F *GenMassFitSigFr, *GenMassFitMean, *GenMassFitFwhm;
  TH1F *RecMassFitSigFr, *RecMassFitMean, *RecMassFitFwhm;
  TH1F *GenMassFitNevtTot, *GenMassFitNsigTot, *GenMassFitNbkgTot;
  TH1F *RecMassFitNevtTot, *RecMassFitNsigTot, *RecMassFitNbkgTot;
  TH1F *GenMassFitNobs, *GenMassFitNfit, *GenMassFitNsig, *GenMassFitNbkg;
  TH1F *RecMassFitNobs, *RecMassFitNfit, *RecMassFitNsig, *RecMassFitNbkg;
  TH1F *GenMassFitKSDistSig, *GenMassFitKSDistBkg;
  TH1F *GenMassFitKSProbSig, *GenMassFitKSProbBkg;
  TH1F *RecMassFitKSDistSig, *RecMassFitKSDistBkg;
  TH1F *RecMassFitKSProbSig, *RecMassFitKSProbBkg;
  TH1F *GenMassFitSc1, *GenMassFitScl, *GenMassFitSc12, *GenMassFitSc2;
  TH1F *GenMassFitSl2;
  TH1F *RecMassFitSc1, *RecMassFitScl, *RecMassFitSc12, *RecMassFitSc2;
  TH1F *RecMassFitSl2;

  // Significances
  double sc1_gen, sc2_gen, sc12_gen, scl_gen, sl2_gen;
  double sc1_rec, sc2_rec, sc12_rec, scl_rec, sl2_rec;
};

#endif // ZP2MUMASSREACH_H
