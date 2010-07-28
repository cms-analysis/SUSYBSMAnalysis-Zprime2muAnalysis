//
// Author: Slava Valuev, UCLA
//

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TPad.h"
#include "TPostScript.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TText.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Functions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/UnbinnedFitter.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerDecision.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Instead of bringing Zprime2muMassReach into the new framework, it
// is easier to keep it around in this way by having the (needed bits
// of the) Zprime2muAnalysis class defined here. The Zprime2muAnalysis
// class is not left in its own file so that no one thinks it is
// available to be or should be used.

class Zprime2muAnalysis : public edm::EDAnalyzer {
 public:
  explicit Zprime2muAnalysis(const edm::ParameterSet&);
  virtual ~Zprime2muAnalysis() {}

  virtual void beginJob() {}
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}

 protected:
  bool doingElectrons;
  bool useGen;
  int lBest;

  edm::Service<TFileService> fs;
  HardInteraction hardInteraction;
  TriggerDecision trigDecision;
  RecLevelHelper recLevelHelper;

  reco::CandidateBaseRefVector allLeptons[MAX_LEVELS];
  reco::CompositeCandidateCollection allDileptons[MAX_LEVELS];

  double resonanceMass(const reco::CompositeCandidate& dil) const;
};

Zprime2muAnalysis::Zprime2muAnalysis(const edm::ParameterSet& config) 
  : doingElectrons(config.getParameter<bool>("doingElectrons")),
    useGen(config.getParameter<bool>("useGen")),
    lBest(config.getParameter<int>("bestRecLevel")),
    hardInteraction(doingElectrons, true)
{
  InitROOT();
  trigDecision.init(config, false);
  recLevelHelper.init(config);
}

void Zprime2muAnalysis::analyze(const edm::Event& event, const edm::EventSetup&) {
  trigDecision.initEvent(event);
  if (useGen) hardInteraction.Fill(event);

  recLevelHelper.initEvent(event);
  for (int irec = 0; irec < MAX_LEVELS; irec++) {
    recLevelHelper.getLeptons(event, irec, allLeptons[irec]);
    recLevelHelper.getDileptons(event, irec, allDileptons[irec]);
  }
}

double Zprime2muAnalysis::resonanceMass(const reco::CompositeCandidate& dil) const {
  // Gen level muons don't have photons matched to them; we can take
  // the mass of the resonance directly from the MC record.
  if (recLevelHelper.recLevel(dil) == lGN) {
    if (hardInteraction.resonance)
      return hardInteraction.resonance->mass();
    else
      return dil.mass();
  }

  const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
  const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);

  reco::Particle::LorentzVector p4 = dil.p4(); // == lep0->p4() + lep1->p4() by construction.

  // If the leptons are electrons, or there were no photons in the
  // cone dR < 0.1, these four-vectors will be 0.
  reco::Particle::LorentzVector ph0 = recLevelHelper.closestPhoton(lep0);
  reco::Particle::LorentzVector ph1 = recLevelHelper.closestPhoton(lep1);

  p4 += ph0;
  // Don't double-count.
  if (ph1 != ph0) p4 += ph1;

  return p4.mass();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
  unsigned eventNum;

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
  static unsigned fit_genevent[NF][MASS_FIT_ARRAY_SIZE];
  static unsigned fit_recevent[NF][MASS_FIT_ARRAY_SIZE];
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

// Number of events in unbinned fits of mass.
int Zprime2muMassReach::nfit_genmass_used[NF] = {0};
int Zprime2muMassReach::nfit_recmass_used[NF] = {0};

// Generated and reconstructed mass, gen. and rec. event number, and sample id.
double Zprime2muMassReach::fit_genmass[NF][MASS_FIT_ARRAY_SIZE] = {{0.}};
double Zprime2muMassReach::fit_recmass[NF][MASS_FIT_ARRAY_SIZE] = {{0.}};
unsigned Zprime2muMassReach::fit_genevent[NF][MASS_FIT_ARRAY_SIZE] = {{0}};
unsigned Zprime2muMassReach::fit_recevent[NF][MASS_FIT_ARRAY_SIZE] = {{0}};
int Zprime2muMassReach::fit_sid[NF] = {0};

const double Zprime2muMassReach::fwhm_over_sigma = 2.*sqrt(2.*log(2.));

Zprime2muMassReach::Zprime2muMassReach(const edm::ParameterSet& config) 
  : Zprime2muAnalysis(config), eventCount(0), fileNum(-1)
{
  // Postscript file with mass reach plots.
  psFile    = config.getUntrackedParameter<std::string>("psFile", "mass_fits.ps");

  kDYEvents       = config.getParameter<bool>("DYEvents");
  kBackgroundFit  = config.getParameter<bool>("BackgroundFit");
  kGenuineEvents  = config.getParameter<bool>("GenuineEvents");
  kGMR            = config.getParameter<bool>("GMR");
  kFixedMass      = config.getParameter<bool>("FixedMass");
  kFixedFWHM      = config.getParameter<bool>("FixedFWHM");
  kSmoothedSample = config.getParameter<bool>("SmoothedSample");
  kRandomSeed     = config.getParameter<bool>("RandomSeed");
  kExpPlots       = config.getParameter<bool>("ExpPlots");
  kBinnedFit      = config.getParameter<bool>("BinnedFit");
  kNexp           = config.getParameter<unsigned int>("Nexp");
  intLumi         = config.getParameter<double>("intLumi");

  fitGenMass = config.getParameter<bool>("fitGenMass");
  fitRecMass = config.getParameter<bool>("fitRecMass");

  // get the parameters specific to the data sample on which we are running
  std::string dataSet = config.getParameter<std::string>("dataSet");
  edm::ParameterSet dataSetConfig =
    config.getParameter<edm::ParameterSet>(dataSet);

  resModel     = dataSetConfig.getParameter<std::string>("resModel");
  resMassId    = dataSetConfig.getParameter<unsigned int>("resMassId");
  nBins        = dataSetConfig.getParameter<unsigned int>("nBins");
  massWin      = dataSetConfig.getParameter<std::vector<double> >("massWin");
  lowerGenMass = dataSetConfig.getParameter<std::vector<double> >("lowerGenMass");
  upperGenMass = dataSetConfig.getParameter<std::vector<double> >("upperGenMass");
  XSec         = dataSetConfig.getParameter<std::vector<double> >("XSec");
  KFactor      = dataSetConfig.getParameter<std::vector<double> >("KFactor");
  effFilter    = dataSetConfig.getParameter<std::vector<double> >("effFilter");
  nGenEvents   = dataSetConfig.getParameter<std::vector<unsigned int> >("nGenEvents");

  bookMassHistos();
  bookMassFitHistos();
}

void Zprime2muMassReach::analyze(const edm::Event& event, 
				 const edm::EventSetup& eSetup) {
  // delegate filling our muon vectors to the parent class
  Zprime2muAnalysis::analyze(event, eSetup);

  eventNum = event.id().event(); // not unique in data but the code that uses this needs reworking anyway...

  // JMTBAD Previously, we loaded three separate files in each job:
  // signal, DY background, and DY auxillary background. We can still
  // do that in CMSSW with PoolSource's fileNames array, but we have
  // no way of getting at the current filename or even which file
  // number we're on! (See the HN thread
  // https://hypernews.cern.ch/HyperNews/CMS/get/edmFramework/500.html
  // .) So, for now we do a hacky trick: we count events and augment the file
  // number when the number of events reaches some pre-defined value.
  // The total number of events in each sample must therefore be known
  // beforehand (and defined in the nGenEvents[fileNum] array).
  //
  // Or: can we set the run # (e.g., = 0,1,2) in the edm root files
  // when we generate them? that would at least deal with the hocus
  // pocus of watching the event number for a reset...
  if (eventCount == 0) {
    fileNum++;
    edm::LogInfo("Zprime2muMassReach") << "\n New sample type found. \n";
  }
  eventCount++;
  if (static_cast<unsigned int>(eventCount) == nGenEvents[fileNum]) eventCount = 0;

  if (fileNum < 0 || fileNum > 2)
    throw cms::Exception("MassReach") << "+++ invalid fileNum! +++\n";

  // fill the dilepton mass dist histos
  dilMassPlots();

  // Skip first ("signal") sample if fit DY (kDYEvents=true), and
  // second (DY closest in mass) sample if fit signal (kDYEvents=false).
  if ((kDYEvents && fileNum != SIGNAL) || (!kDYEvents && fileNum != BKG1))
    fillMassArrays();
}

void Zprime2muMassReach::endJob() {
  Zprime2muAnalysis::endJob();

  scaleMassHistos();
  drawMassHistos();

  // Check whether the fitter set up to analyse this model/mass or not.
  bool avail = false;
  if (resModel == "Zssm" || resModel == "Zpsi" || resModel == "Zeta" ||
      resModel == "Zchi" || resModel == "Zlrm" || resModel == "Zalrm") {
    // 1, 3, and 5 TeV Z's
    if (resMassId == 1000 || resMassId == 1200 || resMassId == 1500 ||
	resMassId == 2000 || resMassId == 3000 || resMassId == 5000) {
      avail = true;
    }
  }
  // Add more models and masses later...  They might work out of the box,
  // but this needs to be checked.

  if (!avail) {
    edm::LogWarning("endJob") 
      << "+++ Sorry, mass reach analysis for resModel = " << resModel
      << " and resMassId = " << resMassId << " is not yet available +++";
    return;
  }

  fitMass();
  drawUnbinnedMassFitsHistos();
}

void Zprime2muMassReach::bookMassHistos() {
  double xmin = massWin[0];
  double xmax = massWin[1];
  for (unsigned int i = 0; i < 3; i++) {
    GenDilMass[i] = fs->make<TH1F>(nameHist("GenDilMass", i),
				   "True Dilepton Mass",
				   nBins, xmin, xmax);
    HltDilMass[i] = fs->make<TH1F>(nameHist("HltDilMass", i),
				   "Reconstructed (GMR) Dilepton Mass",
				   nBins, xmin, xmax);
    OffDilMass[i] = fs->make<TH1F>(nameHist("OffDilMass", i),
				   "Reconstructed (best) Dilepton Mass",
				   nBins, xmin, xmax);
  }
}

void Zprime2muMassReach::dilMassPlots(const bool debug) {
  // Plots only made if there is generator information.
  if (!useGen) return;

  double genm, recm;

  // Fill GenDilMass here if you want ALL generated events; uncomment Fill
  // in if-statement below if you want gen. mass plot only for events
  // with a reconstructed mu+mu- pair.
  if (allDileptons[lGN].size() > 0) {
    genm = resonanceMass(allDileptons[lGN].at(0));
    GenDilMass[fileNum]->Fill(genm);
  }

  // Check to see if event passed trigger cuts
  if (!trigDecision.pass()) return;

  // "Off-line" dileptons in events passing the trigger and quality cuts
  for (int i = 0; i < 2; i++) {
    const reco::CompositeCandidateCollection& dileptons =
      i == 1 ? allDileptons[lBest] : allDileptons[lGR];

    if (dileptons.size() > 0) {
      recm = resonanceMass(dileptons[0]); // highest mass dilepton
      if (allDileptons[lGN].size() > 0) {
	genm = resonanceMass(allDileptons[lGN].at(0));

	// Fill generated mass plot only for those events in which
	// TMR dilepton was found.
	// if (i == 1) GenDilMass->Fill(genm);
	
	// Reconstructed dilepton invariant mass within specific
	  // generated-mass regions (to avoid double-counting).
	if (genm >= lowerGenMass[fileNum] && genm < upperGenMass[fileNum]) {
	  if (i == 0)      HltDilMass[fileNum]->Fill(recm);
	  else if (i == 1) OffDilMass[fileNum]->Fill(recm);
	}
      }
    }
  }
}

void Zprime2muMassReach::fillMassArrays() {
  // Fill arrays for unbinned maximum likelihood fit
  
  // Keep track of how many times we're called.
  // We only fit NF (=2) files.
  static int sample_id = -1, idx = -1;
  if (sample_id != fileNum) {
    if (++idx < NF)
      fit_sid[idx] = sample_id = fileNum;
    else
      throw cms::Exception("fillMassArrays") 
	<< "+++ tried to fit too many files! idx = " << idx << " +++\n";
  }

  // True mass
  if (useGen) {
    for (unsigned idi = 0; idi < allDileptons[lGN].size(); idi++) {
      if (nfit_genmass_used[idx] < MASS_FIT_ARRAY_SIZE) {
	fit_genmass[idx][nfit_genmass_used[idx]] = resonanceMass(allDileptons[lGN].at(idi));
	fit_genevent[idx][nfit_genmass_used[idx]] = eventNum;
	nfit_genmass_used[idx]++;
      }
      else
	edm::LogWarning("fillMassArrays")
	  << "+++ MASS_FIT_ARRAY_SIZE is too small to keep all gen events";
    }
  }

  // Make sure the event passed the trigger
  if (!trigDecision.pass()) return;

  // Reconstructed dileptons
  const reco::CompositeCandidateCollection& dileptons = 
    kGMR ? allDileptons[lGR] : allDileptons[lBest];
  for (unsigned idi = 0; idi < dileptons.size(); idi++) {
    // Check the "off-line" track quality and apply the cuts
    if (nfit_recmass_used[idx] < MASS_FIT_ARRAY_SIZE) {
      double resmass = resonanceMass(dileptons[idi]);
      fit_recmass[idx][nfit_recmass_used[idx]] = resmass;

      // TESTS ONLY!!!
      //double epsil = 0.042;
      //double gsmear = gRandom->Gaus(0., epsil);
      //fit_recmass[idx][nfit_recmass_used[idx]] = resmass*(1 + gsmear);

      fit_recevent[idx][nfit_recmass_used[idx]] = eventNum;
      nfit_recmass_used[idx]++;
    }
    else
      edm::LogWarning("fillMassArrays")
	<< "+++ MASS_FIT_ARRAY_SIZE is too small to keep all rec events";
  }
}

void Zprime2muMassReach::scaleMassHistos() {
  // Normalize appropriate signal and background histograms according to the
  // cross sections and int. luminosity of interest.
  double scale;

  for (unsigned int i = 0; i < 3; i++) {
    if (XSec[i] > 0.)
      scale = XSec[i]*KFactor[i]*intLumi/(nGenEvents[i]/effFilter[i]);
    else
      scale = 1.;

    // To keep number of signal events as is, and normalize backgrounds
    // relative to it.  Useful for tests.
    /* 
    if (i == BKG1)
      scale = (XSec[BKG1]*KFactor[BKG1])/(XSec[SIGNAL]*KFactor[SIGNAL]);
    else if (i == BKG2)
      scale = (XSec[BKG2]*KFactor[BKG2])/(XSec[SIGNAL]*KFactor[SIGNAL]);
    else
      scale = 1.;
    */

    GenDilMass[i]->Scale(scale);
    HltDilMass[i]->Scale(scale);
    OffDilMass[i]->Scale(scale);
  }
}

void Zprime2muMassReach::drawMassHistos() {
  // Plots of Signal + Background histograms super-imposed.

  // Combine contents of histos: total background and the sum of
  // the signal and the total background.
  TH1F* HltDilMassBgTot = (TH1F*)HltDilMass[BKG1]->Clone();
  TH1F* OffDilMassBgTot = (TH1F*)OffDilMass[BKG1]->Clone();
  TH1F* GenDilMassTot   = (TH1F*)GenDilMass[SIGNAL]->Clone();
  TH1F* HltDilMassTot   = (TH1F*)HltDilMass[SIGNAL]->Clone();
  TH1F* OffDilMassTot   = (TH1F*)OffDilMass[SIGNAL]->Clone();

  // Total background
  HltDilMassBgTot->Add(HltDilMass[BKG1], HltDilMass[BKG2], 1., 1.);
  OffDilMassBgTot->Add(OffDilMass[BKG1], OffDilMass[BKG2], 1., 1.);

  // Sum of the signal and backgrounds.
  // Adding Drell-Yan background in the same mass range to signal
  // double-counts the same Drell-Yan background for Z', but not for
  // gravitons (G* generation in PYTHIA does not do interference).
  if (resModel != "G") {
    HltDilMassTot->Add(HltDilMass[SIGNAL], HltDilMass[BKG2], 1., 1.);
    OffDilMassTot->Add(OffDilMass[SIGNAL], OffDilMass[BKG2], 1., 1.);
  }
  else {
    GenDilMassTot->Add(GenDilMass[SIGNAL], GenDilMass[BKG1], 1., 1.);
    HltDilMassTot->Add(HltDilMass[SIGNAL], HltDilMassBgTot, 1., 1.);
    OffDilMassTot->Add(OffDilMass[SIGNAL], OffDilMassBgTot, 1., 1.);
  }

  // Histograms of running integrals.
  TH1F* HltMassSigInt = (TH1F*)HltDilMass[SIGNAL]->Clone();
  TH1F* OffMassSigInt = (TH1F*)OffDilMass[SIGNAL]->Clone();
  TH1F* HltMassBgInt  = (TH1F*)HltDilMassBgTot->Clone();
  TH1F* OffMassBgInt  = (TH1F*)OffDilMassBgTot->Clone();
  TH1F* HltMassTotInt = (TH1F*)HltDilMassTot->Clone();
  TH1F* OffMassTotInt = (TH1F*)OffDilMassTot->Clone();
  int nbins = OffMassTotInt->GetNbinsX();
  for (int ibin = 0; ibin < nbins; ibin++) {
    HltMassSigInt->SetBinContent(ibin, HltDilMass[SIGNAL]->Integral(1, ibin));
    OffMassSigInt->SetBinContent(ibin, OffDilMass[SIGNAL]->Integral(1, ibin));
    HltMassBgInt->SetBinContent(ibin,  HltDilMassBgTot->Integral(1, ibin));
    OffMassBgInt->SetBinContent(ibin,  OffDilMassBgTot->Integral(1, ibin));
    HltMassTotInt->SetBinContent(ibin, HltDilMassTot->Integral(1, ibin));
    OffMassTotInt->SetBinContent(ibin, OffDilMassTot->Integral(1, ibin));
  }

  // Plot histos and save canvas.
  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 700); // portrait
  TPad *pad[10];
  for (int i_page = 0; i_page < 10; i_page++) {
    pad[i_page] = new TPad("","", .05, .05, .95, .93);
  }
  int page = 0;

  std::ostringstream strtitl;
  strtitl << resModel << " (" << resMassId/1000.
	  << " TeV) vs Drell-Yan background";
  TPaveLabel *title;

  std::ostringstream strlumi, strbinw, strmass;
  strlumi << intLumi;
  strbinw << (massWin[1]-massWin[0])/nBins;
  strmass << resMassId;
  std::string ytit = "Events/" + strbinw.str() + " GeV/" + strlumi.str() + " fb-1";
  std::string psfile = resModel + strmass.str() + "_sig_bckg.ps";
  double mass_min = massWin[0], mass_max = massWin[1]; // for prints only.
  TPostScript *ps  = new TPostScript(psfile.c_str(), 111);

  // Generated mass
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, strtitl.str().c_str());
  title->SetFillColor(10);
  title->Draw();
  gStyle->SetOptStat(0);
  gStyle->SetOptLogy(0);
  pad[page]->Draw();
  pad[page]->Divide(1,1);
  std::ostringstream out;
  out
    << "-------------------------------------------------------------\n"
    << " Signal-vs-background plots: "
    << "\n   Generated signal events = " << GenDilMass[SIGNAL]->Integral()
    << "; gen. between " << mass_min << " and " << mass_max << " = "
    << GenDilMass[SIGNAL]->Integral(GenDilMass[SIGNAL]->FindBin(mass_min),
				    GenDilMass[SIGNAL]->FindBin(mass_max)-1)
    << "\n   Closest   bckg.  events = " << GenDilMass[BKG1]->Integral()
    << "; gen. between " << mass_min << " and " << mass_max << " = "
    << GenDilMass[BKG1]->Integral(GenDilMass[BKG1]->FindBin(mass_min),
				  GenDilMass[BKG1]->FindBin(mass_max)-1)
    << "\n   Far       bckg.  events = " << GenDilMass[BKG2]->Integral()
    << "; gen. between " << mass_min << " and " << mass_max << " = "
    << GenDilMass[BKG2]->Integral(GenDilMass[BKG2]->FindBin(mass_min),
				  GenDilMass[BKG2]->FindBin(mass_max)-1)
    << "\n-------------------------------------------------------------";
  LogTrace("drawMassHistos") << out.str();
  //GenDilMassTot->GetXaxis()->SetNdivisions(5);
  GenDilMassTot->GetXaxis()->SetTitle("mu+ mu- mass");
  GenDilMassTot->GetYaxis()->SetTitle(ytit.c_str());
  GenDilMassTot->GetYaxis()->SetTitleOffset(1.2);
  GenDilMassTot->SetLineColor(kRed);
  if (GenDilMass[BKG1]->GetMaximum() > GenDilMassTot->GetMaximum())
    GenDilMassTot->SetMaximum(1.15*GenDilMassTot->GetMaximum());
  pad[page]->cd(1);
  GenDilMassTot->Draw();
  if (resModel == "G") {
    GenDilMass[SIGNAL]->SetLineColor(kGreen);
    GenDilMass[SIGNAL]->Draw("same");
  }
  GenDilMass[BKG1]->SetLineColor(kBlue);
  GenDilMass[BKG1]->Draw("same");
  page++;
  c1->Update();

  // Set ymax for all reconstruction plots
  double ymax_hlt =
    (HltDilMassBgTot->GetMaximum() > HltDilMass[SIGNAL]->GetMaximum()) ?
     HltDilMassBgTot->GetMaximum() : HltDilMass[SIGNAL]->GetMaximum();
  double ymax_off =
    (OffDilMassBgTot->GetMaximum() > OffDilMass[SIGNAL]->GetMaximum()) ?
     OffDilMassBgTot->GetMaximum() : OffDilMass[SIGNAL]->GetMaximum();
  double ymax = (ymax_hlt > ymax_off) ? 1.1*ymax_hlt : 1.1*ymax_off;

  // Intermediate plots for the reconstructed mass: individual components
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptDate(0);
  title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, strtitl.str().c_str());
  title->SetFillColor(10);
  title->Draw();
  pad[page]->Draw();
  pad[page]->Divide(1,2);
  // L3 or GMR muons (depending on the flags).
  out.str("");
  out 
    << "-------------------------------------------------------------\n"
    << "   L3/GMR  signal   events = " << HltDilMass[SIGNAL]->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << HltDilMass[SIGNAL]->Integral(HltDilMass[SIGNAL]->FindBin(mass_min),
				    HltDilMass[SIGNAL]->FindBin(mass_max)-1)
    << "\n   L3/GMR  bckg #1  events = " << HltDilMass[BKG1]->Integral()
    << "\n       made above gen mass = " << HltDilMass[BKG1]->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << HltDilMass[BKG1]->Integral(HltDilMass[BKG1]->FindBin(mass_min),
				  HltDilMass[BKG1]->FindBin(mass_max)-1)
    << "\n   L3/GMR  bckg #2  events = " << HltDilMass[BKG2]->Integral()
    << "\n       made below gen mass = " << HltDilMass[BKG2]->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << HltDilMass[BKG2]->Integral(HltDilMass[BKG2]->FindBin(mass_min),
				  HltDilMass[BKG2]->FindBin(mass_max)-1)
    << "\n   L3/GMR tot. bckg events = " << HltDilMassBgTot->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << HltDilMassBgTot->Integral(HltDilMassBgTot->FindBin(mass_min),
				 HltDilMassBgTot->FindBin(mass_max)-1)
    << "\n   L3/GMR total     events = " << HltDilMassTot->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << HltDilMassTot->Integral(HltDilMassTot->FindBin(mass_min),
			       HltDilMassTot->FindBin(mass_max)-1)
    << "\n-------------------------------------------------------------";
  LogTrace("drawMassPlots") << out.str();

  HltDilMass[SIGNAL]->GetXaxis()->SetTitle("mu+ mu- mass");
  HltDilMass[SIGNAL]->GetYaxis()->SetTitle(ytit.c_str());
  HltDilMass[SIGNAL]->SetMaximum(ymax);
  HltDilMass[SIGNAL]->SetLineColor(kRed);
  pad[page]->cd(1);
  HltDilMass[SIGNAL]->Draw();
  HltDilMass[BKG1]->SetLineColor(kBlue);
  HltDilMass[BKG1]->Draw("same");
  HltDilMass[BKG2]->SetLineColor(kGreen);
  HltDilMass[BKG2]->Draw("same");

  // "Best" offline reconstruction.
  out.str("");
  out
    << "-------------------------------------------------------------\n"
    << "   Optimized  sign. events = " << OffDilMass[SIGNAL]->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << OffDilMass[SIGNAL]->Integral(OffDilMass[SIGNAL]->FindBin(mass_min),
				    OffDilMass[SIGNAL]->FindBin(mass_max)-1)
    << "\n   Optimized bckg #1 events= " << OffDilMass[BKG1]->Integral()
    << "\n       made above gen mass = " << OffDilMass[BKG1]->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << OffDilMass[BKG1]->Integral(OffDilMass[BKG1]->FindBin(mass_min),
				  OffDilMass[BKG1]->FindBin(mass_max)-1)
    << "\n   Optimized bckg #2 events= " << OffDilMass[BKG2]->Integral()
    << "\n       made below gen mass = " << OffDilMass[BKG2]->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << OffDilMass[BKG2]->Integral(OffDilMass[BKG2]->FindBin(mass_min),
				  OffDilMass[BKG2]->FindBin(mass_max)-1)
    << "\n   Optimized total bckg    = " << OffDilMassBgTot->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << OffDilMassBgTot->Integral(OffDilMassBgTot->FindBin(mass_min),
				 OffDilMassBgTot->FindBin(mass_max)-1)
    << "\n   Optimized total events  = " << OffDilMassTot->Integral()
    << "; rec. between " << mass_min << " and " << mass_max << " = "
    << OffDilMassTot->Integral(OffDilMassTot->FindBin(mass_min),
			       OffDilMassTot->FindBin(mass_max)-1)
    << "\n-------------------------------------------------------------";
  LogTrace("drawMassPlots") << out.str();

  OffDilMass[SIGNAL]->GetXaxis()->SetTitle("mu+ mu- mass");
  OffDilMass[SIGNAL]->GetYaxis()->SetTitle(ytit.c_str());
  OffDilMass[SIGNAL]->SetMaximum(ymax);
  OffDilMass[SIGNAL]->SetLineColor(kRed);
  pad[page]->cd(2);
  OffDilMass[SIGNAL]->Draw();
  OffDilMass[BKG1]->SetLineColor(kBlue);
  OffDilMass[BKG1]->Draw("same");
  OffDilMass[BKG2]->SetLineColor(kGreen);
  OffDilMass[BKG2]->Draw("same");
  page++;
  c1->Update();

  // Signal plus background vs total background
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, strtitl.str().c_str());
  title->SetFillColor(10);
  //title->Draw();
  pad[page]->Draw();
  pad[page]->Divide(1,1);
  pad[page]->cd(1);                 
  HltDilMassTot->GetXaxis()->SetTitle("mu+ mu- mass");
  HltDilMassTot->GetYaxis()->SetTitle(ytit.c_str());
  HltDilMassTot->GetYaxis()->SetTitleOffset(1.2);
  HltDilMassTot->SetMaximum(ymax);
  HltDilMassTot->SetLineColor(kRed);  HltDilMassTot->Draw();
  //HltDilMassBgTot->SetLineStyle(2);
  //HltDilMassBgTot->SetFillColor(1); HltDilMassBgTot->SetFillStyle(3013);
  HltDilMassBgTot->SetLineColor(kBlue);  HltDilMassBgTot->Draw("same");
  //HltDilMass[BKG1]->SetLineColor(4); HltDilMass[BKG1]->Draw("same");
  //HltDilMass[BKG2]->SetLineColor(5); HltDilMass[BKG2]->Draw("same");
  /*
    pad[page]->cd(2);
    gPad->SetGrid(1);
    HltMassSigInt->GetXaxis()->SetTitle("mu+ mu- mass");
    HltMassSigInt->GetYaxis()->SetTitleOffset(1.2);
    HltMassSigInt->GetYaxis()->SetTitle("Integrated Events/100 GeV/100 fb-1");
    HltMassSigInt->SetMaximum(300.);
    HltMassSigInt->SetLineColor(2); HltMassSigInt->Draw();
    //HltMassTotInt->SetMaximum(420.);  HltMassTotInt->Draw();
    HltMassBgInt->SetLineColor(4);  HltMassBgInt->Draw("same");
  */
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  //title->Draw();
  pad[page]->Draw();
  pad[page]->Divide(1,1);
  pad[page]->cd(1);
  OffDilMassTot->GetXaxis()->SetTitle("mu+ mu- mass");
  OffDilMassTot->GetYaxis()->SetTitle(ytit.c_str());
  OffDilMassTot->GetYaxis()->SetTitleOffset(1.2);
  OffDilMassTot->SetMaximum(ymax);
  OffDilMassTot->SetLineColor(kRed);  OffDilMassTot->Draw();
  OffDilMassBgTot->SetLineColor(kBlue);  OffDilMassBgTot->Draw("same");
  //OffDilMass[BKG1]->SetLineColor(4); OffDilMass[BKG1]->Draw("same");
  //OffDilMass[BKG2]->SetLineColor(5); OffDilMass[BKG2]->Draw("same");
  /*
    pad[page]->cd(2);
    gPad->SetGrid(1);
    OffMassSigInt->GetXaxis()->SetTitle("mu+ mu- mass");
    OffMassSigInt->GetYaxis()->SetTitleOffset(1.2);
    OffMassSigInt->GetYaxis()->SetTitle("Integrated Events/100 GeV/100 fb-1");
    OffMassSigInt->SetMaximum(300.);
    OffMassSigInt->SetLineColor(2); OffMassSigInt->Draw();
    //OffMassTotInt->SetMaximum(420.);  OffMassTotInt->Draw();
    OffMassBgInt->SetLineColor(4);  OffMassBgInt->Draw("same");
  */
  page++;
  c1->Update();
  ps->Close();

//#define FORTALKS
#ifdef FORTALKS
  TCanvas *c3 = new TCanvas("c3", "", 0, 0, 500, 500);

  // Generated mass
  TPostScript *eps1 = new TPostScript(modmass+"_sb_gen.eps", 113);
  eps1->NewPage();
  c3->Clear();
  c3->cd(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptLogy(0);
  pad[page]->Draw();
  pad[page]->cd(1);
  GenDilMassTot->SetTitle("");
  GenDilMassTot->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass (GeV)");
  GenDilMassTot->GetXaxis()->SetLabelOffset(0.01);
  GenDilMassTot->GetYaxis()->SetTitle(ytit.c_str());
  GenDilMassTot->GetYaxis()->SetTitleOffset(1.15);
  GenDilMassTot->SetLineColor(kBlue);
  GenDilMass[BKG1]->SetLineColor(kBlack);
  GenDilMass[BKG1]->SetFillColor(kGreen);
  if      (resMassId == 1000) GenDilMassTot->SetMaximum(11.0);
  else if (resMassId == 3000) GenDilMassTot->SetMaximum( 7.0);
  else if (resMassId == 5000) GenDilMassTot->SetMaximum( 1.6);
  // Just to make plot look nicer: exclude downward fluctuations of signal
  // at low (background-dominated) masses.
  for (int ibin = 1; ibin <= 10; ibin++) {
    if (GenDilMass[BKG1]->GetBinContent(ibin) >
	GenDilMassTot->GetBinContent(ibin))
      GenDilMassTot->SetBinContent(ibin,GenDilMass[BKG1]->GetBinContent(ibin));
  }
  GenDilMassTot->Draw();
  GenDilMass[BKG1]->Draw("same");  GenDilMassTot->Draw("same");
  if (resModel == "G") {
    GenDilMass[SIGNAL]->SetLineColor(kGreen);
    GenDilMass[SIGNAL]->Draw("same");
  }
  c3->RedrawAxis();
  page++;
  c3->Update();
  eps1->Close();

  TPostScript *eps2 = new TPostScript(modmass+"_sb_rec.eps", 113);
  eps2->NewPage();
  c3->Clear();
  c3->cd(0);
  pad[page]->Draw();
  pad[page]->cd(1);
  OffDilMassTot->SetTitle("");
  OffDilMassTot->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass (GeV)");
  OffDilMassTot->GetXaxis()->SetLabelOffset(0.01);
  OffDilMassTot->GetYaxis()->SetTitle(ytit.c_str());
  OffDilMassTot->GetYaxis()->SetTitleOffset(1.25);
  OffDilMassTot->SetLineColor(kBlue);
  OffDilMassBgTot->SetLineColor(kBlack); OffDilMassBgTot->SetFillColor(kGreen);
  if      (resMassId == 1000) OffDilMassTot->SetMaximum(11.0);
  else if (resMassId == 3000) OffDilMassTot->SetMaximum( 5.0);
  else if (resMassId == 5000) OffDilMassTot->SetMaximum( 1.6);
  for (int ibin = 1; ibin <= 10; ibin++) {
    if (OffDilMassBgTot->GetBinContent(ibin) >
	OffDilMassTot->GetBinContent(ibin))
      OffDilMassTot->SetBinContent(ibin, OffDilMassBgTot->GetBinContent(ibin));
  }
  OffDilMassTot->Draw();
  OffDilMassBgTot->Draw("same");  OffDilMassTot->Draw("same");
  c3->RedrawAxis();
  page++;
  c3->Update();
  eps2->Close();
#endif
}

void Zprime2muMassReach::bookMassFitHistos() {
  // Generated and reconstructed mass histo
  GenMassFitData = fs->make<TH1F>("GenMassFitData", "Generated mass",     100, massWin[0], massWin[1]);
  RecMassFitData = fs->make<TH1F>("RecMassFitData", "Reconstructed mass", 100, massWin[0], massWin[1]);
  GenMassFitData->Sumw2(); // needed for Kolmogorov-Smirnov tests
  RecMassFitData->Sumw2();

  if (kSmoothedSample) {
    RecMassFitDataSmoothed = fs->make<TH1F>("RecMassFitDataSmoothed", "Smoothed rec. mass", 100, massWin[0], massWin[1]);
    RecMassFitDataSmoothed->Sumw2();
  }

  // Signal fraction, mean and FWHM
  GenMassFitSigFr  = fs->make<TH1F>("GenMassFitSigFr", "Fraction of signal events", 100, 0., 1.);
  RecMassFitSigFr  = fs->make<TH1F>("RecMassFitSigFr", "Fraction of signal events", 100, 0., 1.);

  double mean_min = 0, mean_max = 10000.;
  double fwhm_min = 0, fwhm_max =  1000.;
  if (resMassId == 1000) {
    mean_min =  800.; mean_max = 1200.; fwhm_max =  100.;
  }
  else if (resMassId == 1500) {
    mean_min = 1000.; mean_max = 2000.; fwhm_max =  200.;
  }
  else if (resMassId == 2000) {
    mean_min = 1000.; mean_max = 3000.; fwhm_max =  400.;
  }
  else if (resMassId == 3000) {
    mean_min = 2500.; mean_max = 3500.; fwhm_max =  500.;
  }
  else if (resMassId == 5000) {
    mean_min = 3000.; mean_max = 7000.; fwhm_max = 1000.;
  }
  GenMassFitMean = fs->make<TH1F>("GenMassFitMean", "Signal mean", 100, mean_min, mean_max);
  RecMassFitMean = fs->make<TH1F>("RecMassFitMean", "Signal mean", 100, mean_min, mean_max);
  GenMassFitFwhm = fs->make<TH1F>("GenMassFitFwhm", "Signal FWHM", 100, fwhm_min, fwhm_max);
  RecMassFitFwhm = fs->make<TH1F>("RecMassFitFwhm", "Signal FWHM", 100, fwhm_min, fwhm_max);

  // Number of signal and background events
  GenMassFitNevtTot = fs->make<TH1F>("Nevt_all_gen", "Nevt_all", 100, 0., 25.);
  GenMassFitNsigTot = fs->make<TH1F>("Nsig_all_gen", "Nsig_all", 100, 0., 10.);
  GenMassFitNbkgTot = fs->make<TH1F>("Nbkg_all_gen", "Nbkg_all", 100, 0., 10.);

  RecMassFitNevtTot = fs->make<TH1F>("Nevt_all_rec", "Nevt_all", 100, 0., 50.);
  RecMassFitNsigTot = fs->make<TH1F>("Nsig_all_rec", "Nsig_all", 100, 0., 50.);
  RecMassFitNbkgTot = fs->make<TH1F>("Nbkg_all_rec", "Nbkg_all", 100, 0., 50.);

  GenMassFitNobs = fs->make<TH1F>("Nobs_gen", "Nobs", 100, 0., 25.); // observed
  GenMassFitNfit = fs->make<TH1F>("Nfit_gen", "Nfit", 100, 0., 25.); // fitted
  GenMassFitNsig = fs->make<TH1F>("Nsig_gen", "Nsig", 100, 0., 25.);
  GenMassFitNbkg = fs->make<TH1F>("Nbkg_gen", "Nbkg", 100, 0., 10.);

  RecMassFitNobs = fs->make<TH1F>("Nobs_rec", "Nobs", 100, 0., 25.);
  RecMassFitNfit = fs->make<TH1F>("Nfit_rec", "Nfit", 100, 0., 25.);
  RecMassFitNsig = fs->make<TH1F>("Nsig_rec", "Nsig", 100, 0., 25.);
  RecMassFitNbkg = fs->make<TH1F>("Nbkg_rec", "Nbkg", 100, 0., 10.);

  // Kolmogorov-Smirnov goodness-of-fit test
  GenMassFitKSDistSig = fs->make<TH1F>("DistSig_gen", "DistSig", 100, 0., 1.);
  GenMassFitKSDistBkg = fs->make<TH1F>("DistBkg_gen", "DistBkg", 100, 0., 1.);
  RecMassFitKSDistSig = fs->make<TH1F>("DistSig_rec", "DistSig", 100, 0., 1.);
  RecMassFitKSDistBkg = fs->make<TH1F>("DistBkg_rec", "DistBkg", 100, 0., 1.);

  GenMassFitKSProbSig = fs->make<TH1F>("ProbSig_gen", "ProbSig", 100, 0., 1.);
  GenMassFitKSProbBkg = fs->make<TH1F>("ProbBkg_gen", "ProbBkg", 100, 0., 1.);
  RecMassFitKSProbSig = fs->make<TH1F>("ProbSig_rec", "ProbSig", 100, 0., 1.);
  RecMassFitKSProbBkg = fs->make<TH1F>("ProbBkg_rec", "ProbBkg", 100, 0., 1.);

  // Counting and likelihood significances
  GenMassFitSc1  = fs->make<TH1F>("Sc1_gen",  "Sc1",  100, 0., 50.);
  GenMassFitSc2  = fs->make<TH1F>("Sc2_gen",  "Sc2",  100, 0., 25.);
  GenMassFitSc12 = fs->make<TH1F>("Sc12_gen", "Sc12", 100, 0., 25.);
  GenMassFitScl  = fs->make<TH1F>("ScL_gen",  "ScL",  100, 0., 25.);
  GenMassFitSl2  = fs->make<TH1F>("SL_gen",   "SL",   100, 0., 25.);

  double smax = 25.;
  if      (resMassId == 3000) smax = 15.;
  else if (resMassId == 5000) smax = 10.;
  RecMassFitSc1  = fs->make<TH1F>("Sc1_rec",  "Sc1",  100, 0., smax);
  RecMassFitSc2  = fs->make<TH1F>("Sc2_rec",  "Sc2",  100, 0., smax);
  RecMassFitSc12 = fs->make<TH1F>("Sc12_rec", "Sc12", 100, 0., smax);
  RecMassFitScl  = fs->make<TH1F>("ScL_rec",  "ScL",  100, 0., smax);
  RecMassFitSl2  = fs->make<TH1F>("SL_rec",   "SL",   100, 0., smax);
}

void Zprime2muMassReach::fitMass() {
  // Do unbinned maximum likelihood fit on the dilepton mass data
  bool   debug = true;
  int    igen, irec, irec_first[NF] = {0};
  int    ngen[NF], nrec[NF], ngen_tot[NF] = {0}, nrec_tot[NF] = {0};
  double genmass[MASS_FIT_ARRAY_SIZE], recmass[MASS_FIT_ARRAY_SIZE];
  double *genweight = 0, *recweight = 0; // for unit weights
  double w[NF];
  double genwght[MASS_FIT_ARRAY_SIZE], recwght[MASS_FIT_ARRAY_SIZE];
  double logML_gen_sb = -999.0, logML_gen_b = -999.0;
  double logML_rec_sb = -999.0, logML_rec_b = -999.0;

  // Print all events in input arrays
  if (debug) {
    std::ostringstream out;
    for (int idx = 0; idx < NF; idx++) {
      out << "\n" << " List of events for file # " << idx << ":"  << "\n";
      out << " gen. # | gen. mass | rec. # | rec. mass |" << "\n";
      int ifirst = 0;
      for (igen = 0; igen < nfit_genmass_used[idx]; igen++) {
	out << std::setw(7) << igen << " | "
	    << std::setw(9) << fit_genmass[idx][igen] << " |";
	for (irec = ifirst; irec < nfit_recmass_used[idx]; irec++) {
	  if (fit_recevent[idx][irec] == fit_genevent[idx][igen]) {
	    out << std::setw(7) << irec << " | "
		<< std::setw(9) << fit_recmass[idx][irec] << " |";
	    ifirst = irec+1;
	    break;
	  }
	  if (fit_recevent[idx][irec] > fit_genevent[idx][igen]) break;
	}
	out << "\n";
      }
    }
    LogTrace("fitMass") << out.str();
  }

  // Check samples and get weights relative to the main (first) sample
  for (int idx = 0; idx < NF; idx++) {
    //if (fit_sid[idx] < 0. || fit_sid[idx] >= NSAMPLES ||
    if (XSec[fit_sid[idx]] < 0.) {
      edm::LogWarning("fitMass") 
	<< "+++ negative XSec entered for file " << fit_sid[idx]
	<< "... refusing to do fits +++";
      return;
    }
    // XSec is sigma*Br of a sample
    w[idx] =
      (XSec[fit_sid[idx]]*KFactor[fit_sid[idx]]*effFilter[fit_sid[idx]])/
      (XSec[fit_sid[0]]*KFactor[fit_sid[0]]*effFilter[fit_sid[0]]);
    LogTrace("fitMass")
      << "fitMass(): file # " << idx // << " is " << SAMPLES[fit_sid[idx]]
      << "; sigma*Br = " << XSec[fit_sid[idx]] << " fb"
      << "; Kfactor = " << KFactor[fit_sid[idx]]
      << "; effFilter = " << effFilter[fit_sid[idx]]
      << "; weight = " << w[idx];
  }

  // Fill histograms that can be used to generate very large numbers of
  // experiments (not used - except for drawing - if kGenuineEvents is true).
  double gen_mass;
  for (int idx = 0; idx < NF; idx++) {
    for (igen = 0; igen < nfit_genmass_used[idx]; igen++) {
      gen_mass = fit_genmass[idx][igen];
      if (gen_mass >= lowerGenMass[fit_sid[idx]] &&
	  gen_mass <  upperGenMass[fit_sid[idx]]) {
	GenMassFitData->Fill(gen_mass, w[idx]);
      }
    }
    for (irec = 0; irec < nfit_recmass_used[idx]; irec++) {
      gen_mass = fit_genmass[idx][fit_recevent[idx][irec]-1];
      if (gen_mass >= lowerGenMass[fit_sid[idx]] &&
	  gen_mass <  upperGenMass[fit_sid[idx]]) {
	RecMassFitData->Fill(fit_recmass[idx][irec], w[idx]);
      }
    }
  }

  // Smooth DY histogram according to a function with its best-fit parameters
  // if required.  Had been used for 5 TeV samples only.
  TF1 *smoother = 0;
  if (kSmoothedSample && resMassId == 5000) {
    smoother = new TF1("smoother", lorentzianPlusExpbckgNorm,
		       massWin[0], massWin[1], 5);
    smoother->SetParNames("SigFr", "FWHM", "Mean", "BckgSlope", "BckgInt");
    smoother->SetParameters(0,       150.,  5000.,  -0.001598,    5.17234);
  // This histogram below is needed only for checks: we sample events directly
  // from smoother function.
  /* for (int iexp = 0; iexp < 50000; iexp++)
        RecMassFitDataSmoothed->Fill(smoother->GetRandom(), 1.); */
  }

  // Overall selection efficiency.  Needed to set Poisson mean for the number
  // of reconstructed events when events are sampled from RecMassFitData
  // histogram.  Calculated from nfit_genmass_used[0] instead of
  // GenMassFitData->GetSumOfWeights() because in the latter some
  // of the events fall below lowerGenMass (and therefore are not included),
  // whereas nfit_genmass_used[0] contains all events corresponding to
  // a given XSEC.
  //double eff = RecMassFitData->GetSumOfWeights()/nfit_genmass_used[0];
  double eff = nfit_recmass_used[0]/float(nfit_genmass_used[0]);
  std::ostringstream out;
  out
    << " Input to the mass fits: " << "\n"
    << "  generated events in main sample = "     << nfit_genmass_used[0]
    << "; reconstructed events in main sample = " << nfit_recmass_used[0]
    << "; efficiency = " << eff << "\n"
    << "  events generated in mass interval = "
    << GenMassFitData->GetSumOfWeights()
    << "; events reconstructed in mass interval = "
    << RecMassFitData->GetSumOfWeights();
  LogTrace("fitMass") << out.str();

  // Poisson means for numbers of generated and reconstructed events.
  // The latter one is used only when kGenuineEvents is false.
  double pmean_gen =
    XSec[fit_sid[0]]*KFactor[fit_sid[0]]*effFilter[fit_sid[0]]*intLumi;
  double pmean_rec = pmean_gen*eff;

  // Modify random generator seed if required.  The seed is set to the current
  // machine clock.
  if (kRandomSeed) gRandom->SetSeed(0);

  // Loop over kNexp "toy" experiments, prepare arrays with the "data" and fit.
  LogTrace("fitMass") 
    << " Start a loop over " << kNexp << " experiments for an"
    << " int. luminosity = " << intLumi << " (Ngen_mean = " << pmean_gen
    << ")" << " Nrec_mean = " << pmean_rec;

  for (unsigned int iexp = 1; iexp <= kNexp; iexp++) {
    int ngen_sum = 0, nrec_sum = 0;
    /* To track memory leaks in ROOT.  You also need to edit or add the
       two following lines in the file .rootrc:
       Root.MemStat:    1
       Root.ObjectStat: 1
    */
    // gObjectTable->Print();

    if (debug || (iexp % 10 == 0))
      LogTrace("fitMass") << "\n*** Experiment # " << iexp;

    if (kGenuineEvents) {
      // Use genuine simulated events

      for (int idx = 0; idx < NF; idx++) {
	//ngen[idx] = nGenEvents[1]; // if use entire available data samples
	if (idx == 0) ngen[0]   = gRandom->Poisson(pmean_gen);
	else          ngen[idx] = ngen[0];

	// Check whether we have sufficient events to perform the experiment
	if (ngen_tot[idx] + ngen[idx] > nfit_genmass_used[idx]) {
	  bool quit = false;
	  std::ostringstream out;
	  out << "+++ Not enough generated events for experiment # " << iexp
	      << " sample # " << idx << "! +++" << "\n"
	      << "Need " << ngen[idx] << " events, but only "
	      << nfit_genmass_used[idx]-ngen_tot[idx] << " are available.";
	  if (idx == 0) { // signal; exit
	    out << " Exiting fits... ";
	    quit = true;
	  }
	  else { // auxilliary backgrounds; try to continue
	    if (ngen[idx] > nfit_genmass_used[idx]) {
	      out << "\n+++ Not enough generated events for a single"
		  << " experiment. Exiting ...";
	      quit = true;
	    }
	    else {
	      out << " Use the same array from the beginning... ";
	      ngen_tot[idx] = 0;
	      nrec_tot[idx] = 0;
	    }
	  }
	  edm::LogWarning("fitMass") << out.str();
	  if (quit) return;
	}

	// Generated events
	int ngen_used = 0;
	int igen_first = ngen_tot[idx];
	int igen_last  = ngen_tot[idx] + ngen[idx] - 1;
	for (igen = igen_first; igen <= igen_last; igen++) {
	  double gen_mass = fit_genmass[idx][igen];
	  if (gen_mass >= lowerGenMass[fit_sid[idx]] &&
	      gen_mass <  upperGenMass[fit_sid[idx]]) {
	    if (ngen_sum < MASS_FIT_ARRAY_SIZE) {
	      genmass[ngen_sum] = gen_mass;
	      genwght[ngen_sum] = w[idx];
	      ngen_used++;  ngen_sum++;
	    }
	    else {
	      edm::LogWarning("fitMass")
		<< "+++ MASS_FIT_ARRAY_SIZE is too small "
		<< "to keep all events in experiment # " << iexp << "+++";
	      return;
	    }
	  }
	}

	// Reconstructed events
	int nrec_used = 0;
	nrec[idx] = 0;
	for (igen = igen_first; igen <= igen_last; igen++) {
	  for (irec = irec_first[idx]; irec < nfit_recmass_used[idx]; irec++) {
	    if (fit_recevent[idx][irec] == fit_genevent[idx][igen]) {
	      nrec[idx]++;
	      double gen_mass = fit_genmass[idx][igen];
	      if (gen_mass >= lowerGenMass[fit_sid[idx]] &&
		  gen_mass <  upperGenMass[fit_sid[idx]]) {
		if (nrec_sum < MASS_FIT_ARRAY_SIZE) {
		  recmass[nrec_sum] = fit_recmass[idx][irec];
		  recwght[nrec_sum] = w[idx];
		  nrec_used++;  nrec_sum++;
		}
		else {
		  edm::LogWarning("fitMass")
		    << "+++ MASS_FIT_ARRAY_SIZE is too small "
		    << "to keep all events in experiment # " << iexp << "+++";
		  return;
		}
	      }
	      irec_first[idx] = irec+1;
	      break;
	    }
	    if (fit_recevent[idx][irec] > fit_genevent[idx][igen]) break;
	  }
	}

	// Use the following loop to select fixed number of rec. events
	// (for checks)
	/* for (int ind = 0; ind < 100; ind++) {
	   recmass[nrec_sum] = fit_recmass[idx][nrec_tot[idx]+ind];
	   recwght[nrec_sum] = w[idx];
	   nrec[idx]++;
	   nrec_used++;  nrec_sum++;
	   } */

	// Set pointers to filled arrays with weights
	genweight = genwght;
	recweight = recwght;

	if (debug) {
	  std::ostringstream out;
	  out << "idx = " << idx << "; events: gen = " << ngen[idx]
	      << " [" << igen_first << "-" << igen_last << "]"
	      << " used = " << ngen_used << " used total = " << ngen_sum
	      << " read = " << ngen_tot[idx] + ngen[idx] << "\n";
	  int ievt = 0;
	  for (igen = ngen_sum - ngen_used; igen < ngen_sum; igen++) {
	    out << " " << std::setw(7) << genmass[igen];
	    if ((++ievt)%9 == 0 && igen < ngen_sum-1) out << "\n";
	  }
	  out << "\n";
	  out << "idx = " << idx << "; events: rec = " << nrec[idx]
	      << " [" << nrec_tot[idx] << "-" << nrec_tot[idx] + nrec[idx] - 1
	      << "]"
	      << " used = " << nrec_used << " used total = " << nrec_sum
	      << " read = " << nrec_tot[idx] + nrec[idx] << "\n";
	  ievt = 0;
	  for (irec = nrec_sum - nrec_used; irec < nrec_sum; irec++) {
	    out << " " << std::setw(7) << recmass[irec];
	    if ((++ievt)%9 == 0 && irec < nrec_sum-1) out << "\n";
	  }
	  LogTrace("fitMass") << out.str();
	}

	ngen_tot[idx] += ngen[idx];
	nrec_tot[idx] += nrec[idx];
      }

      if (debug) {
	if (fitGenMass) {
	  LogTrace("fitMass") 
	    << " Gen. events: used = " << ngen_sum << "; total (s1+s2) = "
	    << ngen_tot[0] << " + " << ngen_tot[1];
	}
	if (fitRecMass) {
	  LogTrace("fitMass")
	    << " Rec. events: used = " << nrec_sum << "; total (s1+s2) = "
	    << nrec_tot[0] << " + " << nrec_tot[1];
	}
      }
    }
    else {
      // Events distributed according to the histogram contents.
      // Generated events
      if (fitGenMass) {
	ngen_sum = gRandom->Poisson(pmean_gen);
	if (ngen_sum < MASS_FIT_ARRAY_SIZE) {
	  for (igen = 0; igen < ngen_sum; igen++) {
	    genmass[igen] = GenMassFitData->GetRandom();
	  }
	}
	else {
	  edm::LogWarning("fitMass")
	    << "+++ MASS_FIT_ARRAY_SIZE is too small "
	    << "to keep all events in experiment # " << iexp << " +++";
	  return;
	}
	ngen_tot[0] += ngen_sum;
      }

      // Reconstructed events
      if (fitRecMass) {
	nrec_sum = gRandom->Poisson(pmean_rec);
	if (nrec_sum < MASS_FIT_ARRAY_SIZE) {
	  for (irec = 0; irec < nrec_sum; irec++) {
	    if (!kSmoothedSample) recmass[irec] = RecMassFitData->GetRandom();
	    else                  recmass[irec] = smoother->GetRandom();
	  }
	}
	else {
	  edm::LogWarning("fitMass")
	    << "+++ MASS_FIT_ARRAY_SIZE is too small "
	    << "to keep all events in experiment # " << iexp << " +++";
	  return;
	}
	nrec_tot[0] += nrec_sum;
      }
    }

    if (fitGenMass) {
      sc1_gen = sc2_gen = sc12_gen = scl_gen = sl2_gen = 0.;
      /*
      logML_gen_sb = unbinnedMassFit(nfit_genmass_used[idx],
                                     fit_genmass[idx], genweight, "lor_expbg");
      logML_gen_b  = unbinnedMassFit(nfit_genmass_used[idx],
                                     fit_genmass[idx], genweight, "lor_expbg",
                                     true);
      */
      if (ngen_sum > 0) {
	logML_gen_sb = unbinnedMassFit(ngen_sum, genmass, genweight,
				       "lor_expbg");
	logML_gen_b  = unbinnedMassFit(ngen_sum, genmass, genweight,
				       "lor_expbg", true);
	if (logML_gen_sb >= logML_gen_b)
	  sl2_gen = sqrt(2.*(logML_gen_sb - logML_gen_b));
      }
    }
    if (fitRecMass) {
      sc1_rec = sc2_rec = sc12_rec = scl_rec = sl2_rec = 0.;
      /*
      logML_rec_sb = unbinnedMassFit(nfit_recmass_used[idx],
                                     fit_recmass[idx], recweight,
				     "lorgau_expbg");
      logML_rec_b  = unbinnedMassFit(nfit_recmass_used[idx],
                                     fit_recmass[idx], recweight,
				     "lorgau_expbg", true);
      */
      if (nrec_sum > 0) {
	logML_rec_sb = unbinnedMassFit(nrec_sum, recmass, recweight,
				       "lorgau_expbg");
	logML_rec_b  = unbinnedMassFit(nrec_sum, recmass, recweight,
				       "lorgau_expbg", true);
	if (logML_rec_sb > logML_rec_b)
	  sl2_rec = sqrt(2.*(logML_rec_sb - logML_rec_b));
      }
    }

    if (debug) {
      std::ostringstream out;
      if (fitGenMass) {
	out
	  << " logML_gen_sb = " << logML_gen_sb
	  << " logML_gen_b  = " << logML_gen_b
	  << " sl2_gen = "  << sl2_gen  << "\n"
	  << " sc1_gen = "  << sc1_gen  << " sc2_gen = " << sc2_gen
	  << " sc12_gen = " << sc12_gen << " scl_gen = " << scl_gen;
      }
      if (fitRecMass) {
	out
	  << " logML_rec_sb = " << logML_rec_sb
	  << " logML_rec_b  = " << logML_rec_b
	  << " sl2_rec = "  << sl2_rec  << "\n"
	  << " sc1_rec = "  << sc1_rec  << " sc2_rec = " << sc2_rec
	  << " sc12_rec = " << sc12_rec << " scl_rec = " << scl_rec;
      }
      LogTrace("fitMass") << out.str();
    }

    // Fill histograms
    if (fitGenMass) {
      GenMassFitSc1->Fill(sc1_gen);
      GenMassFitSc2->Fill(sc2_gen);
      GenMassFitSc12->Fill(sc12_gen);
      GenMassFitScl->Fill(scl_gen);
      GenMassFitSl2->Fill(sl2_gen);
    }
    if (fitRecMass) {
      RecMassFitSc1->Fill(sc1_rec);
      RecMassFitSc2->Fill(sc2_rec);
      RecMassFitSc12->Fill(sc12_rec);
      RecMassFitScl->Fill(scl_rec);
      RecMassFitSl2->Fill(sl2_rec);
    }
  }

  if (kSmoothedSample) delete smoother;
}

double Zprime2muMassReach::unbinnedMassFit(int nentries,
					   double *mass_data, double *weight,
					   const std::string fittype,
					   const bool bckgonly) {
  // Interface to unbinned mass fitter.

  // Input: number of events, followed by two arrays containing the pointers
  // to the mass and weights.  If weight=0, unit weights will be used; there
  // is no need to fill an array with 1's.
  // Fittype defines the fitting function. The choices are:
  //   "lor"          - simple Lorentzian,
  //   "lor_bg"       - sum of Lorentzian and polynomial background,
  //   "lor_expbg"    - the sum of Lorentzian and exponential background
  //                    with the slope determined from the fit to Drell-Yan,
  //   "lorgau_expbg" - the sum of a BW-Gaussian convolution and exponential
  //                    background with the slope determined from the fit
  //                    to Drell-Yan.
  // Bckgonly is true for the fits with number of signal events set to zero
  //   (background-only fits).

  static bool debug = true;
  double log_ML = 0.;    // log of maximum likelihood of the fit
  TF1 *fmass_unb = 0;

  // Mass window
  double mass_min = massWin[0];
  double mass_max = massWin[1];

  if (fittype == "lor") {
    // JMTBAD fixed for 3 TeV?
    mass_min = 2000.;
    mass_max = 4000.;
    fmass_unb = new TF1("fmass_unb", Lorentzian, mass_min, mass_max, 3);
    fmass_unb->SetParNames("Norm", "FWHM", "Mean");
    fmass_unb->SetParameters(1.,    300.,   2000.);
    fmass_unb->FixParameter(0, 1.); // must fix norm for unbinned fit
  }
  else if (fittype == "lor_bg") {
    fmass_unb = new TF1("fmass_unb", lorentzianPlusBckgNorm,
			mass_min, mass_max, 3);
    fmass_unb->SetParNames("Norm", "FWHM", "Mean");
    fmass_unb->SetParameters(0.5,    300.,  2000.);
  }
  else if (fittype == "lor_expbg" || fittype == "lorgau_expbg") {
    // Slope of the exponential background; fixed from the fits to Drell-Yan
    double bkgslope = 0.;
    if (fittype == "lor_expbg") {
      bkgslope = -0.002828;
      fmass_unb = new TF1("fmass_unb", lorentzianPlusExpbckgNorm,
			  mass_min, mass_max, 5);
      fmass_unb->SetParNames("SigFr", "FWHM", "Mean", "BckgSlope", "BckgInt");
      fmass_unb->SetParameters(0.5,     100.,  2500.,  bkgslope,      5.07);
    }
    else if (fittype == "lorgau_expbg") {
      if (resMassId == 1000 || resMassId == 1500 || resMassId == 2000) {
	bkgslope = -2.00;           // for 1, 1.5 and 2 TeV
      }
      else if (resMassId == 3000) { // for 3 TeV
	bkgslope = -1.993;
      }
      else if (resMassId == 5000) { // for 5 TeV
	bkgslope = -2.064;
      }
      int npar = 5;
// JMTBAD should these different "Extended Maximum Likelihood" options
// be controllable from the config file, or kept as compile-time defines?
#if defined(EML_1) || defined(EML_2)
      npar = 6;
#endif
      fmass_unb = new TF1("fmass_unb", lorengauPlusExpbckgNorm,
      		  mass_min, mass_max, npar);
#ifndef EML_2
      fmass_unb->SetParNames("SigFr", "FWHM", "Mean", "BckgSlope", "BckgInt");
#else
      fmass_unb->SetParNames("Nsig",  "FWHM", "Mean", "BckgSlope", "BckgInt");
#endif
      if (resMassId == 1000)
	fmass_unb->SetParameters(0.5,   100.,   900.,  bkgslope,   1.);
      else if (resMassId == 1500)
	fmass_unb->SetParameters(0.5,   100.,  1500.,  bkgslope,   1.);
      else if (resMassId == 2000)
	fmass_unb->SetParameters(0.5,   100.,  2000.,  bkgslope,   1.);
      else if (resMassId == 3000)
	fmass_unb->SetParameters(0.5,   300.,  2500.,  bkgslope,   11.2);
      else if (resMassId == 5000)
	fmass_unb->SetParameters(0.5,   300.,  4500.,  bkgslope,   1.);
#if defined(EML_1) || defined(EML_2)
      fmass_unb->SetParameter(5, 1.);
#ifdef EML_1
      fmass_unb->SetParName(5, "Norm");
#else
      fmass_unb->SetParName(5, "Nbkg");
      fmass_unb->SetParLimits(5, 0., 1000.);
#endif
#endif
    }

    // Fix the slope of the background, except when fitting DY to get the slope
    if (!kBackgroundFit || (kBackgroundFit && !bckgonly))
      fmass_unb->FixParameter(3, bkgslope);

    // Integral of the background function in a chosen mass range
    TF1 *fbkg = new TF1("fbkg", expBckg, mass_min, mass_max, 3);
    fbkg->SetParameters(1., bkgslope, 1.);
    fmass_unb->FixParameter(4, fbkg->Integral(mass_min, mass_max));
    delete fbkg;
  }
  else if (fittype == "bckg") {
    // Fits to Drell-Yan to get a slope of the exponential background.
    // This is the test mode; the same result can be obtained with two
    // modes above, by fixing the first parameter to zero.
    fmass_unb = new TF1("fmass_unb", expBckgNorm, mass_min, mass_max, 3);
    fmass_unb->SetParNames("Norm", "BckgSlope", "BckgInt");
    fmass_unb->SetParameters(1.,     -0.003,        1.);
    fmass_unb->FixParameter(0,  1.); // fix both norms
    fmass_unb->FixParameter(2,  1.);
  }
  else {
    edm::LogWarning("unbinnedMassFit")
      << "+++ Unknown fit type " << fittype << "! +++";
    return log_ML;
  }

  // Fix parameters or set limits
  if (fittype != "bckg") {
    if (bckgonly) {
      fmass_unb->FixParameter(0, 0.);   // set Ns to zero
      fmass_unb->FixParameter(1, 100.); // only to speed-up the fit
      fmass_unb->SetParLimits(2, mass_min, mass_max);
    }
    else {
#ifndef EML_2
      fmass_unb->SetParLimits(0, 0., 1.);
#else
      fmass_unb->SetParLimits(0, 0., 1000.);
#endif

      // Fix FWHM if required
      if (kFixedFWHM) {
	if      (resMassId == 1000) fmass_unb->FixParameter(1,  30.);
	else if (resMassId == 1500) fmass_unb->FixParameter(1,  45.);
	else if (resMassId == 2000) fmass_unb->FixParameter(1,  60.);
	else if (resMassId == 3000) fmass_unb->FixParameter(1, 100.);
	else if (resMassId == 5000) fmass_unb->FixParameter(1, 150.);
      }
      else {
	if      (resMassId == 1000) fmass_unb->SetParLimits(1,  1.,  200.);
	else if (resMassId == 1500) fmass_unb->SetParLimits(1,  5.,  500.);
	else if (resMassId == 2000) fmass_unb->SetParLimits(1, 10.,  500.);
	else if (resMassId == 3000) fmass_unb->SetParLimits(1, 30., 1000.);
	else if (resMassId == 5000) fmass_unb->SetParLimits(1, 50., 1000.);
      }

      // Fix mean mass if required
      if (kFixedMass) {
	if      (resMassId == 1000) fmass_unb->FixParameter(2, 1000.);
	else if (resMassId == 1500) fmass_unb->FixParameter(2, 1500.);
	else if (resMassId == 2000) fmass_unb->FixParameter(2, 2000.);
	else if (resMassId == 3000) fmass_unb->FixParameter(2, 3000.);
	else if (resMassId == 5000) fmass_unb->FixParameter(2, 4900.);
      }
      else {
	fmass_unb->SetParLimits(2, mass_min, mass_max);
      }
    }
  }

  if (debug) {
    std::ostringstream out;
    out
      << "Fit option: " << fittype << ", " << (bckgonly ? "B" : "S+B") << "\n"
      << "Fit region: " << mass_min << " to " << mass_max << " GeV" << "\n";
    LogTrace("unbinnedMassFit") << out.str();
  }

  // Filter data arrays for the fitter: only use events within the mass
  // interval of interest.  Evt structure contains invariant mass and weight.
  double sumw = 0.;
  int nevents = 0;
  std::vector<evt> events(nentries);
  for (int ievt = 0; ievt < nentries; ievt++) {
    if (mass_data[ievt] >= mass_min && mass_data[ievt] <= mass_max) {
      // event[nevents].mass = fsum->GetRandom(); // to test K-S stat.
      events[nevents].mass = mass_data[ievt];
      if (weight) events[nevents].weight = weight[ievt];
      else        events[nevents].weight = 1.;
      sumw += events[nevents].weight;
      nevents++;
    }
  }
  if (nevents != nentries) events.resize(nevents);

  if (debug) {
    std::ostringstream out;
    out << "Total events in the fit region: " << nevents << "\n";
    for (int ievt = 0; ievt < nevents; ievt++) {
      out << " " << std::setw(7) << events[ievt].mass
	  << "(" << std::setw(4) << std::setprecision(4) << events[ievt].weight
	  << std::setprecision(6) << ")";
      if ((ievt+1)%5 == 0 && ievt < nevents-1) out << "\n";
    }
    LogTrace("unbinnedMassFit") << out.str();
  }

  if (debug) {
    double sum = 0;
    double var[1], par[] = {1., 100., 3000., -0.002828, 5.07};
    double mass_step = 10.;
    for (double x = mass_min; x < mass_max; x += mass_step) {
      var[0] = x;
      if (fittype == "lor")
	sum += Lorentzian(var, par);
      else if (fittype == "lor_bg")
	sum += lorentzianPlusBckgNorm(var, par);
      else if (fittype == "lor_expbg")
	sum += lorentzianPlusExpbckgNorm(var, par);
      else if (fittype == "lorgau_expbg")
	sum += lorengauPlusExpbckgNorm(var, par);
      else if (fittype == "bckg")
	sum += expBckgNorm(var, &par[2]);
    }
    LogTrace("unbinnedMassFit")
      << "Normalization test: norm = " << sum*mass_step;
  }

  // Do it!
  Option_t *opt = (debug) ? "V" : "Q";
  int fit_status = unbinnedFitExec("fmass_unb", opt, events, log_ML);

  // Check normalization
  double fmass_int = fmass_unb->Integral(mass_min, mass_max);
  if (fabs(fmass_int - 1.) > 0.001) {
    edm::LogWarning("unbinnedMassFit")
      << "+++ integral of fmass_unb in the fit region is "
      <<  fmass_int << " +++";
  }

  if (fit_status >= 1) {
    if (debug)
      LogTrace("unbinnedMassFit")
	<< "Fit converged successfully, return status is " << fit_status;

    int    nfpars = fmass_unb->GetNpar();
    double *fpars = fmass_unb->GetParameters();

    if (debug) {
      std::ostringstream out;
      out << " Fitted values of parameters:";
      for (int ifpars = 0; ifpars < nfpars; ifpars++)
	out << "  " << fpars[ifpars];
      LogTrace("unbinnedMassFit") << out.str();
    }

    // Fill histograms
    if (!bckgonly) {
      if (fittype == "lor_expbg") {
	GenMassFitSigFr->Fill(fpars[0]);
	GenMassFitMean->Fill(fpars[2]);
	GenMassFitFwhm->Fill(fpars[1]);
      }
      else if (fittype == "lorgau_expbg") {
	RecMassFitSigFr->Fill(fpars[0]);
	RecMassFitMean->Fill(fpars[2]);
	RecMassFitFwhm->Fill(fpars[1]);
      }
    }

    if (fittype == "lor_expbg" || fittype == "lorgau_expbg") {
      if (debug) {
	LogTrace("unbinnedMassFit")
	  << " Events:          total = " << nevents
	  << "; signal = "     << nevents*fpars[0]
	  << "; background = " << nevents*(1.-fpars[0]);
	if (weight) {
	  LogTrace("unbinnedMassFit")
	    << " Weighted events: total = " << sumw
	    << "; signal = "     << sumw*fpars[0]
	    << "; background = " << sumw*(1.-fpars[0]);
	}
      }

      // Analyze and plot the result
      analyzeUnbinnedMassFits(sumw, events, nfpars, fpars,
			      fittype, bckgonly, mass_min, mass_max);
    }
  }
  else {
    edm::LogWarning("unbinnedMassFit")
      << "+++ fit has not converged: status = " << fit_status
      << ", skipping this experiment +++!";
    log_ML = 0.;
  }
  
  delete fmass_unb;
  return log_ML;
}

int Zprime2muMassReach::unbinnedFitExec(const Char_t *funcname,
					const Option_t *option,
					const std::vector<evt>& events,
					double& log_ML) {
  // Small interface
  bool useW = false;
  const int nevents = events.size();
  const int MAXSIZE = 10000;
  double data[MAXSIZE] = {0.}, wght[MAXSIZE] = {0.};
  double *nulldata2 = 0; // won't need second data array of fitter
  double *nulldata3 = 0; // won't need third  data array of fitter
  double *weight = 0;
  
  if (nevents >= MAXSIZE)
    throw cms::Exception("+++ unbinnedFitExec: MAXSIZE is too small +++\n");
			 
  for (int iev = 0; iev < nevents; iev++) {
    data[iev] = events[iev].mass;
    wght[iev] = events[iev].weight;
    if (!useW && fabs(wght[iev]-1.) > 1.e-6) useW = true;
  }
  
  UnbinnedFitter *fitter = new UnbinnedFitter();
  if (useW) weight = wght;
  int status = fitter->unbinnedFitExec("fmass_unb", option, nevents, data,
				       nulldata2, nulldata3, weight, log_ML);
  delete fitter;
  return status;
}

void Zprime2muMassReach::analyzeUnbinnedMassFits(
           const double sumw,     const std::vector<evt>& events,
           const int npars,       const double *pars,
           const std::string fittype,  const bool bckgonly,
           const double mass_min, const double mass_max) {
  
  static bool debug = true;
  static int nexp_gen_b  = 0, nexp_rec_b  = 0;
  static int nexp_gen_sb = 0, nexp_rec_sb = 0, nexp_unknown = 0;
  TF1 *fsum = 0, *fbkg = 0;
  
  const int nbins = 100;
  double bin_width = (mass_max-mass_min)/nbins;

#if !defined(EML_1) && !defined(EML_2)
  const int npsum = 6; // total number of parameters
#else
  const int npsum = 6+1;
#endif
  const int npbkg = 3; // number of parameters controlling background
  double par_sum[npsum], par_bkg[npbkg];
  double nobs = 0, ntot = 0, nbkg = 0, nsig = 0;
  if (fittype == "lor_expbg" || fittype == "lorgau_expbg") {
    if (npars != npsum-1) {
      edm::LogWarning("analyzeUMF")
	<< "+++ unexpected number of parameters!  Exiting... +++";
      return;
    }

    if (fittype == "lor_expbg")         // generated mass fit
      fsum = new TF1("fsum", lorentzianPlusExpbckgN, mass_min, mass_max,
		     npsum);
    else if (fittype == "lorgau_expbg") // reconstructed mass fit
      fsum = new TF1("fsum", lorengauPlusExpbckgN,   mass_min, mass_max,
		     npsum);
    par_sum[0] = bin_width; // normalization
#if !defined (EML_1) && !defined (EML_2)
    par_sum[0] *= sumw;
#endif
    for (int ip = 1; ip < npsum; ip++) par_sum[ip] = pars[ip-1];
    fsum->SetParameters(par_sum);

    fbkg = new TF1("fbkg", expBckg, mass_min, mass_max, npbkg);
#ifndef EML_2
    par_bkg[0] = (1.-pars[0])*sumw*bin_width; // normalization
#else
    par_bkg[0] = pars[5]*bin_width;
#endif
    par_bkg[1] = pars[3];
    par_bkg[2] = pars[4];
    fbkg->SetParameters(par_bkg);

    if (debug) {
      TF1 *fsig = new TF1("fsig", lorengauPlusExpbckgN, mass_min, mass_max,
			  npsum);
      par_sum[1] = 1.; // fr_bkg = 0
      fsig->SetParameters(par_sum);
      LogTrace("analyzeUMF") 
	<< " TEST: total = " << fsum->Integral(mass_min, mass_max)/bin_width
	<< " signal = "
	<< pars[0]*fsig->Integral(mass_min, mass_max)/bin_width
	<< " bkgtot = "  << fbkg->Integral(mass_min, mass_max)/bin_width
	<< " nevents = " << events.size() << " sumw = " << sumw;
      delete fsig;
    }

    // Kolmogorov-Smirnov test
    double d = -999.0, prob = -999.0;
    ks(events, fsum, mass_min, mass_max, d, prob);
    if (debug)
      LogTrace("analyzeUMF") << " Kolmogorov-Smirnov test: dist = " << d
			     << " prob = " << prob;
    if (!bckgonly) {
      if (fittype == "lor_expbg") {         // generated mass fit
	GenMassFitKSDistSig->Fill(d);
	GenMassFitKSProbSig->Fill(prob);
      }
      else if (fittype == "lorgau_expbg") { // reconstructed mass fit
	RecMassFitKSDistSig->Fill(d);
	RecMassFitKSProbSig->Fill(prob);
      }
    }
    else {
      if (fittype == "lor_expbg") {         // generated mass fit
	GenMassFitKSDistBkg->Fill(d);
	GenMassFitKSProbBkg->Fill(prob);
      }
      else if (fittype == "lorgau_expbg") { // reconstructed mass fit
	RecMassFitKSDistBkg->Fill(d);
	RecMassFitKSProbBkg->Fill(prob);
      }
    }

    if (!bckgonly) {
      // Number of signal and background events in a certain mass range
      // around the peak
      double fwhm      = par_sum[2];
      double mass_mean = par_sum[3];
      double int_sigma = fwhm/fwhm_over_sigma;
      double res_sigma = 0.0, width = 100.0, width_sigmas = 3.0;
      if (fittype == "lor_expbg") {
	width = int_sigma;
	width_sigmas = 5.0;
      }
      else if (fittype == "lorgau_expbg") {
	res_sigma = mass_resolution(&par_sum[3], par_sum); // f(mass)
	width = sqrt(int_sigma*int_sigma + res_sigma*res_sigma);
	width_sigmas = 2.0;
      }
      double mass_low  = mass_mean - width_sigmas*width;
      double mass_high = mass_mean + width_sigmas*width;
      if (mass_low  < mass_min) mass_low  = mass_min;
      if (mass_high > mass_max) mass_high = mass_max;

      // Total, signal, and background events between mass_low and mass_high
      // from the fit
      ntot = fsum->Integral(mass_low, mass_high)/bin_width;
      nbkg = fbkg->Integral(mass_low, mass_high)/bin_width;
      nsig = ntot - nbkg;

      // Number of events *observed* between mass_low and mass_high
      std::vector<evt>::const_iterator p;
      for (p = events.begin(); p != events.end(); p++) {
	if (p->mass > mass_low && p->mass < mass_high) nobs += p->weight;
      }

      if (debug) {
	std::ostringstream out;
	out << " Width:  FWHM = " << fwhm << " (sigma = " << int_sigma;
	if (fittype == "lor_expbg")
	  out << ")";
	else
	  out << "); sigma_Gauss = " << res_sigma;
	out << "; total = " << width << "\n";
	out << " Mass:   mean = " << mass_mean
	    << "; low = " << mass_low << "; high = " << mass_high << "\n";
	out << " Events within +/- " << width_sigmas << " sigma: " << "\n";
	out << "         obs = " << nobs << " total = " << ntot
	    << "; signal = " << nsig << "; background = " << nbkg;
	LogTrace("analyzeUMF") << out.str();
      }

      if (fittype == "lor_expbg") {
	GenMassFitNevtTot->Fill(sumw);
	GenMassFitNsigTot->Fill(sumw*pars[0]);
	GenMassFitNbkgTot->Fill(sumw*(1.-pars[0]));

	GenMassFitNobs->Fill(nobs);
	GenMassFitNfit->Fill(ntot);
	GenMassFitNsig->Fill(nsig);
	GenMassFitNbkg->Fill(nbkg);
      }
      else if (fittype == "lorgau_expbg") {
	RecMassFitNevtTot->Fill(sumw);
	RecMassFitNsigTot->Fill(sumw*pars[0]);
	RecMassFitNbkgTot->Fill(sumw*(1.-pars[0]));

	RecMassFitNobs->Fill(nobs);
	RecMassFitNfit->Fill(ntot);
	RecMassFitNsig->Fill(nsig);
	RecMassFitNbkg->Fill(nbkg);
      }

      // Significance from counting methods
      double sc1 = -999, sc2 = -999, sc12 = -999, scl = -999;
      if (nbkg > 0.)  sc1  = nsig/sqrt(nbkg);
      if (ntot > 0.)  sc2  = nsig/sqrt(ntot);
      if (nbkg >= 0.) sc12 = 2.*(sqrt(ntot) - sqrt(nbkg));
      if (nbkg > 0.) {
	scl = 2.*(ntot*TMath::Log(ntot/nbkg)-nsig);
	if (scl > 0.) scl = sqrt(scl);
	else          scl = 0.;
      }
      if (fittype == "lor_expbg") {
	sc1_gen = sc1;  sc2_gen = sc2;  sc12_gen = sc12;  scl_gen = scl;
      }
      else if (fittype == "lorgau_expbg") {
	sc1_rec = sc1;  sc2_rec = sc2;  sc12_rec = sc12;  scl_rec = scl;
      }
    }
  }

  if (kExpPlots) {
    int   iexp = 0;
    std::string filename, htitle;
    if (fittype == "lor_expbg") {
      if (bckgonly) {filename = "genmass_b";  iexp = ++nexp_gen_b;}
      else          {filename = "genmass_sb"; iexp = ++nexp_gen_sb;}
      htitle = "Generated dilepton mass";
    }
    else if (fittype == "lorgau_expbg") {
      if (bckgonly) {filename = "recmass_b";  iexp = ++nexp_rec_b;}
      else          {filename = "recmass_sb"; iexp = ++nexp_rec_sb;}
      htitle = "Reconstructed dilepton mass";
    }
    else {
      edm::LogWarning("analyzeUMF")
	<< "+++ Warning in analyzeUnbinnedMassFits: fit function is not"
	<< " defined. Only the histograms will be plotted. +++";
      filename = "tmp";  iexp = ++nexp_unknown;
      htitle   = "Dilepton mass";
    }

    if (iexp > 10) { // plot only a few first experiments
      kExpPlots = false;
      return;
    }
    // if (!bckgonly && nbkg > 1.) kExpPlots = true;

    // if (kExpPlots) {
    TH1F *massHisto = fs->make<TH1F>("massHisto", htitle.c_str(), nbins, mass_min, mass_max);
    std::vector<evt>::const_iterator p;
    for (p = events.begin(); p != events.end(); p++) {
      massHisto->Fill(p->mass, p->weight);
    }

    if (!bckgonly)
      LogTrace("analyzeUMF") 
	<< " Binned K-S prob: " << massHisto->KolmogorovTest(RecMassFitData);

    TCanvas *c1 = new TCanvas("c4", "", 0, 0, 500, 700); // for .ps
    // TCanvas *c1 = new TCanvas("c4", "", 0, 0, 500, 500); // for .eps

    std::ostringstream s_iexp;
    s_iexp << "_exp" << iexp;
    system("mkdir -p mass_plots");
    std::string psname = "./mass_plots/" + filename + s_iexp.str() + ".ps";
    TPostScript *ps = new TPostScript(psname.c_str(), 111);

    std::ostringstream strlumi, strbinw;
    strlumi << intLumi;
    strbinw << bin_width;

    gStyle->SetOptStat(1111);
    ps->NewPage();
    // c1->Clear();
    c1->cd(0);
    massHisto->SetTitle("");
    massHisto->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass (GeV)");
    std::string ytit = "Events/" + strbinw.str() + " GeV/"
      + strlumi.str() + " fb^{-1}";
    massHisto->GetYaxis()->SetTitle(ytit.c_str());
    massHisto->GetYaxis()->SetTitleOffset(1.15);
    if (fittype == "lor_expbg")
      massHisto->SetMaximum(1.2*massHisto->GetMaximum());
    massHisto->Draw();
    if (fittype == "lor_expbg" || fittype == "lorgau_expbg") {
      fsum->SetLineStyle(1); fsum->SetLineWidth(2); fsum->SetLineColor(kBlue);
      fbkg->SetLineStyle(2); fbkg->SetLineWidth(2); fbkg->SetLineColor(kGreen);
      fsum->Draw("same");
      fbkg->Draw("same");

      TF1 *fbkg_only= new TF1("fbkg_only", expBckg, mass_min, mass_max, npbkg);
      par_bkg[0] = sumw*bin_width; // normalization
      par_bkg[1] = pars[3];
      par_bkg[2] = pars[4];
      fbkg_only->SetParameters(par_bkg);
      fbkg_only->SetLineStyle(3); fbkg_only->SetLineWidth(2);
      fbkg_only->SetLineColor(kGreen);
      fbkg_only->Draw("same");
      // delete fbkg_only;

      c1->RedrawAxis();
      massHisto->Draw("same"); // to have stat. box on the top of axes
    }
    c1->Update();

    // Binned fit, to compare and check the goodness of fit
    if (kBinnedFit) {
    //if (!bckgonly) {
      if (fittype == "lor_expbg" || fittype == "lorgau_expbg") {
	ps->NewPage();
	c1->cd(0);
	TF1 *fbinned = 0;
	if (fittype == "lor_expbg")         // generated mass fit
	  fbinned = new TF1("fbinned", lorentzianPlusExpbckgN,
			    mass_min, mass_max, npsum);
	else if (fittype == "lorgau_expbg") // reconstructed mass fit
	  fbinned = new TF1("fbinned", lorengauPlusExpbckgN,
			    mass_min, mass_max, npsum);
	fbinned->SetParNames("Norm", "SigFr", "FWHM", "Mean", "BckgSlope",
			     "BckgInt");
	fbinned->SetParameters(par_sum);
	if (bckgonly) {
	  fbinned->FixParameter(1, 0.); // set SigFr to zero
	  fbinned->FixParameter(2, 100.);
	  fbinned->FixParameter(3, 3000.);
	}
	else {
	  fbinned->SetParLimits(1, 0., 1.);
	  fbinned->SetParLimits(2, 50., 1000.);
	  fbinned->SetParLimits(3, mass_min, mass_max);
	}
	// Fix the slope of the background, except when fitting DY to get
	// the slope.
	if (!kBackgroundFit || (kBackgroundFit && !bckgonly))
	  fbinned->FixParameter(4, par_sum[4]); // fixed from Drell-Yan
	fbinned->FixParameter(5, par_sum[5]); // fixed from Drell-Yan
	LogTrace("analyzeUMF")
	  << "Binned ML fit, fit option: " << fittype << ", "
	  << (bckgonly ? "B" : "S+B");
	massHisto->Fit("fbinned", "LVR");
	c1->Update();
	delete fbinned;
      }
    }

    ps->Close();
    delete ps;
    delete massHisto;
    delete c1;
  }
  delete fsum;
  //delete fbkg; //??
}
  
void Zprime2muMassReach::ks(std::vector<evt> data, TF1 *func,
			    const double& mass_min, const double& mass_max,
			    double& d, double& prob) {
  // Given a vector data, and given a user-supplied function of a
  // single variable func which is a cumulative distribution function ranging
  // from 0 (for smallest values of its argument) to 1 (for largest values
  // of its argument), this routine returns the K-S statistic d, and the
  // significance level prob. Small values of prob show that the cumulative
  // distribution function of data is significantly different from func.
  // The array data is modified by being sorted into ascending order.

  double d1, d2, dt, ff, fn, fo = 0.0;
  std::vector<evt>::const_iterator p;

  sort(data.begin(), data.end(), std::less<evt>());

  double en = 0.0, sumw = 0.0;
  for (p = data.begin(); p != data.end(); p++) {
    en += p->weight; // effective number of data points
  }
  d = 0.0;
  double norm = func->Integral(mass_min, mass_max);
  // Loop over the sorted data points.
  for (p = data.begin(); p != data.end(); p++) {
    sumw += p->weight;
    fn = sumw/en;  // Data's c.d.f. after this step.
    ff = func->Integral(mass_min, p->mass)/norm;
    d1 = fabs(fo-ff);
    d2 = fabs(fn-ff);
    dt = (d1 > d2) ? d1 : d2; // Maximum distance.
    if (dt > d) d = dt;
    fo = fn;
  }
  en = sqrt(en);
  prob = probks((en+0.12+0.11/en)*d); // Compute significance.
}

void Zprime2muMassReach::drawUnbinnedMassFitsHistos() {
  // Draw the results of unbinned maximum likelihood mass fits.
  TCanvas *c1 = new TCanvas("c5", "", 0, 0, 500, 700);

  int page = 0;
  std::ostringstream strpage;
  TPaveLabel *title;
  TPad *pad;

  TText t;
  t.SetTextFont(12);
  t.SetTextSize(0.03);

  // Gaussian with a positive mean
  TF1 *fgaus = new TF1("fgaus", "gaus");
  fgaus->SetParLimits(1, 0., 500.);

  gStyle->SetOptStat(111111);
  //gStyle->SetOptFit(0);
  //gStyle->SetOptStat(0);

  std::ostringstream strmass;
  strmass << resMassId;
  std::string psfile = resModel + strmass.str() + "_" + psFile;
  TPostScript *ps = new TPostScript(psfile.c_str(), 111);
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, "MC Mass Spectra");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad = new TPad("","", .05, .05, .95, .93);
  pad->Draw();
  pad->Divide(1,2);
  pad->cd(1);  GenMassFitData->Draw("histo");
  pad->cd(2);  RecMassFitData->Draw("histo");
  /*
    double scf =
      RecMassFitData->GetSumOfWeights()/RecMassFitDataSmoothed->GetSumOfWeights();
    RecMassFitDataSmoothed->Scale(scf);
    RecMassFitDataSmoothed->SetLineColor(kBlue);
    RecMassFitDataSmoothed->Draw("samehisto");
    LogTrace("drawUMFH")
      << " K-S prob: "
      << RecMassFitDataSmoothed->KolmogorovTest(RecMassFitData);
  */
  c1->Update();

  // Generated mass
  if (fitGenMass) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, "Generated Mass Fits");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    gStyle->SetStatW(0.11); // width of stat. box
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(1,3);
    pad->cd(1);  GenMassFitSigFr->Draw();
    if (GenMassFitSigFr->Integral() > 1) GenMassFitSigFr->Fit("gaus","LQ");
    pad->cd(2);  GenMassFitMean->Draw();
    if (GenMassFitMean->Integral() > 1)  GenMassFitMean->Fit("gaus","LQ");
    pad->cd(3);  GenMassFitFwhm->Draw();
    if (GenMassFitFwhm->Integral() > 1)  GenMassFitFwhm->Fit("fgaus","LQ");
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			   "Generated Mass Fits: all events");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(1,3);
    pad->cd(1);  GenMassFitNevtTot->Draw();
    pad->cd(2);  GenMassFitNsigTot->Draw();
    pad->cd(3);  GenMassFitNbkgTot->Draw();
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			   "Generated Mass Fits: signal region");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(2,2);
    pad->cd(1);  GenMassFitNobs->Draw();
    if (GenMassFitNobs->Integral() > 1) GenMassFitNobs->Fit("fgaus","LQ");
    pad->cd(2);  GenMassFitNfit->Draw();
    if (GenMassFitNfit->Integral() > 1) GenMassFitNfit->Fit("fgaus","LQ");
    pad->cd(3);  GenMassFitNsig->Draw();
    if (GenMassFitNsig->Integral() > 1) GenMassFitNsig->Fit("fgaus","LQ");
    pad->cd(4);  GenMassFitNbkg->Draw();
    if (GenMassFitNbkg->Integral() > 1) GenMassFitNbkg->Fit("fgaus","LQ");
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			   "Generated Mass Fits: K-S g.o.f. test");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    gStyle->SetStatW(0.25); // width of stat. box
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(2,2);
    pad->cd(1);  GenMassFitKSDistSig->Draw();
    pad->cd(2);  GenMassFitKSDistBkg->Draw();
    pad->cd(3);  GenMassFitKSProbSig->Draw();
    pad->cd(4);  GenMassFitKSProbBkg->Draw();
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			   "Generated Mass Fits: significance");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(2,3);
    pad->cd(1);  GenMassFitSc1->Draw();
    if (GenMassFitSc1->Integral() > 1)  GenMassFitSc1->Fit("gaus","LQ");
    pad->cd(2);  GenMassFitScl->Draw();
    if (GenMassFitScl->Integral() > 1)  GenMassFitScl->Fit("gaus","LQ");
    pad->cd(3);  GenMassFitSc12->Draw();
    if (GenMassFitSc12->Integral() > 1) GenMassFitSc12->Fit("gaus","LQ");
    pad->cd(4);  GenMassFitSc2->Draw();
    if (GenMassFitSc2->Integral() > 1)  GenMassFitSc2->Fit("gaus","LQ");
    pad->cd(5);  GenMassFitSl2->Draw();
    if (GenMassFitSl2->Integral() > 1)  GenMassFitSl2->Fit("gaus","LQ");
    c1->Update();
  }

  // Reconstructed mass
  if (fitRecMass) {
    //TPostScript *ps = new TPostScript("mass_results.eps", 113);
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, "Reconstructed Mass Fits");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    //TPaveText *stats = (TPaveText*)gPad->GetPrimitive("stats");
    //stats->SetTextFont(82);
    //stats->SetTextSize(0.03);
    //stats->SetFillColor(42);
    gStyle->SetStatW(0.17); // width of stat. box
    RecMassFitSigFr->GetXaxis()->SetLabelSize(0.06);
    RecMassFitSigFr->GetYaxis()->SetLabelSize(0.06);
    RecMassFitMean->GetXaxis()->SetLabelSize(0.06);
    RecMassFitMean->GetYaxis()->SetLabelSize(0.06);
    RecMassFitMean->GetXaxis()->SetTitle("m_{0}, GeV");
    RecMassFitMean->GetXaxis()->SetTitleOffset(1.5);
    RecMassFitMean->GetXaxis()->SetTitleSize(0.06);
    RecMassFitFwhm->GetXaxis()->SetLabelSize(0.06);
    RecMassFitFwhm->GetYaxis()->SetLabelSize(0.06);
    RecMassFitFwhm->GetXaxis()->SetTitleOffset(1.5);
    RecMassFitFwhm->GetXaxis()->SetTitleSize(0.06);
    RecMassFitFwhm->GetXaxis()->SetTitle("FWHM, GeV");
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(1,3);
    pad->cd(1);  RecMassFitSigFr->Draw();
    if (RecMassFitSigFr->Integral() > 1) RecMassFitSigFr->Fit("gaus","LQ");
    pad->cd(2);  RecMassFitMean->Draw();
    if (RecMassFitMean->Integral() > 1)  RecMassFitMean->Fit("gaus","LQ");
    pad->cd(3);  RecMassFitFwhm->Draw();
    if (RecMassFitFwhm->Integral() > 1)  RecMassFitFwhm->Fit("fgaus","LQ");
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			   "Reconstructed Mass Fits: all events");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(1,3);
    pad->cd(1);  RecMassFitNevtTot->Draw();
    pad->cd(2);  RecMassFitNsigTot->Draw();
    pad->cd(3);  RecMassFitNbkgTot->Draw();
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			   "Reconstructed Mass Fits: signal region");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(2,2);
    pad->cd(1);  RecMassFitNobs->Draw();
    if (!kDYEvents && RecMassFitNobs->Integral() > 1)
      RecMassFitNobs->Fit("fgaus","LQ");
    pad->cd(2);  RecMassFitNfit->Draw();
    if (!kDYEvents && RecMassFitNfit->Integral() > 1)
      RecMassFitNfit->Fit("fgaus","LQ");
    pad->cd(3);  RecMassFitNsig->Draw();
    if (!kDYEvents && RecMassFitNsig->Integral() > 1)
      RecMassFitNsig->Fit("fgaus","LQ");
    pad->cd(4);  RecMassFitNbkg->Draw();
    if (!kDYEvents && RecMassFitNbkg->Integral() > 1)
      RecMassFitNbkg->Fit("fgaus","LQ");
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			   "Reconstructed Mass Fits: K-S g.o.f. test");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    gStyle->SetStatW(0.25); // width of stat. box
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(2,2);
    pad->cd(1);  RecMassFitKSDistSig->Draw();
    pad->cd(2);  RecMassFitKSDistBkg->Draw();
    pad->cd(3);  RecMassFitKSProbSig->Draw();
    pad->cd(4);  RecMassFitKSProbBkg->Draw();
    c1->Update();

    // TPostScript *ps = new TPostScript("mass_fits.eps", 113);
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			   "Reconstructed Mass Fits: significance");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    if (kDYEvents) gStyle->SetOptLogy(1);
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(2,3);
    pad->cd(1);  RecMassFitSc1->Draw();
    if (!kDYEvents && RecMassFitSc1->Integral() > 1)
      RecMassFitSc1->Fit("gaus","LQ");
    pad->cd(2);  RecMassFitScl->Draw();
    if (!kDYEvents && RecMassFitScl->Integral() > 1)
      RecMassFitScl->Fit("gaus","LQ");
    pad->cd(3);  RecMassFitSc12->Draw();
    if (!kDYEvents && RecMassFitSc12->Integral() > 1)
      RecMassFitSc12->Fit("gaus","LQ");
    pad->cd(4);  RecMassFitSc2->Draw();
    if (!kDYEvents && RecMassFitSc2->Integral() > 1)
      RecMassFitSc2->Fit("gaus","LQ");
    pad->cd(5);  RecMassFitSl2->Draw();
    if (!kDYEvents && RecMassFitSl2->Integral() > 1)
      RecMassFitSl2->Fit("gaus","LQ");
    if (kDYEvents) gStyle->SetOptLogy(0);
    c1->Update();
  }
  ps->Close();
  delete c1;
  delete fgaus;

  // Selected plots for talks/papers
#ifdef FORTALKS
  TCanvas *c2 = new TCanvas("c2", "", 0, 0, 500, 500);
  if (fitRecMass) {
    TPostScript *eps1 = new TPostScript("m0.eps", 113);
    eps1->NewPage();
    c2->Clear();
    c2->cd(0);
    //RecMassFitMean->SetTitle("1 TeV/c^{2} Z_{#psi}, 0.1 fb^{-1}");
    RecMassFitMean->SetTitle("3 TeV/c^{2} Z_{SSM}, 10 fb^{-1}");
    gStyle->SetStatW(0.17); // width of stat. box
    RecMassFitMean->GetXaxis()->SetLabelSize(0.04); // default
    RecMassFitMean->GetYaxis()->SetLabelSize(0.04);
    RecMassFitMean->GetXaxis()->SetTitle("m_{0} (GeV/c^{2})");
    RecMassFitMean->GetXaxis()->SetTitleOffset(0.95);
    RecMassFitMean->GetXaxis()->SetTitleSize(0.05); // default
    RecMassFitMean->GetXaxis()->SetNdivisions(1005);
    RecMassFitMean->GetYaxis()->SetTitle("Number of experiments");
    RecMassFitMean->GetYaxis()->SetTitleOffset(1.15);
    RecMassFitMean->Draw();
    if (RecMassFitMean->Integral() > 1)  RecMassFitMean->Fit("gaus","LQ");
    c2->Update();
    eps1->Close();

    TPostScript *eps2 = new TPostScript("s.eps", 113);
    eps2->NewPage();
    c2->Clear();
    c2->cd(0);
    if (kDYEvents) gStyle->SetOptLogy(1);
    pad = new TPad("","", .05, .05, .95, .93);
    pad->Draw();
    pad->Divide(2,3);
    pad->cd(1);  RecMassFitSc1->SetTitle("");  RecMassFitSc1->Draw();
    if (!kDYEvents) RecMassFitSc1->Fit("gaus","LQ");
    pad->cd(2);  RecMassFitScl->SetTitle("");  RecMassFitScl->Draw();
    if (!kDYEvents) RecMassFitScl->Fit("gaus","LQ");
    pad->cd(3);  RecMassFitSc12->SetTitle("");  RecMassFitSc12->Draw();
    if (!kDYEvents) RecMassFitSc12->Fit("gaus","LQ");
    pad->cd(4);  RecMassFitSc2->SetTitle("");  RecMassFitSc2->Draw();
    if (!kDYEvents) RecMassFitSc2->Fit("gaus","LQ");
    pad->cd(5);  RecMassFitSl2->SetTitle("");  RecMassFitSl2->Draw();
    if (!kDYEvents) RecMassFitSl2->Fit("gaus","LQ");
    if (kDYEvents) gStyle->SetOptLogy(0);
    c2->Update();
    eps2->Close();
  }
#endif

  // Save a few histos to a file
  // JMTBAD can this go away? we're already saving these
  TFile *fout = new TFile("sign.root", "RECREATE");

  RecMassFitNevtTot->Write();
  RecMassFitNsigTot->Write();
  RecMassFitNbkgTot->Write();
  RecMassFitNfit->Write();
  RecMassFitNsig->Write();
  RecMassFitNbkg->Write();

  RecMassFitSc1->Write();
  RecMassFitScl->Write();
  RecMassFitSc12->Write();
  RecMassFitSc2->Write();
  RecMassFitSl2->Write();
  fout->Close();
  delete fout;
}

DEFINE_FWK_MODULE(Zprime2muMassReach);
