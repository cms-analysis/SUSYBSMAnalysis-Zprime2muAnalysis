//
// Authors: Jason Mumford, Jordan Tucker, Slava Valuev, UCLA
//

#include <fstream>
#include <string>
#include <vector>

#include "TBranch.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFitManager.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
// JMTBAD perhaps a more descriptive name for this file
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Functions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TFMultiD.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/UnbinnedFitter.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAsymmetry.h"

using namespace std;

// fixed parameters concerning how to do the 2D/6D fit
const double       Zprime2muAsymmetry::nSigma = 5.;
const unsigned int Zprime2muAsymmetry::nSmearPoints = 100000;
const unsigned int Zprime2muAsymmetry::nNormPoints = 100000;

Zprime2muAsymmetry::Zprime2muAsymmetry(const edm::ParameterSet& config) 
  : Zprime2muAnalysis(config) {
  // initialize array counters since we don't use automatically-sized
  // arrays yet
  for (int i = 0; i < 3; i++)
    nfit_used[i] = 0;

  // verbosity controls the amount of debugging information printed;
  // levels are defined in an enum
  verbosity = VERBOSITY(config.getUntrackedParameter<int>("verbosity", 0));

  // turn off the fit and only make the histograms -- useful to get the
  // recSigma information from the asymHistos.ps
  noFit = config.getParameter<bool>("noFit");

  // if requested, only evaluate the log-likelihood ratio for each
  // data point
  onlyEvalLLR = config.getParameter<bool>("onlyEvalLLR");

  // choose the type of fit (magic number in the cfg, but
  // is an enum FITTYPE here in the code)
  FITTYPE fitType = FITTYPE(config.getParameter<int>("fitType"));

  // convenient flag for checking to see if we're doing graviton studies
  doingGravFit = fitType >= GRAV;

  // only do so many of the fits (up to 6; from gen+gen to rec+rec)
  numFits = config.getParameter<int>("numFits");

  // limit on how many events to read from the parameterization sample
  maxParamEvents = config.getParameter<int>("maxParamEvents");

  // whether to use the cached parameters (depending on the existence
  // of the cache root file)
  useCachedParams = config.getParameter<bool>("useCachedParams");

  // the name of the root file containing the cached parameters
  paramCacheFile = config.getParameter<string>("paramCacheFile");

  // whether only to calculate the parameterization and then quit
  calcParamsOnly = config.getParameter<bool>("calcParamsOnly");

  // whether to use the on-peak fit window or the off-peak one
  onPeak = config.getParameter<bool>("onPeak");

  // whether bremsstrahlung was switched on in the generator for the
  // parameterization sample
  internalBremOn = config.getParameter<bool>("internalBremOn");

  // whether to fix b for the simple 1-D fits
  fixbIn1DFit = config.getParameter<bool>("fixbIn1DFit");

  // whether to use cos_true in the 2/6-D fits
  useCosTrueInFit = config.getParameter<bool>("useCosTrueInFit");

  // whether to correct the cos_cs values using MC truth
  artificialCosCS = config.getParameter<bool>("artificialCosCS");

  // whether to use the probabilistic mistag correction; always turn
  // it off if doing graviton fit, if we are artificially
  // correcting with MC truth, or if we are fitting to cos_theta_true
  correctMistags = doingGravFit || artificialCosCS || useCosTrueInFit ? false :
    config.getParameter<bool>("correctMistags");

  // whether to use the calculation of the mistag probability omega(y,M)
  // instead of the parameterization obtained from the MC
  calculateMistag = config.getParameter<bool>("calculateMistag");

  // get the parameters specific to the data sample on which we are running
  string dataSet = config.getParameter<string>("dataSet");
  edm::ParameterSet dataSetConfig =
    config.getParameter<edm::ParameterSet>(dataSet);

  outputFileBase = 
    dataSetConfig.getUntrackedParameter<string>("outputFileBase");
  genSampleFiles = dataSetConfig.getParameter<vector<string> >("genSampleFiles");
  peakMass = dataSetConfig.getParameter<double>("peakMass");

  // let the asymFitManager take care of setting mass_type, etc.
  asymFitManager.setConstants(dataSetConfig, onPeak, peakMass);

  // turn on/off the fit prints in AsymFunctions
  asymDebug = verbosity >= VERBOSITY_SIMPLE;

  ostringstream out;
  out << "------------------------------------------------------------\n"
      << "Zprime2muAsymmetry parameter summary:\n"
      << "Using data set " << dataSet << " and writing to files *." 
      << outputFileBase << ".*" << endl;
  if (noFit)
    out << "Not performing any fit or evaluation of the likelihood ratio\n";
  else {
    if (onlyEvalLLR)
      out << "Only evaluating the likelihood ratios for each data point\n";
    out << "Using pdf " << asymFitManager.getPDFName() << " for fits\n";
    out << "Using ";
    if (maxParamEvents < 0) out << "all";
    else out << maxParamEvents;
    out << " events from " << genSampleFiles.size() << " files beginning with "
	<< genSampleFiles[0] << " for parameterizations\n";
  }
  out << asymFitManager << endl;
  out << "Peak mass: " << peakMass << "; accordingly, looking ";
  if (onPeak) out << "on ";
  else out << "off ";
  out << "peak\n";
  if (!correctMistags) out << "NOT ";
  out << "Correcting for mistags in data using probabilistic method\n";
  if (calculateMistag && correctMistags)
    out << "Calculating mistag probability on event-by-event basis\n";
  else
    out << "Obtaining mistag probability from parameterization\n";
  if (!internalBremOn) out << "NOT ";
  out << "Correcting for bremsstrahlung in parameterization sample\n";
  if (fixbIn1DFit) out << "Fixing b = 1.0 in 1-D fits\n";
  if (useCosTrueInFit) out << "Using cos_theta_true in 2/6-D fits\n";
  if (artificialCosCS)
    out << "Artificially correcting mistags using MC truth\n";
  
  out << "------------------------------------------------------------";
  edm::LogInfo("Zprime2muAsymmetry") << out.str();

  bookFrameHistos();
  bookFitHistos();

  bookParamHistos();
}

Zprime2muAsymmetry::~Zprime2muAsymmetry() {
  deleteHistos();
}

void Zprime2muAsymmetry::analyze(const edm::Event& event, 
				 const edm::EventSetup& eSetup) {
  // delegate filling our muon vectors to the parent class
  Zprime2muAnalysis::analyze(event, eSetup);
  eventNum = event.id().event();

  // JMTBAD break out generator-level stuffs from fillFitData
  fillFitData(event);

  fillFrameHistos();
}

void Zprime2muAsymmetry::beginJob(const edm::EventSetup& eSetup) {
  getAsymParams();
  if (calcParamsOnly)
    // JMTBAD this is the wrong way to do this!
    throw cms::Exception("Zprime2muAsymmetry") 
      << "user requested to calculate parameters only; stopping\n";
}

void Zprime2muAsymmetry::endJob() {
  calcFrameAsym();

  drawFrameHistos();
  drawFitHistos();

  dumpFitData();
  
  if (noFit)
    return;

  if (onlyEvalLLR)
    evalLikelihoods();
  else
    fitAsymmetry();
}

void Zprime2muAsymmetry::bookFitHistos() {
  const double pi = TMath::Pi(), twoPi = 2*pi;
  const AsymFitManager& amg = asymFitManager;

  h_cos_theta_true     = fs->make<TH1F>("h_cos_theta_true",     "cos #theta^{*}_{true}",                   50, -1., 1.);
  h_cos_theta_cs       = fs->make<TH1F>("h_cos_theta_cs",       "cos #theta^{*}_{CS}",                     50, -1., 1.);
  h_cos_theta_cs_acc   = fs->make<TH1F>("h_cos_theta_cs_acc",   "cos #theta^{*}_{CS} (in acceptance)",     50, -1., 1.);
  h_cos_theta_cs_fixed = fs->make<TH1F>("h_cos_theta_cs_fixed", "cos #theta^{*}_{CS} (mistags corrected)", 50, -1., 1.);
  h_cos_theta_cs_rec   = fs->make<TH1F>("h_cos_theta_cs_rec",   "cos #theta^{*}_{CS} (rec)",               50, -1., 1.);

  h_b_smass = fs->make<TH1F>("h_b_smass", "Smeared Mass of Backward Events", 50, amg.fit_win(0), amg.fit_win(1));
  h_f_smass = fs->make<TH1F>("h_f_smass", "Smeared Mass of Forward Events",  50, amg.fit_win(0), amg.fit_win(1));
  h_b_mass  = fs->make<TH1F>("h_b_mass",  "Mass of Backward Events",         50, amg.fit_win(0), amg.fit_win(1));
  h_f_mass  = fs->make<TH1F>("h_f_mass",  "Mass of Forward Events",          50, amg.fit_win(0), amg.fit_win(1));

  h_gen_sig[0] = fs->make<TH1F>("h_gen_sig0", "Mass of Forward Events in Signal",  50, amg.fit_win(0), amg.fit_win(1));
  h_gen_sig[1] = fs->make<TH1F>("h_gen_sig1", "Mass of Backward Events in Signal", 50, amg.fit_win(0), amg.fit_win(1));

  h2_rap_cos_d[0]       = fs->make<TH2F>("h2_rap_cos_d0",       "rap vs cos_cs dil (gen)",                25, -1., 1., 25, -3.5, 3.5);
  h2_rap_cos_d[1]       = fs->make<TH2F>("h2_rap_cos_d1",       "rap vs cos_true dil (gen)",              25, -1., 1., 25, -3.5, 3.5);
  h2_rap_cos_d_uncut[0] = fs->make<TH2F>("h2_rap_cos_d_uncut0", "rap vs cos_cs dil, no acc. cut (gen)",   25, -1., 1., 25, -3.5, 3.5);
  h2_rap_cos_d_uncut[1] = fs->make<TH2F>("h2_rap_cos_d_uncut1", "rap vs cos_true dil, no acc. cut (gen)", 25, -1., 1., 25, -3.5, 3.5);
  h2_rap_cos_d_rec      = fs->make<TH2F>("h2_rap_cos_d_rec",    "rap vs cos_cs dil (rec)",                25, -1., 1., 25, -3.5, 3.5);

  mistagProbEvents[0] = fs->make<TH1F>("mistagProbEvents0", "mistag prob, gen. events",              25, 0, 0.5);
  mistagProbEvents[1] = fs->make<TH1F>("mistagProbEvents1", "mistag prob, rec. events",              25, 0, 0.5);
  mistagProbEvents[2] = fs->make<TH1F>("mistagProbEvents2", "mistag prob (calculated), rec. events", 25, 0, 0.5);

  h_genCosNoCut      = fs->make<TH1F>("h_genCosNoCut", "Dil Gen cos_cs (all)", 20, -1., 1.);

  AsymFitHistoGen[0] = fs->make<TH1F>("AsymFitHistoGen0", "Dil Gen Pt",   50,  0., amg.max_pt());
  AsymFitHistoGen[1] = fs->make<TH1F>("AsymFitHistoGen1", "Dil Gen Rap",  50, -amg.max_rap(), amg.max_rap());
  AsymFitHistoGen[2] = fs->make<TH1F>("AsymFitHistoGen2", "Dil Gen Phi",  50,  0., twoPi);
  AsymFitHistoGen[3] = fs->make<TH1F>("AsymFitHistoGen3", "Dil Gen Mass", 50, amg.fit_win(0), amg.fit_win(1));
  AsymFitHistoGen[4] = fs->make<TH1F>("AsymFitHistoGen4", "Gen cos_cs",   50, -1., 1.);
  AsymFitHistoGen[5] = fs->make<TH1F>("AsymFitHistoGen5", "Gen phi_cs",   50,  0., twoPi);

  AsymFitHistoRec[0] = fs->make<TH1F>("AsymFitHistoRec0", "Dil Rec Pt",   50,  0., amg.max_pt());
  AsymFitHistoRec[1] = fs->make<TH1F>("AsymFitHistoRec1", "Dil Rec Rap",  50, -amg.max_rap(), amg.max_rap());
  AsymFitHistoRec[2] = fs->make<TH1F>("AsymFitHistoRec2", "Dil Rec Phi",  50,  0., twoPi);
  AsymFitHistoRec[3] = fs->make<TH1F>("AsymFitHistoRec3", "Dil Rec Mass", 50, amg.fit_win(0), amg.fit_win(1));
  AsymFitHistoRec[4] = fs->make<TH1F>("AsymFitHistoRec4", "Rec cos_cs",   50, -1., 1.);
  AsymFitHistoRec[5] = fs->make<TH1F>("AsymFitHistoRec5", "Rec phi_cs",   50, 0., twoPi);

  AsymFitSmearHisto[0] = fs->make<TH2F>("AsymFitSmearHisto0", "Dil Rec Pt vs Gen Pt",     20, 0., amg.max_pt(), 20, 0., amg.max_pt());
  AsymFitSmearHisto[1] = fs->make<TH2F>("AsymFitSmearHisto1", "Dil Rec Rap vs Gen Rap",   20, -amg.max_rap(), amg.max_rap(), 20, -amg.max_rap(), amg.max_rap());
  AsymFitSmearHisto[2] = fs->make<TH2F>("AsymFitSmearHisto2", "Dil Rec Phi vs Gen Phi",   20, -pi, pi, 20, -pi, pi);
  AsymFitSmearHisto[3] = fs->make<TH2F>("AsymFitSmearHisto3", "Dil Rec Mass vs Gen Mass", 20, amg.fit_win(0), amg.fit_win(1), 20, amg.fit_win(0), amg.fit_win(1));
  AsymFitSmearHisto[4] = fs->make<TH2F>("AsymFitSmearHisto4", "Rec cos_CS vs Gen cos_CS", 20, -1., 1., 20, -1., 1.);
  AsymFitSmearHisto[5] = fs->make<TH2F>("AsymFitSmearHisto5", "Rec phi_CS vs Gen phi_CS", 20, 0., twoPi, 20, 0., twoPi);

  AsymFitSmearHistoDif[0] = fs->make<TH1F>("AsymFitSmearHistoDif1", "Dil Rec Pt - Gen Pt",     50, -200., 200.);
  AsymFitSmearHistoDif[1] = fs->make<TH1F>("AsymFitSmearHistoDif2", "Dil Rec Rap - Gen Rap",   50, -.1, .1);
  AsymFitSmearHistoDif[2] = fs->make<TH1F>("AsymFitSmearHistoDif3", "Dil Rec Phi - Gen Phi",   50, -1., 1.);
  AsymFitSmearHistoDif[3] = fs->make<TH1F>("AsymFitSmearHistoDif4", "Dil Rec Mass - Gen Mass", 50, -300.,300.);
  AsymFitSmearHistoDif[4] = fs->make<TH1F>("AsymFitSmearHistoDif5", "Rec cos_cs - Gen cos_cs", 50, -1.e-3, 1.e-3);
  AsymFitSmearHistoDif[5] = fs->make<TH1F>("AsymFitSmearHistoDif6", "Rec phi_cs - Gen phi_cs", 50, -1., 1.);

  // keep separate histos by type (gg/qqbar) for gravitons
  for (int type = 0; type < 2; type++) {
    AsymFitHistoGenByType[type][0] = fs->make<TH1F>(nameHist("AsymFitHistoGenByType", type, 0).c_str(), "Dil Gen Pt",   50, 0., amg.max_pt());
    AsymFitHistoGenByType[type][1] = fs->make<TH1F>(nameHist("AsymFitHistoGenByType", type, 1).c_str(), "Dil Gen Rap",  50, -amg.max_rap(), amg.max_rap());
    AsymFitHistoGenByType[type][2] = fs->make<TH1F>(nameHist("AsymFitHistoGenByType", type, 2).c_str(), "Dil Gen Phi",  50, 0., twoPi);
    AsymFitHistoGenByType[type][3] = fs->make<TH1F>(nameHist("AsymFitHistoGenByType", type, 3).c_str(), "Dil Gen Mass", 50, amg.fit_win(0), amg.fit_win(1));
    AsymFitHistoGenByType[type][4] = fs->make<TH1F>(nameHist("AsymFitHistoGenByType", type, 4).c_str(), "Gen cos_cs",   50,-1., 1.);
    AsymFitHistoGenByType[type][5] = fs->make<TH1F>(nameHist("AsymFitHistoGenByType", type, 5).c_str(), "Gen phi_cs",   50, 0., twoPi);

    AsymFitHistoRecByType[type][0] = fs->make<TH1F>(nameHist("AsymFitHistoRecByType", type, 0).c_str(), "Dil Rec Pt",   50, 0., amg.max_pt());
    AsymFitHistoRecByType[type][1] = fs->make<TH1F>(nameHist("AsymFitHistoRecByType", type, 1).c_str(), "Dil Rec Rap",  50, -amg.max_rap(), amg.max_rap());
    AsymFitHistoRecByType[type][2] = fs->make<TH1F>(nameHist("AsymFitHistoRecByType", type, 2).c_str(), "Dil Rec Phi",  50, 0., twoPi);
    AsymFitHistoRecByType[type][3] = fs->make<TH1F>(nameHist("AsymFitHistoRecByType", type, 3).c_str(), "Dil Rec Mass", 50, amg.fit_win(0), amg.fit_win(1));
    AsymFitHistoRecByType[type][4] = fs->make<TH1F>(nameHist("AsymFitHistoRecByType", type, 4).c_str(), "Rec cos_cs",   50,-1., 1.);
    AsymFitHistoRecByType[type][5] = fs->make<TH1F>(nameHist("AsymFitHistoRecByType", type, 5).c_str(), "Rec phi_cs",   50, 0., twoPi);
  }

  AsymFitHistoGenSmeared[0] = fs->make<TH1F>("AsymFitHistoGenSmeared0", "Dil Smeared Gen Pt",   50,  0., amg.max_pt());
  AsymFitHistoGenSmeared[1] = fs->make<TH1F>("AsymFitHistoGenSmeared1", "Dil Smeared Gen Rap",  50, -amg.max_rap(), amg.max_rap());
  AsymFitHistoGenSmeared[2] = fs->make<TH1F>("AsymFitHistoGenSmeared2", "Dil Smeared Gen Phi",  50,  0., twoPi);
  AsymFitHistoGenSmeared[3] = fs->make<TH1F>("AsymFitHistoGenSmeared3", "Dil Smeared Gen Mass", 50, amg.fit_win(0), amg.fit_win(1));
  AsymFitHistoGenSmeared[4] = fs->make<TH1F>("AsymFitHistoGenSmeared4", "Gen Smeared cos_cs",   50, -1., 1.);
  AsymFitHistoGenSmeared[5] = fs->make<TH1F>("AsymFitHistoGenSmeared5", "Gen Smeared phi_cs",   50,  0., twoPi);
}

void Zprime2muAsymmetry::bookFrameHistos() {
  // Variable bins size histogram: number of bins and chosen binning.
  const int NBIN = 9;
  const double DMBINS[NBIN] =
    {200., 400., 500., 750., 1000., 1250., 1500., 2000., 3000.};
  /*
  const int NBIN = 33;
  double DMBINS[NBIN];
  for (int i = 0; i < NBIN; i++)
    DMBINS[i] = 40 + 5. * i;
  */
  //  for (int j = 0; j < NUM_REC_LEVELS; j++) {
  for (int j = 0; j < MAX_LEVELS; j++) {
    cosGJ[j][0] = fs->make<TH1F>(nameHist("cosGJ", j, 0).c_str(), "Cos theta in Gottfried-Jackson Frame",        100, -1., 1.);
    cosGJ[j][1] = fs->make<TH1F>(nameHist("cosGJ", j, 1).c_str(), "Cos theta in tagged Gottfried-Jackson Frame", 100, -1., 1.);
    cosCS[j][0] = fs->make<TH1F>(nameHist("cosCS", j, 0).c_str(), "Cos theta in Collins-Soper Frame",            100, -1., 1.);
    cosCS[j][1] = fs->make<TH1F>(nameHist("cosCS", j, 1).c_str(), "Cos theta in analytic Collins-Soper Frame",	 100, -1., 1.);
    cosBoost[j] = fs->make<TH1F>(nameHist("cosBoost", j).c_str(), "Cos theta in Boost Frame",                    100, -1., 1.);
    cosW[j]     = fs->make<TH1F>(nameHist("cosW", j).c_str(),     "Cos theta in Wulz Frame",                     100, -1., 1.);

    if (j == 0 || j == 3)
      rap_vs_cosCS[j] = fs->make<TH2F>(nameHist("rap_vs_cosCS", j).c_str(),
				       "Rap vs Cos theta CS",
				       50, -1., 1.,  50, -4., 4.);

    FMassGJ[j][0] = fs->make<TH1F>(nameHist("FMassGJ", j, 0).c_str(), "F (M), GJ frame",        NBIN-1, DMBINS);
    FMassGJ[j][1] = fs->make<TH1F>(nameHist("FMassGJ", j, 0).c_str(), "F (M), tagged GJ frame", NBIN-1, DMBINS);
    FMassCS[j][0] = fs->make<TH1F>(nameHist("FMassCS", j, 0).c_str(), "F (M), CS frame",        NBIN-1, DMBINS);
    FMassCS[j][1] = fs->make<TH1F>(nameHist("FMassCS", j, 0).c_str(), "F (M), CS anal frame",   NBIN-1, DMBINS);
    FMassBoost[j] = fs->make<TH1F>(nameHist("FMassBoost", j).c_str(), "F (M), Boost frame",     NBIN-1, DMBINS);
    FMassW[j]     = fs->make<TH1F>(nameHist("FMassW", j).c_str(),     "F (M), Wulz frame",      NBIN-1, DMBINS);

    BMassGJ[j][0] = fs->make<TH1F>(nameHist("BMassGJ", j, 0).c_str(), "B (M), GJ frame",        NBIN-1, DMBINS);
    BMassGJ[j][1] = fs->make<TH1F>(nameHist("BMassGJ", j, 1).c_str(), "B (M), tagged GJ frame", NBIN-1, DMBINS);
    BMassCS[j][0] = fs->make<TH1F>(nameHist("BMassCS", j, 0).c_str(), "B (M), CS frame",        NBIN-1, DMBINS);
    BMassCS[j][1] = fs->make<TH1F>(nameHist("BMassCS", j, 1).c_str(), "B (M), CS anal frame",   NBIN-1, DMBINS);
    BMassBoost[j] = fs->make<TH1F>(nameHist("BMassBoost", j).c_str(), "B (M), Boost frame",     NBIN-1, DMBINS);
    BMassW[j]     = fs->make<TH1F>(nameHist("BMassW", j).c_str(),     "B (M), Wulz frame",      NBIN-1, DMBINS);

    AMassGJ[j][0] = fs->make<TH1F>(nameHist("AMassGJ", j, 0).c_str(), "A (M), GJ frame",        NBIN-1, DMBINS);
    AMassGJ[j][1] = fs->make<TH1F>(nameHist("AMassGJ", j, 1).c_str(), "A (M), tagged GJ frame", NBIN-1, DMBINS);
    AMassCS[j][0] = fs->make<TH1F>(nameHist("AMassCS", j, 0).c_str(), "A (M), CS frame",        NBIN-1, DMBINS);
    AMassCS[j][1] = fs->make<TH1F>(nameHist("AMassCS", j, 1).c_str(), "A (M), analyt CS frame", NBIN-1, DMBINS);
    AMassBoost[j] = fs->make<TH1F>(nameHist("AMassBoost", j).c_str(), "A (M), Boost frame",     NBIN-1, DMBINS);
    AMassW[j]     = fs->make<TH1F>(nameHist("AMassW", j).c_str(),     "A (M), Wulz frame",      NBIN-1, DMBINS);

    FRapGJ[j][0] = fs->make<TH1F>(nameHist("FRapGJ", j, 0).c_str(), "F (Rap), GJ frame",        20, -2.5, 2.5);
    FRapGJ[j][1] = fs->make<TH1F>(nameHist("FRapGJ", j, 1).c_str(), "F (Rap), tagged GJ frame", 20, -2.5, 2.5);
    FRapCS[j][0] = fs->make<TH1F>(nameHist("FRapCS", j, 0).c_str(), "F (Rap), CS frame",        20, -2.5, 2.5);
    FRapCS[j][1] = fs->make<TH1F>(nameHist("FRapCS", j, 1).c_str(), "F (Rap), CS anal frame",   20, -2.5, 2.5);
    FRapBoost[j] = fs->make<TH1F>(nameHist("FRapBoost", j).c_str(), "F (Rap), Boost frame",     20, -2.5, 2.5);
    FRapW[j]     = fs->make<TH1F>(nameHist("FRapW", j).c_str(),     "F (Rap), Wulz frame",      20, -2.5, 2.5);

    BRapGJ[j][0] = fs->make<TH1F>(nameHist("BRapGJ", j, 0).c_str(), "B (Rap), GJ frame",        20, -2.5, 2.5);
    BRapGJ[j][1] = fs->make<TH1F>(nameHist("BRapGJ", j, 1).c_str(), "B (Rap), tagged GJ frame", 20, -2.5, 2.5);
    BRapCS[j][0] = fs->make<TH1F>(nameHist("BRapCS", j, 0).c_str(), "B (Rap), CS frame",        20, -2.5, 2.5);
    BRapCS[j][1] = fs->make<TH1F>(nameHist("BRapCS", j, 1).c_str(), "B (Rap), CS anal frame",   20, -2.5, 2.5);
    BRapBoost[j] = fs->make<TH1F>(nameHist("BRapBoost", j).c_str(), "B (Rap), Boost frame",     20, -2.5, 2.5);
    BRapW[j]     = fs->make<TH1F>(nameHist("BRapW", j).c_str(),     "B (Rap), Wulz frame",      20, -2.5, 2.5);

    ARapGJ[j][0] = fs->make<TH1F>(nameHist("ARapGJ", j, 0).c_str(), "A (y), GJ frame",        20, -2.5, 2.5);
    ARapGJ[j][1] = fs->make<TH1F>(nameHist("ARapGJ", j, 1).c_str(), "A (y), tagged GJ frame", 20, -2.5, 2.5);
    ARapCS[j][0] = fs->make<TH1F>(nameHist("ARapCS", j, 0).c_str(), "A (y), CS frame",        20, -2.5, 2.5);
    ARapCS[j][1] = fs->make<TH1F>(nameHist("ARapCS", j, 1).c_str(), "A (y), analyt CS frame", 20, -2.5, 2.5);
    ARapBoost[j] = fs->make<TH1F>(nameHist("ARapBoost", j).c_str(), "A (y), Boost frame",     20, -2.5, 2.5);
    ARapW[j]     = fs->make<TH1F>(nameHist("ARapW", j).c_str(),     "A (y), Wulz frame",      20, -2.5, 2.5);

    FPseudGJ[j]    = fs->make<TH1F>(nameHist("FPseudGJ", j, 0).c_str(), "F (Pseud), tagged GJ frame", 50, -6., 6.);
    FPseudCS[j]    = fs->make<TH1F>(nameHist("FPseudCS", j, 1).c_str(), "F (Pseud), CS frame",        50, -6., 6.);
    FPseudBoost[j] = fs->make<TH1F>(nameHist("FPseudBoost", j).c_str(), "F (Pseud), Boost frame",     50, -6., 6.);
    FPseudW[j]     = fs->make<TH1F>(nameHist("FPseudW", j).c_str(),     "F (Pseud), Wulz frame",      50, -6., 6.);

    BPseudGJ[j]    = fs->make<TH1F>(nameHist("BPseudGJ", j, 0).c_str(), "B (Pseud), tagged GJ frame", 50, -6., 6.);
    BPseudCS[j]    = fs->make<TH1F>(nameHist("BPseudCS", j, 1).c_str(), "B (Pseud), CS frame",        50, -6., 6.);
    BPseudBoost[j] = fs->make<TH1F>(nameHist("BPseudBoost", j).c_str(), "B (Pseud), Boost frame",     50, -6., 6.);
    BPseudW[j]     = fs->make<TH1F>(nameHist("BPseudW", j).c_str(),     "B (Pseud), Wulz frame",      50, -6., 6.);

    FMBoostCut[j][0]= fs->make<TH1F>(nameHist("FMBoostCut", j, 0).c_str(), "F(mass), boost, fabs(rap)<0.4", NBIN-1, DMBINS);
    FMBoostCut[j][1]= fs->make<TH1F>(nameHist("FMBoostCut", j, 1).c_str(), "F(mass), boost, 0.4<fabs(rap)<0.8", NBIN-1, DMBINS);
    FMBoostCut[j][2]= fs->make<TH1F>(nameHist("FMBoostCut", j, 2).c_str(), "F(mass), boost, 0.8<fabs(rap)<2.4", NBIN-1, DMBINS);
    FMBoostCut[j][3]= fs->make<TH1F>(nameHist("FMBoostCut", j, 3).c_str(), "F(mass), boost, fabs(pseudorap)<0.4", NBIN-1, DMBINS);
    FMBoostCut[j][4]= fs->make<TH1F>(nameHist("FMBoostCut", j, 4).c_str(), "F(mass), boost, 0.4<fabs(pseudorap)<0.8", NBIN-1, DMBINS);
    FMBoostCut[j][5]= fs->make<TH1F>(nameHist("FMBoostCut", j, 5).c_str(), "F(mass), boost, 0.8<fabs(pseudorap)<2.4", NBIN-1, DMBINS);

    BMBoostCut[j][0]= fs->make<TH1F>(nameHist("BMBoostCut", j, 0).c_str(), "B(mass), boost, fabs(rap)<0.4", NBIN-1, DMBINS);
    BMBoostCut[j][1]= fs->make<TH1F>(nameHist("BMBoostCut", j, 1).c_str(), "B(mass), boost, 0.4<fabs(rap)<0.8", NBIN-1, DMBINS);
    BMBoostCut[j][2]= fs->make<TH1F>(nameHist("BMBoostCut", j, 2).c_str(), "B(mass), boost, 0.8<fabs(rap)<2.4", NBIN-1, DMBINS);
    BMBoostCut[j][3]= fs->make<TH1F>(nameHist("BMBoostCut", j, 3).c_str(), "B(mass), boost, fabs(pseudorap)<0.4", NBIN-1, DMBINS);
    BMBoostCut[j][4]= fs->make<TH1F>(nameHist("BMBoostCut", j, 4).c_str(), "B(mass), boost, 0.4<fabs(pseudorap)<0.8", NBIN-1, DMBINS);
    BMBoostCut[j][5]= fs->make<TH1F>(nameHist("BMBoostCut", j, 5).c_str(), "B(mass), boost, 0.8<fabs(pseudorap)<2.4", NBIN-1, DMBINS);

    AsymMBoostCut[j][0]= fs->make<TH1F>(nameHist("AsymMBoostCut", j, 0).c_str(), "A(mass), boost, abs(rap)<0.4", NBIN-1, DMBINS);
    AsymMBoostCut[j][1]= fs->make<TH1F>(nameHist("AsymMBoostCut", j, 1).c_str(), "A(mass), boost, 0.4<abs(rap)<0.8", NBIN-1, DMBINS);
    AsymMBoostCut[j][2]= fs->make<TH1F>(nameHist("AsymMBoostCut", j, 2).c_str(), "A(mass), boost, 0.8<abs(rap)<2.4", NBIN-1, DMBINS);
    AsymMBoostCut[j][3]= fs->make<TH1F>(nameHist("AsymMBoostCut", j, 3).c_str(), "A(mass), boost, abs(pseudorap)<0.4", NBIN-1, DMBINS);
    AsymMBoostCut[j][4]= fs->make<TH1F>(nameHist("AsymMBoostCut", j, 4).c_str(), "A(mass), boost, 0.4<abs(pseudorap)<0.8", NBIN-1, DMBINS);
    AsymMBoostCut[j][5]= fs->make<TH1F>(nameHist("AsymMBoostCut", j, 5).c_str(), "A(mass), boost, 0.8<abs(pseudorap)<2.4", NBIN-1, DMBINS);
  }

  // Resolution histograms
  cosCSRes[0] = fs->make<TH1F>("cosCSRes0", "L1 cos CS - Gen cos CS", 100, -0.5,   0.5);
  cosCSRes[1] = fs->make<TH1F>("cosCSRes1", "L2 cos CS - Gen cos CS", 100, -0.25,  0.25);
  cosCSRes[2] = fs->make<TH1F>("cosCSRes2", "L3 cos CS - Gen cos CS", 100, -0.005, 0.005);
  for (int j = 3; j < MAX_LEVELS; j++) {
    string title = recLevelHelper.levelName(j) + " cos CS - Gen cos CS";
    cosCSRes[j] = fs->make<TH1F>(nameHist("cosCSRes", j).c_str(), title.c_str(), 100, -0.005, 0.005);
  }
  cosCS3_diffsq_vs_cosCS0 = fs->make<TProfile>("cosCS3_diffsq_vs_cosCS0",
					       "L3 cos theta CS diffsq vs Gen cos theta CS",
					       20, -1., 1., 0., 0.0625);
  rap3_vs_rap0 = fs->make<TH2F>("rap3_vs_rap0","L3 Rap vs Gen Rap", 50, -4., 4., 50, -4., 4.);
}

void Zprime2muAsymmetry::makeFakeData(int nEvents, double A_FB, double b) {
  // To see whether fitter is working properly, sometimes it is helpful
  // to use artificially created data taken from distributions.

  gRandom->SetSeed(12191982);

  // Chosen sets of mistag and rapidity parameters
  //mistag_pars[0] = -.462;   mistag_pars[1] = .108; mistag_pars[2] = 2.15;
  //mistag_pars[3] = .000762; mistag_pars[4] = .281; mistag_pars[5] = 14.1;
  //rap_pars[0] = 4.94e+3; rap_pars[1] = -3.82; rap_pars[2] = 6.1;
  //rap_pars[3] = -1.11;   rap_pars[4] =  6.3e4;

  // asymetry and rapiditiy pdfs
  TF1* f_cos = new TF1("f_cos", asym_3_PDF, -1., 1., 3);
  TF1* f_rap = new TF1("f_rap", yDistRTC, -3.5, 3.5);
  TF1* f_mass = new TF1("f_mass", massDist,
			asymFitManager.fit_win(0), asymFitManager.fit_win(1));

  // Set desired asymmetry values
  f_cos->SetParameters(1., A_FB, b);

  for (int i = 0; i < nEvents; i++) {
    // Get random rapidity and cos_true values, and store in arrays
    double rap = f_rap->GetRandom();
    fake_rap[i] = rap;
    double cos_true = f_cos->GetRandom();
    fake_cos_true[i] = cos_true;
    double mass = f_mass->GetRandom();

    // Get probability for mistag from rapidity and cos_true
    double w_rap = mistagProbVsRap(rap, mass);
    double w = mistagProb(rap, cos_true, mass);

    // Generate random number, if random number is less than total mistag
    // probability take cos_cs to have opposite sign from cos_true.
    double rand = gRandom->Rndm();
    double cos_cs = cos_true;
    if (rand < w) {
      cos_cs *= -1.;
      fake_mistag_cs[i] = 1;
    }
    else
      fake_mistag_cs[i] = 0;
    fake_cos_cs[i] = cos_cs;

    // Also if random number is less than mistag from rapidity, store
    // set mistag_true value to 1
    fake_mistag_true[i] = rand < w_rap ? 1 : 0;
  }
}

void Zprime2muAsymmetry::smearGenData(const AsymFitData& data) {
  // Smear generated data with gaussians and store to be fit

  // alias for typing
  const AsymFitManager& afm = asymFitManager;

  // First smear the mass... If the data still lies within the mass window
  // then proceed to smear and store cos and rapidity
  double smear_mass = gRandom->Gaus(data.mass, afm.recSigmaMass());

  if (smear_mass >= afm.fit_win(0) && smear_mass <= afm.fit_win(1)) {
    double smear_cos = gRandom->Gaus(data.cos_cs, afm.recSigmaCosCS());
    // JMTBAD cos_true smear uses resolution for cos_cs
    double smear_cos_true = gRandom->Gaus(data.cos_true, afm.recSigmaCosCS());
    double smear_rap = gRandom->Gaus(data.rapidity, afm.recSigmaRap());
    double smear_pt = gRandom->Gaus(data.pT, afm.recSigmaPt());
    double smear_phi = gRandom->Gaus(data.phi, afm.recSigmaPhi());
    double smear_phi_cs = gRandom->Gaus(data.phi_cs, afm.recSigmaPhiCS());

    // Clamp smeared values appropriately
    if (smear_cos > 1.)
      smear_cos =  2. - smear_cos;
    else if (smear_cos < -1.)
      smear_cos = -2. - smear_cos;

    if (smear_pt < 0)
      smear_pt *= -1;

    const double twopi = 2*TMath::Pi();
    if (smear_phi < 0)
      smear_phi += twopi;
    else if (smear_phi > twopi)
      smear_phi -= twopi; 

    if (smear_phi_cs < 0)
      smear_phi_cs += twopi;
    else if (smear_phi_cs > twopi)
      smear_phi_cs -= twopi;

    // Store data in arrays and histograms
    cos_theta_cs_data[2][nfit_used[2]] = smear_cos;
    rap_dil_data[2][nfit_used[2]] = smear_rap;
    pt_dil_data[2][nfit_used[2]] = smear_pt;
    phi_dil_data[2][nfit_used[2]] = smear_phi;
    mass_dil_data[2][nfit_used[2]] = smear_mass;
    phi_cs_data[2][nfit_used[2]] = smear_phi_cs;
    nfit_used[2]++;

    AsymFitHistoGenSmeared[0]->Fill(smear_pt);
    AsymFitHistoGenSmeared[1]->Fill(smear_rap);
    AsymFitHistoGenSmeared[2]->Fill(smear_phi);
    AsymFitHistoGenSmeared[3]->Fill(smear_mass);
    AsymFitHistoGenSmeared[4]->Fill(smear_cos);
    AsymFitHistoGenSmeared[5]->Fill(smear_phi_cs);

    // Artificially smear of mass and cos_true to compare with number
    // extracted from fit to reconstructed data (this is for a
    // counting-type measurement for forward-backward asymmetry;
    // different than the 2D fit done with smeared data)
    if (smear_cos_true < 0)
      h_b_smass->Fill(smear_mass);
    else
      h_f_smass->Fill(smear_mass);
  
    // also store un-smeared values
    if (data.cos_true < 0)
      h_b_mass->Fill(data.mass);
    else
      h_f_mass->Fill(data.mass);
  }
}

void Zprime2muAsymmetry::fillFitData(const edm::Event& event) {
  // the debug dumps are the outputs of the calcCosTheta*/etc. functions
  bool debug = verbosity >= VERBOSITY_TOOMUCH; 

  double gen_cos_cs = 0., rec_cos_cs = 0., gen_phi_cs = 0., rec_phi_cs = 0.;
  unsigned int n_dil = bestDileptons.size();
  unsigned int n_gen = allDileptons[lgen].size();

  int* type = new int[n_gen];
  for (unsigned i = 0; i < n_gen; i++) type[i] = -1;

  if (!reconstructedOnly) {
    edm::Handle<reco::CandidateCollection> genParticles;
    event.getByLabel("genParticleCandidates", genParticles);

    AsymFitData data;
    // JMTBAD redundant with some calculations below
    if (!computeFitQuantities(*genParticles, eventNum, data)) {
      // if finding the Z'/etc failed, skip this event
      edm::LogWarning("fillFitData") << "could not compute fit quantities!";
      return;
    }

    // remove mistag correction for gravitons (where there are only
    // terms even in cos_cs anyway)
    if (doingGravFit && data.pL < 0)
      angDist.push_back(-data.cos_cs);
    else
      angDist.push_back(data.cos_cs);

    // mass distributions for forward and backward events separately for
    // the entire spectrum, not just in the reconstructed window
    h_gen_sig[data.cos_true >= 0 ? 0 : 1]->Fill(data.mass);

    if (data.mass > asymFitManager.fit_win(0) &&
	data.mass < asymFitManager.fit_win(1)) {
      // 1D histograms with cos_true and cos_cs to be fit with binned fit.
      h_cos_theta_true->Fill(data.cos_true);
      h_cos_theta_cs->Fill(data.cos_cs); 
      //if (data.mistag_true == 0) h_cos_theta_cs_fixed->Fill(data.cos_cs);
      if (data.mistag_cs == 0)
	h_cos_theta_cs_fixed->Fill(data.cos_cs);
      else
	h_cos_theta_cs_fixed->Fill(-data.cos_cs);
    
      h2_rap_cos_d_uncut[0]->Fill(data.cos_cs, data.rapidity);
      h2_rap_cos_d_uncut[1]->Fill(data.cos_true, data.rapidity);
      
      // See if dimuon is within detector acceptance.
      if (data.cut_status == NOTCUT) {
	//1D histo for fit with data containing mistags and acceptance
	h_cos_theta_cs_acc->Fill(data.cos_cs); 
      
	// Store data used for doing resolution-smeared type fits
	// JMTBAD currently do not do the fits
	smearGenData(data);

	// Artificially correct for mistags using MC truth
	//if (data.mistag_cs != 0 && correctMistags) data.cos_cs *= -1;
	//if (data.mistag_true == 1 && correctMistags) data.cos_cs *= -1.;

	// These lines can be used to test mistag probability (start with
	// cos_true and introduce mistag artificially (either with mistag
	// definition or from sampling from mistag as a function of 
	// rapidity distribution.
	//double rndm = gRandom->Rndm();
	//double w = mistagProb(data.rapidity, data.cos_true);
	//double w = mistagProbVsRap(data.rapidity);
	//if (rndm < w) data.cos_cs *= -1.;

	h2_rap_cos_d[0]->Fill(data.cos_cs, data.rapidity);
	h2_rap_cos_d[1]->Fill(data.cos_true, data.rapidity);

	double w = mistagProb(data.rapidity, data.cos_cs, data.mass);
	mistagProbEvents[0]->Fill(w);
      }
    }

    // Loop over all generated dimuons
    for (unsigned int i_dil = 0; i_dil < n_gen; i_dil++) {
      const reco::Candidate& gen_dil = allDileptons[lgen][i_dil];
      const reco::CandidateBaseRef& gen_mum = 
	dileptonDaughterByCharge(gen_dil, -1);
      const reco::CandidateBaseRef& gen_mup = 
	dileptonDaughterByCharge(gen_dil, +1);

      gen_cos_cs = calcCosThetaCSAnal(gen_mum->pz(), gen_mum->energy(), 
				      gen_mup->pz(), gen_mup->energy(), 
				      gen_dil.pt(), gen_dil.pz(),
				      gen_dil.mass(), debug);
      gen_phi_cs = calcPhiCSAnal(gen_mum->px(), gen_mum->py(), gen_mup->px(),
				 gen_mup->py(), gen_dil.pt(), gen_dil.eta(), 
				 gen_dil.phi(), gen_dil.mass(), debug);

      int grannyID = abs(grandmotherId(gen_mum));
      if (grannyID > 0 && grannyID <= 6) type[i_dil] = 0;
      else if (grannyID == 21) type[i_dil] = 1;

      double phi = gen_dil.phi();
      if (phi < 0.)
	phi += 2*TMath::Pi();

      if (type[i_dil] >= 0) {
	AsymFitHistoGenByType[type[i_dil]][0]->Fill(gen_dil.pt());
	AsymFitHistoGenByType[type[i_dil]][1]->Fill(gen_dil.rapidity());
	AsymFitHistoGenByType[type[i_dil]][2]->Fill(phi);
	AsymFitHistoGenByType[type[i_dil]][3]->Fill(gen_dil.mass());
	AsymFitHistoGenByType[type[i_dil]][4]->Fill(gen_cos_cs);
	AsymFitHistoGenByType[type[i_dil]][5]->Fill(gen_phi_cs); 
      }

      // Check to see if generated dimuons lie within acceptance cut and
      // generated window
      if (gen_mum->eta() > MUM_ETA_LIM[0] && gen_mum->eta() < MUM_ETA_LIM[1] &&
	  gen_mup->eta() > MUP_ETA_LIM[0] && gen_mup->eta() < MUP_ETA_LIM[1]) {
	if (gen_dil.mass() > asymFitManager.fit_win(0) && 
	    gen_dil.mass() < asymFitManager.fit_win(1)) {
	  if (nfit_used[0] == FIT_ARRAY_SIZE - 1)
	    throw cms::Exception("Zprime2muAsymmetry")
	      << "data arrays not large enough! nfit_used[0]="
	      << nfit_used[0] << " FIT_ARRAY_SIZE=" << FIT_ARRAY_SIZE << endl;

	  // Store generated quantities in histograms and arrays used in fit.
	  pt_dil_data[0][nfit_used[0]] = gen_dil.pt();
	  phi_dil_data[0][nfit_used[0]] = phi;
	  mass_dil_data[0][nfit_used[0]] = gen_dil.mass();
	  rap_dil_data[0][nfit_used[0]] = gen_dil.rapidity();
	  if (useCosTrueInFit)
	    cos_theta_cs_data[0][nfit_used[0]] = data.cos_true;
	  else {
	    if (artificialCosCS && data.mistag_cs)
	      gen_cos_cs *= -1;
	    cos_theta_cs_data[0][nfit_used[0]] = gen_cos_cs;
	  }
	  phi_cs_data[0][nfit_used[0]] = gen_phi_cs;
	  AsymFitHistoGen[0]->Fill(pt_dil_data[0][nfit_used[0]]);
	  AsymFitHistoGen[1]->Fill(rap_dil_data[0][nfit_used[0]]);
	  AsymFitHistoGen[2]->Fill(phi_dil_data[0][nfit_used[0]]);
	  AsymFitHistoGen[3]->Fill(mass_dil_data[0][nfit_used[0]]);
	  AsymFitHistoGen[4]->Fill(cos_theta_cs_data[0][nfit_used[0]]);
	  AsymFitHistoGen[5]->Fill(phi_cs_data[0][nfit_used[0]]);
	  // keep track of which data point is from qqbar or gg
	  if (type[i_dil] == 0) {
	    qqbar_weights[0][nfit_used[0]] = 1.0;
	    gg_weights[0][nfit_used[0]] = 0.0;
	  }
	  else if (type[i_dil] == 1) {
	    qqbar_weights[0][nfit_used[0]] = 0.0;
	    gg_weights[0][nfit_used[0]] = 1.0;
	  }
	  nfit_used[0]++;
	}

	// if number of generated and reconstructed dimuons are the same,
	// fill histos with resolutions between reconstructed and generated
	// dimuon information.  This will be used for obtaining sigmas used
	// in convolutions for asymmetry fits.
	if (n_dil == n_gen) {      
	  const reco::Candidate& rec_dil = bestDileptons[i_dil];
	  const reco::CandidateBaseRef& rec_mum = 
	    dileptonDaughterByCharge(rec_dil, -1);
	  const reco::CandidateBaseRef& rec_mup = 
	    dileptonDaughterByCharge(rec_dil, +1);
	  rec_cos_cs = calcCosThetaCSAnal(rec_mum->pz(), rec_mum->energy(), 
					  rec_mup->pz(), rec_mup->energy(), 
					  rec_dil.pt(), rec_dil.pz(),
					  rec_dil.mass(), debug);
	  rec_phi_cs = calcPhiCSAnal(rec_mum->px(), rec_mum->py(),
				     rec_mup->px(), rec_mup->py(),
				     rec_dil.pt(), rec_dil.eta(),
				     rec_dil.phi(), rec_dil.mass(), debug);
	  AsymFitSmearHisto[0]->Fill(gen_dil.pt(), rec_dil.pt());
	  AsymFitSmearHistoDif[0]->Fill(rec_dil.pt() - gen_dil.pt());
	  AsymFitSmearHisto[1]->Fill(gen_dil.rapidity(), rec_dil.rapidity());
	  AsymFitSmearHistoDif[1]->Fill(rec_dil.rapidity() - gen_dil.rapidity());
	  AsymFitSmearHisto[2]->Fill(gen_dil.phi(), rec_dil.phi());
	  AsymFitSmearHistoDif[2]->Fill(rec_dil.phi() - gen_dil.phi());
	  AsymFitSmearHisto[3]->Fill(gen_dil.mass(), rec_dil.mass());
	  AsymFitSmearHistoDif[3]->Fill(rec_dil.mass() - gen_dil.mass());
	  AsymFitSmearHisto[4]->Fill(gen_cos_cs, rec_cos_cs);
	  AsymFitSmearHistoDif[4]->Fill(rec_cos_cs - gen_cos_cs);
	  AsymFitSmearHisto[5]->Fill(gen_phi_cs, rec_phi_cs);
	  AsymFitSmearHistoDif[5]->Fill(rec_phi_cs - gen_phi_cs);
	}
      }
      if (gen_dil.mass() > asymFitManager.fit_win(0) &&
	  gen_dil.mass() < asymFitManager.fit_win(1)) {
	gen_cos_cs = calcCosThetaCSAnal(gen_mum->pz(), gen_mum->energy(), 
					gen_mup->pz(), gen_mup->energy(), 
					gen_dil.pt(), gen_dil.pz(),
					gen_dil.mass(), debug);
      
	// Use this histogram for value of cos_true to be compared with 
	// with reconstructed.
	h_genCosNoCut->Fill(gen_cos_cs);
      }
    }
  }

  // Now loop over all reconstructed dimuons and store those to be fitted
  // (which lie inside desired reconstructed window).
  for (unsigned int i_dil = 0; i_dil < n_dil; i_dil++) {
    const reco::Candidate& rec_dil = bestDileptons[i_dil];
    if (rec_dil.mass() > asymFitManager.fit_win(0) &&
	rec_dil.mass() < asymFitManager.fit_win(1)) {
      if (nfit_used[0] == FIT_ARRAY_SIZE - 1)
	throw cms::Exception("Zprime2muAsymmetry")
	  << "data arrays not large enough! nfit_used[1]="
	  << nfit_used[1] << " FIT_ARRAY_SIZE=" << FIT_ARRAY_SIZE << endl;
      
      const reco::CandidateBaseRef& rec_mum = 
	dileptonDaughterByCharge(rec_dil, -1);
      const reco::CandidateBaseRef& rec_mup = 
	dileptonDaughterByCharge(rec_dil, +1);
      
      rec_cos_cs = calcCosThetaCSAnal(rec_mum->pz(), rec_mum->energy(), 
				      rec_mup->pz(), rec_mup->energy(), 
				      rec_dil.pt(), rec_dil.pz(),
				      rec_dil.mass(), debug);
      rec_phi_cs = calcPhiCSAnal(rec_mum->px(), rec_mum->py(), rec_mup->px(),
				 rec_mup->py(), rec_dil.pt(), rec_dil.eta(), 
				 rec_dil.phi(), rec_dil.mass(), debug);
      pt_dil_data[1][nfit_used[1]] = rec_dil.pt();
      phi_dil_data[1][nfit_used[1]] = rec_dil.phi();
      if (phi_dil_data[1][nfit_used[1]] < 0.) 
	phi_dil_data[1][nfit_used[1]] += 2.*TMath::Pi();
      mass_dil_data[1][nfit_used[1]] = rec_dil.mass();
      double rec_rap = rap_dil_data[1][nfit_used[1]] = rec_dil.rapidity();
      cos_theta_cs_data[1][nfit_used[1]] = rec_cos_cs;
      phi_cs_data[1][nfit_used[1]] = rec_phi_cs;
      AsymFitHistoRec[0]->Fill(pt_dil_data[1][nfit_used[1]]);
      AsymFitHistoRec[1]->Fill(rap_dil_data[1][nfit_used[1]]);
      AsymFitHistoRec[2]->Fill(phi_dil_data[1][nfit_used[1]]);
      AsymFitHistoRec[3]->Fill(mass_dil_data[1][nfit_used[1]]);
      AsymFitHistoRec[4]->Fill(cos_theta_cs_data[1][nfit_used[1]]);
      AsymFitHistoRec[5]->Fill(phi_cs_data[1][nfit_used[1]]);

      h_cos_theta_cs_rec->Fill(rec_cos_cs);
      h2_rap_cos_d_rec->Fill(rec_cos_cs, rec_rap);

      double w = mistagProb(rec_rap, rec_cos_cs, rec_dil.mass());
      mistagProbEvents[1]->Fill(w);

      // if we get exactly one generated dimuon and one reconstructed
      // dimuon, we can assume that the grandma of the reconstructed
      // one is the same as the generated one
      if (n_gen >= 1 && (n_gen != n_dil))
	edm::LogWarning("Zprime2muAsymmetry")
	  << "don't know how to get grandmaType for sure";

      if (type >= 0 && i_dil < n_gen) {
	AsymFitHistoRecByType[type[i_dil]][0]->Fill(pt_dil_data[1][nfit_used[1]]);
	AsymFitHistoRecByType[type[i_dil]][1]->Fill(rap_dil_data[1][nfit_used[1]]);
	AsymFitHistoRecByType[type[i_dil]][2]->Fill(phi_dil_data[1][nfit_used[1]]);
	AsymFitHistoRecByType[type[i_dil]][3]->Fill(mass_dil_data[1][nfit_used[1]]);
	AsymFitHistoRecByType[type[i_dil]][4]->Fill(cos_theta_cs_data[1][nfit_used[1]]);
	AsymFitHistoRecByType[type[i_dil]][5]->Fill(phi_cs_data[1][nfit_used[1]]); 
      }
      // keep track of which data point is from qqbar or gg
      if (type[i_dil] == 0) {
	qqbar_weights[1][nfit_used[1]] = 1.0;
	gg_weights[1][nfit_used[1]] = 0.0;
      }
      else if (type[i_dil] == 1) {
	qqbar_weights[1][nfit_used[1]] = 0.0;
	gg_weights[1][nfit_used[1]] = 1.0;
      }

      nfit_used[1]++;
    }
  }

  delete[] type;
}

void Zprime2muAsymmetry::fillFrameHistos() {
  // Calculate cosine theta^* in various frames.
  bool debug = verbosity >= VERBOSITY_TOOMUCH;

  double cos_gj, cos_gj_tag, cos_cs_anal;
  const unsigned int MAX_DILEPTONS = 2;
  double cos_cs[MAX_LEVELS][MAX_DILEPTONS];
  double cos_boost, cos_wulz;

  TLorentzVector pp1(0., 0., 7000., 7000.);
  TLorentzVector pp2(0., 0.,-7000., 7000.);
  TLorentzVector pb_star, pt_star, pmum_star, pmup_star;
  TVector3 temp3;

  //for (int i_rec = 0; i_rec < NUM_REC_LEVELS; i_rec++) {
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    // Trigger decision
    if (i_rec < NUM_REC_LEVELS && (cutTrig[i_rec] && !passTrigger(i_rec)))
      break;

    //Look for an opposite-sign dilepton at this level of reconstruction.
    for (unsigned i_dil = 0; i_dil < allDileptons[i_rec].size(); i_dil++) {
      const reco::Candidate& dil = allDileptons[i_rec][i_dil];
      const reco::CandidateBaseRef& mum = dileptonDaughterByCharge(dil, -1);
      const reco::CandidateBaseRef& mup = dileptonDaughterByCharge(dil, +1);
      
      TLorentzVector vmum, vmup, vdil;
      vmum.SetPtEtaPhiM(mum->pt(), mum->eta(), mum->phi(), leptonMass);
      vmup.SetPtEtaPhiM(mup->pt(), mup->eta(), mup->phi(), leptonMass);
      vdil.SetPtEtaPhiM(dil.pt(),  dil.eta(),  dil.phi(),  dil.mass());

      if (debug) {
	ostringstream out;
	TLorentzVector tempV;
	out << recLevelHelper.levelName(i_rec) << " values:" << endl;
	out << "#|Charge |   Eta   |   Phi   |    P    |"
	    << "    Pt   |    Pz   |   Rap   |  Mass  " << endl;
	out << "-------------------------------------------------"
	    << "------------------------------" << endl;
	for (int i_part = 0; i_part < 3; i_part++) {
	  if (i_part == 0)      tempV = vmum;
	  else if (i_part == 1) tempV = vmup;
	  else if (i_part == 2) tempV = vdil;
	  out << i_part << "|  ";
	  if (i_part == 0)  out << setw(4) << mum->charge() << "  | ";
	  else if (i_part == 1) out << setw(4) << mup->charge() << "  | ";
	  else out << "     | ";
	  out << endl << tempV;
	  LogTrace("fillFrameHistos") << out.str();
	}
      }

      // Do this to prevent dividing by zero further on.
      if (vmup.Pz() != 0. && vmum.Pz() != 0.) {
	double mass = vdil.M();
	double rap  = vdil.Rapidity();
	double eta  = vdil.PseudoRapidity();

	//-----------------------------------------------------------------
	// Compute cosine theta^* in Gottfried-Jackson frame.
	// This frame is defined as the frame with the z^* axis parallel 
	// to the beam momentum.  This makes the pT of the proton > 0 in
	// the Z' CMS frame.  This is a good frame for p-pbar collisions
	// but a bad one for pp collisions because we don't know which
	// of the protons had the quark that interacted.

	// Lorentz boost mu- and beam to dilepton frame.
	pmum_star = LorentzBoost(vdil, vmum);
	pb_star = LorentzBoost(vdil, pp1);

	// Do Lorentz TX alternatively using Root routines and compare -Bob
	TVector3 beta;
	beta = -vdil.BoostVector();  //beta = (p_x/E, p_y/E, p_z/E)
	TLorentzVector pb_alt = pp1;
	pb_alt.Boost(beta);
	for (int i=0; i<4; i++) {
	  if (fabs( (pb_star[i]-pb_alt[i])/(fabs(pb_star[i])+1.e-6)) > 1.e-5) {
	    edm::LogWarning("fillFrameHistos") 
	      << "LTX diff for i=" << i << " pabstar = " 
	      << pb_star[i] << " pb_alt = " << pb_alt[i] << " diff = "
	      << pb_star[i]-pb_alt[i] << " relative " 
	      << fabs( (pb_star[i]-pb_alt[i])/(fabs(pb_star[i])+1.e-6));
	  }
	}

	// cosine Gottfried Jackson calculation.
	cos_gj = (pb_star.Vect()*pmum_star.Vect())/
	  ((pb_star.Vect().Mag())*(pmum_star).Vect().Mag());

	if (debug) 
	  LogTrace("fillFrameHistos")
	    << "Cosine theta^* in Gottfried-Jackson frame: " << cos_gj;
	// fill Cosin theta histogram and forward and backward mass histos.
	cosGJ[i_rec][0]->Fill(cos_gj);
	if (cos_gj > 0) {
	  FMassGJ[i_rec][0]->Fill(mass);
	  FRapGJ[i_rec][0]->Fill(rap);
	}
	else {
	  BMassGJ[i_rec][0]->Fill(mass);
	  BRapGJ[i_rec][0]->Fill(rap);
	}

	// Now adopt this frame to the pp collisions, by tagging the "beam"
	// proton using the direction of the dilepton system w.r.t the 
	// pp collision axis.

	if (pp1.Pz()*vdil.Pz() > 0.)
	  pb_star = LorentzBoost(vdil, pp1);
	else
	  pb_star = LorentzBoost(vdil, pp2);
	cos_gj_tag = (pb_star.Vect()*pmum_star.Vect())/
	  ((pb_star.Vect().Mag())*(pmum_star).Vect().Mag());
	if (debug)
	  LogTrace("fillFrameHistos")
	    << "Cos theta^* in tagged Gottfried-Jackson frame: " 
	    << cos_gj_tag;

	// Fill forward and backward mass, rapidity, and eta histograms
	cosGJ[i_rec][1]->Fill(cos_gj_tag);
	if (cos_gj_tag > 0) {
	  FMassGJ[i_rec][1]->Fill(mass);
	  FRapGJ[i_rec][1]->Fill(rap);
	  FPseudGJ[i_rec]->Fill(eta);
	}
	else {
	  BMassGJ[i_rec][1]->Fill(mass);
	  BRapGJ[i_rec][1]->Fill(rap);
	  BPseudGJ[i_rec]->Fill(eta);
	}

	//---------------------------------------------------------------
	// Compute cosine theta^* in Collins-Soper frame.  The z^* axis
	// is chosen so that it bisects the angle between the pBeam and
	// -pTarget.  This is the best frame in p-pbar collisions since
	// it reduces the uncertainty introduced by the fact that p and
	// p_bar are not parallel in the dilepton rest frame, and the
	// quark directions are not the same as the p and p_bar directions
	// Choose the pBeam using the sign of the dilepton system.
	//
	// Identify the beam as the p that is in the same direction as 
	// the dilepton.  Then boost the beam and target into the dilepton
	// CMS frame.
	if (pp1.Pz()*vdil.Pz() > 0.) {
	  pb_star = LorentzBoost(vdil, pp1);
	  pt_star = LorentzBoost(vdil, pp2);
	}
	else {
	  pb_star = LorentzBoost(vdil, pp2);
	  pt_star = LorentzBoost(vdil, pp1);
	}

	// Lengthy calculation for Cosine theta
	temp3 = (pb_star.Vect()*(1/pb_star.Vect().Mag()))-
	  (pt_star.Vect()*(1/pt_star.Vect().Mag()));
	cos_cs[i_rec][i_dil] = (temp3*pmum_star.Vect())/
	  (temp3.Mag()*pmum_star.Vect().Mag());

	if (debug)
	  LogTrace("fillFrameHistos")
	    << "Cosine theta^* in Collins-Soper frame: " 
	    << cos_cs[i_rec][i_dil];

	// Fill the cosine theta, mass, rap, and eta histograms.
	cosCS[i_rec][0]->Fill(cos_cs[i_rec][i_dil]);
	if (i_rec == 0 || i_rec == 3)
	  rap_vs_cosCS[i_rec]->Fill(cos_cs[i_rec][i_dil], rap);

	// A few resolution plots
	if (i_rec > 0 && allDileptons[lgen].size() == allDileptons[i_rec].size()) {
	  cosCSRes[i_rec-1]->Fill(cos_cs[i_rec][i_dil]-cos_cs[0][i_dil]);
	  if (i_rec == 3) {
	    cosCS3_diffsq_vs_cosCS0->Fill(cos_cs[0][i_dil],
                          (cos_cs[i_rec][i_dil]-cos_cs[0][i_dil])*
                          (cos_cs[i_rec][i_dil]-cos_cs[0][i_dil]));
	    rap3_vs_rap0->Fill(allDileptons[lgen][i_dil].rapidity(), rap);
	  }
	}
	if (cos_cs[i_rec][i_dil] > 0.) {
	  FMassCS[i_rec][0]->Fill(mass);
	  FRapCS[i_rec][0]->Fill(rap);
	  FPseudCS[i_rec]->Fill(eta);
	}
	else {
	  BMassCS[i_rec][0]->Fill(mass);
	  BRapCS[i_rec][0]->Fill(rap);
	  BPseudCS[i_rec]->Fill(eta);
	}

	// Now compute the analytic expression.  Why do we
	// get exactly the same answer as for the exact formula?
	double pplus_mum  = (1./sqrt(2.))*(vmum.E() + vmum.Pz());
	double pplus_mup  = (1./sqrt(2.))*(vmup.E() + vmup.Pz());
	double pminus_mum = (1./sqrt(2.))*(vmum.E() - vmum.Pz());
	double pminus_mup = (1./sqrt(2.))*(vmup.E() - vmup.Pz());
	double term_muon = (pplus_mum*pminus_mup)-(pminus_mum*pplus_mup);
	cos_cs_anal = 2.*term_muon/
	  (mass*(sqrt((mass*mass)+(vdil.Pt()*vdil.Pt()))));
	if (vdil.Pz() != 0.)
	  cos_cs_anal = (cos_cs_anal*vdil.Pz())/abs(vdil.Pz());
	if (debug)
	  LogTrace("fillFrameHistos")
	    << "Cosine theta^* in Collins-Soper analytic frame: "
	    << cos_cs_anal;

	// Fill the cosine theta, mass, rapidity, and eta histos.
	cosCS[i_rec][1]->Fill(cos_cs_anal);
	if (cos_cs_anal > 0) {
	  FMassCS[i_rec][1]->Fill(mass);
	  FRapCS[i_rec][1]->Fill(rap);
	}
	else {
	  BMassCS[i_rec][1]->Fill(mass);
	  BRapCS[i_rec][1]->Fill(rap);
	}

	//-----------------------------------------------------
	// Now calculate in Baur-boost frame
	// This frame take the quark direction as the boost direction
	// of the dilepton system.  This approach is described in
	// M. Dittmar, Phys. Rev. D55 (1997) 161;  U. Baur et al.,
	// hep-ph/9707301.
	cos_boost = (vdil.Vect()*pmum_star.Vect())/
	  (vdil.Vect().Mag()*pmum_star.Vect().Mag());
	if (debug)
	  LogTrace("fillFrameHistos")
	    << "Cosine theta^* in Baur-boost frame: " << cos_boost;

	// Again fill the histos.
	cosBoost[i_rec]->Fill(cos_boost);
	if (cos_boost > 0) {
	  FMassBoost[i_rec]->Fill(mass);
	  FRapBoost[i_rec]->Fill(rap);
	  FPseudBoost[i_rec]->Fill(eta);
	}
	else {
	  BMassBoost[i_rec]->Fill(mass);
	  BRapBoost[i_rec]->Fill(rap);
	  BPseudBoost[i_rec]->Fill(eta);
	}

	//See how this approximation works in various intervals of rapidity.
	if (cos_boost > 0.) {
	  if ((abs(rap)) < 0.4)
	    FMBoostCut[i_rec][0]->Fill(mass);
	  else if (fabs(rap) < 0.8)
	    FMBoostCut[i_rec][1]->Fill(mass);
	  else if (fabs(rap) < 2.4)
	    FMBoostCut[i_rec][2]->Fill(mass);

	  if (fabs(eta) < 0.4)
	    FMBoostCut[i_rec][3]->Fill(mass);
	  else if (fabs(eta) < 0.8)
	    FMBoostCut[i_rec][4]->Fill(mass);
	  else if (fabs(eta) < 2.4)
	    FMBoostCut[i_rec][5]->Fill(mass);
	}
	else {
	  if (fabs(rap) < 0.4)
	    BMBoostCut[i_rec][0]->Fill(mass);
	  else if (fabs(rap) < 0.8)
	    BMBoostCut[i_rec][1]->Fill(mass);
	  else if (fabs(rap) < 2.4)
	    BMBoostCut[i_rec][2]->Fill(mass);

	  if (fabs(eta) < 0.4)
	    BMBoostCut[i_rec][3]->Fill(mass);
	  else if (fabs(eta) < 0.8)
	    BMBoostCut[i_rec][4]->Fill(mass);
	  else if (fabs(eta) < 2.4)
	    BMBoostCut[i_rec][5]->Fill(mass);
	}

	//----------------------------------------------------
	// Wulz Frame
	// This is a variant of Baur frame as described in CMS-TN/93-107
	// Cosine theta^* is measured between mu+ and Z' for negative
	// rapidity values, and between mu- and Z' for positive rapidity.
	// The result is an odd function of rapidity and when integrated
	// over the full rapidity interval, gives zero.
	if (rap > 0.) {
	  cos_wulz = (vdil.Vect()*pmum_star.Vect())/
	    (vdil.Vect().Mag()*pmum_star.Vect().Mag());
	}
	else {
	  pmup_star = LorentzBoost(vdil, vmup);
	  cos_wulz = (vdil.Vect()*pmup_star.Vect())/
	    (vdil.Vect().Mag()*pmup_star.Vect().Mag());
	}
	if (debug)
	  LogTrace("fillFrameHistos")
	    << "Cosine theta^* in Wulz frame: " << cos_wulz;
	cosW[i_rec]->Fill(cos_wulz);
	if (cos_wulz > 0) {
	  FMassW[i_rec]->Fill(mass);
	  FRapW[i_rec]->Fill(rap);
	  FPseudW[i_rec]->Fill(eta);
	}
	else {
	  BMassW[i_rec]->Fill(mass);
	  BRapW[i_rec]->Fill(rap);
	  BPseudW[i_rec]->Fill(eta);  
	}
      }
    }
  }
}

double Zprime2muAsymmetry::calcAFBError(double f, double b) {
  return (f > 0 && b > 0) ? 2*f*b/(f+b)/(f+b)*sqrt(1/f+1/b) : 1;
}

void Zprime2muAsymmetry::calcAsymmetry(double f, double b,
				       double& A_FB, double& e_A_FB) {
  // Calculate A_FB and its error from f(orward) and b(ackward) counts
  // and return by reference.
  if (f+b > 0)
    A_FB = (f-b)/(f+b);
  else
    A_FB = 0;

  e_A_FB = calcAFBError(f,b);
}

void Zprime2muAsymmetry::calcAsymmetry(const TH1F* h_cos,
				       double& A_FB, double& e_A_FB) {
  // With h_cos a histogram of cos_theta from -1 to 1, calculate the
  // asymmetry using the integral from -1 to 0 and and the integral
  // from 0 to 1
  int nBins = h_cos->GetNbinsX();
  int nBack = int(h_cos->Integral(1, int(nBins/2.)));
  int nFor = int(h_cos->Integral((int(nBins/2.) + 1), nBins));
  calcAsymmetry(nFor, nBack, A_FB, e_A_FB);
}

void Zprime2muAsymmetry::calcAsymmetry(const TH1F* IdF, const TH1F* IdB,
				       TH1F* IdA, fstream& out) {
  // This routine takes a forward mass histogram (IdF), and a 
  // backward mass histogram (IdB), and creates an asymmetry 
  // histogram (IdA) with (IdF-IdB)/(IdF+IdB).

  // Get the number of bins in each histogram and check they are the same.
  int nBinsF = IdF->GetNbinsX();
  int nBinsB = IdB->GetNbinsX();
  if (nBinsF != nBinsB) {
    out << "error: nBinsF " << nBinsF << " != nBinsB " << nBinsB << endl;
    return;
  }

  // A_FB = (F-B)/(F+B)
  TH1F* temp1 = (TH1F*)IdF->Clone();  temp1->SetNameTitle("temp1", "F-B");
  TH1F* temp2 = (TH1F*)IdF->Clone();  temp2->SetNameTitle("temp2", "F+B");
  temp1->Add(IdF, IdB, 1., -1.);
  temp2->Add(IdF, IdB, 1.,  1.);
  IdA->Divide(temp1, temp2, 1., 1.);

  // Sigma(A_FB) = 2FB*Sqrt((SIGMA(F)/F)**2 + (SIGMA(B)/B)**2)/((F+B)**2)
  Stat_t f_bin, b_bin, sigma_bin;
  for (int ibin = 1; ibin <= nBinsF; ibin++) {
    f_bin = IdF->GetBinContent(ibin);
    b_bin = IdB->GetBinContent(ibin);
    sigma_bin = calcAFBError(f_bin, b_bin);
    IdA->SetBinError(ibin, sigma_bin);
  }

  // Find the asymmetry from the total number of forward and backward events.
  double F_total = IdF->GetEntries();
  double B_total = IdB->GetEntries();
  double A_FB, sigma_AS;
  calcAsymmetry(F_total, B_total, A_FB, sigma_AS);

  const int NBIN = 9;
  if (nBinsF == NBIN-1) {
    out << setprecision(5);
    out << "F = " << setw(5) << F_total << " B = " << setw(5) << B_total
	 << " A_FB = " << setw(5) << A_FB
	 << " +/- "    << setw(5) << sigma_AS << endl;
  }
}

void Zprime2muAsymmetry::calcFrameAsym() {
  string dumpfile = "frameAsym." + outputFileBase + ".txt";
  fstream out(dumpfile.c_str(), ios::out);
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    out << "Rec level: " << i_rec << endl;
    out << "Asymmetry in GJ frame:        ";
    calcAsymmetry(FMassGJ[i_rec][0], BMassGJ[i_rec][0],
		  AMassGJ[i_rec][0], out);
    calcAsymmetry(FRapGJ[i_rec][0],  BRapGJ[i_rec][0],  
		  ARapGJ[i_rec][0], out);
    out << "Asymmetry in tagged GJ frame: ";
    calcAsymmetry(FMassGJ[i_rec][1], BMassGJ[i_rec][1], 
		  AMassGJ[i_rec][1], out);
    calcAsymmetry(FRapGJ[i_rec][1],  BRapGJ[i_rec][1],  
		  ARapGJ[i_rec][1], out);
    out << "Asymmetry in CS frame:        ";
    calcAsymmetry(FMassCS[i_rec][0], BMassCS[i_rec][0], 
		  AMassCS[i_rec][0], out);
    calcAsymmetry(FRapCS[i_rec][0],  BRapCS[i_rec][0],
		  ARapCS[i_rec][0], out);
    out << "Asymmetry in analyt CS frame: ";
    calcAsymmetry(FMassCS[i_rec][1], BMassCS[i_rec][1],
		  AMassCS[i_rec][1], out);
    calcAsymmetry(FRapCS[i_rec][1],  BRapCS[i_rec][1],
		  ARapCS[i_rec][1], out);
    out << "Asymmetry in Boost frame:     ";
    calcAsymmetry(FMassBoost[i_rec], BMassBoost[i_rec],
		  AMassBoost[i_rec], out);
    calcAsymmetry(FRapBoost[i_rec],  BRapBoost[i_rec],
		  ARapBoost[i_rec], out);
    out << "Asymmetry in Wulz frame:      ";
    calcAsymmetry(FMassW[i_rec], BMassW[i_rec], AMassW[i_rec], out);
    calcAsymmetry(FRapW[i_rec], BRapW[i_rec], ARapW[i_rec], out);
    
    for (int i = 0; i < 6; i++) {
      out << "Asymmetry in MBCut[" << i << "] frame:  ";
      calcAsymmetry(FMBoostCut[i_rec][i], BMBoostCut[i_rec][i],
		    AsymMBoostCut[i_rec][i], out);
    }
  }
  out.close();
}

void Zprime2muAsymmetry::dumpFitData() {
  string fn = "dumpFitData." + outputFileBase + ".txt";
  fstream f(fn.c_str(), ios::out);

  for (int i = 0; i < 2; i++) {
    f << "#i=" << i << " nfit_used[i]=" << nfit_used[i] << endl;
    for (int j = 0; j < nfit_used[i]; j++)
      f << cos_theta_cs_data[i][j] << '\t' << mass_dil_data[i][j] << '\t'
	<< rap_dil_data[i][j] << '\t' << pt_dil_data[i][j] << '\t'
	<< phi_dil_data[i][j] << '\t' << phi_cs_data[i][j] << endl;
  }

  f.close();

  // dump the generated pre-acceptance-cut cos_cs distribution
  fn = "angDist." + outputFileBase + ".txt";
  fstream f2(fn.c_str(), ios::out);
  f2 << "# angular distribution\n"
     << "# format: cos_theta_cs \t mistag?\n";
  for (unsigned int i = 0; i < angDist.size(); i++)
    f2 << angDist[i] << endl;
  f2.close();
}

void Zprime2muAsymmetry::bookParamHistos() {
  nMistagBins = 20;
  TH1::AddDirectory(false);
  const AsymFitManager& amg = asymFitManager;
  // Histograms here are not managed by TFileService, but by ourselves
  // (since we save some of them to the parameter cache file later).
  h_pt_dil               = new TH1F("3", "pT^2 dilepton",     100, 0., 3.e6);
  h_rap_dil[0]           = new TH1F("4", "abs(rap) dilepton", 50, 0., 3.5);
  h_rap_dil[1]           = new TH1F("5", "rap dilepton",    50, -3.5, 3.5);
  h_mass_dil[0] = new TH1F("6", "mass dil, full range",
			   100, amg.gen_win(0) - 50., amg.gen_win(1) + 50);
  h_mass_dil[1] = new TH1F("7", "mass dil, fit range",
			   100, amg.fit_win(0), amg.fit_win(1));
  h_phi_cs       = new TH1F("9", "Collins-Soper Phi", 100, 0., TMath::Pi());
  h_rap_mistag   = new TH1F("10","rap dil, mistag=1", 50, 0., 3.5);
  h_rap_nomistag = new TH1F("11","rap dil, mistag=0", 50, 0., 3.5);
  h_cos_mistag   = new TH1F("10","cos#theta_{CS}, mistag^{1}#neq mistag^{2} ", 
			    50, 0., 1.);
  h_cos_cs       = new TH1F("11","cos#theta_{CS}", 50, 0., 1.);
  h_cos_mistag_prob = new TH1F("10/4","cos#theta_{CS} mistag frac",
			       50, 0., 1.);
  h_cos_mistag_prob->Sumw2();
  h2_rap_cos_mistag = new TH2F("12", "rap vs cos dil, mistag=1", 
			       nMistagBins, 0., 1., nMistagBins, 0., 3.5);
  h2_rap_cos_nomistag = new TH2F("13", "rap vs cos dil, mistag=0", 
				 nMistagBins, 0., 1., nMistagBins, 0., 3.5);
  h2_rap_cos_p = new TH2F("14", "rap vs cos dil", nMistagBins, 
			  0., 1., nMistagBins, 0., 3.5);
  h_mistag[0][0] = new TH1F("mist01", "rap all ", 50, 0., 3.5);
  h_mistag[0][1] = new TH1F("mist02", "rap mist==1", 50, 0., 3.5);
  h_mistag[0][2] = new TH1F("mist03", "rap frac mist", 50, 0., 3.5);
  h_mistag[0][2]->Sumw2();
  h_mistag[1][0] = new TH1F("mist11", "pT all ", 50, 0., 20000.);
  h_mistag[1][1] = new TH1F("mist12", "pT mist==1", 50, 0., 20000.);
  h_mistag[1][2] = new TH1F("mist13", "pT frac mist", 50, 0., 20000.);
  h_mistag[1][2]->Sumw2();
  h_mistag[2][0] = new TH1F("mist21", "mass all ",
			    50, amg.fit_win(0), amg.fit_win(1));
  h_mistag[2][1] = new TH1F("mist22", "mass mist==1", 
			    50, amg.fit_win(0), amg.fit_win(1));
  h_mistag[2][2] = new TH1F("mist23", "mass frac mist",
			    50, amg.fit_win(0), amg.fit_win(1));
  h_mistag[2][2]->Sumw2();
  h_mistag[3][0] = new TH1F("mist31", "cos_true all ", 50, 0., 1.);
  h_mistag[3][1] = new TH1F("mist32", "cos_true mist==1", 50, 0., 1.);
  h_mistag[3][2] = new TH1F("mist33", "cos_true frac mist", 50, 0., 1.);
  h_mistag[3][2]->Sumw2();
  h_mistag[4][0] = new TH1F("mist41", "dil pL all ",      50, 0., 5000.);
  h_mistag[4][1] = new TH1F("mist42", "dil pL mist==1",   50, 0., 5000.);
  h_mistag[4][2] = new TH1F("mist43", "dil pL frac mist", 50, 0., 5000.);
  h_mistag[4][2]->Sumw2();
  h_mistag[5][0] = new TH1F("mist41", "q pL all ",      50, 0., 1000.);
  h_mistag[5][1] = new TH1F("mist42", "q pL mist==1",   50, 0., 1000.);
  h_mistag[5][2] = new TH1F("mist43", "q pL frac mist", 50, 0., 1000.);
  h_mistag[5][2]->Sumw2();
  h2_mistag[0] = new TH2F("mist41", "q dil pL all", 
			50, 0., 5000., 50, 0., 5000.);
  h2_mistag[1] = new TH2F("mist42", "q dil pL mist==1",   
			50, 0., 5000., 50, 0., 5000.);
  h2_mistag[2] = new TH2F("mist41", "q dil pL all (reduced range)", 
			25, 0., 500., 25, 0., 5000.);
  h2_mistag[3] = new TH2F("mist42", "q dil pL mist==1 (reduced range)",   
			25, 0., 500., 25, 0., 5000.);
  h2_mistagProb = new TH2F("h2_mistagProb",
			   "Mistag prob cos rap", nMistagBins, 
			   0., 1., nMistagBins, 0., 3.5);
  h2_pTrap = new TH2F("h2_pTrap", "Rap vs pT", nMistagBins, 0., 700., 
		      nMistagBins, 0., 3.5);
  h2_pTrap_mistag = new TH2F("h2_pTrap", "Rap vs pT, mistag == 1", nMistagBins,
			     0., 700., nMistagBins, 0., 3.5);
  h_rap_mistag_prob = new TH1F("h_rap_mistag_prob",
			       "Mistag prob rap", 50, 0., 3.5);
  h_cos_true_mistag_prob = new TH1F("h_cos_true_mistag_prob",
				    "Mistag prob cos", 50, 0., 1.);
  h_pL_mistag_prob = new TH1F("h_pL_mistag_prob",
			      "Mistag prob pL", 50, 0., 5000.);
  h2_pL_mistag_prob = new TH2F("h2_pL_mistag_prob",
			       "Mistag prob dil pL vs quark pL", 
			       25, 0., 500., 25, 0., 5000.);
  h2_cos_cs_vs_true = new TH2F("h2_cos_cs_vs_true",
			       "cos #theta_{true} vs cos #theta_{CS}",
			       50, -1., 1., 50, -1., 1.);
  h2_pTrap_mistag_prob = new TH2F("h2_pTrap_mistag_prob",
				  "Mistag prob rap vs pT",
				  nMistagBins, 0., 700., nMistagBins, 0., 3.5);
}

void Zprime2muAsymmetry::fillParamHistos(bool fakeData) {
  // JMTBAD will we want to change these ever?
  const bool PARAM_MISTAG = true;
  const bool PARAM_MASS = true;
  const bool PARAM_PT = true;
  const bool PARAM_RAP = true;
  const bool PARAM_PHICS = true;
  const bool IGNORE_MASS_RANGE = false;

  // Open the root files containing the extra generated sample, from which
  // we extract the mistag parameterization
  fwlite::ChainEvent ev(genSampleFiles);

  int nentries = maxParamEvents > 0 ? maxParamEvents : ev.size();

  edm::LogInfo("Zprime2muAsymmetry")
    << "Loading sample parameters for asym fit using " << nentries
    << " events from " << genSampleFiles.size() << " files beginning with "
    << genSampleFiles[0];

  int jentry = 0;
  for (ev.toBegin(); !ev.atEnd() && jentry < nentries; ++ev, ++jentry) {
    fwlite::Handle<reco::CandidateCollection> genParticles;
    genParticles.getByLabel(ev, "genParticleCandidates");

    if (verbosity >= VERBOSITY_SIMPLE && jentry % 1000 == 0)
      LogTrace("fillParamHistos") << "fillParamHistos: " << jentry
				  << " events processed";

    // Remove effect from internal brem
    AsymFitData data;
    if (!computeFitQuantities(genParticles.ref(), jentry, data))
      continue; // if finding the Z'/etc failed, skip this event

    /*
    // For some studies, it is useful to replace real data with data
    // artificially created with any parameterization desired
    if (fakeData) {
      data.rapidity    = fake_rap[jentry];
      data.cos_true    = fake_cos_true[jentry];
      data.cos_cs      = fake_cos_cs[jentry];
      data.mistag_true = fake_mistag_true[jentry];
      data.mistag_cs = fake_mistag_cs[jentry];
    }
    */

    // Store full mass spectrum of sample
    if (PARAM_MASS)
      h_mass_dil[0]->Fill(data.mass); 

    // Fill events over window where sample was generated
    if ((data.mass > asymFitManager.fit_win(0) &&
	 data.mass < asymFitManager.fit_win(1)) || 
	IGNORE_MASS_RANGE) {

      // Fill mistag histos
      if (PARAM_MISTAG) {
	if (data.mistag_true == 0)
	  h_rap_nomistag->Fill(fabs(data.rapidity));
	else {
	  h_rap_mistag->Fill(fabs(data.rapidity));
	  h_mistag[0][1]->Fill(fabs(data.rapidity));
	  h_mistag[1][1]->Fill(data.pT*data.pT);
	  h_mistag[2][1]->Fill(data.mass);
	  h_mistag[3][1]->Fill(fabs(data.cos_true));
	  h_mistag[4][1]->Fill(fabs(data.pL));
	  h_mistag[5][1]->Fill(fabs(data.qpL));
	  h2_mistag[1]->Fill(fabs(data.qpL), fabs(data.pL));
	  h2_mistag[3]->Fill(fabs(data.qpL), fabs(data.pL));
	  h2_pTrap_mistag->Fill(fabs(data.pT*data.pT), fabs(data.rapidity));
	}
	h_mistag[0][0]->Fill(fabs(data.rapidity)); 
	h_mistag[1][0]->Fill(data.pT*data.pT);
	h_mistag[2][0]->Fill(data.mass);
	h_mistag[3][0]->Fill(fabs(data.cos_true));
	h_mistag[4][0]->Fill(fabs(data.pL));
	h_mistag[5][0]->Fill(fabs(data.qpL));
	h2_mistag[0]->Fill(fabs(data.qpL), fabs(data.pL));
	h2_mistag[2]->Fill(fabs(data.qpL), fabs(data.pL));
	h2_pTrap->Fill(fabs(data.pT*data.pT), fabs(data.rapidity));

	// 2D mistag probability histos for cos_cs and rapidity
	if (data.mistag_cs == 0)
	  h2_rap_cos_nomistag->Fill(fabs(data.cos_cs), fabs(data.rapidity)); 
	else
	  h2_rap_cos_mistag->Fill(fabs(data.cos_cs), fabs(data.rapidity)); 
	h2_rap_cos_p->Fill(fabs(data.cos_cs), fabs(data.rapidity));
	
	// These histos give small correction from fact of using CS cos 
	// instead of true
	if (data.mistag_true != data.mistag_cs)
	  h_cos_mistag->Fill(fabs(data.cos_cs)); 
	h_cos_cs->Fill(fabs(data.cos_cs));

	// If you want to do convolution between CS cos and true
	if (data.mistag_cs == data.mistag_true) {
	  if (data.mistag_true == 0)
	    h2_cos_cs_vs_true->Fill(data.cos_true, data.cos_cs);
	  else 
	    h2_cos_cs_vs_true->Fill(data.cos_true, -data.cos_cs);
	}
      }

      // Fill mass histos, one histo with all events in generated window,
      // one histo for events in reconstructed mass window, and another 
      // for events within pythia signal window.
	  
      if (PARAM_MASS)
	h_mass_dil[1]->Fill(data.mass);

      // Parametrization for pT is done with pT^2.
      if (PARAM_PT)
	h_pt_dil->Fill(data.pT*data.pT);

      // Fill histo with rapidity
      if (PARAM_RAP) {
	h_rap_dil[0]->Fill(fabs(data.rapidity));
	h_rap_dil[1]->Fill(data.rapidity);
      }

      // Phi Collins-Soper is filled from 0. to Pi. (since probability
      // between Pi and 2Pi is same as from 0 to Pi). 
      if (PARAM_PHICS) {
	if (data.phi_cs > TMath::Pi()) 
	  h_phi_cs->Fill(data.phi_cs-TMath::Pi());
	else
	  h_phi_cs->Fill(data.phi_cs);
      }
    }
  }
  
  if (maxParamEvents < 0) maxParamEvents = jentry;
  LogTrace("fillParamHistos") << "fillParamHistos: " << jentry
			      << " events processed";
}

// This routine calculates the quantities to be used in the multiple
// dimension fit. Variables used in the fit are stored in the
// structure "data" returned by reference.  The flag "internalBremOn",
// is for switching on and off the effect of internal bremsstrahlung.
bool Zprime2muAsymmetry::computeFitQuantities(const reco::CandidateCollection& genParticles, 
					      int eventNum,
					      AsymFitData& data) {
  static const bool debug = verbosity >= VERBOSITY_TOOMUCH;

  InteractionParticles ip;
  if (!storeInteractionParticles(genParticles, eventNum, ip))
    return false;

  // Copy the four-vectors into TLorentzVectors, since our code uses
  // those already
  TLorentzVector v_dil, v_mum, v_mup, v_muq;
  TLorentzVector v_my_dil, v_my_mum, v_my_mup;

  v_muq.SetPxPyPzE(ip.genQuark->p4().x(), ip.genQuark->p4().y(), ip.genQuark->p4().z(),
		   ip.genQuark->p4().e());
  if (internalBremOn) {
    v_my_mum.SetPxPyPzE(ip.genLepMinus->p4().x(), ip.genLepMinus->p4().y(),
			ip.genLepMinus->p4().z(), ip.genLepMinus->p4().e());
    v_my_mup.SetPxPyPzE(ip.genLepPlus->p4().x(), ip.genLepPlus->p4().y(),
			ip.genLepPlus->p4().z(), ip.genLepPlus->p4().e());
  } 
  else {
    v_my_mum.SetPxPyPzE(ip.genLepMinusNoIB->p4().x(), ip.genLepMinusNoIB->p4().y(),
			ip.genLepMinusNoIB->p4().z(), ip.genLepMinusNoIB->p4().e());
    v_my_mup.SetPxPyPzE(ip.genLepPlusNoIB->p4().x(), ip.genLepPlusNoIB->p4().y(),
			ip.genLepPlusNoIB->p4().z(), ip.genLepPlusNoIB->p4().e());
  }

  // The 4-vector for the dimuon
  v_my_dil = v_my_mum + v_my_mup;

  // Store pt, rap, phi, mass, phi_cs, and cos_cs in a structure
  data.cos_true = calcCosThetaTrue(v_muq, v_my_mum, v_my_dil, debug);
  data.pT = v_my_dil.Pt();
  data.rapidity = v_my_dil.Rapidity();
  data.phi = v_my_dil.Phi();
  data.mass = v_my_dil.M(); 
  data.qpL = v_muq.Pz();
  data.pL = v_my_dil.Pz();

  // Cosine theta^* between mu- and the quark in the Z' CMS frame ("true"
  // cos theta^*).
  // JMTBAD calculated by calcCosThetaTrue above already
  //math::XYZTLorentzVector mum_star = LorentzBoost(ip.genMom.res, ip.genMom.mum);
  //math::XYZTLorentzVector quark_star = LorentzBoost(ip.genMom.res, ip.genMom.quark);
  //float gen_cos_true = cosTheta(mum_star, quark_star);
  
  data.cos_cs = calcCosThetaCSAnal(v_my_mum.Pz(), v_my_mum.E(),
				   v_my_mup.Pz(), v_my_mup.E(),
				   v_my_dil.Pt(), v_my_dil.Pz(), 
				   data.mass, debug);
  data.phi_cs = calcPhiCSAnal(v_my_mum.Px(), v_my_mum.Py(),
			      v_my_mup.Px(), v_my_mup.Py(),
			      v_my_dil.Pt(), v_my_dil.Eta(), v_my_dil.Phi(),
			      data.mass, debug);

  // we will later assume the p_z of the dilepton is the p_z of the
  // quark we've mistagged this if the quark was actually going the
  // other direction.  translating from kir_anal.F, we looked at the
  // dilepton formed from adding the muon final-state four-vectors
  // (i.e., we ignored small effect of any brem):
  // data.mistag_true = (v_my_dil.Pz()*ip.genQuark->p4().z() < 0) ? 1 : 0;
  // but, since we have the true resonance four-vector:
  data.mistag_true = (ip.genResonance->p4().z()*ip.genQuark->p4().z() < 0) ? 1 : 0;
  // alternatively, we've mistagged if the signs of cos_cs and
  // cos_true are different
  data.mistag_cs = (data.cos_cs/data.cos_true > 0) ? 0 : 1;

  // Do some extra calculations to determine if event passed acceptance.
  calc4Vectors(data, v_dil, v_mum, v_mup, debug);
  data.cut_status = diRapAccept(v_dil, v_mum, v_mup);

  return true;
}

void Zprime2muAsymmetry::getAsymParams() {
  TFile* paramFile = 0; 
    
  if (useCachedParams) {
    // first try to get the already calculated params from the file
    paramFile = new TFile(paramCacheFile.c_str(), "read");
    if (paramFile->IsOpen()) {
      // JMTBAD from here, assume that all the reading of the file works
      edm::LogInfo("getAsymParams") << "Using cached parameterizations from "
				    << paramCacheFile;
      TArrayD* arr;
      paramFile->GetObject("mistag_pars", arr);
      for (int i = 0; i < 6; i++) mistag_pars[i] = arr->At(i);
      paramFile->GetObject("mass_pars", arr);
      for (int i = 0; i < 7; i++) mass_pars[i] = arr->At(i);
      paramFile->GetObject("rap_pars", arr);
      for (int i = 0; i < 5; i++) rap_pars[i] = arr->At(i);
      paramFile->GetObject("pt_pars", arr);
      for (int i = 0; i < 5; i++) pt_pars[i] = arr->At(i);
      paramFile->GetObject("phi_cs_pars", arr);
      for (int i = 0; i < 5; i++) phi_cs_pars[i] = arr->At(i);
      TArrayI* iarr;
      paramFile->GetObject("nMistagBins", iarr);
      nMistagBins = iarr->At(0);
      // JMTBAD having to delete the booked histograms
      delete h2_mistagProb;
      h2_mistagProb = (TH2F*)paramFile->Get("h2_mistagProb");
      h2_mistagProb->SetDirectory(0);
      delete h_rap_mistag_prob;
      h_rap_mistag_prob = (TH1F*)paramFile->Get("h_rap_mistag_prob");
      h_rap_mistag_prob->SetDirectory(0);
      delete h_cos_true_mistag_prob;
      h_cos_true_mistag_prob = (TH1F*)paramFile->Get("h_cos_true_mistag_prob");
      h_cos_true_mistag_prob->SetDirectory(0);
      delete h_pL_mistag_prob;
      h_pL_mistag_prob = (TH1F*)paramFile->Get("h_pL_mistag_prob");
      h_pL_mistag_prob->SetDirectory(0);
      delete h2_pL_mistag_prob;
      h2_pL_mistag_prob = (TH2F*)paramFile->Get("h2_pL_mistag_prob");
      h2_pL_mistag_prob->SetDirectory(0);
      delete h2_cos_cs_vs_true;
      h2_cos_cs_vs_true = (TH2F*)paramFile->Get("h2_cos_cs_vs_true");
      h2_cos_cs_vs_true->SetDirectory(0);
      delete h2_pTrap_mistag_prob;
      h2_pTrap_mistag_prob = (TH2F*)paramFile->Get("h2_pTrap_mistag_prob");
      h2_pTrap_mistag_prob->SetDirectory(0);
      paramFile->Close();
      delete paramFile;
      return;
    }
  }
  
  fillParamHistos(false);

  TCanvas *c1 = new TCanvas("c1", "", 0, 1, 500, 700);
  c1->SetTicks();

  string filename = "fitParams.";
  if (internalBremOn) filename += "ib.";
  if (onPeak)
    filename += "on";
  else
    filename += "off";
  filename += "_peak." + outputFileBase + ".ps";

  TPostScript* ps = new TPostScript(filename.c_str(), 111);

  const int NUM_PAGES = 15;
  TPad *pad[NUM_PAGES];
  for (int i_page=0; i_page<NUM_PAGES; i_page++)
    pad[i_page] = new TPad("","", .05, .05, .95, .93);

  ostringstream page_print;
  int page = 0;

  TLatex ttl;

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  ttl.DrawLatex(.4, .95, "y of #mu^{+}#mu^{-} (mistag^{1})");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  TH1F* h_temp_rap = (TH1F*) h_rap_mistag->Clone();
  h_temp_rap->Add(h_rap_nomistag);
  h_temp_rap->Draw();
  pad[page]->cd(2); 
  h_rap_mistag->Draw();
  pad[page]->cd(3); 
  h_rap_nomistag->Draw();
  h_rap_mistag_prob->Sumw2();
  h_rap_mistag_prob->Divide(h_rap_mistag, h_temp_rap, 1., 1.);
  pad[page]->cd(4);  h_rap_mistag_prob->Draw();
  
  const double MISTAG_LIM = 3.0;

  //Fit rapidity efficiency to 0.5 + ax + bx^2
  TF1 *f_rapmis=new TF1("f_rapmis","0.5+[0]*x+[1]*x*x", 0. , MISTAG_LIM);
  f_rapmis->SetParameters(0., -0.2);
  f_rapmis->SetParNames("p1","p2");
  // JMTBAD here we use cout explicitly so that it will go to the same
  // place as MINUIT's output in order to annotate it
  cout << "\n#### Fitting quadratic f(x=0)=0.5 to rapidity mistag prob" << endl;
  h_rap_mistag_prob->Fit("f_rapmis","VIR");
  double min_rap = f_rapmis->GetMinimumX(0., MISTAG_LIM);
  cout << "\n#### Minimum of function " << f_rapmis->GetMinimum(0., MISTAG_LIM)
       << " occurs at rap = " << min_rap << endl;
  mistag_pars[0] = f_rapmis->GetParameter(0);
  mistag_pars[1] = f_rapmis->GetParameter(1);
  mistag_pars[2] = min_rap;
  delete f_rapmis;
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  ttl.DrawLatex(.4, .95, 
	  "cos#theta^{*}_{CS} of #mu^{+}#mu^{-} (mistag^{1} #neq mistag^{2})");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  h_cos_mistag->Draw();
  pad[page]->cd(2); 
  h_cos_cs->Draw();
  pad[page]->cd(3); 
  h_cos_mistag_prob->Divide(h_cos_mistag, h_cos_cs, 1., 1.);
  h_cos_mistag_prob->Draw();  
  pad[page]->cd(4); 
  gPad->SetLogy(1);
  h_cos_mistag_prob->Draw();  
  TF1 *f_cosmis=new TF1("f_cosmis","[0]+[1]*exp(-[2]*x)", 0. , 1.);
  f_cosmis->SetParameters(.25, .25, 10.);
  f_cosmis->SetParNames("p0","p1","p2");
  f_cosmis->SetParLimits(0, 0., 1.);
  f_cosmis->SetParLimits(1, 0., 1.);
  f_cosmis->SetParLimits(2, 0., 50.);
  cout << "\n#### Fitting falling exp. to cos mistag prob" << endl;
  h_cos_mistag_prob->Fit("f_cosmis","VIR");
  mistag_pars[3] = f_cosmis->GetParameter(0);
  mistag_pars[4] = f_cosmis->GetParameter(1);
  mistag_pars[5] = f_cosmis->GetParameter(2);
  delete f_cosmis;
  page++;
  c1->Update();
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "y vs cos#theta^{*}_{CS} (mistag^{2})");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  gPad->SetPhi(210); h2_rap_cos_mistag->Draw("lego2");
  pad[page]->cd(2);
  gPad->SetPhi(210); h2_rap_cos_nomistag->Draw("lego2");
  pad[page]->cd(3);
  gPad->SetPhi(210); h2_rap_cos_p->Draw("lego2");
  h2_mistagProb->Divide(h2_rap_cos_mistag, h2_rap_cos_p, 1., 1.);
  pad[page]->cd(4);  
  gPad->SetPhi(210); h2_mistagProb->Draw("lego2");
  page++;
  c1->Update();  
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);
  ttl.DrawLatex(.4, .95, "Mass of #mu^{+}#mu^{-}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  TF1 *f_mass;
  int nPars = 3;
  double par_norm = 7000.; // initial values of normalization constant,
  double par_mean = peakMass; // mean, and
  double par_fwhm = 2.*sqrt(par_mean);      // fwhm

  const AsymFitManager& amg = asymFitManager;

  if (amg.mass_type() == MASS_EXP) {
    cout << "\n#### fitting falling exp to dil mass\n";
    f_mass = new TF1("f_mass", expBckg,
		     amg.fit_win(0), amg.fit_win(1), nPars);
    f_mass->SetParNames("Norm",  "Slope", "Integral");
    f_mass->SetParameters(par_norm, 0., 1.);
    //f_mass->SetParLimits(0, 100., 1.e9);
    //f_mass->SetParLimits(1, 0., -10.);
    f_mass->FixParameter(2, 1.);
  }
  else if (amg.mass_type() == MASS_LOR) {
    cout << "\n#### fitting Lorentzian to dil mass\n";
    f_mass = new TF1("f_mass", Lorentzian, 
		     amg.fit_win(0), amg.fit_win(1), nPars);
    f_mass->SetParNames("Norm", "FWHM", "Mean");
    f_mass->SetParameters(par_norm, par_fwhm, par_mean);
  }
  else if (amg.mass_type() == MASS_LOREXP) {
    nPars = 6;
    cout << "\n#### fitting Lorentzian plus exp background to dil mass\n";
    f_mass = new TF1("f_mass", lorentzianPlusExpbckg, 
		     amg.fit_win(0), amg.fit_win(1), nPars);
    f_mass->SetParNames("NormSign", "FWHM", "Mean", "NormBckg",
			   "SlopeBckg", "IntBckg");
    f_mass->SetParameters(par_norm, par_fwhm, par_mean, 1000., -0.01, 1.);
    // set these limits
    f_mass->SetParLimits(0, 0, 1e9);
    f_mass->SetParLimits(1, 0., 1000.);
    f_mass->FixParameter(2, par_mean);
    f_mass->SetParLimits(3, 0., 1e9);
    f_mass->FixParameter(5, 1.);
  }
  else
    throw cms::Exception("Zprime2muAsymmetry")
      << amg.mass_type() << " is not a known mass fit type!\n";

  h_mass_dil[0]->Draw(); //
  pad[page]->cd(2); h_mass_dil[1]->Fit("f_mass", "VIL", "", amg.fit_win(0), amg.fit_win(1));;
  double mass_norm = f_mass->Integral(amg.fit_win(0), amg.fit_win(1));
  for (int i = 0; i < nPars; i++) mass_pars[i] = f_mass->GetParameter(i);
  mass_pars[nPars] = mass_norm;
  delete f_mass;
  page++;
  c1->Update();
    
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "Fits to y of #mu^{+}#mu^{-}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  cout << "\n#### Fitting revised thermalized cylinder model to rapidity h_rap_dil\n";
  TF1* f_rap_rtc = new TF1("f_rap_rtc", "[0]*(tanh(([1]*x)+[2]+[3])-tanh(([1]*x)-[2]+[3])+tanh(([1]*x)+[2]-[3])-tanh(([1]*x)-[2]-[3]))", 0., 4.0);
  f_rap_rtc->SetParameters(50., 1., 1., 1.);
  f_rap_rtc->SetParNames("p0", "p1", "p2", "p3");
  h_rap_dil[0]->Fit("f_rap_rtc", "VL", "", 0., 3.5);
  for (int i = 0; i < 4; i++) rap_pars[i] = f_rap_rtc->GetParameter(i);
  rap_pars[4] = 2.*f_rap_rtc->Integral(0., 3.5);
  delete f_rap_rtc;
  page++;
  c1->Update();
    
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "pT^{2} of #mu^{+}#mu^{-}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  gPad->SetLogy(1);
  h_pt_dil->Draw();
  h_pt_dil->SetName("2 term exponential fit"); h_pt_dil->Draw();
  cout << "\n#### Fitting 2 term exponential to dilepton pT^2" << endl;
  //TF1 *f_2exp = new TF1("f_2exp", 
  //"[0]*exp(-[1]*sqrt(x))*(x<5.e4)+[2]*exp(-[3]*sqrt(x))*(x>5.e4)", 0., 3.e6);
  TF1 *f_2exp = new TF1("f_2exp", 
       		"[0]*exp(-[1]*sqrt(x))+[2]*exp(-[3]*sqrt(x))", 0., 3.e6);
  f_2exp->SetParameters(1.3e6, 3.9e-2, 9.5e3, 1.2e-2);
  f_2exp->SetParNames("p0","p1","p2","p3");
  f_2exp->SetParLimits(0, 0., 1.e7);
  f_2exp->SetParLimits(1, 0., 1.);
  f_2exp->SetParLimits(2, 0., 1.e7);
  f_2exp->SetParLimits(3, 0., 1.);
  h_pt_dil->Fit("f_2exp", "VIL", "", 0., 3.e6);
  for (int i = 0; i < 4; i++) pt_pars[i] = f_2exp->GetParameter(i);
  pt_pars[4] = f_2exp->Integral(0., 3.e6);
  delete f_2exp;
  page++;
  c1->Update();
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "#phi^{*}_{CS} of #mu^{+}#mu^{-}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  gStyle->SetOptStat(0000);  h_phi_cs->Draw();
  h_phi_cs->SetName("x^y fit"); h_phi_cs->Draw();
  cout << "\n#### Fitting to Collins-Soper phi" << endl;
  TF1 *f_phics = new TF1("f_phics", "[0]+[1]*((x-[2])^8)", 0., 3.14);
  f_phics->SetParNames("p0","p1","p2");
  f_phics->SetParameter(0, 500.);
  f_phics->SetParameter(1, 1.);
  // set these limits
  f_phics->SetParLimits(0, 0, 1e6);
  f_phics->SetParLimits(1, 0, 1e6);
  f_phics->FixParameter(2, 1.57);
  h_phi_cs->Fit("f_phics", "VIL", "", 0., 3.14);
  if (internalBremOn) {
    for (int i = 0; i < 3; i++) phi_cs_pars[i] = f_phics->GetParameter(i);
    phi_cs_pars[3] = 8.;
    phi_cs_pars[4] = f_phics->Integral(0., 3.14);
  }
  else
    for (int i = 0; i < 5; i++) phi_cs_pars[i] = -999.;
  delete f_phics;
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "cos #theta_{true} vs cos #theta_{CS}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->cd(0);
  h2_cos_cs_vs_true->Draw("lego2");
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  ttl.DrawLatex(.3, .95, "Mistag Probability for Various Quantities");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(3, 6);
  for (int i = 0; i < 6; i++) {
    h_mistag[i][2]->Divide(h_mistag[i][1], h_mistag[i][0], 1., 1., "B");
    pad[page]->cd(3*i+1);
    h_mistag[i][0]->Draw();
    pad[page]->cd(3*i+2);
    h_mistag[i][1]->Draw();
    pad[page]->cd(3*i+3);
    h_mistag[i][2]->Draw("hist E");
  }

  // Calculate histogram which can be used to obtain mistag probability
  h_cos_true_mistag_prob->Divide(h_mistag[3][1], h_mistag[3][0], 1., 1., "B");
  h_pL_mistag_prob->Divide(h_mistag[4][1], h_mistag[4][0], 1., 1., "B");
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  ttl.DrawLatex(.3, .95, "2D Mistag Probability (quark and dil pL)");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 3);
  for (int i = 0; i < 4; i++) {
    h2_mistag[i]->GetXaxis()->SetTitle("quark pL");
    h2_mistag[i]->GetYaxis()->SetTitle("dilepton pL");
    h2_mistag[i]->GetYaxis()->SetTitleOffset(1.6);
    pad[page]->cd(i+1); h2_mistag[i]->Draw();
  }
  h2_pL_mistag_prob->Sumw2();
  h2_pL_mistag_prob->Divide(h2_mistag[3], h2_mistag[2], 1., 1., "B");
  pad[page]->cd(5); 
  h2_pL_mistag_prob->GetXaxis()->SetTitle("quark pL");
  h2_pL_mistag_prob->GetXaxis()->SetTitleOffset(1.4);
  h2_pL_mistag_prob->GetYaxis()->SetTitle("dilepton pL");
  h2_pL_mistag_prob->GetYaxis()->SetTitleOffset(2.0);
  h2_pL_mistag_prob->GetZaxis()->SetTitle("mistag probability");
  h2_pL_mistag_prob->GetZaxis()->SetTitleOffset(1.5);
  gPad->SetPhi(210); h2_pL_mistag_prob->Draw("lego2");
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  ttl.DrawLatex(.3, .95, "2D Mistag Probability dil pT and rap");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 2);
  pad[page]->cd(1); h2_pTrap_mistag->Draw();
  pad[page]->cd(2); h2_pTrap->Draw();
  pad[page]->cd(3);
  h2_pTrap_mistag_prob->Sumw2();
  h2_pTrap_mistag_prob->Divide(h2_pTrap_mistag, h2_pTrap, 1., 1., "B");
  h2_pTrap_mistag_prob->GetXaxis()->SetTitle("pT*pT");
  h2_pTrap_mistag_prob->GetXaxis()->SetTitleOffset(1.4);
  h2_pTrap_mistag_prob->GetYaxis()->SetTitle("Y");
  h2_pTrap_mistag_prob->GetYaxis()->SetTitleOffset(2.0);
  h2_pTrap_mistag_prob->GetZaxis()->SetTitle("mistag probability");
  h2_pTrap_mistag_prob->GetZaxis()->SetTitleOffset(1.5);
  gPad->SetPhi(210); h2_pTrap_mistag_prob->Draw("lego2");
  page++;
  c1->Update();

  // now write out all the calculated parameters and associated histos
  // to the cache file
  paramFile = new TFile(paramCacheFile.c_str(), "recreate");
  if (paramFile->IsOpen()) {
    TArrayD arr;
    arr.Set(6, mistag_pars);
    paramFile->WriteObject(&arr, "mistag_pars");
    arr.Set(7, mass_pars);
    paramFile->WriteObject(&arr, "mass_pars");
    arr.Set(5, rap_pars);
    paramFile->WriteObject(&arr, "rap_pars");
    arr.Set(5, pt_pars);
    paramFile->WriteObject(&arr, "pt_pars");
    arr.Set(5, phi_cs_pars);
    paramFile->WriteObject(&arr, "phi_cs_pars");
    TArrayI iarr(1);
    iarr.AddAt(nMistagBins, 0);
    paramFile->WriteObject(&iarr, "nMistagBins");

    h2_mistagProb->Write();
    h2_mistagProb->SetDirectory(0);
    h_rap_mistag_prob->Write();
    h_rap_mistag_prob->SetDirectory(0);
    h_cos_true_mistag_prob->Write(); 
    h_cos_true_mistag_prob->SetDirectory(0);
    h_pL_mistag_prob->Write();
    h_pL_mistag_prob->SetDirectory(0);
    h2_pL_mistag_prob->Write();
    h2_pL_mistag_prob->SetDirectory(0);
    h2_cos_cs_vs_true->Write();
    h2_cos_cs_vs_true->SetDirectory(0);
    h2_pTrap_mistag_prob->Write();
    h2_pTrap_mistag_prob->SetDirectory(0);

    paramFile->Close();
  }
  delete paramFile;

  ps->Close(); 

  delete ps;
  delete c1;
}

// Initialize and return a new TFMultiD object which interfaces the 
// executive asymmetry fit functions
// caller takes ownership of the pointer (i.e. caller must delete)
TFMultiD* Zprime2muAsymmetry::getNewFitFcn(int fitType) {
  TFMultiD* f_recmd = 0;

  if (fitType == 0) {
    f_recmd = new TFMultiD("f_recmd", execAsym2D, 4, 4);
    f_recmd->SetParameters(1., .5, .9, nNormPoints);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints");
  }
  else if (fitType == 1) {
    f_recmd = new TFMultiD("f_recmd", execAsym6D, 6, 4);
    f_recmd->SetParameters(1., .5, .9, nNormPoints);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints");
  }
  else if (fitType == 2) {
    f_recmd = new TFMultiD("f_recmd", execAsym2D, 4, 4);
    f_recmd->SetParameters(1., .5, .9, nNormPoints);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints");
  }
  else if (fitType == 3) {
    f_recmd = new TFMultiD("f_recmd", execAsym6D, 6, 4);
    f_recmd->SetParameters(1., .5, .9, nNormPoints);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints");
  }
  else if (fitType == 4) {
    f_recmd = new TFMultiD("f_recmd", recAsym2D, 4, 6);
    f_recmd->SetParameters(1., .5, .9, nNormPoints, nSmearPoints, nSigma);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->FixParameter(4, nSmearPoints);
    f_recmd->FixParameter(5, nSigma);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints", 
			 "nSmearPoints", "nSigma");
  }
  else if (fitType == 5) {
    f_recmd = new TFMultiD("f_recmd", recAsym6D, 6, 6);
    f_recmd->SetParameters(1., .5, .9, nNormPoints, nSmearPoints, nSigma);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->FixParameter(4, nSmearPoints);
    f_recmd->FixParameter(5, nSigma);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints", 
			 "nSmearPoints", "nSigma");
  }
  else
    throw cms::Exception("getNewFitFcn")
      << "Invalid fit type number: " << fitType << endl;

  return f_recmd;
}

void Zprime2muAsymmetry::evalLikelihoods() {
#if 0
  // include DY bkgnd in calculating fractions
  const bool includeBgnd = true;
  // which fit function to use to get the likelihood;
  //   default is recAsym2D
  const int fitTypeLR = 5; 

  // JMTBAD pass down from GenKineAna-equiv
  const double numEvents = 1;
  const double sigmas[3][2] = {{0,0},{0,0},{0,0}};
  const double ggqqcounts[3][2] = {{0,0},{0,0},{0,0}};
  double coeffs[3][2] = {
    {0., 1.},
    {0., 0.},
    {0., 0.}
  };
  // 1.5TeV :
  //    double f_qq = 0.221, f_gg = 0.779;
  //    double sigma_grav = 1.0289E-03;
  //    double sigma_dy = 1.7382E-04;
  int MASS = 0; //(MASS_INDEX-11)%3; // 0 = 1 TeV, 1 = 3 TeV, 2 = 5 TeV
  double f_gg = ggqqcounts[MASS][0]/numEvents;
  double f_qq = ggqqcounts[MASS][1]/numEvents;
  double sigma_grav = sigmas[MASS][0];
  double sigma_dy = sigmas[MASS][1];
  double frac_grav = sigma_grav/(sigma_grav+sigma_dy);

  double x_qq, x_gg, x_b;
  if (includeBgnd) {
    x_qq = f_qq*frac_grav;
    x_gg = f_gg*frac_grav;
    x_b = 1-frac_grav;
  }
  else {
    x_qq = f_qq;
    x_gg = f_gg;
    x_b = 0.0;
  }
  //    coeffs[1][0] = -3*x_qq+x_b; // coeff of cos^2
  //    coeffs[1][1] = 4*x_qq-x_gg; // coeff of cos^4  
  double C = 8.0/(3+2*(x_qq+x_gg));
  coeffs[1][0] = 3.0*C/8.0*(x_b - 5.0*x_qq);
  coeffs[1][1] = 5.0*C/8.0*(4.0*x_qq - x_gg);
  coeffs[2][0] = 3.0/(1+frac_grav)*x_b;
  coeffs[2][1] = 0.0;

  int i_rec = fitTypeLR > 1 ? 1 : 0;
  int n2fit = nfit_used[i_rec];

  double *L[3];
  for (int lx = 0; lx < 3; lx++)
    L[lx] = new double[n2fit];
    
  TFMultiD* f_recmd = getNewFitFcn(fitTypeLR);
  UnbinnedFitter *unbfitter = new UnbinnedFitter();

  // outer loop is pdf: faster since it only has to renormalize with
  // each new pdf
  for (int pdf = 0; pdf < 3; pdf++) {
    if (pdf < 2)
      asymFitManager.setPDF(AsymFitManager::PDFTYPE(pdf));
    f_recmd->FixParameter(1, coeffs[pdf][0]);
    f_recmd->FixParameter(2, coeffs[pdf][1]);
	
    for (int cnt = 0; cnt < n2fit; cnt++) {
      // make new arrays, each for just the current point so we can
      // reuse the unbinned fit code (which expects its data to be in
      // a vector of arrays)
      double cos_theta_cs[1], rap_dil[1], pt_dil[1], phi_dil[1], 
	mass_dil[1], phi_cs[1];
      cos_theta_cs[0] = fabs(cos_theta_cs_data[i_rec][cnt]);
      rap_dil[0] = rap_dil_data[i_rec][cnt];
      pt_dil[0] = pt_dil_data[i_rec][cnt];
      phi_dil[0] = phi_dil_data[i_rec][cnt];
      mass_dil[0] = mass_dil_data[i_rec][cnt];
      phi_cs[0] = phi_cs_data[i_rec][cnt];
      
      vector<double *> fit_data;
      fit_data.push_back(cos_theta_cs);
      fit_data.push_back(rap_dil);
      // if we're doing a 6D fit, append the other 4 variables
      if (fitTypeLR == 1 || fitTypeLR == 3 || fitTypeLR == 5) {
	fit_data.push_back(pt_dil);
	fit_data.push_back(phi_dil);
	fit_data.push_back(mass_dil);
	fit_data.push_back(phi_cs);
      }
      
      double log_ML;
      int cov_status = unbfitter->unbinnedFitExec("f_recmd", "c",
						  1, fit_data, 
						  0, log_ML);
      double neg2logML = -2.0*log_ML;
      L[pdf][cnt] = neg2logML;

      if (verbosity >= VERBOSITY_TOOMUCH)
	LogTrace("evalLikelihoods")
	  << "pdftype=" << pdf << " -2.0*log_ML=" << neg2logML
	  << ", cov_status = " <<  cov_status << endl;
      if (L[pdf][cnt] > 60)
	edm::LogWarning("Zprime2muAsymmetry")
	  << "Warning in unbinned fit: L=" << L[pdf][cnt]
	  << " cos_theta_cs=" << cos_theta_cs[0]
	  << " rap_dil=" << rap_dil[0] << " pt_dil=" << pt_dil[0]
	  << " phi_dil=" << phi_dil[0] << " mass_dil=" << mass_dil[0]
	  << " phi_cs=" << phi_cs[0] << endl;
    }
  }

  const string fit_type[6] = {
    "gen+exec2D", "gen+exec6D", "rec+exec2D",
    "rec+exec6D", "rec+rec2D", "rec+rec6D"
  };

  fstream LRfile;
  string filename = "likelihoods." + outputFileBase + "." +
    fit_type[fitTypeLR] + ".txt";
  LRfile.open(filename.c_str(), ios::out);

  LRfile << "#n2fit=" << n2fit << endl;
  if (includeBgnd)
    LRfile << "#DY background included\n";
  LRfile << "#x_qq=" << x_qq << " x_gg=" << x_gg << " x_b=" << x_b << endl;
  LRfile << "#spin2: a=" << coeffs[1][0] << " b=" << coeffs[1][1] << endl;
  LRfile << "#spin0: a=" << coeffs[2][0] << " b=" << coeffs[2][1] << endl;
  LRfile << "#format:\n";
  LRfile << "#spin1\tspin2\tspin0\n";
  int skipcnt = 0;
  for (int cnt = 0; cnt < n2fit; cnt++) {
    if (L[0][cnt] > 60. && L[1][cnt] > 60.) {
      skipcnt++;
      continue;
    }
    for (int lx = 0; lx < 3; lx++)
      LRfile << setprecision(10) << L[lx][cnt] << '\t';
    LRfile << setprecision(10) << mass_dil_data[i_rec][cnt] << '\t';
    LRfile << setprecision(10) << cos_theta_cs_data[i_rec][cnt] << '\t';
    LRfile << endl;
  }
  LRfile << "#skipcnt=" << skipcnt << endl;
  LRfile.close();

  for (int lx = 0; lx < 3; lx++)
    delete L[lx];
  delete f_recmd;
  delete unbfitter;
#endif
}

void Zprime2muAsymmetry::fitAsymmetry() {
  // Do the asymmetry fit using different combinations:
  // generated/reconstructed data, smearing on/off

  string filename = "recAsymFit." + string(onPeak ? "on" : "off") +
    "_peak." + outputFileBase + ".txt";
  fstream outfile;
  outfile.open(filename.c_str(), ios::out);

  outfile << " Asym Fit Results" << endl;
  outfile << "---------------------------------------------\n";
  outfile << "Parameterization " << (internalBremOn ? "after" : "before")
	  << " internal bremsstrahlung.\n";
  outfile << "Signal window is (" << asymFitManager.fit_win(0) << ", "
	  << asymFitManager.fit_win(1) << ")\n";
  outfile << "Acceptance in diRapAccept: "
	  << MUM_ETA_LIM[0] << " < eta mu- < " << MUM_ETA_LIM[1] << "; "
	  << MUP_ETA_LIM[0] << " < eta mu+ < " << MUP_ETA_LIM[1] << endl;
  outfile << "Mistag correction in asym2D/6D is " 
	  << (correctMistags ? "on" : "off") << endl;
  outfile << "Parameterization from " << maxParamEvents << " events from "
	  << genSampleFiles.size() << " files beginning with "
	  << genSampleFiles[0] << endl;

  if (useMistagHist)
    outfile << "2D histogram used for mistag probability\n";
  else {
    outfile << endl;
    outfile << setw(15) << "mistag_pars = ";
    for (int i = 0; i < 6; i++)
      outfile << setw(9) << setprecision(3) << mistag_pars[i] << " "; 
  }
  outfile << endl << setw(15) << "rap_pars = ";
  for (int i = 0; i < 5; i++)
    outfile << setw(9) << setprecision(3) << rap_pars[i] << " "; 
  outfile << endl << setw(15) << "pt_pars = ";
  for (int i = 0; i < 5; i++)
    outfile << setw(9) << setprecision(3) << pt_pars[i] << " "; 
  outfile << endl << setw(15) << "phi_cs_pars = ";
  for (int i = 0; i < 5; i++)
    outfile << setw(9) << setprecision(3) << phi_cs_pars[i] << " "; 
  outfile << endl << setw(15) << "mass_pars = ";
  for (int i = 0; i < 7; i++)
    outfile << setw(9) << setprecision(3) << mass_pars[i] << " "; 
  outfile << endl << setw(15) << "sigma = ";
  for (int i = 0; i < 6; i++)
    outfile << setw(9) << setprecision(3) 
	    << asymFitManager.rec_sigma(i) << " ";

  outfile << "\n\n" << nSmearPoints << " points used in smear\n";
  outfile << nSigma << " sigma used in smear\n";
  outfile << nNormPoints << " points used in normalization\n\n";

  // Count forward and backward generated events, f&b smeared generated events,
  // and f&b generated events in signal area.
  double b = h_b_mass->Integral(), f = h_f_mass->Integral(),
    bs = h_b_smass->Integral(), fs = h_f_smass->Integral(),
    f_sig = h_gen_sig[0]->Integral(), b_sig = h_gen_sig[1]->Integral();

  // Calculate Asymmetry for all, smeared and signal
  double asym, sigma, asyms, sigmas, asym_sig, sigma_sig;
  calcAsymmetry(f, b, asym, sigma);
  calcAsymmetry(fs, bs, asyms, sigmas);
  calcAsymmetry(f_sig, b_sig, asym_sig, sigma_sig);

  outfile << "countSmearAsym" << endl
	  << "  Forward/Backward events in signal = " << f_sig 
	  << "/" << b_sig << endl
	  << "  Forward events before/after mass smear = " << f << "/" << fs 
	  << "\n  Backward events before/after mass smear = " << b << "/" << bs
	  << "\n  Asymmetry in signal = " << asym_sig << "+/-" 
	  << sigma_sig << endl
	  << "  Asymmetry in reconstructed window before smear = " << endl
	  << "\t" << asym << "+/-" << sigma << endl
	  << "  Asymmetry in reconstructed window after smear = " 
	  << asyms << "+/-" << sigmas << endl;

  //Get asymmetry for generated data (no trigger cuts) using counting method
  double genAfb, e_genAfb;
  calcAsymmetry(h_genCosNoCut, genAfb, e_genAfb);
  outfile << "A_FB from counting generated data = " << genAfb 
	  << " +/- " << e_genAfb << endl;

  // and also from a 1D fit
  TF1* f_gen = new TF1("f_gen", asym_3_PDF, -1., 1., 3);
  int nEntries = int(h_genCosNoCut->GetEntries());
  int nBins = h_genCosNoCut->GetNbinsX();
  f_gen->SetParameters(1., .6, .9);
  f_gen->FixParameter(0, 2.*nEntries/nBins);
  f_gen->SetParNames("Norm", "A_fb", "b");
  h_genCosNoCut->Fit(f_gen, "ILNVR");
  
  double par, eplus, eminus, eparab, globcc;
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();
  outfile << "From fit to gen cos_cs histo: \n";
  par = fitter->GetParameter(1);
  fitter->GetErrors(1, eplus, eminus, eparab, globcc);
  outfile << "A_FB = " << setw(8) << setprecision(3) << par 
	  << " +/- " << eparab << endl;
  par = fitter->GetParameter(2);
  fitter->GetErrors(2, eplus, eminus, eparab, globcc);
  outfile << "   b = " << setw(8) << setprecision(3) << par 
	  << " +/- " << eparab << endl << endl;

  // Get asymmetry for generated data with acceptance cuts (same as above?)
  calcAsymmetry(h_cos_theta_cs_acc, genAfb, e_genAfb);
  outfile << "A_FB from counting generated data with acceptance cuts = "
	  << genAfb << " +/- " << e_genAfb << endl << endl;

  // also the reconstructed 1D afb
  double recAfb, e_recAfb;
  calcAsymmetry(h_cos_theta_cs_rec, recAfb, e_recAfb);
  outfile << "A_FB from counting reconstructed data = "
	  << recAfb << " +/- " << e_recAfb << endl << endl;

  nEntries = int(h_cos_theta_cs_rec->GetEntries());
  nBins = h_cos_theta_cs_rec->GetNbinsX();
  f_gen->SetParameters(1., .6, .9);
  f_gen->FixParameter(0, 2.*nEntries/nBins);
  h_cos_theta_cs_rec->Fit(f_gen, "ILNVR");
  
  fitter = TVirtualFitter::GetFitter();
  outfile << "From fit to rec cos_cs histo: \n";
  par = fitter->GetParameter(1);
  fitter->GetErrors(1, eplus, eminus, eparab, globcc);
  outfile << "A_FB = " << setw(8) << setprecision(3) << par 
	  << " +/- " << eparab << endl;
  par = fitter->GetParameter(2);
  fitter->GetErrors(2, eplus, eminus, eparab, globcc);
  outfile << "   b = " << setw(8) << setprecision(3) << par 
	  << " +/- " << eparab << endl << endl;

  // loop over the 6 high-level fitting functions

  // labels for the different fits
  const string fit_type[6] = {
    "Gen Fit with execAsym2D", "Gen Fit with execAsym6D", 
    "Rec Fit with execAsym2D", "Rec Fit with execAsym6D",
    "Rec Fit with recAsym2D",  "Rec Fit with recAsym6D"
  };

  vector<double *> fit_data;

  for (int i = 0; i < numFits; i++) {
    TFMultiD* f_recmd = getNewFitFcn(i);

    ostringstream fit_announce;
    fit_announce << fit_type[i];
    if (fitType == GRAV_QQBAR)
      fit_announce << " (qqbar)";
    else if (fitType == GRAV_GG)
      fit_announce << " (gg)";
    fit_announce << " using " << asymFitManager.getPDFName();
    if (verbosity >= VERBOSITY_SIMPLE)
      edm::LogVerbatim("fitAsymmetry") << endl << fit_announce.str() << endl;

    fit_data.clear();
    int i_select = 0;
    int i_rec = 0;
    if (i == 1 || i == 3 || i == 5) { i_select = 1; }
    if (i > 1) { i_rec = 1; }
    int n2fit = nfit_used[i_rec];

    // turn off some of the events depending on whether they came
    // from q/qbar or gg collisions to fit the pdfs seperately
    double* weights = 0;
    if (fitType == GRAV_QQBAR)
      weights = qqbar_weights[i_rec];
    else if (fitType == GRAV_GG)
      weights = gg_weights[i_rec];

    if (i_select == 0) {
      fit_data.push_back(cos_theta_cs_data[i_rec]);
      fit_data.push_back(rap_dil_data[i_rec]);
    }
    else {
      fit_data.push_back(cos_theta_cs_data[i_rec]);
      fit_data.push_back(rap_dil_data[i_rec]);
      fit_data.push_back(pt_dil_data[i_rec]);
      fit_data.push_back(phi_dil_data[i_rec]);
      fit_data.push_back(mass_dil_data[i_rec]);
      fit_data.push_back(phi_cs_data[i_rec]);
    } 

    // Do the fit
    UnbinnedFitter *unbfitter = new UnbinnedFitter();

    // fix the "additional normalization" constant
    f_recmd->FixParameter(0,  1.);

    // either the x (A_fb) term for the asymmetry fits
    // or the x^2 term for the graviton fits
    // can be negative, so set the appropriate limit
    f_recmd->SetParLimits(1, -5., 5.);

    // the x^4 term for gravitons can be negative (gg)
    // or positive (qqbar) depending on the relative
    // fraction of each, so don't clamp that term
    // positive when doing gravitons
    if (fitType < GRAV) {
      asymFitManager.setPDF(AsymFitManager::ASYM);
      f_recmd->SetParLimits(2,  0., 5.);
    }
    else if (fitType == GRAV_THEORY) {
      asymFitManager.setPDF(AsymFitManager::GRAVTH);
      // now parameter 2 is unused; the only free parameter
      // in the fit is "s", the relative ratio of the qqbar
      // contribution to the gg one
      f_recmd->FixParameter(2, 1.); 
    }
    else
      asymFitManager.setPDF(AsymFitManager::GRAV);
      
    double logML;
    int cov_status = unbfitter->unbinnedFitExec("f_recmd", "", n2fit, fit_data, weights, logML);
      
    double* covmat = unbfitter->getFitter()->GetCovarianceMatrix();
    int npar = ((TFMultiD*)gROOT->GetFunction("f_recmd"))->GetNpar();

    ostringstream covMatrix;
    // minuit reorganizes the order of the parameters to put the unfixed
    // at the end, so A_fb and b are now the first two pars (apparently)
    for (int mi = 0; mi < 2; mi++) {
      for (int mj = 0; mj < 2; mj++)
	covMatrix << setw(10) << covmat[mi+npar*mj] << "\t";
      covMatrix << endl;
    }

    // correlation coefficient between A_fb and b
    double rho = covmat[1]/sqrt(covmat[0]*covmat[1+npar]);
    
    // Save results to external file
    outfile << fit_announce.str() << endl
	    << n2fit << " events in fit" << endl;
    if (weights) {
      double wsum = 0;
      for (int wndx = 0; wndx < n2fit; wndx++) wsum += weights[wndx];
      outfile << "Using weights, sum of which is " << wsum << endl;
    }
    fitter = TVirtualFitter::GetFitter();
    par = fitter->GetParameter(1);
    fitter->GetErrors(1, eplus, eminus, eparab, globcc);
    outfile << "A_FB = " << setw(8) << setprecision(3) << par 
	    << " +/- " << eparab << endl;
    par = fitter->GetParameter(2);
    fitter->GetErrors(2, eplus, eminus, eparab, globcc);
    outfile << "   b = " << setw(8) << setprecision(3) << par 
	    << " +/- " << eparab << endl
	    << "-2*log_ML = " << setprecision(6) << -2.*logML << endl
	    << "covariance matrix:\n" << covMatrix.str()
	    << "correlation coefficient: " << rho << endl
	    << "covariance matrix status: " << cov_status << "\n\n"; 
    outfile.flush();

    delete f_recmd;
    delete unbfitter;
  }

  delete f_gen;

  outfile.close();
}

void Zprime2muAsymmetry::drawFitHistos() {
  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 700);
  string filename = "fitHistos." + outputFileBase + ".ps";
  TPostScript *ps = new TPostScript(filename.c_str(), 111);

  const int NUM_PAGES = 20;
  TPad *pad[NUM_PAGES];
  for (int i_page = 0; i_page <= NUM_PAGES; i_page++)
    pad[i_page] = new TPad("","", .05, .05, .95, .93);

  int page = 0;
  ostringstream strpage;
  string tit;
  TPaveLabel *title;
  TText t;

  ps->NewPage();
  c1->Clear(); 
  c1->cd(0); 
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "Generated cos #theta_{CS}, (no trigger cuts)");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(1,2);
  gStyle->SetOptStat(111111);
  pad[page]->cd(1); 
  h_genCosNoCut->SetStats(false);
  h_genCosNoCut->Draw();
  pad[page]->cd(2);
  TH1F* h_genCopy = (TH1F*) h_genCosNoCut->Clone("h_genCopy");
  h_genCopy->SetStats(true);
  TF1* f_genAsym = new TF1("f_genAsym", asym_3_PDF, -1., 1., 3);
  int nEntries = int(h_genCosNoCut->GetEntries());
  int nBins = int(h_genCosNoCut->GetNbinsX());
  f_genAsym->SetParameters(1., .6, .9);
  f_genAsym->FixParameter(0, 2.*nEntries/nBins);
  f_genAsym->SetParNames("Norm", "A_fb", "b");
  h_genCopy->Fit(f_genAsym, "ILQR");
  c1->Update();
  delete f_genAsym;

  ps->NewPage();
  c1->Clear(); 
  c1->cd(0); 
  title = new TPaveLabel(0.1,0.94,0.9,0.98,"Generated PDFs of quantities");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  for (int i = 0; i < 6; i++) {
    pad[page]->cd(i+1);  
    if (i == 0) gPad->SetLogy(1);
    AsymFitHistoGen[i]->Draw(); 
  }
  c1->Update();

  ps->NewPage();
  c1->Clear(); 
  c1->cd(0); 
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "Generated PDFs of quantities, smeared");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  for (int i = 0; i < 6; i++) {
    pad[page]->cd(i+1);  
    if (i == 0) gPad->SetLogy(1);
    AsymFitHistoGenSmeared[i]->Draw(); 
  }
  c1->Update();

  ps->NewPage();
  c1->Clear(); 
  c1->cd(0); 
  title = new TPaveLabel(0.1,0.94,0.9,0.98,"Reconstructed PDFs of quantities");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  for (int i = 0; i < 6; i++) {
    pad[page]->cd(i+1);
    if (i == 0) gPad->SetLogy(1);
    AsymFitHistoRec[i]->Draw(); 
  }
  c1->Update();

  // draw the by type pages, but only bother to do so if we're doing
  // gravitons (Z' all come from qqbar, of course)
  if (doingGravFit) {
    ps->NewPage();
    c1->Clear(); 
    c1->cd(0); 
    title = new TPaveLabel(0.1,0.94,0.9,0.98,"Generated PDFs of quantities (qqbar)");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    for (int i = 0; i < 6; i++) {
      pad[page]->cd(i+1);  
      if (i == 0) gPad->SetLogy(1);
      AsymFitHistoGenByType[0][i]->Draw(); 
    }
    c1->Update();
  
    ps->NewPage();
    c1->Clear(); 
    c1->cd(0); 
    title = new TPaveLabel(0.1,0.94,0.9,0.98,"Reconstructed PDFs of quantities (qqbar)");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    for (int i = 0; i < 6; i++) {
      pad[page]->cd(i+1);  
      if (i == 0) gPad->SetLogy(1);
      AsymFitHistoRecByType[0][i]->Draw(); 
    }
    c1->Update();

    ps->NewPage();
    c1->Clear(); 
    c1->cd(0); 
    title = new TPaveLabel(0.1,0.94,0.9,0.98,"Generated PDFs of quantities (gg)");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    for (int i = 0; i < 6; i++) {
      pad[page]->cd(i+1);  
      if (i == 0) gPad->SetLogy(1);
      AsymFitHistoGenByType[1][i]->Draw(); 
    }
    c1->Update();
  
    ps->NewPage();
    c1->Clear(); 
    c1->cd(0); 
    title = new TPaveLabel(0.1,0.94,0.9,0.98,"Reconstructed PDFs of quantities (gg)");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    for (int i = 0; i < 6; i++) {
      pad[page]->cd(i+1);  
      if (i == 0) gPad->SetLogy(1);
      AsymFitHistoRecByType[1][i]->Draw(); 
    }
    c1->Update();

    ps->NewPage();
    c1->Clear(); 
    c1->cd(0); 
    title = new TPaveLabel(0.1,0.94,0.9,0.98,"Generated PDFs of quantities by type");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(1,2);
    pad[page]->cd(1);
    AsymFitHistoGenByType[0][4]->Draw(); 
    pad[page]->cd(2);
    AsymFitHistoGenByType[1][4]->Draw(); 
    c1->Update();

    ps->NewPage();
    c1->Clear(); 
    c1->cd(0); 
    title = new TPaveLabel(0.1,0.94,0.9,0.98,"Reconstructed PDFs of quantities by type");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(1,2);
    pad[page]->cd(1);
    AsymFitHistoRecByType[0][4]->Draw(); 
    pad[page]->cd(2);
    AsymFitHistoRecByType[1][4]->Draw(); 
    c1->Update();

    ps->NewPage();
    c1->Clear(); 
    c1->cd(0); 
    title = new TPaveLabel(0.1,0.94,0.9,0.98,"Generated rapidity (by type)");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(1,3);
    pad[page]->cd(1);
    AsymFitHistoGenByType[0][1]->Draw();
    pad[page]->cd(2);
    AsymFitHistoGenByType[1][1]->Draw();
    pad[page]->cd(3);
    TH1F *hist = (TH1F*)AsymFitHistoGenByType[0][1]->Clone("hist");
    hist->SetTitle("y_qqbar/y_gg");
    hist->Divide(AsymFitHistoGenByType[1][1]);
    hist->Draw();
    c1->Update();
  }

  ps->NewPage();
  c1->Clear(); 
  c1->cd(0); 
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "Resolution plots for asymmetry fit");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  for (int i = 0; i < 3; i++) {
    pad[page]->cd(2*i+1);  
    AsymFitSmearHisto[i]->Draw();
    pad[page]->cd(2*i+2);  
    AsymFitSmearHistoDif[i]->Draw();  AsymFitSmearHistoDif[i]->Fit("gaus","Q");
  }
  c1->Update();

  ps->NewPage();
  c1->Clear(); 
  c1->cd(0); 
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "Resolution plots for asymmetry fit");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  for (int i = 0; i < 3; i++) {
    pad[page]->cd(2*i+1);  
    AsymFitSmearHisto[3+i]->Draw();
    pad[page]->cd(2*i+2);  
    AsymFitSmearHistoDif[3+i]->Draw();  
    AsymFitSmearHistoDif[3+i]->Fit("gaus","Q");
  }
  c1->Update();

  c1->SetTicks();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "cos #theta_{CS} vs rapidity, no acceptance cut");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  gPad->SetPhi(120); h2_rap_cos_d_uncut[0]->Draw("lego2");
  pad[page]->cd(2); 
  gPad->SetPhi(120); h2_rap_cos_d_uncut[1]->Draw("lego2");
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "cos #theta_{CS} vs rapidity, acceptance cut applied");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  gPad->SetPhi(120); h2_rap_cos_d[0]->Draw("lego2");
  pad[page]->cd(2); 
  gPad->SetPhi(120); h2_rap_cos_d[1]->Draw("lego2");
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "cos #theta_{CS} vs rapidity, acceptance cut effects");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(1, 2);
  TH1F* hdiff[2];
  for (int j = 0; j < 2; j++) {
    hdiff[j] = (TH1F*)h2_rap_cos_d_uncut[j]->Clone();
    hdiff[j]->SetTitle("rap vs cos dil, all - accepted");
    hdiff[j]->Add(h2_rap_cos_d[j], -1);
    pad[page]->cd(j+1); 
    gPad->SetPhi(120); hdiff[j]->Draw("lego2");
  }
  c1->Update();
  for (int j = 0; j < 2; j++)
    delete hdiff[j];

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "cos #theta_{CS} vs rapidity, reconstructed events");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(1, 2);
  pad[page]->cd(1);
  gPad->SetPhi(120);
  h2_rap_cos_d_rec->Draw("lego2");
  c1->Update();
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "Mass Forward/Backward Events with Smear");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(2, 3);
  pad[page]->cd(1);  h_b_mass->Draw();
  pad[page]->cd(2);  h_b_smass->Draw();
  pad[page]->cd(3);  h_f_mass->Draw();
  pad[page]->cd(4);  h_f_smass->Draw();
  TH1F *h_smass = (TH1F*)h_b_smass->Clone();
  TH1F *h_mass = (TH1F*)h_b_mass->Clone();
  h_smass->Add(h_f_smass);  h_mass->Add(h_f_mass);
  pad[page]->cd(5);  h_mass->Draw();
  h_mass->SetTitle("Mass of All Events");
  pad[page]->cd(6);  h_smass->Draw();
  h_smass->SetTitle("Smeared Mass of All Events");
  c1->Update();

  delete h_smass;
  delete h_mass;

  /*
  if (fitSmearedData) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title = new TPaveLabel(0.1,0.94,0.9,0.98, "Smeared Events");
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->cd(0);
    pad[page]->Divide(2, 3);
    pad[page]->cd(1); h_data[3]->Draw();
    pad[page]->cd(3); h_data[1]->Draw();
    pad[page]->cd(5); h_data[4]->Draw();
    for (Int_t i = 0; i < 3; i++){
      pad[page]->cd(2*i+2);
      h_smear_data[i]->Draw();
    }
    c1->Update();
  }
  */

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(1111);
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "cos#theta (no acceptance required)");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  TF1* f_cos = new TF1("f_cos", asym_3_PDF, -1., 1., 3);
  f_cos->SetParameters(1., .5, .9);
  f_cos->SetLineWidth(2);
  //f_cos->FixParameter(0,  1.);
  if (doingGravFit) {
    f_cos->SetParLimits(1, -5., 5.);
    f_cos->SetParLimits(2, -5., 5.);
  }
  else {
    f_cos->SetParLimits(1, -2., 2.);
    f_cos->SetParLimits(2,  0., 5.);
  }
  if (fixbIn1DFit)
    f_cos->FixParameter(2, 1);
  f_cos->SetParNames("Norm","AFB","b");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  gStyle->SetStatX(.5);
  gStyle->SetStatY(.85);
  gStyle->SetTitleY(1.);
  // figure out a common scale
  double m,max = 0;
  if ((m = h_cos_theta_true->GetMaximum()) > max) max = m;
  if ((m = h_cos_theta_cs->GetMaximum()) > max) max = m;
  if ((m = h_cos_theta_cs_fixed->GetMaximum()) > max) max = m;
  if ((m = h_cos_theta_cs_acc->GetMaximum()) > max) max = m;
  h_cos_theta_true->SetMaximum(1.1*max);
  h_cos_theta_true->Fit(f_cos, "ILVER","", -1., 1.);
  h_cos_theta_true->GetXaxis()->SetTitle("cos #theta^{*}");
  h_cos_theta_true->GetYaxis()->SetTitle("Entries");
  h_cos_theta_true->GetYaxis()->SetTitleOffset(1.6);
  TLatex tex;
  tex.DrawLatex(.7, 2700., "(a)");
  pad[page]->cd(2);  
  h_cos_theta_cs->SetMaximum(1.1*max);
  h_cos_theta_cs->Fit(f_cos, "ILVER","", -1., 1.);
  h_cos_theta_cs->GetXaxis()->SetTitle("cos #theta^{*}");
  h_cos_theta_cs->GetYaxis()->SetTitle("Entries");
  h_cos_theta_cs->GetYaxis()->SetTitleOffset(1.6);
  tex.DrawLatex(.7, 2700., "(b)");
  pad[page]->cd(3);  
  h_cos_theta_cs_fixed->SetMaximum(1.1*max);
  h_cos_theta_cs_fixed->Fit(f_cos, "ILVER","", -1., 1.);
  h_cos_theta_cs_fixed->GetXaxis()->SetTitle("cos #theta^{*}");
  h_cos_theta_cs_fixed->GetYaxis()->SetTitle("Entries");
  h_cos_theta_cs_fixed->GetYaxis()->SetTitleOffset(1.6);
  tex.DrawLatex(.7, 2700., "(c)");
  pad[page]->cd(4);  
  h_cos_theta_cs_acc->SetMaximum(1.1*max);
  h_cos_theta_cs_acc->Fit(f_cos, "ILVER","", -1., 1.);
  h_cos_theta_cs_acc->GetXaxis()->SetTitle("cos #theta^{*}");
  h_cos_theta_cs_acc->GetYaxis()->SetTitle("Entries");
  h_cos_theta_cs_acc->GetYaxis()->SetTitleOffset(1.6);
  tex.DrawLatex(.7, 2700., "(d)");
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title = new TPaveLabel(0.1,0.94,0.9,0.98, "cos#theta_{CS} (rec)");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  //f_cos->FixParameter(0,  1.);
  h_cos_theta_cs_rec->Fit(f_cos, "ILVER", "", -.8, .8);
  c1->Update();

  delete f_cos;
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  title = new TPaveLabel(0.1,0.94,0.9,0.98,
			 "Mistag probability for accepted events");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(1, 2);
  for (int j = 0; j < 2; j++) {
    pad[page]->cd(j+1);
    mistagProbEvents[j]->Draw();
  }
  c1->Update();

  ps->Close();

  delete c1;
  delete ps;
}

void Zprime2muAsymmetry::fitCosCS(TH1F* cos_hist, int params) {
  // Function which fits histograms of x = cos(theta) to quadratic in x.
  // Inputs: pointer to histogram to be fit, and number of params (2 or 3)
  // 2-parameter fit fixes x^0 and x^2 terms to have same coefficient
  //   (as in typical S.M.)
  // 3-parameter fits adds tunable factor "b" to x^2 term.

  //Give each fit a unique name based on how many calls to this function
  static int ncalls=0; //number of calls
  ostringstream fitname;
  fitname << "asymfit" << ++ncalls;

  string title = cos_hist->GetTitle();

  //Get total contents of histo. From documentation on Integral: 
  //By default the integral is computed as the sum of bin contents in the range.
  double histo_integral = cos_hist->Integral(); //Returns sum of bins contents
  if (histo_integral == 0.) {
    edm::LogWarning("fitCosCS")
      << "no events in histogram \"" << title << '"';
    return;
  }

  //Get bin width
  double bin_width = cos_hist->GetBinWidth(1);
  if (bin_width == 0.) {
    edm::LogWarning("fitCosCS") << "bin width zero for \"" << title << '"';
    return;
  }

  //Estimate normalization of fitted function.
  //The fitted function is dN/d(bin) = dN/dx * dx/d(bin) = dN/dx * (bin width)
  //If function given to fitter is normalized to unity, 
  //and if N events on histo, then norm of fitted function should
  //be N * (bin width)
  double norm = histo_integral*bin_width;

  if (verbosity >= VERBOSITY_TOOMUCH)
    LogTrace("fitCosCS")
      << "fitHistos: " << title << " integral = " << histo_integral
      << " bin width = " << bin_width
      << " norm = " << norm;

  TF1* f1;
  if (params == 3) {
    f1 = new TF1(fitname.str().c_str(), asym_3_PDF, -1, 1, 3);
    //Add 1 to norm just to keep fitter from completely cheating
    f1->SetParameters(norm+1, .3, 1.);
    f1->SetParNames("Norm","AFB","b");
  }
  else { //default to 2
    f1 = new TF1(fitname.str().c_str(), asym_2_PDF, -1, 1, 2);
    f1->SetParameters(norm+1, .3);
    f1->SetParNames("Norm","AFB");
  }

  //I = integrate function over bin instead of using center
  //L = use likelihood rather than chi-sqaure
  //V = verbose (can delete or change to Q)
  //E = error estimation using MINOS
  //R = use range specified in function range
  cos_hist->Fit(fitname.str().c_str(), "ILER");

  //Check for consistency: compare integral of fitted function
  // to sum of contents of histo.
  double function_integral = f1->Integral(-1.,1.);
  if (verbosity >= VERBOSITY_TOOMUCH)
    LogTrace("fitCosCS")
      << "Integral of function -1 to 1: " << function_integral
      << "\nIntegral of function / bin width: " << function_integral/bin_width
      << "\nSum of bins of histo: " <<  histo_integral;
  delete f1;
}

void Zprime2muAsymmetry::drawFrameHistos() {
  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 700);
  string fname = "diffFrameAsym." + outputFileBase + ".ps";
  TPostScript *ps = new TPostScript(fname.c_str(), 111);

  const int NUM_PAGES = 80;
  TPad *pad[NUM_PAGES];
  for (int i_page = 0; i_page <= NUM_PAGES; i_page++) {
    pad[i_page] = new TPad("","", .05, .05, .95, .93);
  }

  int page = 0;
  ostringstream strpage;
  string tit;
  TPaveLabel *title;
  TText t;
  gStyle->SetOptStat(1110);

  // Draw histos for Gottfried-Jackson frame
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    if (i_rec == 1 || i_rec == 2) continue; // skip L1 and L2
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = recLevelHelper.levelName(i_rec) + " Gottfried-Jackson Frame";
    title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,2);
    pad[page]->cd(1); cosGJ[i_rec][0]->Draw(); fitCosCS(cosGJ[i_rec][0],3);
    pad[page]->cd(2); cosGJ[i_rec][1]->Draw(); fitCosCS(cosGJ[i_rec][1],3);
    pad[page]->cd(3); AMassGJ[i_rec][0]->Draw();
    pad[page]->cd(4); AMassGJ[i_rec][1]->Draw();
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,2);
    pad[page]->cd(1);  ARapGJ[i_rec][0]->Draw();
    pad[page]->cd(2);  ARapGJ[i_rec][1]->Draw();
    pad[page]->cd(3);  FPseudGJ[i_rec]->Draw();
    pad[page]->cd(4);  BPseudGJ[i_rec]->Draw();
    c1->Update();
  }

  // Collins-Soper frame
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    if (i_rec == 1 || i_rec == 2) continue; // skip L1 and L2
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = recLevelHelper.levelName(i_rec) + " Collins-Soper Frame";
    title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    pad[page]->cd(1); cosCS[i_rec][0]->Draw(); fitCosCS(cosCS[i_rec][0],3);
    pad[page]->cd(3); AMassCS[i_rec][0]->Draw();
    pad[page]->cd(4); ARapCS[i_rec][0]->Draw();
    pad[page]->cd(5); FPseudCS[i_rec]->Draw();
    pad[page]->cd(6); BPseudCS[i_rec]->Draw();
    c1->Update();
  }

  // Collins-Soper frame, special scatter plots and fits
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  tit = " Collins-Soper Frame Y vs cos, Gen and L3";
  title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(1,2);
  pad[page]->cd(1); rap_vs_cosCS[0]->Draw();
  pad[page]->cd(2); rap_vs_cosCS[3]->Draw();
  c1->Update();

  // Resolution of cos(theta) CS
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  tit = " Cos theta CS Resolution";
  title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,2);
  // pad[page]->cd(1);  rap3_vs_rap0->Draw();
  pad[page]->cd(1);  cosCSRes[0]->Draw();  cosCSRes[0]->Fit("gaus","Q");
  pad[page]->cd(2);  cosCSRes[1]->Draw();  cosCSRes[1]->Fit("gaus","Q");
  pad[page]->cd(3);  cosCSRes[2]->Draw();  cosCSRes[2]->Fit("gaus","Q");
  Stat_t f_bin;
  int nbins = cosCS3_diffsq_vs_cosCS0->GetNbinsX();
  TH1F* cosCS3_diffsq_sqrt
    = fs->make<TH1F>("cosCS3_diffsq_sqrt",
		     "Sqrt(Var(L3-Gen cos theta CS)) vs Gen cos theta CS", nbins,
		     cosCS3_diffsq_vs_cosCS0->GetXaxis()->GetXmin(),
		     cosCS3_diffsq_vs_cosCS0->GetXaxis()->GetXmax());
  for (int ibin = 1; ibin <= nbins; ibin++) {
    f_bin = cosCS3_diffsq_vs_cosCS0->GetBinContent(ibin);
    if (f_bin > 0.) {f_bin = sqrt(f_bin);}
    else            {f_bin = 0.;}
    cosCS3_diffsq_sqrt->SetBinContent(ibin, f_bin);
  }
  pad[page]->cd(4);  cosCS3_diffsq_sqrt->Draw();
  c1->Update();
  delete cosCS3_diffsq_sqrt;

  // Quark direction is chosen as the boost direction of the dilepton system
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    if (i_rec == 1 || i_rec == 2) continue; // skip L1 and L2
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = recLevelHelper.levelName(i_rec) + " Baur-Boost Frame";
    title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    pad[page]->cd(1); cosBoost[i_rec]->Draw(); fitCosCS(cosBoost[i_rec],3);
    pad[page]->cd(3); AMassBoost[i_rec]->Draw();
    pad[page]->cd(4); ARapBoost[i_rec]->Draw();
    pad[page]->cd(5); FPseudBoost[i_rec]->Draw();
    pad[page]->cd(6); BPseudBoost[i_rec]->Draw();
    c1->Update();

    // Drawing the histos for the Asymmetry cuts
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    for (int j = 0; j < 3; j++) {
      pad[page]->cd(2*j+1);  AsymMBoostCut[i_rec][j]->Draw();
      pad[page]->cd(2*j+2);  AsymMBoostCut[i_rec][j+3]->Draw();
    }
    c1->Update();
  }

  // Wulz frame
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    if (i_rec == 1 || i_rec == 2) continue; // skip L1 and L2
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = recLevelHelper.levelName(i_rec) + " Wulz Frame";
    title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    pad[page]->cd(1);  cosW[i_rec]->Draw();  fitCosCS(cosW[i_rec],3);
    pad[page]->cd(3);  AMassW[i_rec]->Draw();
    pad[page]->cd(4);  ARapW[i_rec]->Draw();
    pad[page]->cd(5);  FPseudW[i_rec]->Draw();
    pad[page]->cd(6);  BPseudW[i_rec]->Draw();
    c1->Update();
  }

  gStyle->SetOptStat(111111);
  ps->Close();

  delete c1;
  delete ps;
}

void Zprime2muAsymmetry::deleteHistos() {
  // ParamHistos
  delete h_pt_dil;
  for (int i = 0; i < 2; i++) {
    delete h_rap_dil[i];
    delete h_mass_dil[i];
  }
  delete h_phi_cs;
  delete h_rap_mistag;
  delete h_rap_nomistag;
  delete h_cos_mistag;
  delete h_cos_cs;
  delete h_cos_mistag_prob;
  delete h2_rap_cos_mistag;
  delete h2_rap_cos_nomistag;
  delete h2_rap_cos_p;
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 3; j++)
      delete h_mistag[i][j];
  for (int i = 0; i < 4; i++)
    delete h2_mistag[i];
  delete h2_mistagProb;
  delete h2_pTrap;
  delete h2_pTrap_mistag;
  delete h_rap_mistag_prob;
  delete h_cos_true_mistag_prob;
  delete h_pL_mistag_prob; 
  delete h2_pL_mistag_prob;
  delete h2_cos_cs_vs_true;
  delete h2_pTrap_mistag_prob;
}
