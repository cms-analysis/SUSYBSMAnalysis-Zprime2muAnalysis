//
// Authors: Bob Cousins, Jason Mumford, Jordan Tucker, Slava Valuev, UCLA
//

#include <fstream>

#include "TBackCompFitter.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFitData.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymmetryHelpers.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Functions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TFMultiD.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/UnbinnedFitter.h"

const int FIT_ARRAY_SIZE = 50000; // max size of arrays for unbinned fits. 

class Zprime2muAsymmetry : public edm::EDAnalyzer {
public:
  explicit Zprime2muAsymmetry(const edm::ParameterSet&);

  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

private:
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

  const edm::InputTag gen_particle_src;
  const edm::InputTag gen_dilepton_src;
  const edm::InputTag dilepton_src;

  // Whether we have MC truth available.
  const bool useGen;
  // Whether to turn off the fit and only make the histograms; useful
  // to get the recSigma information from fitHistos.ps.
  const bool noFit;
  // How many of the fits to do (up to 6; from gen+gen to rec+rec).
  const int numFits;
  // Whether bremsstrahlung was switched on in the generator for the
  // parameterization sample.
  const bool internalBremOn;
  // Whether to fix b for the simple 1D fits.
  const bool fixbIn1DFit;
  // Whether to use cos_true in the 2/6D fits instead of cos_cs, when
  // MC truth is available.
  const bool useCosTrueInFit;
  // Whether to correct the cos_cs values for mistag using MC truth.
  const bool artificialCosCS;
  // Fixed parameters concerning how to do the 2D/6D fit -- how many
  // sigma to use in the convolution, and how many points to use in
  // the integrals.
  const double nSigma;
  const unsigned nNormPoints;
  const unsigned nSmearPoints;

  TH1F* h_genCosNoCut;
  TH2F* AsymFitSmearHisto[6];
  TH1F* AsymFitHistoGen[6];
  TH1F* AsymFitHistoRec[6];
  TH1F* AsymFitSmearHistoDif[6];
  TH1F* AsymFitHistoGenSmeared[6];

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

  // arrays for generating fake data from distributions
  double fake_cos_true[FIT_ARRAY_SIZE];
  double fake_cos_cs[FIT_ARRAY_SIZE];
  double fake_rap[FIT_ARRAY_SIZE];
  int fake_mistag_true[FIT_ARRAY_SIZE];
  int fake_mistag_cs[FIT_ARRAY_SIZE];
};

using namespace std;

Zprime2muAsymmetry::Zprime2muAsymmetry(const edm::ParameterSet& cfg)
  : gen_particle_src(cfg.getParameter<edm::InputTag>("gen_particle_src")),
    gen_dilepton_src(cfg.getParameter<edm::InputTag>("gen_dilepton_src")),
    dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    useGen(cfg.getParameter<bool>("useGen")),
    noFit(cfg.getParameter<bool>("noFit")),
    numFits(cfg.getParameter<int>("numFits")),
    internalBremOn(cfg.getParameter<bool>("internalBremOn")),
    fixbIn1DFit(cfg.getParameter<bool>("fixbIn1DFit")),
    useCosTrueInFit(cfg.getParameter<bool>("useCosTrueInFit")),
    artificialCosCS(cfg.getParameter<bool>("artificialCosCS")),
    nSigma(5.),
    nNormPoints(100000),
    nSmearPoints(100000)
{
  // Initialize array counters since we don't use automatically-sized
  // arrays yet.
  for (int i = 0; i < 3; i++)
    nfit_used[i] = 0;

  // Let the asymFitManager take care of setting mass_type, etc.
  asymFitManager.setConstants(cfg);
  asymFitManager.loadParametrization(cfg.getParameter<std::string>("paramCacheFile").c_str());

  // Whether to use the probabilistic mistag correction; always turn
  // it off if doing graviton fit, if we are artificially correcting
  // with MC truth, or if we are fitting to cos_true.
  if (asymFitManager.correct_mistags() && (artificialCosCS || useCosTrueInFit))
    throw cms::Exception("Zprime2muAsymmetry") << "Shouldn't be applying mistag correction when one of artificialCosCS || useCosTrueInFit is set.\n";

  ostringstream out;
  out << "------------------------------------------------------------\n"
      << "Zprime2muAsymmetry parameter summary:\n";
  if (noFit)
    out << "Not performing any fit\n";
  out << asymFitManager << endl;

  if (fixbIn1DFit) out << "Fixing b = 1.0 in 1-D fits\n";
  if (useCosTrueInFit) out << "Using cos_theta_true in 2/6-D fits\n";
  if (artificialCosCS)
    out << "Artificially correcting mistags using MC truth\n";
  
  out << "------------------------------------------------------------";
  edm::LogInfo("Zprime2muAsymmetry") << out.str();

  InitROOT();
  bookFitHistos();
}

void Zprime2muAsymmetry::analyze(const edm::Event& event, const edm::EventSetup&) {
  fillFitData(event);
}

void Zprime2muAsymmetry::endJob() {
  drawFitHistos();
  dumpFitData();
  if (!noFit)
    fitAsymmetry();
}

void Zprime2muAsymmetry::bookFitHistos() {
  const double pi = TMath::Pi(), twoPi = 2*pi;
  const AsymFitManager& amg = asymFitManager;
  edm::Service<TFileService> fs;

  h_cos_theta_true     = fs->make<TH1F>("h_cos_theta_true",     "cos #theta^{*}_{true}",                   50, -1., 1.);
  h_cos_theta_cs       = fs->make<TH1F>("h_cos_theta_cs",       "cos #theta^{*}_{CS}",                     50, -1., 1.);
  h_cos_theta_cs_acc   = fs->make<TH1F>("h_cos_theta_cs_acc",   "cos #theta^{*}_{CS} (in acceptance)",     50, -1., 1.);
  h_cos_theta_cs_fixed = fs->make<TH1F>("h_cos_theta_cs_fixed", "cos #theta^{*}_{CS} (mistags corrected)", 50, -1., 1.);
  h_cos_theta_cs_acc_fixed   = fs->make<TH1F>("h_cos_theta_cs_acc_fixed",   "cos #theta^{*}_{CS} (mistags corrected, in acceptance)",     50, -1., 1.);
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

  AsymFitSmearHistoDif[0] = fs->make<TH1F>("AsymFitSmearHistoDif1", "Dil Rec Pt - Gen Pt",     50, -50., 50.);
  AsymFitSmearHistoDif[1] = fs->make<TH1F>("AsymFitSmearHistoDif2", "Dil Rec Rap - Gen Rap",   50, -.06, .06);
  AsymFitSmearHistoDif[2] = fs->make<TH1F>("AsymFitSmearHistoDif3", "Dil Rec Phi - Gen Phi",   50, -1., 1.);
  AsymFitSmearHistoDif[3] = fs->make<TH1F>("AsymFitSmearHistoDif4", "Dil Rec Mass - Gen Mass", 50, -100.,100.);
  AsymFitSmearHistoDif[4] = fs->make<TH1F>("AsymFitSmearHistoDif5", "Rec cos_cs - Gen cos_cs", 50, -1.e-3, 1.e-3);
  AsymFitSmearHistoDif[5] = fs->make<TH1F>("AsymFitSmearHistoDif6", "Rec phi_cs - Gen phi_cs", 50, -1., 1.);

  AsymFitHistoGenSmeared[0] = fs->make<TH1F>("AsymFitHistoGenSmeared0", "Dil Smeared Gen Pt",   50,  0., amg.max_pt());
  AsymFitHistoGenSmeared[1] = fs->make<TH1F>("AsymFitHistoGenSmeared1", "Dil Smeared Gen Rap",  50, -amg.max_rap(), amg.max_rap());
  AsymFitHistoGenSmeared[2] = fs->make<TH1F>("AsymFitHistoGenSmeared2", "Dil Smeared Gen Phi",  50,  0., twoPi);
  AsymFitHistoGenSmeared[3] = fs->make<TH1F>("AsymFitHistoGenSmeared3", "Dil Smeared Gen Mass", 50, amg.fit_win(0), amg.fit_win(1));
  AsymFitHistoGenSmeared[4] = fs->make<TH1F>("AsymFitHistoGenSmeared4", "Gen Smeared cos_cs",   50, -1., 1.);
  AsymFitHistoGenSmeared[5] = fs->make<TH1F>("AsymFitHistoGenSmeared5", "Gen Smeared phi_cs",   50,  0., twoPi);
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
  double gen_cos_cs = 0., rec_cos_cs = 0., gen_phi_cs = 0., rec_phi_cs = 0.;
  
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);
  const unsigned int n_dil = dileptons->size();

  if (useGen) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    event.getByLabel(gen_particle_src, genParticles);

    edm::Handle<reco::CompositeCandidateCollection> genDileptons;
    event.getByLabel(gen_dilepton_src, genDileptons);
    const unsigned int n_gen = genDileptons->size();

    AsymFitData data;
    // JMTBAD redundant with some calculations below
    if (!computeFitQuantities(*genParticles, asymFitManager.doing_electrons(), internalBremOn, data)) {
      // if finding the Z'/etc failed, skip this event
      edm::LogWarning("fillFitData") << "could not compute fit quantities!";
      return;
    }

    // mass distributions for forward and backward events separately for
    // the entire spectrum, not just in the reconstructed window
    h_gen_sig[data.cos_true >= 0 ? 0 : 1]->Fill(data.mass);

    if (data.mass > asymFitManager.fit_win(0) && data.mass < asymFitManager.fit_win(1)) {
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
	h_cos_theta_cs_acc_fixed->Fill(data.cos_cs * (data.mistag_cs ? -1 : 1));
      
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
      const reco::CompositeCandidate& gen_dil = genDileptons->at(i_dil);
      printf("run %u event %u gen dimuon pt %f eta %f phi %f mass %f\n", event.id().run(), event.id().event(), gen_dil.pt(), gen_dil.eta(), gen_dil.phi(), gen_dil.mass());
      const reco::Candidate* gen_mum = gen_dil.daughter(0)->charge() == -1 ? gen_dil.daughter(0) : gen_dil.daughter(1);
      const reco::Candidate* gen_mup = gen_dil.daughter(0)->charge() ==  1 ? gen_dil.daughter(0) : gen_dil.daughter(1);
      printf("  mu- status %i + %i\n", gen_mum->status(), gen_mup->status());

      gen_cos_cs = calcCosThetaCSAnal(gen_mum->pz(), gen_mum->energy(), 
				      gen_mup->pz(), gen_mup->energy(), 
				      gen_dil.pt(), gen_dil.pz(), gen_dil.mass());
      gen_phi_cs = calcPhiCSAnal(gen_mum->px(), gen_mum->py(),
				 gen_mup->px(), gen_mup->py(),
				 gen_dil.pt(), gen_dil.eta(), gen_dil.phi(), gen_dil.mass());

      double phi = gen_dil.phi();
      if (phi < 0.)
	phi += 2*TMath::Pi();

      // Check to see if generated dimuons lie within acceptance cut and
      // generated window
      if (gen_mum->eta() > asymFitManager.mum_eta_lim_lo() && gen_mum->eta() < asymFitManager.mum_eta_lim_hi() &&
	  gen_mup->eta() > asymFitManager.mup_eta_lim_lo() && gen_mup->eta() < asymFitManager.mup_eta_lim_hi()) {
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
	  nfit_used[0]++;
	}

	// if number of generated and reconstructed dimuons are the same,
	// fill histos with resolutions between reconstructed and generated
	// dimuon information.  This will be used for obtaining sigmas used
	// in convolutions for asymmetry fits.
	if (n_dil == n_gen) {      
	  const pat::CompositeCandidate& rec_dil = dileptons->at(i_dil);
	  const reco::CandidateBaseRef& rec_mum = dileptonDaughterByCharge(rec_dil, -1);
	  const reco::CandidateBaseRef& rec_mup = dileptonDaughterByCharge(rec_dil, +1);
	  rec_cos_cs = calcCosThetaCSAnal(rec_mum->pz(), rec_mum->energy(), 
					  rec_mup->pz(), rec_mup->energy(), 
					  rec_dil.pt(), rec_dil.pz(), rec_dil.mass());
	  rec_phi_cs = calcPhiCSAnal(rec_mum->px(), rec_mum->py(),
				     rec_mup->px(), rec_mup->py(),
				     rec_dil.pt(), rec_dil.eta(), rec_dil.phi(), rec_dil.mass());
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
					gen_dil.pt(), gen_dil.pz(), gen_dil.mass());
      
	// Use this histogram for value of cos_true to be compared with 
	// with reconstructed.
	h_genCosNoCut->Fill(gen_cos_cs);
      }
    }
  }

  // Now loop over all reconstructed dimuons and store those to be fitted
  // (which lie inside desired reconstructed window).
  for (unsigned int i_dil = 0; i_dil < n_dil; i_dil++) {
    const pat::CompositeCandidate& rec_dil = dileptons->at(i_dil);
    if (rec_dil.mass() > asymFitManager.fit_win(0) && rec_dil.mass() < asymFitManager.fit_win(1)) {
      if (nfit_used[0] == FIT_ARRAY_SIZE - 1)
	throw cms::Exception("Zprime2muAsymmetry") << "data arrays not large enough! nfit_used[1]=" << nfit_used[1] << " FIT_ARRAY_SIZE=" << FIT_ARRAY_SIZE << endl;
      
      const reco::CandidateBaseRef& rec_mum = dileptonDaughterByCharge(rec_dil, -1);
      const reco::CandidateBaseRef& rec_mup = dileptonDaughterByCharge(rec_dil, +1);
      
      rec_cos_cs = calcCosThetaCSAnal(rec_mum->pz(), rec_mum->energy(), 
				      rec_mup->pz(), rec_mup->energy(), 
				      rec_dil.pt(), rec_dil.pz(), rec_dil.mass());
      rec_phi_cs = calcPhiCSAnal(rec_mum->px(), rec_mum->py(),
				 rec_mup->px(), rec_mup->py(),
				 rec_dil.pt(), rec_dil.eta(), rec_dil.phi(), rec_dil.mass());
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

      nfit_used[1]++;
    }
  }
}

void Zprime2muAsymmetry::dumpFitData() {
  double cos_cs, rapidity, mass, pT, phi, phi_cs;
  edm::Service<TFileService> fs;

  for (int i = 0; i < 2; ++i) {
    TTree* t = fs->make<TTree>(i == 0 ? "fitData_gen" : "fitData_rec", "");
    t->Branch("cos_cs", &cos_cs,   "cos_cs/D");
    t->Branch("rap",    &rapidity, "rap/D");
    t->Branch("mass",   &mass,     "mass/D");
    t->Branch("pt",     &pT,       "pt/D");
    t->Branch("phi",    &phi,      "phi/D");
    t->Branch("phi_cs", &phi_cs,   "phi_cs/D");

    for (int j = 0; j < nfit_used[i]; ++j) {
      cos_cs   = cos_theta_cs_data[i][j];
      rapidity = rap_dil_data     [i][j];
      mass     = mass_dil_data    [i][j];
      pT       = pt_dil_data      [i][j];
      phi      = phi_dil_data     [i][j];
      phi_cs   = phi_cs_data      [i][j];
      
      t->Fill();
    }
  }
}

void Zprime2muAsymmetry::countAsymmetry(const int which, double& Afb, double& eAfb) {
  double F = 0, B = 0, dil_avg = 0;
  for (int i = 0; i < nfit_used[which]; ++i) { 
    double cos = cos_theta_cs_data[which][i]; 
    double rap = rap_dil_data[which][i]; 
    double omega = asymFitManager.h_rap_mistag_prob->GetBinContent(asymFitManager.h_rap_mistag_prob->FindBin(fabs(rap)));
    if (cos > 0)
      F += 1;
    else
      B += 1;

    dil_avg += 1 - 2*omega;
  }

  double N = F+B;
  dil_avg /= N;

  Afb = (F-B)/N/dil_avg;
  eAfb = 2/N/dil_avg*sqrt(F*B/N);
}
    
void Zprime2muAsymmetry::asymmetryByEventWeighting(const int which, double& Afb, double& eAfb) {
  // Following A. Bodek, arXiv:0911.2850, section 5, rapidity bins only.
  double A1 = 0, A2 = 0, B1 = 0, B2 = 0, VA1 = 0, VA2 = 0, VB1 = 0, VB2 = 0;
  for (int i = 0; i < nfit_used[which]; ++i) {
    double cos = cos_theta_cs_data[which][i];
    double rap = rap_dil_data[which][i];

    double omega = asymFitManager.h_rap_mistag_prob->GetBinContent(asymFitManager.h_rap_mistag_prob->FindBin(fabs(rap)));
    double kB = 1-2*omega;
    double kA = kB*kB;
    
    if (cos > 0) {
      A1 += kA;
      B1 += kB;
      VA1 += kA*kA;
      VB1 += kB*kB;
    }
    else {
      A2 += kA;
      B2 += kB;
      VA2 += kA*kA;
      VB2 += kB*kB;
    }
  }

  double cross = A2*B1 + A1*B2;
  double E12 = VB1/B1/B1;
  double E22 = VB2/B2/B2;
  double A = A1 + A2;

  Afb = (B1-B2)/(A1+A2);
  eAfb = sqrt(E12+E22)*cross/A/A;
}
    
void Zprime2muAsymmetry::asymmetryByEventWeightingWithAngularInfo(const int which, double& Afb, double& eAfb, bool use_h, bool bin_x) {
  // Following A. Bodek, arXiv:0911.2850, section 7.
  double A1 = 0, A2 = 0, B1 = 0, B2 = 0, VA1 = 0, VA2 = 0, VB1 = 0, VB2 = 0;
  for (int i = 0; i < nfit_used[which]; ++i) {
    double cos = cos_theta_cs_data[which][i];
    double rap = rap_dil_data[which][i];
    double pt = pt_dil_data[which][i];
    double mass = mass_dil_data[which][i];

    double omega = asymFitManager.h_rap_mistag_prob->GetBinContent(asymFitManager.h_rap_mistag_prob->FindBin(fabs(rap)));
    double kB = 1-2*omega;
    double kA = kB*kB;
    
    double x = fabs(cos);
    if (bin_x) // use 10 bins in x.
      x = int(x*10)/10.;
    
    double h = 0;
    if (use_h) {
      h = pt/mass;
      h = h*h;
      h = 0.5 * h/(1+h) * (1 - 3*x*x);
    }

    double zd = 1 + x*x + h;
    double z1 = 0.5*x*x/zd/zd/zd;
    double z2 = 0.5*x/zd/zd;

    if (cos > 0) {
      A1 += z1*kA;
      B1 += z2*kB;
      VA1 += z1*z1*kA*kA;
      VB1 += z2*z2*kB*kB;
    }
    else {
      A2 += z1*kA;
      B2 += z2*kB;
      VA2 += z1*z1*kA*kA;
      VB2 += z2*z2*kB*kB;
    }
  }

  double cross = A2*B1 + A1*B2;
  double E12 = VB1/B1/B1;
  double E22 = VB2/B2/B2;
  double A = A1 + A2;

  Afb = 0.375*(B1-B2)/(A1+A2);
  eAfb = 0.375*sqrt(E12+E22)*cross/A/A;
}

// Initialize and return a new TFMultiD object which interfaces the 
// executive asymmetry fit functions
// caller takes ownership of the pointer (i.e. caller must delete)
TFMultiD* Zprime2muAsymmetry::getNewFitFcn(int which) {
  TFMultiD* f_recmd = 0;

  if (which == 0) {
    f_recmd = new TFMultiD("f_recmd", execAsym2D, 4, 4);
    f_recmd->SetParameters(1., .5, .9, nNormPoints);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints");
  }
  else if (which == 1) {
    f_recmd = new TFMultiD("f_recmd", execAsym6D, 6, 4);
    f_recmd->SetParameters(1., .5, .9, nNormPoints);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints");
  }
  else if (which == 2) {
    f_recmd = new TFMultiD("f_recmd", execAsym2D, 4, 4);
    f_recmd->SetParameters(1., .5, .9, nNormPoints);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints");
  }
  else if (which == 3) {
    f_recmd = new TFMultiD("f_recmd", execAsym6D, 6, 4);
    f_recmd->SetParameters(1., .5, .9, nNormPoints);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints");
  }
  else if (which == 4) {
    f_recmd = new TFMultiD("f_recmd", recAsym2D, 4, 6);
    f_recmd->SetParameters(1., .5, .9, nNormPoints, nSmearPoints, nSigma);
    f_recmd->FixParameter(3, nNormPoints);
    f_recmd->FixParameter(4, nSmearPoints);
    f_recmd->FixParameter(5, nSigma);
    f_recmd->SetParNames("Norm", "A_fb", "b", "nNormPoints", 
			 "nSmearPoints", "nSigma");
  }
  else if (which == 5) {
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
      << "Invalid fit type number: " << which << endl;

  return f_recmd;
}

void Zprime2muAsymmetry::fitAsymmetry() {
  // Do the asymmetry fit using different combinations:
  // generated/reconstructed data, smearing on/off

  std::ofstream outfile("recAsymFit.txt", std::ios::trunc);

  outfile << " Asym Fit Results" << endl;
  outfile << "---------------------------------------------\n";
  outfile << "Parameterization " << (internalBremOn ? "after" : "before")
	  << " internal bremsstrahlung.\n";
  outfile << "Signal window is (" << asymFitManager.fit_win(0) << ", "
	  << asymFitManager.fit_win(1) << ")\n";
  outfile << "Acceptance in diRapAccept: "
	  << asymFitManager.mum_eta_lim_lo() << " < eta mu- < " << asymFitManager.mum_eta_lim_hi() << "; "
	  << asymFitManager.mup_eta_lim_lo() << " < eta mu+ < " << asymFitManager.mup_eta_lim_hi() << endl;

  outfile << "AsymFitManager:\n" << asymFitManager << "\n";

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

  enum { gen, rec };
  outfile << "gen level:\n";
  countAsymmetry(gen, recAfb, e_recAfb);
  outfile << "counting with avg dilution:                                    " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeighting(gen, recAfb, e_recAfb);
  outfile << "event weighting:                                               " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeightingWithAngularInfo(gen, recAfb, e_recAfb);
  outfile << "event weighting with angular info:                             " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeightingWithAngularInfo(gen, recAfb, e_recAfb, false);
  outfile << "event weighting with angular info, don't use h(theta):         " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeightingWithAngularInfo(gen, recAfb, e_recAfb, true, false);
  outfile << "event weighting with angular info, don't bin cos_theta:        " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeightingWithAngularInfo(gen, recAfb, e_recAfb, false, false);
  outfile << "event weighting with angular info, don't use h, don't bin cos: " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  outfile << "\nrec level:\n";
  countAsymmetry(rec, recAfb, e_recAfb);
  outfile << "counting with avg dilution:                                    " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeighting(rec, recAfb, e_recAfb);
  outfile << "event weighting:                                               " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeightingWithAngularInfo(rec, recAfb, e_recAfb);
  outfile << "event weighting with angular info:                             " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeightingWithAngularInfo(rec, recAfb, e_recAfb, false);
  outfile << "event weighting with angular info, don't use h(theta):         " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeightingWithAngularInfo(rec, recAfb, e_recAfb, true, false);
  outfile << "event weighting with angular info, don't bin cos_theta:        " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  asymmetryByEventWeightingWithAngularInfo(rec, recAfb, e_recAfb, false, false);
  outfile << "event weighting with angular info, don't use h, don't bin cos: " << setprecision(3) << setw(5) << recAfb << " +/- " << e_recAfb << "\n";
  outfile << "\n";

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

    fit_data.clear();
    int i_select = 0;
    int i_rec = 0;
    if (i == 1 || i == 3 || i == 5) { i_select = 1; }
    if (i > 1) { i_rec = 1; }
    int n2fit = nfit_used[i_rec];

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

    f_recmd->FixParameter(0,  1.); // fix the "additional normalization" constant
    f_recmd->SetParLimits(1, -5., 5.); // A_FB
    f_recmd->SetParLimits(2,  0., 5.); // b
      
    double logML;
    int cov_status = unbfitter->unbinnedFitExec("f_recmd", "", n2fit, fit_data, 0, logML);
      
    TBackCompFitter* jmt_tbcf = dynamic_cast<TBackCompFitter*>(unbfitter->getFitter());
    
    // covariance matrix between afb and b, in order {V_afb,afb, V_afb,b, V_b,b }
    double cov[3] = {-999, -999, -999};
    double* covmat = unbfitter->getFitter()->GetCovarianceMatrix();

    if (covmat) {
      if (jmt_tbcf == 0) {
	int npar = ((TFMultiD*)gROOT->GetFunction("f_recmd"))->GetNpar();
	// minuit reorganizes the order of the parameters to put the unfixed
	// at the end, so A_fb and b are now the first two pars (apparently)
	cov[0] = covmat[0];
	cov[1] = covmat[1];
	cov[2] = covmat[npar+1];
      }
      else {
	// TBackCompFitter's returned covariance matrix doesn't contain
	// the fixed parameters, so it is a nfreepar^2 = 4 length array.
	cov[0] = covmat[0];
	cov[1] = covmat[1];
	cov[2] = covmat[3];
      }
    }

    // correlation coefficient between A_fb and b
    double rho = cov[1]/sqrt(cov[0]*cov[2]);
    
    // Save results to external file
    outfile << fit_type[i] << endl
	    << n2fit << " events in fit" << endl;
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
	    << "covariance matrix:\n"
	    << setw(15) << cov[0] << setw(15) << cov[1] << endl
 	    << setw(30)                       << cov[2] << endl
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
  string filename = "fitHistos.ps";
  TPostScript *ps = new TPostScript(filename.c_str(), 111);

  const int NUM_PAGES = 23;
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
  f_cos->SetParLimits(1, -2., 2.);
  f_cos->SetParLimits(2,  0., 5.);
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
  title = new TPaveLabel(0.1,0.94,0.9,0.98, "acceptance vs cos_theta_cs");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->cd(0);
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  h_cos_theta_cs_acc_fixed->Fit(f_cos, "ILVER", "", -.8, .8);
  pad[page]->cd(2);
  TH1F* hacc = (TH1F*)h_cos_theta_cs_acc_fixed->Clone("hacc");
  hacc->SetTitle("acceptance");
  hacc->SetStats(0);
  hacc->SetMinimum(0);
  hacc->SetMaximum(1.05);
  hacc->Divide(h_cos_theta_cs_fixed);
  hacc->Draw("hist");
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

DEFINE_FWK_MODULE(Zprime2muAsymmetry);
