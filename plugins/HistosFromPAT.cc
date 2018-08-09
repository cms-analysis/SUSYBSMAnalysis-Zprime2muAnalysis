// JMTBAD make this fill lepton histos from all leptons + those that
// make it into dileptons always, not just depending on flag.

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"///
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"

class Zprime2muHistosFromPAT : public edm::EDAnalyzer {
 public:
  explicit Zprime2muHistosFromPAT(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void getBSandPV(const edm::Event&);
  void fillBasicLeptonHistos(const reco::CandidateBaseRef&);
  void fillOfflineMuonHistos(const pat::Muon*);
  void fillOfflineElectronHistos(const pat::Electron*);
  void fillLeptonHistos(const reco::CandidateBaseRef&);
  void fillLeptonHistos(const edm::View<reco::Candidate>&);
  void fillLeptonHistosFromDileptons(const pat::CompositeCandidateCollection&);
  void fillDileptonHistos(const pat::CompositeCandidate&, const edm::Event&, double);
  void fillDileptonHistos(const pat::CompositeCandidateCollection&, const edm::Event&, double);
  double getSmearedMass(const pat::CompositeCandidate&, double);

  edm::InputTag lepton_src;
  edm::InputTag dilepton_src;
  const bool leptonsFromDileptons;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  const bool use_bs_and_pv;

  struct debug_tree_t {
    unsigned run;
    unsigned lumi;
    unsigned event;
    float mass;
    short id;
  };

  debug_tree_t dbg_t;
  TTree* dbg_tree;

  // mmm bare ptrs
  const reco::BeamSpot* beamspot;
  const reco::Vertex*   vertex;
    bool _usePrescaleWeight;
    int  _prescaleWeight;
    
    double eventWeight;///
    bool _useMadgraphWeight;///
    double _madgraphWeight;///
    bool _usekFactor;
    double _kFactor;
    double _kFactor_bb;
    double _kFactor_be;
    double _scaleUncert = 0.01;

  TH1F* NBeamSpot;
  TH1F* NVertices;
  TH1F* NLeptons;
  TH2F* LeptonSigns;
  TH1F* LeptonEta;
  TH1F* LeptonRap;
  TH1F* LeptonPhi;
  TH1F* LeptonPt;
  TH1F* LeptonPz;
  TH1F* LeptonP;
  TProfile* LeptonPVsEta;
  TProfile* LeptonPtVsEta;
  TH1F* IsoSumPt;
  TH1F* RelIsoSumPt;
  TH1F* IsoEcal;   
  TH1F* IsoHcal;   
  TH1F* CombIso;
  TH1F* RelCombIso;
  TH1F* CombIsoNoECAL;
  TH1F* RelCombIsoNoECAL;
  TH1F* IsoNTracks;
  TH1F* IsoNJets;  
  TH1F* NPxHits;
  TH1F* NStHits;
  TH1F* NTkHits;
  TH1F* NMuHits;
  TH1F* NHits;
  TH1F* NInvalidHits;
  TH1F* NPxLayers;
  TH1F* NStLayers;
  TH1F* NTkLayers;
  TH1F* Chi2dof;
  TH1F* TrackD0BS;
  TH1F* TrackDZBS;
  TH1F* TrackD0PV;
  TH1F* TrackDZPV;
  TH1F* NDileptons;
  TH1F* DileptonEta;
  TH1F* DileptonRap;
  TH1F* DileptonPhi;
  TH1F* DileptonPt;
  TH1F* DileptonPz;
  TH1F* DileptonP;
  TProfile* DileptonPVsEta;
  TProfile* DileptonPtVsEta;
  TH1F* ChiDilepton;
  TH1F* CosThetaStarDilepton;
  TH1F* DileptonMass;
  TH1F* DileptonMass_bb;
  TH1F* DileptonMass_be;
  TH1F* DileptonMass_CSPos;
  TH1F* DileptonMass_bb_CSPos;
  TH1F* DileptonMass_be_CSPos;
  TH1F* DileptonMass_CSNeg;
  TH1F* DileptonMass_bb_CSNeg;
  TH1F* DileptonMass_be_CSNeg;
  TH1F* DileptonMassWeight;
  TH1F* DileptonWithPhotonsMass;
  TH1F* DileptonDeltaPt;
  TH1F* DileptonDeltaP;
  TH2F* DimuonMuonPtErrors;
  TH1F* DimuonMuonPtErrOverPt;
  TH1F* DimuonMuonPtErrOverPtM200;
  TH1F* DimuonMuonPtErrOverPtM500;
  TH2F* DileptonDaughterIds;
  TH1F* DileptonDaughterDeltaR;
  TH1F* DileptonDaughterDeltaPhi;
  TH1F* DimuonMassVertexConstrained;
  TH1F* DimuonMassVertexConstrained_bb;
  TH1F* DimuonMassVertexConstrained_be;
  TH1F* DimuonMassVertexConstrained_CSPos;
  TH1F* DimuonMassVertexConstrained_bb_CSPos;
  TH1F* DimuonMassVertexConstrained_be_CSPos;
  TH1F* DimuonMassVertexConstrained_CSNeg;
  TH1F* DimuonMassVertexConstrained_bb_CSNeg;
  TH1F* DimuonMassVertexConstrained_be_CSNeg;
  TH1F* DimuonMassVertexConstrainedSmear;
  TH1F* DimuonMassVertexConstrainedSmear_bb;
  TH1F* DimuonMassVertexConstrainedSmear_be;
  TH1F* DimuonMassVertexConstrainedSmear_CSPos;
  TH1F* DimuonMassVertexConstrainedSmear_bb_CSPos;
  TH1F* DimuonMassVertexConstrainedSmear_be_CSPos;
  TH1F* DimuonMassVertexConstrainedSmear_CSNeg;
  TH1F* DimuonMassVertexConstrainedSmear_bb_CSNeg;
  TH1F* DimuonMassVertexConstrainedSmear_be_CSNeg;
  TH1F* DimuonMassVertexConstrainedScaleUp;
  TH1F* DimuonMassVertexConstrainedScaleUp_bb;
  TH1F* DimuonMassVertexConstrainedScaleUp_be;
  TH1F* DimuonMassVertexConstrainedScaleUp_CSPos;
  TH1F* DimuonMassVertexConstrainedScaleUp_bb_CSPos;
  TH1F* DimuonMassVertexConstrainedScaleUp_be_CSPos;
  TH1F* DimuonMassVertexConstrainedScaleUp_CSNeg;
  TH1F* DimuonMassVertexConstrainedScaleUp_bb_CSNeg;
  TH1F* DimuonMassVertexConstrainedScaleUp_be_CSNeg;
  TH1F* DimuonMassVertexConstrainedScaleDown;
  TH1F* DimuonMassVertexConstrainedScaleDown_bb;
  TH1F* DimuonMassVertexConstrainedScaleDown_be;
  TH1F* DimuonMassVertexConstrainedScaleDown_CSPos;
  TH1F* DimuonMassVertexConstrainedScaleDown_bb_CSPos;
  TH1F* DimuonMassVertexConstrainedScaleDown_be_CSPos;
  TH1F* DimuonMassVertexConstrainedScaleDown_CSNeg;
  TH1F* DimuonMassVertexConstrainedScaleDown_bb_CSNeg;
  TH1F* DimuonMassVertexConstrainedScaleDown_be_CSNeg;

  TH1F* DimuonMassVertexConstrainedWeight;
  TH1F* DimuonMassVtxConstrainedLog;
  TH1F* DimuonMassVtxConstrainedLog_bb;
  TH1F* DimuonMassVtxConstrainedLog_be;
  TH1F* DimuonMassVtxConstrainedLogWeight;
  TH2F* DimuonMassConstrainedVsUn;
  TH2F* DimuonMassVertexConstrainedError;
    //special
    TH1F* DimuonMassVtx_chi2;
    TH1F* DimuonMassVtx_prob;
    //weight
    TH1F* WeightMadGraph;///
    TH1F *kFactorGraph;
    TH1F *kFactorGraph_bb;
    TH1F *kFactorGraph_be;
    
	const bool fill_gen_info;
	HardInteraction* hardInteraction;
	
};

Zprime2muHistosFromPAT::Zprime2muHistosFromPAT(const edm::ParameterSet& cfg)
  : lepton_src(cfg.getParameter<edm::InputTag>("lepton_src")),
    dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    leptonsFromDileptons(cfg.getParameter<bool>("leptonsFromDileptons")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    dbg_tree(0),
    beamspot(0),
    vertex(0),
    _usePrescaleWeight(cfg.getUntrackedParameter<bool>("usePrescaleWeight",false)),
    _prescaleWeight(1),
    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),///
    _madgraphWeight(1.),///
    _usekFactor(cfg.getParameter<bool>("usekFactor")),
    _kFactor(1.),
    _kFactor_bb(1.),
    _kFactor_be(1.),
    fill_gen_info(cfg.existsAs<edm::ParameterSet>("hardInteraction")),
    hardInteraction(fill_gen_info ? new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")) : 0)
{

  consumes<reco::CandidateView>(lepton_src);
  consumes<pat::CompositeCandidateCollection>(dilepton_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));
  if (fill_gen_info) consumes<std::vector<reco::GenParticle>>(hardInteraction->src);
 


  std::string title_prefix = cfg.getUntrackedParameter<std::string>("titlePrefix", "");
  if (title_prefix.size() && title_prefix[title_prefix.size()-1] != ' ')
    title_prefix += " ";
  const TString titlePrefix(title_prefix.c_str());

  edm::Service<TFileService> fs;

  TH1::SetDefaultSumw2(true);

  if (cfg.getUntrackedParameter<bool>("debug", false)) {
    dbg_tree = fs->make<TTree>("t", "");
    dbg_tree->Branch("tt", &dbg_t, "run/i:lumi:event:mass/F:id/S");
  }
 
  // Whole-event things.
  NBeamSpot = fs->make<TH1F>("NBeamSpot", titlePrefix + "# beamspots/event",  2, 0,  2);
  NVertices = fs->make<TH1F>("NVertices", titlePrefix + "# vertices/event",  40, 0, 40);

  // Basic kinematics.

  // Lepton multiplicity.
  NLeptons = fs->make<TH1F>("NLeptons", titlePrefix + "# leptons/event", 10, 0, 10);

  // Opposite/like-sign counts.
  LeptonSigns = fs->make<TH2F>("LeptonSigns", titlePrefix + "lepton sign combinations", 6, 0, 6, 13, -6, 7);
  LeptonSigns->GetXaxis()->SetTitle("# leptons");
  LeptonSigns->GetYaxis()->SetTitle("total charge");

  // Lepton eta, y, phi.
  LeptonEta = fs->make<TH1F>("LeptonEta", titlePrefix + "#eta", 100, -5, 5);
  LeptonRap = fs->make<TH1F>("LeptonRap", titlePrefix + "y",    100, -5, 5);
  LeptonPhi = fs->make<TH1F>("LeptonPhi", titlePrefix + "#phi", 100, -TMath::Pi(), TMath::Pi());

  // Lepton momenta: p, p_T, p_z.
  LeptonPt = fs->make<TH1F>("LeptonPt", titlePrefix + "pT", 5000, 0, 5000);
  LeptonPz = fs->make<TH1F>("LeptonPz", titlePrefix + "pz", 5000, 0, 5000);
  LeptonP  = fs->make<TH1F>("LeptonP",  titlePrefix + "p",  5000, 0, 5000);

  // Lepton momenta versus pseudorapidity.
  LeptonPVsEta  = fs->make<TProfile>("LeptonPVsEta",   titlePrefix + "p vs. #eta",  100, -6, 6);
  LeptonPtVsEta = fs->make<TProfile>("LeptonPtVsEta",  titlePrefix + "pT vs. #eta", 100, -6, 6);

  // Muon specific histos.

  // Delta R < 0.3 isolation variables.
  IsoSumPt         = fs->make<TH1F>("IsoSumPt",         titlePrefix + "Iso. (#Delta R < 0.3) #Sigma pT",           50, 0, 50);
  RelIsoSumPt      = fs->make<TH1F>("RelIsoSumPt",      titlePrefix + "Iso. (#Delta R < 0.3) #Sigma pT / tk. pT",  50, 0, 1);
  IsoEcal          = fs->make<TH1F>("IsoEcal",          titlePrefix + "Iso. (#Delta R < 0.3) ECAL",                50, 0, 50);
  IsoHcal          = fs->make<TH1F>("IsoHcal",          titlePrefix + "Iso. (#Delta R < 0.3) HCAL",                50, 0, 50);
  CombIso          = fs->make<TH1F>("CombIso",          titlePrefix + "Iso. (#Delta R < 0.3), combined",           50, 0, 50);
  RelCombIso       = fs->make<TH1F>("RelCombIso",       titlePrefix + "Iso. (#Delta R < 0.3), combined / tk. pT",  50, 0, 1);
  CombIsoNoECAL    = fs->make<TH1F>("CombIsoNoECAL",    titlePrefix + "Iso. (#Delta R < 0.3), combined (no ECAL)", 50, 0, 50);
  RelCombIsoNoECAL = fs->make<TH1F>("RelCombIsoNoECAL", titlePrefix + "Iso. (#Delta R < 0.3), combined (no ECAL), relative", 50, 0, 50);
  IsoNTracks       = fs->make<TH1F>("IsoNTracks",       titlePrefix + "Iso. (#Delta R < 0.3) nTracks",             10, 0, 10);
  IsoNJets         = fs->make<TH1F>("IsoNJets",         titlePrefix + "Iso. (#Delta R < 0.3) nJets",               10, 0, 10);

  // Track hit counts.
  NPxHits = fs->make<TH1F>("NPxHits", titlePrefix + "# pixel hits",    8, 0,  8);
  NStHits = fs->make<TH1F>("NStHits", titlePrefix + "# strip hits",   30, 0, 30);
  NTkHits = fs->make<TH1F>("NTkHits", titlePrefix + "# tracker hits", 40, 0, 40);
  NMuHits = fs->make<TH1F>("NMuHits", titlePrefix + "# muon hits",    55, 0, 55);

  NHits        = fs->make<TH1F>("NHits",        titlePrefix + "# hits",         78, 0, 78);
  NInvalidHits = fs->make<TH1F>("NInvalidHits", titlePrefix + "# invalid hits", 78, 0, 78);

  NPxLayers = fs->make<TH1F>("NPxLayers", titlePrefix + "# pixel layers",    8, 0,  8);
  NStLayers = fs->make<TH1F>("NStLayers", titlePrefix + "# strip layers",   15, 0, 15);
  NTkLayers = fs->make<TH1F>("NTkLayers", titlePrefix + "# tracker layers", 20, 0, 20);

  // Other track variables.
  Chi2dof = fs->make<TH1F>("Chi2dof", titlePrefix + "#chi^{2}/dof", 100, 0, 10);
  TrackD0BS = fs->make<TH1F>("TrackD0BS", titlePrefix + "|d0 wrt BS|", 100, 0, 0.2);
  TrackDZBS = fs->make<TH1F>("TrackDZBS", titlePrefix + "|dz wrt BS|", 100, 0, 20);
  TrackD0PV = fs->make<TH1F>("TrackD0PV", titlePrefix + "|d0 wrt PV|", 100, 0, 0.2);
  TrackDZPV = fs->make<TH1F>("TrackDZPV", titlePrefix + "|dz wrt PV|", 100, 0, 20);

  // Electron specific histos (none yet).

  // Dilepton quantities.

  // Dilepton multiplicity.
  NDileptons = fs->make<TH1F>("NDileptons", "# dileptons/event" + titlePrefix, 10, 0, 10);

  // Dilepton eta, y, phi.
  DileptonEta = fs->make<TH1F>("DileptonEta", titlePrefix + "dil. #eta", 100, -5,  5);
  DileptonRap = fs->make<TH1F>("DileptonRap", titlePrefix + "dil. y",    100, -5,  5);
  DileptonPhi = fs->make<TH1F>("DileptonPhi", titlePrefix + "dil. #phi", 100, -TMath::Pi(), TMath::Pi());

  // Dilepton momenta: p, p_T, p_z.
  DileptonPt = fs->make<TH1F>("DileptonPt", titlePrefix + "dil. pT", 5000, 0, 5000);
  DileptonPz = fs->make<TH1F>("DileptonPz", titlePrefix + "dil. pz", 5000, 0, 5000);
  DileptonP  = fs->make<TH1F>("DileptonP",  titlePrefix + "dil. p",  5000, 0, 5000);
  
  // Dilepton momenta versus pseudorapidity.
  DileptonPVsEta  = fs->make<TProfile>("DileptonPVsEta",  titlePrefix + "dil. p vs. #eta",  100, -6, 6);
  DileptonPtVsEta = fs->make<TProfile>("DileptonPtVsEta", titlePrefix + "dil. pT vs. #eta", 100, -6, 6);
 
   // Dilepton chi a la dijet
   //
  ChiDilepton            = fs->make<TH1F>("ChiDilepton",            titlePrefix + "dil. chi", 100, 0, 20);
  CosThetaStarDilepton   = fs->make<TH1F>("CosThetaStarDilepton",            titlePrefix + "dil. cos theta star", 100, -1, 1);

  // Dilepton invariant mass.
  DileptonMass            = fs->make<TH1F>("DileptonMass",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DileptonMass_bb         = fs->make<TH1F>("DileptonMass_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DileptonMass_be         = fs->make<TH1F>("DileptonMass_be",            titlePrefix + "dil. mass barrel-endcaps and endcaps-endcaps", 20000, 0, 20000);
  DileptonMass_CSPos            = fs->make<TH1F>("DileptonMass_CSPos",            titlePrefix + "dil. mass for positive cos theta star", 20000, 0, 20000);
  DileptonMass_bb_CSPos         = fs->make<TH1F>("DileptonMass_bb_CSPos",            titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DileptonMass_be_CSPos         = fs->make<TH1F>("DileptonMass_be_CSPos",            titlePrefix + "dil. mass barrel-endcaps and endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DileptonMass_CSNeg            = fs->make<TH1F>("DileptonMass_CSNeg",            titlePrefix + "dil. mass for negative cos theta star", 20000, 0, 20000);
  DileptonMass_bb_CSNeg         = fs->make<TH1F>("DileptonMass_bb_CSNeg",            titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DileptonMass_be_CSNeg         = fs->make<TH1F>("DileptonMass_be_CSNeg",            titlePrefix + "dil. mass barrel-endcaps and endcaps-endcaps for negative cos theta star", 20000, 0, 20000);
  DileptonMassWeight      = fs->make<TH1F>("DileptonMassWeight",      titlePrefix + "dil. mass", 20000, 0, 20000);
  DileptonWithPhotonsMass = fs->make<TH1F>("DileptonWithPhotonsMass", titlePrefix + "res. mass", 20000, 0, 20000);
  
  // Plots comparing the daughter lepton momenta.
  DileptonDeltaPt = fs->make<TH1F>("DileptonDeltaPt",  titlePrefix + "dil. |pT^{1}| - |pT^{2}|", 100, -100, 100);
  DileptonDeltaP  = fs->make<TH1F>("DileptonDeltaP",   titlePrefix + "dil. |p^{1}| - |p^{2}|",   100, -500, 500);

  // pT errors of daughter muons
  DimuonMuonPtErrors        = fs->make<TH2F>("DimuonMuonPtErrors",        titlePrefix + "dil. #sigma_{pT}^{1} v. #sigma_{pT}^{2}", 100, 0, 100, 100, 0, 100);
  DimuonMuonPtErrOverPt     = fs->make<TH1F>("DimuonMuonPtErrOverPt",     titlePrefix + "muon #sigma_{pT}/pT",              200, 0., 1.);
  DimuonMuonPtErrOverPtM200 = fs->make<TH1F>("DimuonMuonPtErrOverPtM200", titlePrefix + "muon #sigma_{pT}/pT, M > 200 GeV", 200, 0., 1.);
  DimuonMuonPtErrOverPtM500 = fs->make<TH1F>("DimuonMuonPtErrOverPtM500", titlePrefix + "muon #sigma_{pT}/pT, M > 500 GeV", 200, 0., 1.);

  // More daughter-related info.
  DileptonDaughterIds = fs->make<TH2F>("DileptonDaughterIds", "", 27, -13, 14, 27, -13, 14);
  DileptonDaughterDeltaR = fs->make<TH1F>("DileptonDaughterDeltaR", "", 100, 0, 5);
  DileptonDaughterDeltaPhi = fs->make<TH1F>("DileptonDaughterDeltaPhi", "", 100, 0, 3.15);

  // Dimuons have a vertex-constrained fit: some associated histograms.
  DimuonMassVertexConstrained = fs->make<TH1F>("DimuonMassVertexConstrained", titlePrefix + "dimu. vertex-constrained mass", 20000, 0, 20000);
  DimuonMassVertexConstrained_bb = fs->make<TH1F>("DimuonMassVertexConstrained_bb", titlePrefix + "dimu. vertex-constrained mass barrel-barrel", 20000, 0, 20000);
  DimuonMassVertexConstrained_be = fs->make<TH1F>("DimuonMassVertexConstrained_be", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps", 20000, 0, 20000);
  DimuonMassVertexConstrained_CSPos = fs->make<TH1F>("DimuonMassVertexConstrained_CSPos", titlePrefix + "dimu. vertex-constrained mass for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_bb_CSPos = fs->make<TH1F>("DimuonMassVertexConstrained_bb_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_be_CSPos = fs->make<TH1F>("DimuonMassVertexConstrained_be_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrained_CSNeg", titlePrefix + "dimu. vertex-constrained mass for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_bb_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrained_bb_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_be_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrained_be_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedWeight = fs->make<TH1F>("DimuonMassVertexConstrainedWeight", titlePrefix + "dimu. vertex-constrained mass", 20000, 0, 20000);
 
  DimuonMassVertexConstrainedSmear = fs->make<TH1F>("DimuonMassVertexConstrainedSmear ", titlePrefix + "dimu. vertex-constrained mass", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_bb = fs->make<TH1F>("DimuonMassVertexConstrainedSmear _bb", titlePrefix + "dimu. vertex-constrained mass barrel-barrel", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_be = fs->make<TH1F>("DimuonMassVertexConstrainedSmear _be", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedSmear _CSPos", titlePrefix + "dimu. vertex-constrained mass for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_bb_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedSmear _bb_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_be_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedSmear _be_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedSmear _CSNeg", titlePrefix + "dimu. vertex-constrained mass for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_bb_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedSmear _bb_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_be_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedSmear _be_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for negative cos theta star", 20000, 0, 20000);
  
  DimuonMassVertexConstrainedScaleUp = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp", titlePrefix + "dimu. vertex-constrained mass", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleUp_bb = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp_bb", titlePrefix + "dimu. vertex-constrained mass barrel-barrel", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleUp_be = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp_be", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleUp_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp_CSPos", titlePrefix + "dimu. vertex-constrained mass for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleUp_bb_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp_bb_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleUp_be_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp_be_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleUp_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp_CSNeg", titlePrefix + "dimu. vertex-constrained mass for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleUp_bb_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp_bb_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleUp_be_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedScaleUp_be_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for negative cos theta star", 20000, 0, 20000);
  
  DimuonMassVertexConstrainedScaleDown = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown", titlePrefix + "dimu. vertex-constrained mass", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleDown_bb = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown_bb", titlePrefix + "dimu. vertex-constrained mass barrel-barrel", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleDown_be = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown_be", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleDown_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown_CSPos", titlePrefix + "dimu. vertex-constrained mass for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleDown_bb_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown_bb_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleDown_be_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown_be_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleDown_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown_CSNeg", titlePrefix + "dimu. vertex-constrained mass for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleDown_bb_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown_bb_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedScaleDown_be_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedScaleDown_be_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for negative cos theta star", 20000, 0, 20000);
  // Mass plot in bins of log(mass)
  const int    NMBINS = 100;
  const double MMIN = 60., MMAX = 2100.;
  double logMbins[NMBINS+1];
  for (int ibin = 0; ibin <= NMBINS; ibin++)
    logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
  DimuonMassVtxConstrainedLog = fs->make<TH1F>("DimuonMassVtxConstrainedLog", titlePrefix + "dimu vtx-constrained mass in log bins", NMBINS, logMbins);
  DimuonMassVtxConstrainedLog_bb = fs->make<TH1F>("DimuonMassVtxConstrainedLog_bb", titlePrefix + "dimu vtx-constrained mass in log bins barrel-barrel", NMBINS, logMbins);
  DimuonMassVtxConstrainedLog_be = fs->make<TH1F>("DimuonMassVtxConstrainedLog_be", titlePrefix + "dimu vtx-constrained mass in log bins barrel-endcaps", NMBINS, logMbins);
  DimuonMassVtxConstrainedLogWeight = fs->make<TH1F>("DimuonMassVtxConstrainedLogWeight", titlePrefix + "dimu vtx-constrained mass in log bins", NMBINS, logMbins);
  DimuonMassConstrainedVsUn = fs->make<TH2F>("DimuonMassConstrainedVsUn", titlePrefix + "dimu. vertex-constrained vs. non-constrained mass", 200, 0, 3000, 200, 0, 3000);
  DimuonMassVertexConstrainedError = fs->make<TH2F>("DimuonMassVertexConstrainedError", titlePrefix + "dimu. vertex-constrained mass error vs. mass", 100, 0, 3000, 100, 0, 400);
    //special
    DimuonMassVtx_chi2 = fs->make<TH1F>("DimuonMassVtx_chi2", titlePrefix + "dimu. vertex #chi^{2}/dof", 300, 0, 30);
    DimuonMassVtx_prob = fs->make<TH1F>("DimuonMassVtx_prob", titlePrefix + "dimu. vertex probability", 100, 0, 1);
    
     //weight
     WeightMadGraph = fs->make<TH1F>("weightperevent", titlePrefix + "weight per event", 4, -2,2);
     kFactorGraph = fs->make<TH1F>("kFactorperevent", titlePrefix + "kFactor per event", 50, 0.4,1.4);
     kFactorGraph_bb = fs->make<TH1F>("kFactorperevent_bb", titlePrefix + "kFactor per event bb", 50, 0.4,1.4);
     kFactorGraph_be = fs->make<TH1F>("kFactorperevent_be", titlePrefix + "kFactor per event be", 50, 0.4,1.4);
}


double Zprime2muHistosFromPAT::getSmearedMass(const pat::CompositeCandidate& dil, double gM){

    double a = 0.;
    double b = 0.;
    double c = 0.;
    double d = 0.;
    double e = 0.;

    double mass = dil.userFloat("vertexM");

   //sigma BB
   if (dil.daughter(0)->eta()<=1.2 && dil.daughter(1)->eta()<=1.2 && dil.daughter(0)->eta()>=-1.2 && dil.daughter(1)->eta()>=-1.2){

   	a=0.00701;
  	b=3.32e-05;
   	c=-1.29e-08;
   	d=2.73e-12;
   	e=-2.05e-16;
   
    }
    else{
   //sigma BE
	a=0.0124;
   	b=3.75e-05;
   	c=-1.52e-08;
   	d=3.44e-12;
   	e=-2.85e-16;

    }
    double res = a + b*gM + c*gM*gM + d*pow(gM,3) + e*pow(gM,4);

    double extraSmear = res*0.15;

	
    TRandom3 *rand = new TRandom3(0);
    return mass*rand->Gaus(1,extraSmear);


}


void Zprime2muHistosFromPAT::getBSandPV(const edm::Event& event) {
  // We store these as bare pointers. Should find better way, but
  // don't want to pass them around everywhere...
  edm::Handle<reco::BeamSpot> hbs;
  event.getByLabel(beamspot_src, hbs);
  beamspot = hbs.isValid() ? &*hbs : 0; // nice and fragile
  NBeamSpot->Fill(beamspot != 0);

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel(vertex_src, vertices);
  vertex = 0;
  int vertex_count = 0;
  for (reco::VertexCollection::const_iterator it = vertices->begin(), ite = vertices->end(); it != ite; ++it) {
    if (it->ndof() > 4 && fabs(it->z()) <= 24 && fabs(it->position().rho()) <= 2) {
      if (vertex == 0)
	vertex = &*it;
      ++vertex_count;
    }
  }
  NVertices->Fill(vertex_count, _madgraphWeight*_kFactor);
}

void Zprime2muHistosFromPAT::fillBasicLeptonHistos(const reco::CandidateBaseRef& lep) {
  LeptonEta->Fill(lep->eta(), _madgraphWeight*_kFactor);
  LeptonRap->Fill(lep->rapidity(), _madgraphWeight*_kFactor);
  LeptonPhi->Fill(lep->phi(), _madgraphWeight*_kFactor);

  LeptonPt->Fill(lep->pt(), _madgraphWeight*_kFactor);
  LeptonPz->Fill(fabs(lep->pz()), _madgraphWeight*_kFactor);
  LeptonP ->Fill(lep->p(), _madgraphWeight*_kFactor);

  LeptonPtVsEta->Fill(lep->eta(), lep->pt(), _madgraphWeight*_kFactor);
  LeptonPVsEta ->Fill(lep->eta(), lep->p(), _madgraphWeight*_kFactor);
}

void Zprime2muHistosFromPAT::fillOfflineMuonHistos(const pat::Muon* mu) {
  const reco::MuonIsolation& iso = mu->isolationR03();
  IsoSumPt   ->Fill(iso.sumPt, _madgraphWeight*_kFactor);
  RelIsoSumPt->Fill(iso.sumPt / mu->innerTrack()->pt(), _madgraphWeight*_kFactor);
  IsoEcal    ->Fill(iso.emEt, _madgraphWeight*_kFactor);
  IsoHcal    ->Fill(iso.hadEt + iso.hoEt, _madgraphWeight*_kFactor);
  CombIso    ->Fill( iso.sumPt + iso.emEt + iso.hadEt + iso.hoEt, _madgraphWeight*_kFactor);
  RelCombIso ->Fill((iso.sumPt + iso.emEt + iso.hadEt + iso.hoEt) / mu->innerTrack()->pt(), _madgraphWeight*_kFactor);
  IsoNTracks ->Fill(iso.nTracks, _madgraphWeight*_kFactor);
  IsoNJets   ->Fill(iso.nJets, _madgraphWeight*_kFactor);

  CombIsoNoECAL   ->Fill( iso.sumPt + iso.hadEt + iso.hoEt, _madgraphWeight*_kFactor);
  RelCombIsoNoECAL->Fill((iso.sumPt + iso.hadEt + iso.hoEt) / mu->innerTrack()->pt(), _madgraphWeight*_kFactor);

  const reco::TrackRef track = patmuon::getPickedTrack(*mu);
  if (track.isAvailable()) {
    Chi2dof->Fill(track->normalizedChi2(), _madgraphWeight*_kFactor);

    if (beamspot != 0) {
      TrackD0BS->Fill(fabs(track->dxy(beamspot->position())), _madgraphWeight*_kFactor);
      TrackDZBS->Fill(fabs(track->dz (beamspot->position())), _madgraphWeight*_kFactor);
    }

    if (vertex != 0) {
      TrackD0PV->Fill(fabs(track->dxy(vertex->position())), _madgraphWeight*_kFactor);
      TrackDZPV->Fill(fabs(track->dz (vertex->position())), _madgraphWeight*_kFactor);
    }

    const reco::HitPattern& hp = track->hitPattern();
    NPxHits->Fill(hp.numberOfValidPixelHits(), _madgraphWeight*_kFactor);
    NStHits->Fill(hp.numberOfValidStripHits(), _madgraphWeight*_kFactor);
    NTkHits->Fill(hp.numberOfValidTrackerHits(), _madgraphWeight*_kFactor);
    NMuHits->Fill(hp.numberOfValidMuonHits(), _madgraphWeight*_kFactor);

    NHits->Fill(hp.numberOfValidHits(), _madgraphWeight*_kFactor);
    NInvalidHits->Fill(hp.numberOfAllHits(reco::HitPattern::TRACK_HITS) - hp.numberOfValidHits(), _madgraphWeight*_kFactor);
    //NInvalidHits->Fill(hp.numberOfAllHits() - hp.numberOfValidHits());
    
    NPxLayers->Fill(hp.pixelLayersWithMeasurement(), _madgraphWeight*_kFactor);
    NStLayers->Fill(hp.stripLayersWithMeasurement(), _madgraphWeight*_kFactor);
    NTkLayers->Fill(hp.trackerLayersWithMeasurement(), _madgraphWeight*_kFactor);
  }
}

void Zprime2muHistosFromPAT::fillOfflineElectronHistos(const pat::Electron* lep) {
  // Can add electron quantities here.
}

void Zprime2muHistosFromPAT::fillLeptonHistos(const reco::CandidateBaseRef& lep) {
  fillBasicLeptonHistos(lep);
  
  const pat::Muon* muon = toConcretePtr<pat::Muon>(lep);
  if (muon) fillOfflineMuonHistos(muon);
  
  const pat::Electron* electron = toConcretePtr<pat::Electron>(lep);
  if (electron) fillOfflineElectronHistos(electron);
}

void Zprime2muHistosFromPAT::fillLeptonHistos(const edm::View<reco::Candidate>& leptons) {
  NLeptons->Fill(leptons.size(), _madgraphWeight*_kFactor);

 // JMTBAD this should use leptonsPassingCuts or whatever
  int total_q = 0;
  for (size_t i = 0; i < leptons.size(); ++i) {
    total_q += leptons[i].charge();
    fillLeptonHistos(leptons.refAt(i));
  }

  LeptonSigns->Fill(leptons.size(), total_q);
}

void Zprime2muHistosFromPAT::fillLeptonHistosFromDileptons(const pat::CompositeCandidateCollection& dileptons) {
  int nleptons = 0;
  int total_q = 0;

  pat::CompositeCandidateCollection::const_iterator dil = dileptons.begin(), dile = dileptons.end();
  for ( ; dil != dile; ++dil)
    for (size_t i = 0; i < dil->numberOfDaughters(); ++i) {
      // JMTBAD if photons ever become daughters of the
      // CompositeCandidate, need to protect against this here
      fillLeptonHistos(dileptonDaughter(*dil, i));
      
      ++nleptons;
      total_q += dileptonDaughter(*dil, i)->charge();
    }

  // These become sanity checks.
  NLeptons->Fill(nleptons, _madgraphWeight*_kFactor);
  LeptonSigns->Fill(nleptons, total_q, _madgraphWeight*_kFactor);
}

void Zprime2muHistosFromPAT::fillDileptonHistos(const pat::CompositeCandidate& dil, const edm::Event& event, double gM) {

    
	kFactorGraph->Fill(_kFactor);
	kFactorGraph_bb->Fill(_kFactor_bb);
	kFactorGraph_be->Fill(_kFactor_be);
	
  if (dbg_tree) {
    dbg_t.mass = dil.mass();
    dbg_t.id = dil.daughter(0)->pdgId() + dil.daughter(1)->pdgId();
    dbg_tree->Fill();
  }
  DileptonEta->Fill(dil.eta(), _madgraphWeight*_kFactor);
  DileptonRap->Fill(dil.rapidity(), _madgraphWeight*_kFactor);
  DileptonPhi->Fill(dil.phi(), _madgraphWeight*_kFactor);

  DileptonPt->Fill(dil.pt(), _madgraphWeight*_kFactor);
  DileptonPz->Fill(fabs(dil.pz()), _madgraphWeight*_kFactor);
  DileptonP ->Fill(dil.p(), _madgraphWeight*_kFactor);

  DileptonPtVsEta->Fill(dil.eta(), dil.pt(), _madgraphWeight*_kFactor);
  DileptonPVsEta ->Fill(dil.eta(), dil.p(), _madgraphWeight*_kFactor);

  DileptonMass->Fill(dil.mass(), _madgraphWeight*_kFactor);
  DileptonMassWeight->Fill(dil.mass(),_prescaleWeight*_madgraphWeight*_kFactor);//?
  DileptonWithPhotonsMass->Fill(resonanceP4(dil).mass(), _madgraphWeight*_kFactor);

  const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
  const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);
  double cos_cs = -999.;
  if (lep0.isNonnull() && lep1.isNonnull()) {
    DileptonDeltaPt->Fill(fabs(lep0->pt()) - fabs(lep1->pt()), _madgraphWeight*_kFactor);
    DileptonDeltaP ->Fill(fabs(lep0->p())  - fabs(lep1->p()), _madgraphWeight*_kFactor);

     cos_cs = calcCosThetaCSAnal(lep0->pz(), lep0->energy(), lep1->pz(), lep1->energy(), dil.pt(), dil.pz(), dil.mass());
     CosThetaStarDilepton->Fill(cos_cs);
     //ChiDilepton->Fill((1+fabs(cos_cs))/(1-fabs(cos_cs)));
     ChiDilepton->Fill(exp(std::abs(lep0->p4().Rapidity()-lep1->p4().Rapidity())));
     if (cos_cs >= 0) DileptonMass_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor);
     else DileptonMass_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor);
    const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
    const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
    if (mu0 && mu1) {
      const reco::Track* tk0 = patmuon::getPickedTrack(*mu0).get();
      const reco::Track* tk1 = patmuon::getPickedTrack(*mu1).get();
      if (tk0 && tk1) {
	DimuonMuonPtErrors->Fill(ptError(tk0), ptError(tk1), _madgraphWeight*_kFactor);
	DimuonMuonPtErrOverPt->Fill(ptError(tk0)/tk0->pt(), _madgraphWeight*_kFactor);
	DimuonMuonPtErrOverPt->Fill(ptError(tk1)/tk1->pt(), _madgraphWeight*_kFactor);
	float mass = -999.;
	// Use mass calculated with the vertex constraint when available
	if (dil.hasUserFloat("vertexM"))
	  mass = dil.userFloat("vertexM");
	else{
	  edm::LogWarning("fillDileptonHistos") << "+++ Mass calculated without vertex constraint! +++";
	  mass = dil.mass();
	}
	if (mass > 200.) {
	  DimuonMuonPtErrOverPtM200->Fill(ptError(tk0)/tk0->pt(), _madgraphWeight*_kFactor);
	  DimuonMuonPtErrOverPtM200->Fill(ptError(tk1)/tk1->pt(), _madgraphWeight*_kFactor);
	}
	if (mass > 500.) {
	  DimuonMuonPtErrOverPtM500->Fill(ptError(tk0)/tk0->pt(), _madgraphWeight*_kFactor);
	  DimuonMuonPtErrOverPtM500->Fill(ptError(tk1)/tk1->pt(), _madgraphWeight*_kFactor);
	}
      }
    }
  }

  DileptonDaughterIds->Fill(dil.daughter(0)->pdgId(), dil.daughter(1)->pdgId(), _madgraphWeight*_kFactor);

  DileptonDaughterDeltaR->Fill(reco::deltaR(*dil.daughter(0), *dil.daughter(1)), _madgraphWeight*_kFactor);
  DileptonDaughterDeltaPhi->Fill(reco::deltaPhi(dil.daughter(0)->phi(), dil.daughter(1)->phi()), _madgraphWeight*_kFactor);

  if (dil.hasUserFloat("vertexM") && dil.hasUserFloat("vertexMError")) {
    float vertex_mass = dil.userFloat("vertexM");
    float vertex_mass_err = dil.userFloat("vertexMError");
      //std::cout<<" filling mass "<<vertex_mass<<std::endl;
    float smearedMass = getSmearedMass(dil,gM);
    DimuonMassVertexConstrained->Fill(vertex_mass, _madgraphWeight*_kFactor);
    DimuonMassVertexConstrainedSmear->Fill(smearedMass, _madgraphWeight*_kFactor);
    DimuonMassVertexConstrainedScaleUp->Fill(vertex_mass*(1+_scaleUncert), _madgraphWeight*_kFactor);
    DimuonMassVertexConstrainedScaleDown->Fill(vertex_mass*(1-_scaleUncert), _madgraphWeight*_kFactor);
    if (cos_cs > -998.){
     if (cos_cs >= 0) DimuonMassVertexConstrained_CSPos->Fill(vertex_mass, _madgraphWeight*_kFactor);
     else DimuonMassVertexConstrained_CSNeg->Fill(vertex_mass, _madgraphWeight*_kFactor);
    }
    DimuonMassVtxConstrainedLog->Fill(vertex_mass, _madgraphWeight*_kFactor);
    DimuonMassConstrainedVsUn->Fill(dil.mass(), vertex_mass, _madgraphWeight*_kFactor);
    DimuonMassVertexConstrainedError->Fill(vertex_mass, vertex_mass_err, _madgraphWeight*_kFactor);
    DimuonMassVertexConstrainedWeight->Fill(vertex_mass,_prescaleWeight*_madgraphWeight*_kFactor);
    DimuonMassVtxConstrainedLogWeight->Fill(vertex_mass,_prescaleWeight*_madgraphWeight*_kFactor);

  
    // plot per categories
  if (dil.daughter(0)->eta()<=1.2 && dil.daughter(1)->eta()<=1.2 && dil.daughter(0)->eta()>=-1.2 && dil.daughter(1)->eta()>=-1.2){
        DimuonMassVertexConstrained_bb->Fill(vertex_mass,_madgraphWeight*_kFactor_bb);
        DimuonMassVtxConstrainedLog_bb->Fill(vertex_mass, _madgraphWeight*_kFactor_bb);
        DimuonMassVertexConstrainedSmear_bb->Fill(smearedMass, _madgraphWeight*_kFactor);
        DimuonMassVertexConstrainedScaleUp_bb->Fill(vertex_mass*(1+_scaleUncert), _madgraphWeight*_kFactor);
        DimuonMassVertexConstrainedScaleDown_bb->Fill(vertex_mass*(1-_scaleUncert), _madgraphWeight*_kFactor);
       DileptonMass_bb->Fill(dil.mass(), _madgraphWeight*_kFactor_bb);
        if (cos_cs > -998.){
    		if (cos_cs >= 0){
			 DimuonMassVertexConstrained_bb_CSPos->Fill(vertex_mass, _madgraphWeight*_kFactor);
			 DimuonMassVertexConstrainedSmear_bb_CSPos->Fill(smearedMass, _madgraphWeight*_kFactor);
    			 DimuonMassVertexConstrainedScaleUp_bb_CSPos->Fill(vertex_mass*(1+_scaleUncert), _madgraphWeight*_kFactor);
    			 DimuonMassVertexConstrainedScaleDown_bb_CSPos->Fill(vertex_mass*(1-_scaleUncert), _madgraphWeight*_kFactor);
		}
     		else{
			 DimuonMassVertexConstrained_bb_CSNeg->Fill(vertex_mass, _madgraphWeight*_kFactor);
		         DimuonMassVertexConstrainedSmear_bb_CSNeg->Fill(smearedMass, _madgraphWeight*_kFactor);
    			 DimuonMassVertexConstrainedScaleUp_bb_CSNeg->Fill(vertex_mass*(1+_scaleUncert), _madgraphWeight*_kFactor);
    			 DimuonMassVertexConstrainedScaleDown_bb_CSNeg->Fill(vertex_mass*(1-_scaleUncert), _madgraphWeight*_kFactor);
		}
   		if (cos_cs >= 0) DileptonMass_bb_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor);
     		else DileptonMass_bb_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor);
   	 }



}

 if (dil.daughter(0)->eta()<-1.2 || dil.daughter(1)->eta()<-1.2 || dil.daughter(0)->eta()>1.2 || dil.daughter(1)->eta()>1.2){
        DimuonMassVertexConstrained_be->Fill(vertex_mass,_madgraphWeight*_kFactor_be);
        DimuonMassVtxConstrainedLog_be->Fill(vertex_mass, _madgraphWeight*_kFactor_be);
        DimuonMassVertexConstrainedSmear_be->Fill(smearedMass, _madgraphWeight*_kFactor);
        DimuonMassVertexConstrainedScaleUp_be->Fill(vertex_mass*(1+_scaleUncert), _madgraphWeight*_kFactor);
        DimuonMassVertexConstrainedScaleDown_be->Fill(vertex_mass*(1-_scaleUncert), _madgraphWeight*_kFactor);
       DileptonMass_be->Fill(dil.mass(), _madgraphWeight*_kFactor_be);
        if (cos_cs > -998.){
   		if (cos_cs >= 0){
			DimuonMassVertexConstrained_be_CSPos->Fill(vertex_mass, _madgraphWeight*_kFactor);
		        DimuonMassVertexConstrainedSmear_be_CSPos->Fill(smearedMass, _madgraphWeight*_kFactor);
                        DimuonMassVertexConstrainedScaleUp_be_CSPos->Fill(vertex_mass*(1+_scaleUncert), _madgraphWeight*_kFactor);
   			DimuonMassVertexConstrainedScaleDown_be_CSPos->Fill(vertex_mass*(1-_scaleUncert), _madgraphWeight*_kFactor);
		}
     		else{
			 DimuonMassVertexConstrained_be_CSNeg->Fill(vertex_mass, _madgraphWeight*_kFactor);
		         DimuonMassVertexConstrainedSmear_be_CSNeg->Fill(smearedMass, _madgraphWeight*_kFactor);
        		 DimuonMassVertexConstrainedScaleUp_be_CSNeg->Fill(vertex_mass*(1+_scaleUncert), _madgraphWeight*_kFactor);
    			 DimuonMassVertexConstrainedScaleDown_be_CSNeg->Fill(vertex_mass*(1-_scaleUncert), _madgraphWeight*_kFactor);
		}
    		if (cos_cs >= 0) DileptonMass_be_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor);
     		else DileptonMass_be_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor);
   	 }



}
	
    // special
    float vertex_chi2 = dil.userFloat("vertex_chi2");
      DimuonMassVtx_chi2->Fill(vertex_chi2, _madgraphWeight*_kFactor);
    if (vertex_chi2 > 0 ) {
        float vertex_ndof = dil.userFloat("vertex_ndof");
        float vertex_chi2_noNormalized = vertex_chi2*vertex_ndof;
        DimuonMassVtx_prob->Fill(TMath::Prob(vertex_chi2_noNormalized, vertex_ndof), _madgraphWeight*_kFactor);}

  }
}

void Zprime2muHistosFromPAT::fillDileptonHistos(const pat::CompositeCandidateCollection& dileptons, const edm::Event& event, double gM) {
  NDileptons->Fill(dileptons.size(), _madgraphWeight*_kFactor);

  pat::CompositeCandidateCollection::const_iterator dil = dileptons.begin(), dile = dileptons.end();
  for ( ; dil != dile; ++dil)
    fillDileptonHistos(*dil, event, gM);
}

void Zprime2muHistosFromPAT::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  if (dbg_tree) {
    memset(&dbg_t, 0, sizeof(debug_tree_t));
    dbg_t.run = event.id().run();
    dbg_t.lumi = event.luminosityBlock();
    dbg_t.event = event.id().event();
  }


//  edm::Handle<int> hltPrescale;
//  edm::Handle<int> l1Prescale;

//  event.getByLabel(edm::InputTag("getPrescales","HLTPrescale","Zprime2muAnalysis"), hltPrescale);
//  event.getByLabel(edm::InputTag("getPrescales","L1Prescale","Zprime2muAnalysis"), l1Prescale);
    if (_usePrescaleWeight) {
        edm::Handle<int> totalPrescale;
        event.getByLabel(edm::InputTag("getPrescales","TotalPrescale","Zprime2muAnalysis"), totalPrescale);
        _prescaleWeight = *totalPrescale;
    }
//    std::cout<<*hltPrescale<<std::endl;
//    std::cout<<l1Prescale<<std::endl;
//    std::cout<<totalPrescale<<std::endl;
    if (_useMadgraphWeight) {
        eventWeight = 1.;
	_madgraphWeight = 1.; 
        edm::Handle<GenEventInfoProduct> gen_ev_info;
        event.getByLabel(edm::InputTag("generator"), gen_ev_info);
	if (gen_ev_info.isValid()){
        	eventWeight = gen_ev_info->weight();
        	_madgraphWeight = ( eventWeight > 0 ) ? 1 : -1;
	}
        WeightMadGraph->Fill(_madgraphWeight);
    }
    

  if (use_bs_and_pv)
    getBSandPV(event);

  edm::Handle<edm::View<reco::Candidate> > leptons;
  event.getByLabel(lepton_src, leptons);

  if (!leptons.isValid())
    edm::LogWarning("LeptonHandleInvalid") << "tried to get " << lepton_src << " with edm::Handle<edm::View<reco::Candidate> > and failed!";
  else {
    if (!leptonsFromDileptons)
      fillLeptonHistos(*leptons);
  }

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);
  
  double gM = -1; 
  if (!dileptons.isValid())
    edm::LogWarning("DileptonHandleInvalid") << "tried to get " << dilepton_src << " and failed!";
  else {
    if (leptonsFromDileptons)
        if (_usekFactor){
    	hardInteraction->Fill(event);
    	if (fill_gen_info) {
// 			if(hardInteraction->IsValid()){
			if(hardInteraction->IsValidForRes()){
    			gM = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).mass();
// 	    		_kFactor = 50.;
// 				if(60 < gM && gM < 120){
// 					_kFactor = 0.9782909298;
// 					_kFactor_bb = 0.9563389736;
// 					_kFactor_be = 0.9986054970;
// 				}
// 				if(120 < gM && gM < 150){
// 					_kFactor = 1.2401820110;
// 					_kFactor_bb = 1.2306943240;
// 					_kFactor_be = 1.2476210179;
// 				}
// 				if(150 < gM && gM < 200){
// 					_kFactor = 1.1119959905;
// 					_kFactor_bb = 1.1166729493;
// 					_kFactor_be = 1.1089891649;
// 				}
// 				if(200 < gM && gM < 300){
// 					_kFactor = 1.0711027694;
// 					_kFactor_bb = 1.0433290922;
// 					_kFactor_be = 1.0867974776;
// 				}
// 				if(300 < gM && gM < 400){
// 					_kFactor = 1.0588399982;
// 					_kFactor_bb = 1.0168447336;
// 					_kFactor_be = 1.0810349902;
// 				}
// 				if(gM > 400){
// 			 		   	_kFactor = 1.039 - 0.0001313 * gM + 4.733e-08 * pow(gM,2) - 7.385e-12 * pow(gM,3);
// 			 		   	_kFactor_bb = 1.012 - 9.968e-5 * gM + 3.321e-08 * pow(gM,2) - 5.694e-12 * pow(gM,3);
// 			 		   	_kFactor_be = 1.056 - 0.0001537 * gM + 6.071e-08 * pow(gM,2) - 9.093e-12 * pow(gM,3);
// 			 	}
				if(gM < 150){
					_kFactor = 1;
					_kFactor_bb = 1;
					_kFactor_be = 1;
				}
				if(gM > 150){
			 		   	_kFactor = 1.053 - 0.0001552 * gM + 5.661e-08 * pow(gM,2) - 8.382e-12 * pow(gM,3);
			 		   	_kFactor_bb = 1.032 - 0.000138 * gM + 4.827e-08 * pow(gM,2) - 7.321e-12 * pow(gM,3);
			 		   	_kFactor_be = 1.064 - 0.0001674 * gM + 6.599e-08 * pow(gM,2) - 9.657e-12 * pow(gM,3);
			 	}
			 	
	    		//std::cout<<"----------------------------------------------------------- GEN MASS = " <<hardInteraction->resonance->mass()<<std::endl;
    			//std::cout<<"------------------------------------------------------------------ kFactor = "<<_kFactor<<" --- BB = "<<_kFactor_bb<<" --- BE = "<<_kFactor_be<<std::endl;
    		} // hardInter
    		else
    			std::cout<<"problems"<<std::endl;
    		} //gen_info
    	} //kFactor
        
      fillLeptonHistosFromDileptons(*dileptons);
    
    fillDileptonHistos(*dileptons, event, gM);
  }
}

DEFINE_FWK_MODULE(Zprime2muHistosFromPAT);
