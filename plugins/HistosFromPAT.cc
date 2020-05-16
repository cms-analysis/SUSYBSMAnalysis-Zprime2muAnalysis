// JMTBAD make this fill lepton histos from all leptons + those that
// make it into dileptons always, not just depending on flag.

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
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
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PUUtilities1617.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PUUtilities18.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/triggerTurnOnElectrons.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"///
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/LRWeightProducer.h"
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
  double ScaleUncert(double, bool, int);
  double getRecoWeight(double, bool, int);
  double getSmearedMass(const pat::CompositeCandidate&, double, int, bool);
  double HEEPSF(double, int);
  double turnOn(double, double, int);
  double L1TurnOn(double, double, int);

  edm::InputTag lepton_src;
  edm::InputTag dilepton_src;
  const bool leptonsFromDileptons;
  const bool doElectrons;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  edm::InputTag pu_src;
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
    double _LRWeight = 1.;///
    bool _usekFactor;
    bool _useTTBarWeight;
    double _kFactor;
    double _kFactor_bb;
    double _kFactor_be;
    double _eleMCFac;
    double _scaleUncertBB = 0.01;
    double _scaleUncertBE = 0.03;
    double _scaleUncertEleBB = 0.02;
    double _scaleUncertEleBE = 0.01;
    int _nTrueInt = 0;
    double _puWeight = 1.0;
    double _puWeight_scaleUp = 1.0;
    double _puWeight_scaleDown = 1.0;
    double _prefireWeight = 1.0;
    double _prefireWeightUp = 1.0;
    double _prefireWeightDown = 1.0;
  TH1F* NBeamSpot;
  TH1F* NVertices;
  TH1F* NVerticesUnweighted;
  TH1F* NTrueInteractions;
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

  TH1F* DielectronMass;
  TH1F* DielectronMass_bbbe;
  TH1F* DielectronMass_bb;
  TH1F* DielectronMass_be;
  TH1F* DielectronMass_ee;
  TH1F* DielectronMass_gen_bb;
  TH1F* DielectronMass_gen_be;
  TH1F* DielectronMass_gen_ee;
  TH1F* DielectronMass_CSPos;
  TH1F* DielectronMass_bb_CSPos;
  TH1F* DielectronMass_bbbe_CSPos;
  TH1F* DielectronMass_be_CSPos;
  TH1F* DielectronMass_ee_CSPos;
  TH1F* DielectronMass_CSNeg;
  TH1F* DielectronMass_bb_CSNeg;
  TH1F* DielectronMass_bbbe_CSNeg;
  TH1F* DielectronMass_be_CSNeg;
  TH1F* DielectronMass_ee_CSNeg;

  TH1F* DielectronMassScaleUp;
  TH1F* DielectronMassScaleUp_bbbe;
  TH1F* DielectronMassScaleUp_bb;
  TH1F* DielectronMassScaleUp_be;
  TH1F* DielectronMassScaleUp_ee;
  TH1F* DielectronMassScaleUp_CSPos;
  TH1F* DielectronMassScaleUp_bb_CSPos;
  TH1F* DielectronMassScaleUp_bbbe_CSPos;
  TH1F* DielectronMassScaleUp_be_CSPos;
  TH1F* DielectronMassScaleUp_ee_CSPos;
  TH1F* DielectronMassScaleUp_CSNeg;
  TH1F* DielectronMassScaleUp_bb_CSNeg;
  TH1F* DielectronMassScaleUp_bbbe_CSNeg;
  TH1F* DielectronMassScaleUp_be_CSNeg;
  TH1F* DielectronMassScaleUp_ee_CSNeg;

  TH1F* DielectronMassScaleDown;
  TH1F* DielectronMassScaleDown_bbbe;
  TH1F* DielectronMassScaleDown_bb;
  TH1F* DielectronMassScaleDown_be;
  TH1F* DielectronMassScaleDown_ee;
  TH1F* DielectronMassScaleDown_CSPos;
  TH1F* DielectronMassScaleDown_bb_CSPos;
  TH1F* DielectronMassScaleDown_bbbe_CSPos;
  TH1F* DielectronMassScaleDown_be_CSPos;
  TH1F* DielectronMassScaleDown_ee_CSPos;
  TH1F* DielectronMassScaleDown_CSNeg;
  TH1F* DielectronMassScaleDown_bb_CSNeg;
  TH1F* DielectronMassScaleDown_bbbe_CSNeg;
  TH1F* DielectronMassScaleDown_be_CSNeg;
  TH1F* DielectronMassScaleDown_ee_CSNeg;

  TH1F* DielectronMassPrefireUp;
  TH1F* DielectronMassPrefireUp_bbbe;
  TH1F* DielectronMassPrefireUp_bb;
  TH1F* DielectronMassPrefireUp_be;
  TH1F* DielectronMassPrefireUp_ee;
  TH1F* DielectronMassPrefireUp_CSPos;
  TH1F* DielectronMassPrefireUp_bb_CSPos;
  TH1F* DielectronMassPrefireUp_bbbe_CSPos;
  TH1F* DielectronMassPrefireUp_be_CSPos;
  TH1F* DielectronMassPrefireUp_ee_CSPos;
  TH1F* DielectronMassPrefireUp_CSNeg;
  TH1F* DielectronMassPrefireUp_bb_CSNeg;
  TH1F* DielectronMassPrefireUp_bbbe_CSNeg;
  TH1F* DielectronMassPrefireUp_be_CSNeg;
  TH1F* DielectronMassPrefireUp_ee_CSNeg;

  TH1F* DielectronMassPrefireDown;
  TH1F* DielectronMassPrefireDown_bbbe;
  TH1F* DielectronMassPrefireDown_bb;
  TH1F* DielectronMassPrefireDown_be;
  TH1F* DielectronMassPrefireDown_ee;
  TH1F* DielectronMassPrefireDown_CSPos;
  TH1F* DielectronMassPrefireDown_bb_CSPos;
  TH1F* DielectronMassPrefireDown_bbbe_CSPos;
  TH1F* DielectronMassPrefireDown_be_CSPos;
  TH1F* DielectronMassPrefireDown_ee_CSPos;
  TH1F* DielectronMassPrefireDown_CSNeg;
  TH1F* DielectronMassPrefireDown_bb_CSNeg;
  TH1F* DielectronMassPrefireDown_bbbe_CSNeg;
  TH1F* DielectronMassPrefireDown_be_CSNeg;
  TH1F* DielectronMassPrefireDown_ee_CSNeg;

  TH1F* DielectronMassPUScaleUp;
  TH1F* DielectronMassPUScaleUp_bbbe;
  TH1F* DielectronMassPUScaleUp_bb;
  TH1F* DielectronMassPUScaleUp_be;
  TH1F* DielectronMassPUScaleUp_ee;
  TH1F* DielectronMassPUScaleUp_CSPos;
  TH1F* DielectronMassPUScaleUp_bb_CSPos;
  TH1F* DielectronMassPUScaleUp_bbbe_CSPos;
  TH1F* DielectronMassPUScaleUp_be_CSPos;
  TH1F* DielectronMassPUScaleUp_ee_CSPos;
  TH1F* DielectronMassPUScaleUp_CSNeg;
  TH1F* DielectronMassPUScaleUp_bb_CSNeg;
  TH1F* DielectronMassPUScaleUp_bbbe_CSNeg;
  TH1F* DielectronMassPUScaleUp_be_CSNeg;
  TH1F* DielectronMassPUScaleUp_ee_CSNeg;

  TH1F* DielectronMassPUScaleDown;
  TH1F* DielectronMassPUScaleDown_bbbe;
  TH1F* DielectronMassPUScaleDown_bb;
  TH1F* DielectronMassPUScaleDown_be;
  TH1F* DielectronMassPUScaleDown_ee;
  TH1F* DielectronMassPUScaleDown_CSPos;
  TH1F* DielectronMassPUScaleDown_bb_CSPos;
  TH1F* DielectronMassPUScaleDown_bbbe_CSPos;
  TH1F* DielectronMassPUScaleDown_be_CSPos;
  TH1F* DielectronMassPUScaleDown_ee_CSPos;
  TH1F* DielectronMassPUScaleDown_CSNeg;
  TH1F* DielectronMassPUScaleDown_bb_CSNeg;
  TH1F* DielectronMassPUScaleDown_bbbe_CSNeg;
  TH1F* DielectronMassPUScaleDown_be_CSNeg;
  TH1F* DielectronMassPUScaleDown_ee_CSNeg;

  TH2F* DielectronMassVsCS;
  TH2F* DielectronMassVsCS_bbbe;
  TH2F* DielectronMassVsCS_bb;
  TH2F* DielectronMassVsCS_be;
  TH2F* DielectronMassVsCS_ee;
  TH2F* DielectronResponse_bb;
  TH2F* DielectronResponse_be;
  TH2F* DielectronResponse_ee;
  TH2F* DielectronResponseMassScaleUp_bb;
  TH2F* DielectronResponseMassScaleUp_be;
  TH2F* DielectronResponseMassScaleUp_ee;
  TH2F* DielectronResponseMassScaleDown_bb;
  TH2F* DielectronResponseMassScaleDown_be;
  TH2F* DielectronResponseMassScaleDown_ee;

  TH1F* GenMass;
  
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
  TH1F* DimuonMassVertexConstrained_gen_bb;
  TH1F* DimuonMassVertexConstrained_gen_be;
  TH1F* DimuonMassVertexConstrained_CSPos;
  TH1F* DimuonMassVertexConstrained_bb_CSPos;
  TH1F* DimuonMassVertexConstrained_be_CSPos;
  TH1F* DimuonMassVertexConstrained_CSNeg;
  TH1F* DimuonMassVertexConstrained_bb_CSNeg;
  TH1F* DimuonMassVertexConstrained_be_CSNeg;
  TH1F* DimuonMassVertexConstrainedMuonID;
  TH1F* DimuonMassVertexConstrainedMuonID_bb;
  TH1F* DimuonMassVertexConstrainedMuonID_be;
  TH1F* DimuonMassVertexConstrainedMuonID_CSPos;
  TH1F* DimuonMassVertexConstrainedMuonID_bb_CSPos;
  TH1F* DimuonMassVertexConstrainedMuonID_be_CSPos;
  TH1F* DimuonMassVertexConstrainedMuonID_CSNeg;
  TH1F* DimuonMassVertexConstrainedMuonID_bb_CSNeg;
  TH1F* DimuonMassVertexConstrainedMuonID_be_CSNeg;
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

  TH2F* DimuonMassVertexConstrainedVsCS;
  TH2F* DimuonMassVertexConstrainedVsCS_bb;
  TH2F* DimuonMassVertexConstrainedVsCS_be;
  TH2F* DimuonResponse_bb;
  TH2F* DimuonResponse_be;
  TH2F* DimuonResponseMassScaleUp_bb;
  TH2F* DimuonResponseMassScaleUp_be;
  TH2F* DimuonResponseMassScaleDown_bb;
  TH2F* DimuonResponseMassScaleDown_be;
  TH2F* DimuonResponseSmear_bb;
  TH2F* DimuonResponseSmear_be;
  TH2F* DimuonResponseMassScaleUpSmear_bb;
  TH2F* DimuonResponseMassScaleUpSmear_be;
  TH2F* DimuonResponseMassScaleDownSmear_bb;
  TH2F* DimuonResponseMassScaleDownSmear_be;


  TH1F* DimuonMassVertexConstrainedWeight;
  TH1F* DimuonMassVtxConstrainedLog;
  TH1F* DimuonMassVtxConstrainedLog_bb;
  TH1F* DimuonMassVtxConstrainedLog_be;
  TH1F* DimuonMassVtxConstrainedLog_gen_bb;
  TH1F* DimuonMassVtxConstrainedLog_gen_be;
  TH1F* DimuonMassVtxConstrainedLogWeight;
  TH2F* DimuonMassConstrainedVsUn;
  TH2F* DimuonMassVertexConstrainedError;
    //special
    TH1F* DimuonMassVtx_chi2;
    TH1F* DimuonMassVtx_prob;
    //weight
    TH1F* WeightMadGraph;///
    TH1F* WeightLR;///
    TH1F *kFactorGraph;
    TH1F *kFactorGraph_bb;
    TH1F *kFactorGraph_be;
    
	const bool fill_gen_info;
	HardInteraction* hardInteraction;
  	std::vector<std::string> pu_info;  
	int year_info;
 	LRWeightProducer lrWeightProducer;
	edm::EDGetTokenT< double > prefweight_token;
	edm::EDGetTokenT< double > prefweightup_token;
	edm::EDGetTokenT< double > prefweightdown_token;
        edm::EDGetTokenT<LHEEventProduct>LHEEventToken_;
};

Zprime2muHistosFromPAT::Zprime2muHistosFromPAT(const edm::ParameterSet& cfg)
  : lepton_src(cfg.getParameter<edm::InputTag>("lepton_src")),
    dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    leptonsFromDileptons(cfg.getParameter<bool>("leptonsFromDileptons")),
    doElectrons(cfg.getParameter<bool>("doElectrons")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    pu_src(cfg.getParameter<edm::InputTag>("pu_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    dbg_tree(0),
    beamspot(0),
    vertex(0),
    _usePrescaleWeight(cfg.getUntrackedParameter<bool>("usePrescaleWeight",false)),
    _prescaleWeight(1),
    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),///
    _madgraphWeight(1.),///
    _usekFactor(cfg.getParameter<bool>("usekFactor")),
    _useTTBarWeight(cfg.getParameter<bool>("useTTBarWeight")),
    _kFactor(1.),
    _kFactor_bb(1.),
    _kFactor_be(1.),
    fill_gen_info(cfg.existsAs<edm::ParameterSet>("hardInteraction")),
    hardInteraction(fill_gen_info ? new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")) : 0),
    pu_info(cfg.getParameter<std::vector<std::string>>("pu_weights")),
    year_info(cfg.getParameter<int>("year")),
    lrWeightProducer(cfg.getParameter<edm::ParameterSet>("lrWeightProducer")),
    LHEEventToken_(consumes<LHEEventProduct>(cfg.getParameter<edm::InputTag>("LHEInfo")))
{

  consumes<reco::CandidateView>(lepton_src);
  consumes<pat::CompositeCandidateCollection>(dilepton_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);
  consumes<std::vector<PileupSummaryInfo>>(pu_src);
  consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  if (fill_gen_info) consumes<std::vector<reco::GenParticle>>(hardInteraction->src);
  if ((year_info==2016 || year_info==2017) && doElectrons) {
	prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
	prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
	prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
  }

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
  NVertices = fs->make<TH1F>("NVertices", titlePrefix + "# vertices/event",  120, 0, 120);
  NVerticesUnweighted = fs->make<TH1F>("NVerticesUnweighted", titlePrefix + "# vertices/event",  120, 0, 120);
  NTrueInteractions = fs->make<TH1F>("NTrueInteractiosn", titlePrefix + "# true interactions/event",  120, 0, 120);

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
  DielectronMass            = fs->make<TH1F>("DielectronMass",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMass_bbbe       = fs->make<TH1F>("DielectronMass_bbbe",       titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMass_bb         = fs->make<TH1F>("DielectronMass_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DielectronMass_be         = fs->make<TH1F>("DielectronMass_be",            titlePrefix + "dil. mass barrel-endcaps", 20000, 0, 20000);
  DielectronMass_ee         = fs->make<TH1F>("DielectronMass_ee",            titlePrefix + "dil. mass endcaps-endcaps", 20000, 0, 20000);
  DielectronMass_gen_bb         = fs->make<TH1F>("DielectronMass_gen_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DielectronMass_gen_be         = fs->make<TH1F>("DielectronMass_gen_be",            titlePrefix + "dil. mass barrel-endcaps", 20000, 0, 20000);
  DielectronMass_gen_ee         = fs->make<TH1F>("DielectronMass_gen_ee",            titlePrefix + "dil. mass endcaps-endcaps", 20000, 0, 20000);
  DielectronMass_CSPos            = fs->make<TH1F>("DielectronMass_CSPos",            titlePrefix + "dil. mass for positive cos theta star", 20000, 0, 20000);
  DielectronMass_bb_CSPos         = fs->make<TH1F>("DielectronMass_bb_CSPos",            titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMass_bbbe_CSPos       = fs->make<TH1F>("DielectronMass_bbbe_CSPos",          titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMass_be_CSPos         = fs->make<TH1F>("DielectronMass_be_CSPos",            titlePrefix + "dil. mass barrel-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMass_ee_CSPos         = fs->make<TH1F>("DielectronMass_ee_CSPos",            titlePrefix + "dil. mass endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMass_CSNeg            = fs->make<TH1F>("DielectronMass_CSNeg",            titlePrefix + "dil. mass for negative cos theta star", 20000, 0, 20000);
  DielectronMass_bb_CSNeg         = fs->make<TH1F>("DielectronMass_bb_CSNeg",            titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMass_bbbe_CSNeg       = fs->make<TH1F>("DielectronMass_bbbe_CSNeg",          titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMass_be_CSNeg         = fs->make<TH1F>("DielectronMass_be_CSNeg",            titlePrefix + "dil. mass barrel-endcaps for negative cos theta star", 20000, 0, 20000);
  DielectronMass_ee_CSNeg         = fs->make<TH1F>("DielectronMass_ee_CSNeg",            titlePrefix + "dil. mass endcaps-endcaps for negative cos theta star", 20000, 0, 20000);

  DielectronMassScaleUp            = fs->make<TH1F>("DielectronMassScaleUp",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassScaleUp_bbbe       = fs->make<TH1F>("DielectronMassScaleUp_bbbe",       titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassScaleUp_bb         = fs->make<TH1F>("DielectronMassScaleUp_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DielectronMassScaleUp_be         = fs->make<TH1F>("DielectronMassScaleUp_be",            titlePrefix + "dil. mass barrel-endcaps", 20000, 0, 20000);
  DielectronMassScaleUp_ee         = fs->make<TH1F>("DielectronMassScaleUp_ee",            titlePrefix + "dil. mass endcaps-endcaps", 20000, 0, 20000);
  DielectronMassScaleUp_CSPos            = fs->make<TH1F>("DielectronMassScaleUp_CSPos",            titlePrefix + "dil. mass for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_bb_CSPos         = fs->make<TH1F>("DielectronMassScaleUp_bb_CSPos",            titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_bbbe_CSPos       = fs->make<TH1F>("DielectronMassScaleUp_bbbe_CSPos",          titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_be_CSPos         = fs->make<TH1F>("DielectronMassScaleUp_be_CSPos",            titlePrefix + "dil. mass barrel-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_ee_CSPos         = fs->make<TH1F>("DielectronMassScaleUp_ee_CSPos",            titlePrefix + "dil. mass endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_CSNeg            = fs->make<TH1F>("DielectronMassScaleUp_CSNeg",            titlePrefix + "dil. mass for negative cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_bb_CSNeg         = fs->make<TH1F>("DielectronMassScaleUp_bb_CSNeg",            titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_bbbe_CSNeg       = fs->make<TH1F>("DielectronMassScaleUp_bbbe_CSNeg",          titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_be_CSNeg         = fs->make<TH1F>("DielectronMassScaleUp_be_CSNeg",            titlePrefix + "dil. mass barrel-endcaps for negative cos theta star", 20000, 0, 20000);
  DielectronMassScaleUp_ee_CSNeg         = fs->make<TH1F>("DielectronMassScaleUp_ee_CSNeg",            titlePrefix + "dil. mass endcaps-endcaps for negative cos theta star", 20000, 0, 20000);

  DielectronMassScaleDown            = fs->make<TH1F>("DielectronMassScaleDown",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassScaleDown_bbbe       = fs->make<TH1F>("DielectronMassScaleDown_bbbe",       titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassScaleDown_bb         = fs->make<TH1F>("DielectronMassScaleDown_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DielectronMassScaleDown_be         = fs->make<TH1F>("DielectronMassScaleDown_be",            titlePrefix + "dil. mass barrel-endcaps", 20000, 0, 20000);
  DielectronMassScaleDown_ee         = fs->make<TH1F>("DielectronMassScaleDown_ee",            titlePrefix + "dil. mass endcaps-endcaps", 20000, 0, 20000);
  DielectronMassScaleDown_CSPos            = fs->make<TH1F>("DielectronMassScaleDown_CSPos",            titlePrefix + "dil. mass for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_bb_CSPos         = fs->make<TH1F>("DielectronMassScaleDown_bb_CSPos",            titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_bbbe_CSPos       = fs->make<TH1F>("DielectronMassScaleDown_bbbe_CSPos",          titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_be_CSPos         = fs->make<TH1F>("DielectronMassScaleDown_be_CSPos",            titlePrefix + "dil. mass barrel-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_ee_CSPos         = fs->make<TH1F>("DielectronMassScaleDown_ee_CSPos",            titlePrefix + "dil. mass endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_CSNeg            = fs->make<TH1F>("DielectronMassScaleDown_CSNeg",            titlePrefix + "dil. mass for negative cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_bb_CSNeg         = fs->make<TH1F>("DielectronMassScaleDown_bb_CSNeg",            titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_bbbe_CSNeg       = fs->make<TH1F>("DielectronMassScaleDown_bbbe_CSNeg",          titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_be_CSNeg         = fs->make<TH1F>("DielectronMassScaleDown_be_CSNeg",            titlePrefix + "dil. mass barrel-endcaps for negative cos theta star", 20000, 0, 20000);
  DielectronMassScaleDown_ee_CSNeg         = fs->make<TH1F>("DielectronMassScaleDown_ee_CSNeg",            titlePrefix + "dil. mass endcaps-endcaps for negative cos theta star", 20000, 0, 20000);

  DielectronMassPrefireUp            = fs->make<TH1F>("DielectronMassPrefireUp",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassPrefireUp_bbbe       = fs->make<TH1F>("DielectronMassPrefireUp_bbbe",       titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassPrefireUp_bb         = fs->make<TH1F>("DielectronMassPrefireUp_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DielectronMassPrefireUp_be         = fs->make<TH1F>("DielectronMassPrefireUp_be",            titlePrefix + "dil. mass barrel-endcaps", 20000, 0, 20000);
  DielectronMassPrefireUp_ee         = fs->make<TH1F>("DielectronMassPrefireUp_ee",            titlePrefix + "dil. mass endcaps-endcaps", 20000, 0, 20000);
  DielectronMassPrefireUp_CSPos            = fs->make<TH1F>("DielectronMassPrefireUp_CSPos",            titlePrefix + "dil. mass for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_bb_CSPos         = fs->make<TH1F>("DielectronMassPrefireUp_bb_CSPos",            titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_bbbe_CSPos       = fs->make<TH1F>("DielectronMassPrefireUp_bbbe_CSPos",          titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_be_CSPos         = fs->make<TH1F>("DielectronMassPrefireUp_be_CSPos",            titlePrefix + "dil. mass barrel-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_ee_CSPos         = fs->make<TH1F>("DielectronMassPrefireUp_ee_CSPos",            titlePrefix + "dil. mass endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_CSNeg            = fs->make<TH1F>("DielectronMassPrefireUp_CSNeg",            titlePrefix + "dil. mass for negative cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_bb_CSNeg         = fs->make<TH1F>("DielectronMassPrefireUp_bb_CSNeg",            titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_bbbe_CSNeg       = fs->make<TH1F>("DielectronMassPrefireUp_bbbe_CSNeg",          titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_be_CSNeg         = fs->make<TH1F>("DielectronMassPrefireUp_be_CSNeg",            titlePrefix + "dil. mass barrel-endcaps for negative cos theta star", 20000, 0, 20000);
  DielectronMassPrefireUp_ee_CSNeg         = fs->make<TH1F>("DielectronMassPrefireUp_ee_CSNeg",            titlePrefix + "dil. mass endcaps-endcaps for negative cos theta star", 20000, 0, 20000);

  DielectronMassPrefireDown            = fs->make<TH1F>("DielectronMassPrefireDown",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassPrefireDown_bbbe       = fs->make<TH1F>("DielectronMassPrefireDown_bbbe",       titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassPrefireDown_bb         = fs->make<TH1F>("DielectronMassPrefireDown_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DielectronMassPrefireDown_be         = fs->make<TH1F>("DielectronMassPrefireDown_be",            titlePrefix + "dil. mass barrel-endcaps", 20000, 0, 20000);
  DielectronMassPrefireDown_ee         = fs->make<TH1F>("DielectronMassPrefireDown_ee",            titlePrefix + "dil. mass endcaps-endcaps", 20000, 0, 20000);
  DielectronMassPrefireDown_CSPos            = fs->make<TH1F>("DielectronMassPrefireDown_CSPos",            titlePrefix + "dil. mass for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_bb_CSPos         = fs->make<TH1F>("DielectronMassPrefireDown_bb_CSPos",            titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_bbbe_CSPos       = fs->make<TH1F>("DielectronMassPrefireDown_bbbe_CSPos",          titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_be_CSPos         = fs->make<TH1F>("DielectronMassPrefireDown_be_CSPos",            titlePrefix + "dil. mass barrel-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_ee_CSPos         = fs->make<TH1F>("DielectronMassPrefireDown_ee_CSPos",            titlePrefix + "dil. mass endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_CSNeg            = fs->make<TH1F>("DielectronMassPrefireDown_CSNeg",            titlePrefix + "dil. mass for negative cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_bb_CSNeg         = fs->make<TH1F>("DielectronMassPrefireDown_bb_CSNeg",            titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_bbbe_CSNeg       = fs->make<TH1F>("DielectronMassPrefireDown_bbbe_CSNeg",          titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_be_CSNeg         = fs->make<TH1F>("DielectronMassPrefireDown_be_CSNeg",            titlePrefix + "dil. mass barrel-endcaps for negative cos theta star", 20000, 0, 20000);
  DielectronMassPrefireDown_ee_CSNeg         = fs->make<TH1F>("DielectronMassPrefireDown_ee_CSNeg",            titlePrefix + "dil. mass endcaps-endcaps for negative cos theta star", 20000, 0, 20000);

  DielectronMassPUScaleUp            = fs->make<TH1F>("DielectronMassPUScaleUp",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassPUScaleUp_bbbe       = fs->make<TH1F>("DielectronMassPUScaleUp_bbbe",       titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassPUScaleUp_bb         = fs->make<TH1F>("DielectronMassPUScaleUp_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DielectronMassPUScaleUp_be         = fs->make<TH1F>("DielectronMassPUScaleUp_be",            titlePrefix + "dil. mass barrel-endcaps", 20000, 0, 20000);
  DielectronMassPUScaleUp_ee         = fs->make<TH1F>("DielectronMassPUScaleUp_ee",            titlePrefix + "dil. mass endcaps-endcaps", 20000, 0, 20000);
  DielectronMassPUScaleUp_CSPos            = fs->make<TH1F>("DielectronMassPUScaleUp_CSPos",            titlePrefix + "dil. mass for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_bb_CSPos         = fs->make<TH1F>("DielectronMassPUScaleUp_bb_CSPos",            titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_bbbe_CSPos       = fs->make<TH1F>("DielectronMassPUScaleUp_bbbe_CSPos",          titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_be_CSPos         = fs->make<TH1F>("DielectronMassPUScaleUp_be_CSPos",            titlePrefix + "dil. mass barrel-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_ee_CSPos         = fs->make<TH1F>("DielectronMassPUScaleUp_ee_CSPos",            titlePrefix + "dil. mass endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_CSNeg            = fs->make<TH1F>("DielectronMassPUScaleUp_CSNeg",            titlePrefix + "dil. mass for negative cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_bb_CSNeg         = fs->make<TH1F>("DielectronMassPUScaleUp_bb_CSNeg",            titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_bbbe_CSNeg       = fs->make<TH1F>("DielectronMassPUScaleUp_bbbe_CSNeg",          titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_be_CSNeg         = fs->make<TH1F>("DielectronMassPUScaleUp_be_CSNeg",            titlePrefix + "dil. mass barrel-endcaps for negative cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleUp_ee_CSNeg         = fs->make<TH1F>("DielectronMassPUScaleUp_ee_CSNeg",            titlePrefix + "dil. mass endcaps-endcaps for negative cos theta star", 20000, 0, 20000);

  DielectronMassPUScaleDown            = fs->make<TH1F>("DielectronMassPUScaleDown",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassPUScaleDown_bbbe       = fs->make<TH1F>("DielectronMassPUScaleDown_bbbe",       titlePrefix + "dil. mass", 20000, 0, 20000);
  DielectronMassPUScaleDown_bb         = fs->make<TH1F>("DielectronMassPUScaleDown_bb",            titlePrefix + "dil. mass barrel-barrel", 20000, 0, 20000);
  DielectronMassPUScaleDown_be         = fs->make<TH1F>("DielectronMassPUScaleDown_be",            titlePrefix + "dil. mass barrel-endcaps", 20000, 0, 20000);
  DielectronMassPUScaleDown_ee         = fs->make<TH1F>("DielectronMassPUScaleDown_ee",            titlePrefix + "dil. mass endcaps-endcaps", 20000, 0, 20000);
  DielectronMassPUScaleDown_CSPos            = fs->make<TH1F>("DielectronMassPUScaleDown_CSPos",            titlePrefix + "dil. mass for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_bb_CSPos         = fs->make<TH1F>("DielectronMassPUScaleDown_bb_CSPos",            titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_bbbe_CSPos       = fs->make<TH1F>("DielectronMassPUScaleDown_bbbe_CSPos",          titlePrefix + "dil. mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_be_CSPos         = fs->make<TH1F>("DielectronMassPUScaleDown_be_CSPos",            titlePrefix + "dil. mass barrel-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_ee_CSPos         = fs->make<TH1F>("DielectronMassPUScaleDown_ee_CSPos",            titlePrefix + "dil. mass endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_CSNeg            = fs->make<TH1F>("DielectronMassPUScaleDown_CSNeg",            titlePrefix + "dil. mass for negative cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_bb_CSNeg         = fs->make<TH1F>("DielectronMassPUScaleDown_bb_CSNeg",            titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_bbbe_CSNeg       = fs->make<TH1F>("DielectronMassPUScaleDown_bbbe_CSNeg",          titlePrefix + "dil. mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_be_CSNeg         = fs->make<TH1F>("DielectronMassPUScaleDown_be_CSNeg",            titlePrefix + "dil. mass barrel-endcaps for negative cos theta star", 20000, 0, 20000);
  DielectronMassPUScaleDown_ee_CSNeg         = fs->make<TH1F>("DielectronMassPUScaleDown_ee_CSNeg",            titlePrefix + "dil. mass endcaps-endcaps for negative cos theta star", 20000, 0, 20000);

  DielectronMassVsCS            = fs->make<TH2F>("DielectronMassVsCS",            titlePrefix + "dil. mass", 200, 0, 20000,100,-1,1);
  DielectronMassVsCS_bbbe       = fs->make<TH2F>("DielectronMassVsCS_bbbe",       titlePrefix + "dil. mass", 200, 0, 20000,100,-1,1);
  DielectronMassVsCS_bb         = fs->make<TH2F>("DielectronMassVsCS_bb",            titlePrefix + "dil. mass barrel-barrel", 200, 0, 20000,100,-1,1);
  DielectronMassVsCS_be         = fs->make<TH2F>("DielectronMassVsCS_be",            titlePrefix + "dil. mass barrel-endcaps", 200, 0, 20000,100,-1,1);
  DielectronMassVsCS_ee         = fs->make<TH2F>("DielectronMassVsCS_ee",            titlePrefix + "dil. mass endcaps-endcaps", 200, 0, 20000,100,-1,1);
  DielectronResponse_bb         = fs->make<TH2F>("DielectronResponse_bb",            titlePrefix + "dil. mass barrel-barrel", 2000, 0, 20000,2000,0,20000);
  DielectronResponse_be         = fs->make<TH2F>("DielectronResponse_be",            titlePrefix + "dil. mass barrel-endcaps", 2000, 0, 20000,2000,0,20000);
  DielectronResponse_ee         = fs->make<TH2F>("DielectronResponse_ee",            titlePrefix + "dil. mass endcaps-endcaps", 2000, 0, 20000,2000,0,20000);
  DielectronResponseMassScaleUp_bb         = fs->make<TH2F>("DielectronResponseMassScaleUp_bb",            titlePrefix + "dil. mass barrel-barrel", 2000, 0, 20000,2000,0,20000);
  DielectronResponseMassScaleUp_be         = fs->make<TH2F>("DielectronResponseMassScaleUp_be",            titlePrefix + "dil. mass barrel-endcaps", 2000, 0, 20000,2000,0,20000);
  DielectronResponseMassScaleUp_ee         = fs->make<TH2F>("DielectronResponseMassScaleUp_ee",            titlePrefix + "dil. mass endcaps-endcaps", 2000, 0, 20000,2000,0,20000);
 DielectronResponseMassScaleDown_bb         = fs->make<TH2F>("DielectronResponseMassScaleDown_bb",            titlePrefix + "dil. mass barrel-barrel", 2000, 0, 20000,2000,0,20000);
  DielectronResponseMassScaleDown_be         = fs->make<TH2F>("DielectronResponseMassScaleDown_be",            titlePrefix + "dil. mass barrel-endcaps", 2000, 0, 20000,2000,0,20000);
  DielectronResponseMassScaleDown_ee         = fs->make<TH2F>("DielectronResponseMassScaleDown_ee",            titlePrefix + "dil. mass endcaps-endcaps", 2000, 0, 20000,2000,0,20000);
  GenMass                 = fs->make<TH1F>("GenMass",            titlePrefix + "dil. mass", 20000, -10 , 20000);
  
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
  DimuonMassVertexConstrained_gen_bb = fs->make<TH1F>("DimuonMassVertexConstrained_gen_bb", titlePrefix + "dimu. vertex-constrained mass barrel-barrel", 20000, 0, 20000);
  DimuonMassVertexConstrained_gen_be = fs->make<TH1F>("DimuonMassVertexConstrained_gen_be", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps", 20000, 0, 20000);
  DimuonMassVertexConstrained_CSPos = fs->make<TH1F>("DimuonMassVertexConstrained_CSPos", titlePrefix + "dimu. vertex-constrained mass for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_bb_CSPos = fs->make<TH1F>("DimuonMassVertexConstrained_bb_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_be_CSPos = fs->make<TH1F>("DimuonMassVertexConstrained_be_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrained_CSNeg", titlePrefix + "dimu. vertex-constrained mass for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_bb_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrained_bb_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrained_be_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrained_be_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedWeight = fs->make<TH1F>("DimuonMassVertexConstrainedWeight", titlePrefix + "dimu. vertex-constrained mass", 20000, 0, 20000);
 
  DimuonMassVertexConstrainedMuonID = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID", titlePrefix + "dimu. vertex-constrained mass", 20000, 0, 20000);
  DimuonMassVertexConstrainedMuonID_bb = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID_bb", titlePrefix + "dimu. vertex-constrained mass barrel-barrel", 20000, 0, 20000);
  DimuonMassVertexConstrainedMuonID_be = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID_be", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps", 20000, 0, 20000);
  DimuonMassVertexConstrainedMuonID_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID_CSPos", titlePrefix + "dimu. vertex-constrained mass for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedMuonID_bb_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID_bb_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedMuonID_be_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID_be_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedMuonID_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID_CSNeg", titlePrefix + "dimu. vertex-constrained mass for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedMuonID_bb_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID_bb_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedMuonID_be_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedMuonID_be_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for negative cos theta star", 20000, 0, 20000);
 
  DimuonMassVertexConstrainedSmear = fs->make<TH1F>("DimuonMassVertexConstrainedSmear", titlePrefix + "dimu. vertex-constrained mass", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_bb = fs->make<TH1F>("DimuonMassVertexConstrainedSmear_bb", titlePrefix + "dimu. vertex-constrained mass barrel-barrel", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_be = fs->make<TH1F>("DimuonMassVertexConstrainedSmear_be", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedSmear_CSPos", titlePrefix + "dimu. vertex-constrained mass for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_bb_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedSmear_bb_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_be_CSPos = fs->make<TH1F>("DimuonMassVertexConstrainedSmear_be_CSPos", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for positive cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedSmear_CSNeg", titlePrefix + "dimu. vertex-constrained mass for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_bb_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedSmear_bb_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-barrel for negative cos theta star", 20000, 0, 20000);
  DimuonMassVertexConstrainedSmear_be_CSNeg = fs->make<TH1F>("DimuonMassVertexConstrainedSmear_be_CSNeg", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps for negative cos theta star", 20000, 0, 20000);
  
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
  // This binning selection gives binwidths of ~30 GeV at 3 TeV
  // It can be rebinned and clipped as necessary after the fact

  DimuonMassVertexConstrainedVsCS = fs->make<TH2F>("DimuonMassVertexConstrainedVsCS", titlePrefix + "dimu. vertex-constrained mass", 200, 0, 20000,100,-1,1);
  DimuonMassVertexConstrainedVsCS_bb = fs->make<TH2F>("DimuonMassVertexConstrainedVsCS_bb", titlePrefix + "dimu. vertex-constrained mass barrel-barrel", 200, 0, 20000,100,-1,1);
  DimuonMassVertexConstrainedVsCS_be = fs->make<TH2F>("DimuonMassVertexConstrainedVsCS_be", titlePrefix + "dimu. vertex-constrained mass barrel-endcaps and endcaps-endcaps", 200, 0, 20000,100,-1,1);
  DimuonResponse_bb = fs->make<TH2F>("DimuonResponse_bb", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
  DimuonResponse_be = fs->make<TH2F>("DimuonResponse_be", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
  DimuonResponseMassScaleUp_bb = fs->make<TH2F>("DimuonResponseMassScaleUp_bb", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
  DimuonResponseMassScaleUp_be = fs->make<TH2F>("DimuonResponseMassScaleUp_be", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000); 
  DimuonResponseMassScaleDown_bb = fs->make<TH2F>("DimuonResponseMassScaleDown_bb", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
  DimuonResponseMassScaleDown_be = fs->make<TH2F>("DimuonResponseMassScaleDown_be", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
 DimuonResponseMassScaleUpSmear_bb = fs->make<TH2F>("DimuonResponseMassScaleUpSmear_bb", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
 DimuonResponseMassScaleUpSmear_be = fs->make<TH2F>("DimuonResponseMassScaleUpSmear_be", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
 DimuonResponseMassScaleDownSmear_bb = fs->make<TH2F>("DimuonResponseMassScaleDownSmear_bb", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
 DimuonResponseMassScaleDownSmear_be = fs->make<TH2F>("DimuonResponseMassScaleDownSmear_be", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
 DimuonResponseSmear_bb = fs->make<TH2F>("DimuonResponseSmear_bb", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);
 DimuonResponseSmear_be = fs->make<TH2F>("DimuonResponseSmear_be", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000,2000,0,20000);

  const int    NMBINS = 500;
  const double MMIN = 50., MMAX = 10000.;
  double logMbins[NMBINS+1];
  for (int ibin = 0; ibin <= NMBINS; ibin++)
    logMbins[ibin] = pow(10,(log10(MMIN) + (log10(MMAX)-log10(MMIN))*ibin/NMBINS));
  DimuonMassVtxConstrainedLog = fs->make<TH1F>("DimuonMassVtxConstrainedLog", titlePrefix + "dimu vtx-constrained mass in log bins", NMBINS, logMbins);
  DimuonMassVtxConstrainedLog_bb = fs->make<TH1F>("DimuonMassVtxConstrainedLog_bb", titlePrefix + "dimu vtx-constrained mass in log bins barrel-barrel", NMBINS, logMbins);
  DimuonMassVtxConstrainedLog_be = fs->make<TH1F>("DimuonMassVtxConstrainedLog_be", titlePrefix + "dimu vtx-constrained mass in log bins barrel-endcaps", NMBINS, logMbins);
 DimuonMassVtxConstrainedLog_gen_bb = fs->make<TH1F>("DimuonMassVtxConstrainedLog_gen_bb", titlePrefix + "dimu vtx-constrained mass in log bins barrel-barrel", NMBINS, logMbins);
  DimuonMassVtxConstrainedLog_gen_be = fs->make<TH1F>("DimuonMassVtxConstrainedLog_gen_be", titlePrefix + "dimu vtx-constrained mass in log bins barrel-endcaps", NMBINS, logMbins);
  DimuonMassVtxConstrainedLogWeight = fs->make<TH1F>("DimuonMassVtxConstrainedLogWeight", titlePrefix + "dimu vtx-constrained mass in log bins", NMBINS, logMbins);
  DimuonMassConstrainedVsUn = fs->make<TH2F>("DimuonMassConstrainedVsUn", titlePrefix + "dimu. vertex-constrained vs. non-constrained mass", 200, 0, 3000, 200, 0, 3000);
  DimuonMassVertexConstrainedError = fs->make<TH2F>("DimuonMassVertexConstrainedError", titlePrefix + "dimu. vertex-constrained mass error vs. mass", 100, 0, 3000, 100, 0, 400);
    //special
    DimuonMassVtx_chi2 = fs->make<TH1F>("DimuonMassVtx_chi2", titlePrefix + "dimu. vertex #chi^{2}/dof", 300, 0, 30);
    DimuonMassVtx_prob = fs->make<TH1F>("DimuonMassVtx_prob", titlePrefix + "dimu. vertex probability", 100, 0, 1);
    
     //weight
     WeightMadGraph = fs->make<TH1F>("weightperevent", titlePrefix + "weight per event", 4, -2,2);
     WeightLR = fs->make<TH1F>("LR weights for CI", titlePrefix + "LR/RL weights", 200, -2,2);
     kFactorGraph = fs->make<TH1F>("kFactorperevent", titlePrefix + "kFactor per event", 50, 0.4,1.4);
     kFactorGraph_bb = fs->make<TH1F>("kFactorperevent_bb", titlePrefix + "kFactor per event bb", 50, 0.4,1.4);
     kFactorGraph_be = fs->make<TH1F>("kFactorperevent_be", titlePrefix + "kFactor per event be", 50, 0.4,1.4);
}

double Zprime2muHistosFromPAT::ScaleUncert(double mass, bool isBB, int year){

	if (year == 2016){
		if (isBB) return 0.01;
		else return 0.03;
	}
	else if (year == 2017){
		if (isBB) return 1. - (1.00037 -1.70128e-06*mass + 3.09916e-09*pow(mass,2) -1.32999e-12*pow(mass,3) +  2.99434e-16*pow(mass,4) + -2.28065e-20*pow(mass,5));
		else return 1. - (1.00263 - 1.04029e-05*mass + 8.9214e-09*pow(mass,2) -3.4176e-12*pow(mass,3) +  6.07934e-16*pow(mass,4)  -3.73738e-20*pow(mass,5));
	}
	else if (year == 2018){
		if (isBB) return 1. - (0.999032 + 3.36979e-06*mass -3.4122e-09*pow(mass,2) + 1.62541e-12*pow(mass,3)  - 3.12864e-16*pow(mass,4) + 2.18417e-20*pow(mass,5));
		else return 1. - (1.00051 - 2.21167e-06*mass + 2.21347e-09*pow(mass,2) -7.72178e-13*pow(mass,3) +  1.28101e-16*pow(mass,4)  - 8.32675e-21*pow(mass,5));
	}
	return 1.0;

}

double Zprime2muHistosFromPAT::getRecoWeight(double mass, bool isBB, int year_info){

	if (isBB) return 0.99;

	else{
		double a = 0;
		double eff_a = 0;
		double b = 0;
		double eff_b = 0;
		double c = 0;
		double eff_c = 0;
		double d = 0;
		double eff_d = 0;
		double e = 0;
		double eff_e = 0;
		double f = 0;

		double eff_default = 0;
		double eff_syst = 0;

		if (year_info == 2017){

			if (mass <= 450){
				a =  13.39;
				b =  6.696;
				c = -4.855e+06;
				d = -7.431e+06;
				e = -108.8;
				f = -1.138;
				eff_default = a - b * TMath::Exp( -( (mass - c) / d) ) + e * pow(mass,f);
			}
			else{
				eff_a     =  0.3148;
				eff_b     =  0.04447;
				eff_c     =  1.42;
				eff_d     = -5108.;
				eff_e     =  713.5;
				eff_default = eff_a + eff_b * pow(mass,eff_c) * TMath::Exp(- ((mass - eff_d ) / eff_e) );
			}
			if (mass <= 450){
				a =  1.33901e+01;
				b =  6.69687e+00;
				c = -4.85589e+06;
				d = -7.43036e+06;
				e = -1.14263e+02;
				f = -1.15028e+00;
				eff_syst= a - b * TMath::Exp( -( (mass - c) / d) ) + e * pow(mass,f);
			}
			else{
				eff_a     =  3.07958e-01;
				eff_b     =  4.63280e-02;
				eff_c     =  1.35632e+00;
				eff_d     = -5.00475e+03;
				eff_e     =  7.38088e+02;
				eff_syst =  eff_a + eff_b * pow(mass,eff_c) * TMath::Exp(- ((mass - eff_d ) / eff_e) );
			}

			return eff_syst/eff_default;


		}

		else{

			if (mass <= 450){
				a =  13.39;
				b =  6.696;
				c = -4.855e+06;
				d = -7.431e+06;
				e = -108.8;
				f = -1.138;
				eff_default = a - b * TMath::Exp( -( (mass - c) / d) ) + e * pow(mass,f);
			}
			else{
				eff_a     =  0.3148;
				eff_b     =  0.04447;
				eff_c     =  1.42;
				eff_d     = -5108.;
				eff_e     =  713.5;
				eff_default = eff_a + eff_b * pow(mass,eff_c) * TMath::Exp(- ((mass - eff_d ) / eff_e) );
			}
			if (mass <= 450){
				a =  1.33901e+01;
				b =  6.69687e+00;
				c = -4.85589e+06;
				d = -7.43036e+06;
				e = -1.14263e+02;
				f = -1.15028e+00;
				eff_syst = a - b * TMath::Exp( -( (mass - c) / d) ) + e * pow(mass,f);
			}
			else{
				eff_a     =  3.07958e-01;
				eff_b     =  4.63280e-02;
				eff_c     =  1.35632e+00;
				eff_d     = -5.00475e+03;
				eff_e     =  7.38088e+02;
				eff_syst =  eff_a + eff_b * pow(mass,eff_c) * TMath::Exp(- ((mass - eff_d ) / eff_e) );

			}

		        return eff_syst/eff_default;


		}


	}

	


}

double Zprime2muHistosFromPAT::HEEPSF(double eta,  int year){
	double result = 1;
	if (year == 2017){
		if (fabs(eta) < 1.4442){
			result = 0.968;
		}
		else if (fabs(eta) < 2.5){
			result = 0.969;
		}
	}
	if (year == 2016){
		if (fabs(eta) < 1.4442){
			result = 0.972;
		}
		else if (fabs(eta) < 2.5){
			result = 0.982;
		}
	}
	return result;

}

double Zprime2muHistosFromPAT::L1TurnOn(double eta, double et, int year){

	double result = 1;
	double P0 = 0;
	double P1 = 0;
	double P2 = 0;
	double P3 = 0;
	double P4 = 0;
	double P5 = 0;
	if (year == 2017){
		if (fabs(eta) < 1.4442){
			P0 = 0.745;
			P1 = 35.3;
			P2 = 3.33;
			P3 = 0.584;
			P4 = 115.;
			P5 = 169.;
		}
		else if (fabs(eta) < 2.5){
			P0 = 0.875;
			P1 = 35.2;
			P2 = 4.45;
			P3 = 0.086;
			P4 = 41.1;
			P5 = 11.5;
		}
		result = std::min(1.,0.5*P0*(1 + TMath::Erf((et-P1)/(pow(2,0.5)*P2))) + 0.5*P3*(1 + TMath::Erf( (et-P4)/(pow(2,0.5)*P5))));
	}
/*	if (year == 2016){
		if (fabs(eta) < 0.79){
			P0 = 0.737;
			P1 = 30.4;
			P2 = 1.51;
			P3 = 0.243;
			P4 = 38.2;
			P5 = 4.31;
		}
		else if (fabs(eta) < 1.10){
			P0 = 0.796;
			P1 = 30.2;
			P2 = 1.69;
			P3 = 0.185;
			P4 = 39.2;
			P5 = 3.47;
		}
		else if (fabs(eta) < 1.442){
			P0 = 0.748;
			P1 = 29.9;
			P2 = 1.8;
			P3 = 0.236;
			P4 = 38;
			P5 = 4.36;
		}
		else if (fabs(eta) < 1.70){
			P0 = 0.704;
			P1 = 28.9;
			P2 = 2.55;
			P3 = 0.287;
			P4 = 36.2;
			P5 = 5.04;
		}
		else if (fabs(eta) < 2.10){
			P0 = 0.501;
			P1 = 30.1;
			P2 = 2.27;
			P3 = 0.483;
			P4 = 34.5;
			P5 = 6.08;
		}
		else if (fabs(eta) < 2.5){
			P0 = 0.388;
			P1 = 30.2;
			P2 = 2.06;
			P3 = 0.544;
			P4 = 34.3;
			P5 = 6.91;
		}
		result = std::min(1.,0.5*P0*(1 + TMath::Erf((et-P1)/(pow(2,0.5)*P2))) + 0.5*P3*(1 + TMath::Erf( (et-P4)/(pow(2,0.5)*P5))));
	}*/
	return result;

}

double Zprime2muHistosFromPAT::turnOn(double eta, double et, int year){

	double result = 1.;
	double P0 = 0;
	double P1 = 0;
	double P2 = 0;
	double P3 = 0;
	double P4 = 0;
	double P5 = 0;
	if (year == 2016){
		if (fabs(eta) < 0.79){
			P0 = 0.099;
			P1 = 32.902;
			P2 = 1.921;
			P3 = 0.900;
			P4 = 33.034;
			P5 = 0.625;
		}
		else if (fabs(eta) < 1.1){
			P0 = 0.868;
			P1 = 33.229;
			P2 = 0.706;
			P3 = 0.132;
			P4 = 33.328;
			P5 = 1.777;
		}
		else if (fabs(eta) < 1.4442){
			P0 = 0.231;
			P1 = 33.331;
			P2 = 1.534;
			P3 = 0.769;
			P4 = 33.347;
			P5 = 0.718;
		}
		else if (fabs(eta) < 1.7){
			P0 = 0.189;
			P1 = 32.638;
			P2 = 2.063;
			P3 = 0.808;
			P4 = 33.047;
			P5 = 0.844;
		}
		else if (fabs(eta) < 2.1){
			P0 = 0.362;
			P1 = 33.510;
			P2 = 1.669;
			P3 = 0.637;
			P4 = 33.264;
			P5 = 0.861;
		}
		else if (fabs(eta) < 2.5){
			P0 = 0.536;
			P1 = 34.688;
			P2 = 1.771;
			P3 = 0.462;
			P4 = 34.155;
			P5 = 1.048;
		}
	}

	if (year == 2017){
		if (fabs(eta) < 0.79){
			P0 = 0.8617;
			P1 = 33.84;
			P2 = 0.4828;
			P3 = 0.1381;
			P4 = 34.6;
			P5 = 1.758;
		}
		else if (fabs(eta) < 1.1){
			P0 = 0.9257;
			P1 = 34.04;
			P2 = 0.6257;
			P3 = 0.07419;
			P4 = 34.24;
			P5 = 2.363;
		}
		else if (fabs(eta) < 1.4442){
			P0 = 0.9283;
			P1 = 34.36;
			P2 = 0.805;
			P3 = 0.07154;
			P4 = 33.73;
			P5 = 2.604;
		}
		else if (fabs(eta) < 1.7){
			P0 = 0.525;
			P1 = 34.38;
			P2 = 0.7307;
			P3 = 0.4732;
			P4 = 35.5;
			P5 = 2.053;
		}
		else if (fabs(eta) < 2.1){
			P0 = 0.4136;
			P1 = 35.88;
			P2 = 1.956;
			P3 = 0.5853;
			P4 = 34.66;
			P5 = 0.7886;
		}
		else if (fabs(eta) < 2.5){
			P0 = 0.5594;
			P1 = 34.44;
			P2 = 1.025;
			P3 = 0.4383;
			P4 = 36.53;
			P5 = 2.355;
		}
	}
	if (year == 2018){
		if (fabs(eta) < 0.79){
			P0 = 0.9667;
			P1 = 25.72;
			P2 = 0.4265;
			P3 = 0.03306;
			P4 = 25.13;
			P5 = 2.485;
		}
		else if (fabs(eta) < 1.1){
			P0 = 0.9518;
			P1 = 25.73;
			P2 = 0.5064;
			P3 = 0.04795;
			P4 = 24.98;
			P5 = 2.502;
		}
		else if (fabs(eta) < 1.4442){
			P0 = 0.9415;
			P1 = 26.05;
			P2 = 0.6502;
			P3 = 0.05819;
			P4 = 25.36;
			P5 = 2.36;
		}
		else if (fabs(eta) < 1.7){
			P0 = 0.8774;
			P1 = 25.92;
			P2 = 0.7452;
			P3 = 0.1219;
			P4 = 25.67;
			P5 = 2.891;
		}
		else if (fabs(eta) < 2.1){
			P0 = 0.8905;
			P1 = 26.05;
			P2 = 0.701;
			P3 = 0.1092;
			P4 = 26.12;
			P5 = 2.077;
		}
		else if (fabs(eta) < 2.5){
			P0 = 0.6658;
			P1 = 26.4;
			P2 = 0.8168;
			P3 = 0.3339;
			P4 = 27.1;
			P5 = 1.692;
		}
	}

	result = std::min(1.,0.5*P0*(1 + TMath::Erf((et-P1)/(pow(2,0.5)*P2))) + 0.5*P3*(1 + TMath::Erf( (et-P4)/(pow(2,0.5)*P5))));
	return result;
}	

double Zprime2muHistosFromPAT::getSmearedMass(const pat::CompositeCandidate& dil, double gM, int year, bool syst){

    double a = 0.;
    double b = 0.;
    double c = 0.;
    double d = 0.;
    double e = 0.;

    double mass = dil.userFloat("vertexM");

   if (year == 2016){
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
    }
   if (year == 2017){
	   //sigma BB
	   if (dil.daughter(0)->eta()<=1.2 && dil.daughter(1)->eta()<=1.2 && dil.daughter(0)->eta()>=-1.2 && dil.daughter(1)->eta()>=-1.2){
		a=0.00606;
		b=3.41e-05;
		c=-1.33e-08;
		d=2.39e-12;
		e=-1.5e-16;
	    }
	    else{
	   //sigma BE
		a=0.0108;
		b=3.25e-05;
		c=-1.18e-08;
		d=2.11e-12;
		e=-1.35e-16;
	    }
    }
   if (year == 2018){
	   //sigma BB
	   if (dil.daughter(0)->eta()<=1.2 && dil.daughter(1)->eta()<=1.2 && dil.daughter(0)->eta()>=-1.2 && dil.daughter(1)->eta()>=-1.2){
		a=0.00608;
		b=3.42e-05;
		c=-1.34e-08;
		d=2.4e-12;
		e=-1.5e-16;
	    }
	    else{
	   //sigma BE
		a=0.0135;
		b=2.83e-05;
		c=-9.71e-09;
		d=1.71e-12;
		e=-1.09e-16;
	    }
    }
    double res = a + b*gM + c*gM*gM + d*pow(gM,3) + e*pow(gM,4);

    double extraSmear = res*0.567;
    double extraSmearSyst = res*0.567;
    if (syst && (year == 2017 || year ==2018)) extraSmearSyst = res*0.42098099719583537;

	
    TRandom3 *rand = new TRandom3(0);
    double result = mass;
    if (!(fabs(dil.daughter(0)->eta())<=1.2 && fabs(dil.daughter(1)->eta())<=1.2)){
	result = mass*rand->Gaus(1,extraSmear);
    	if (syst) result = result*rand->Gaus(1,extraSmearSyst);
    }
    delete rand;
    return result;
}


void Zprime2muHistosFromPAT::getBSandPV(const edm::Event& event) {
  // We store these as bare pointers. Should find better way, but
  // don't want to pass them around everywhere...
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  event.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;

  if(!(event.isRealData())){
  	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	
		int BX = PVI->getBunchCrossing();
	
		if(BX == 0) {
			_nTrueInt = PVI->getTrueNumInteractions(); 
	  		NTrueInteractions->Fill(_nTrueInt);
	  		continue;
		}
      	}
  }
  if (size(pu_info)> 0){
	if (year_info == 2016){
		_puWeight = PU_2016::MC_pileup_weight(_nTrueInt,0,"all");
		_puWeight_scaleUp = PU_2016::MC_pileup_weight(_nTrueInt,1,"all");	
		_puWeight_scaleDown = PU_2016::MC_pileup_weight(_nTrueInt,2,"all");
	}
	else if (year_info == 2017){
		_puWeight = PU_2017::MC_pileup_weight(_nTrueInt,pu_info[0],pu_info[1]);
		_puWeight_scaleUp = PU_2017::MC_pileup_weight(_nTrueInt,pu_info[0],pu_info[1]+TString("_scaleUp"));	
		_puWeight_scaleDown = PU_2017::MC_pileup_weight(_nTrueInt,pu_info[0],pu_info[1]+TString("_scaleDown"));
	}
	else if (year_info == 2018){
		_puWeight = PU_2018::MC_pileup_weight(_nTrueInt,pu_info[0],pu_info[1]);
		_puWeight_scaleUp = PU_2018::MC_pileup_weight(_nTrueInt,pu_info[0],pu_info[1]+TString("_scaleUp"));	
		_puWeight_scaleDown = PU_2018::MC_pileup_weight(_nTrueInt,pu_info[0],pu_info[1]+TString("_scaleDown"));
	}
  }
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
  NVertices->Fill(vertex_count, _madgraphWeight*_kFactor*_puWeight);
  NVerticesUnweighted->Fill(vertex_count, _madgraphWeight*_kFactor);


}

void Zprime2muHistosFromPAT::fillBasicLeptonHistos(const reco::CandidateBaseRef& lep) {
  LeptonEta->Fill(lep->eta(), _madgraphWeight*_kFactor*_puWeight);
  LeptonRap->Fill(lep->rapidity(), _madgraphWeight*_kFactor*_puWeight);
  LeptonPhi->Fill(lep->phi(), _madgraphWeight*_kFactor*_puWeight);

  LeptonPt->Fill(lep->pt(), _madgraphWeight*_kFactor*_puWeight);
  LeptonPz->Fill(fabs(lep->pz()), _madgraphWeight*_kFactor*_puWeight);
  LeptonP ->Fill(lep->p(), _madgraphWeight*_kFactor*_puWeight);

  LeptonPtVsEta->Fill(lep->eta(), lep->pt(), _madgraphWeight*_kFactor*_puWeight);
  LeptonPVsEta ->Fill(lep->eta(), lep->p(), _madgraphWeight*_kFactor*_puWeight);
}

void Zprime2muHistosFromPAT::fillOfflineMuonHistos(const pat::Muon* mu) {
  const reco::MuonIsolation& iso = mu->isolationR03();
  IsoSumPt   ->Fill(iso.sumPt, _madgraphWeight*_kFactor*_puWeight);
  RelIsoSumPt->Fill(iso.sumPt / mu->innerTrack()->pt(), _madgraphWeight*_kFactor*_puWeight);
  IsoEcal    ->Fill(iso.emEt, _madgraphWeight*_kFactor*_puWeight);
  IsoHcal    ->Fill(iso.hadEt + iso.hoEt, _madgraphWeight*_kFactor*_puWeight);
  CombIso    ->Fill( iso.sumPt + iso.emEt + iso.hadEt + iso.hoEt, _madgraphWeight*_kFactor*_puWeight);
  RelCombIso ->Fill((iso.sumPt + iso.emEt + iso.hadEt + iso.hoEt) / mu->innerTrack()->pt(), _madgraphWeight*_kFactor*_puWeight);
  IsoNTracks ->Fill(iso.nTracks, _madgraphWeight*_kFactor*_puWeight);
  IsoNJets   ->Fill(iso.nJets, _madgraphWeight*_kFactor*_puWeight);

  CombIsoNoECAL   ->Fill( iso.sumPt + iso.hadEt + iso.hoEt, _madgraphWeight*_kFactor*_puWeight);
  RelCombIsoNoECAL->Fill((iso.sumPt + iso.hadEt + iso.hoEt) / mu->innerTrack()->pt(), _madgraphWeight*_kFactor*_puWeight);

  const reco::TrackRef track = patmuon::getPickedTrack(*mu);
  if (track.isAvailable()) {
    Chi2dof->Fill(track->normalizedChi2(), _madgraphWeight*_kFactor*_puWeight);

    if (beamspot != 0) {
      TrackD0BS->Fill(fabs(track->dxy(beamspot->position())), _madgraphWeight*_kFactor*_puWeight);
      TrackDZBS->Fill(fabs(track->dz (beamspot->position())), _madgraphWeight*_kFactor*_puWeight);
    }

    if (vertex != 0) {
      TrackD0PV->Fill(fabs(track->dxy(vertex->position())), _madgraphWeight*_kFactor*_puWeight);
      TrackDZPV->Fill(fabs(track->dz (vertex->position())), _madgraphWeight*_kFactor*_puWeight);
    }

    const reco::HitPattern& hp = track->hitPattern();
    NPxHits->Fill(hp.numberOfValidPixelHits(), _madgraphWeight*_kFactor*_puWeight);
    NStHits->Fill(hp.numberOfValidStripHits(), _madgraphWeight*_kFactor*_puWeight);
    NTkHits->Fill(hp.numberOfValidTrackerHits(), _madgraphWeight*_kFactor*_puWeight);
    NMuHits->Fill(hp.numberOfValidMuonHits(), _madgraphWeight*_kFactor*_puWeight);

    NHits->Fill(hp.numberOfValidHits(), _madgraphWeight*_kFactor*_puWeight);
    NInvalidHits->Fill(hp.numberOfAllHits(reco::HitPattern::TRACK_HITS) - hp.numberOfValidHits(), _madgraphWeight*_kFactor*_puWeight);
    //NInvalidHits->Fill(hp.numberOfAllHits() - hp.numberOfValidHits());
    
    NPxLayers->Fill(hp.pixelLayersWithMeasurement(), _madgraphWeight*_kFactor*_puWeight);
    NStLayers->Fill(hp.stripLayersWithMeasurement(), _madgraphWeight*_kFactor*_puWeight);
    NTkLayers->Fill(hp.trackerLayersWithMeasurement(), _madgraphWeight*_kFactor*_puWeight);
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
  NLeptons->Fill(nleptons, _madgraphWeight*_kFactor*_puWeight);
  LeptonSigns->Fill(nleptons, total_q, _madgraphWeight*_kFactor*_puWeight);
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
  DileptonEta->Fill(dil.eta(), _madgraphWeight*_kFactor*_puWeight);
  DileptonRap->Fill(dil.rapidity(), _madgraphWeight*_kFactor*_puWeight);
  DileptonPhi->Fill(dil.phi(), _madgraphWeight*_kFactor*_puWeight);

  DileptonPt->Fill(dil.pt(), _madgraphWeight*_kFactor*_puWeight);
  DileptonPz->Fill(fabs(dil.pz()), _madgraphWeight*_kFactor*_puWeight);
  DileptonP ->Fill(dil.p(), _madgraphWeight*_kFactor*_puWeight);

  DileptonPtVsEta->Fill(dil.eta(), dil.pt(), _madgraphWeight*_kFactor*_puWeight);
  DileptonPVsEta ->Fill(dil.eta(), dil.p(), _madgraphWeight*_kFactor*_puWeight);

  DileptonMass->Fill(dil.mass(), _madgraphWeight*_kFactor*_puWeight);
  DileptonMassWeight->Fill(dil.mass(),_prescaleWeight*_madgraphWeight*_kFactor*_puWeight);//?
  DileptonWithPhotonsMass->Fill(resonanceP4(dil).mass(), _madgraphWeight*_kFactor*_puWeight);

  const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
  const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);
  double cos_cs = -999.;
  if (lep0.isNonnull() && lep1.isNonnull()) {
    DileptonDeltaPt->Fill(fabs(lep0->pt()) - fabs(lep1->pt()), _madgraphWeight*_kFactor*_puWeight);
    DileptonDeltaP ->Fill(fabs(lep0->p())  - fabs(lep1->p()), _madgraphWeight*_kFactor*_puWeight);
     if (lep0->charge()*lep1->charge() == -1){
	     if (lep0->charge() == -1) cos_cs = calcCosThetaCSAnal(lep0->pz(), lep0->energy(), lep1->pz(), lep1->energy(), dil.pt(), dil.pz(), dil.mass());
	     else cos_cs = calcCosThetaCSAnal(lep1->pz(), lep1->energy(), lep0->pz(), lep0->energy(), dil.pt(), dil.pz(), dil.mass());
     }
     else{
	if (fabs(lep0->eta()) < fabs(lep1->eta())){
      		if (lep0->charge() == -1) cos_cs = calcCosThetaCSAnal(lep0->pz(), lep0->energy(), lep1->pz(), lep1->energy(), dil.pt(), dil.pz(), dil.mass());
	     	else cos_cs = calcCosThetaCSAnal(lep1->pz(), lep1->energy(), lep0->pz(), lep0->energy(), dil.pt(), dil.pz(), dil.mass());
	}
	else{ 
   		if (lep1->charge() == -1) cos_cs = calcCosThetaCSAnal(lep1->pz(), lep1->energy(), lep0->pz(), lep0->energy(), dil.pt(), dil.pz(), dil.mass());
	     	else cos_cs = calcCosThetaCSAnal(lep0->pz(), lep0->energy(), lep1->pz(), lep1->energy(), dil.pt(), dil.pz(), dil.mass());
	}
     }
     CosThetaStarDilepton->Fill(cos_cs);
     //ChiDilepton->Fill((1+fabs(cos_cs))/(1-fabs(cos_cs)));
     ChiDilepton->Fill(exp(std::abs(lep0->p4().Rapidity()-lep1->p4().Rapidity())));
     if (cos_cs >= 0) DileptonMass_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_puWeight);
     else DileptonMass_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_puWeight);
    const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
    const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
    if (mu0 && mu1) {
      const reco::Track* tk0 = patmuon::getPickedTrack(*mu0).get();
      const reco::Track* tk1 = patmuon::getPickedTrack(*mu1).get();
      if (tk0 && tk1) {
	DimuonMuonPtErrors->Fill(ptError(tk0), ptError(tk1), _madgraphWeight*_kFactor*_puWeight);
	DimuonMuonPtErrOverPt->Fill(ptError(tk0)/tk0->pt(), _madgraphWeight*_kFactor*_puWeight);
	DimuonMuonPtErrOverPt->Fill(ptError(tk1)/tk1->pt(), _madgraphWeight*_kFactor*_puWeight);
	float mass = -999.;
	// Use mass calculated with the vertex constraint when available
	if (dil.hasUserFloat("vertexM"))
	  mass = dil.userFloat("vertexM");
	else{
	  edm::LogWarning("fillDileptonHistos") << "+++ Mass calculated without vertex constraint! +++";
	  mass = dil.mass();
	}
	if (mass > 200.) {
	  DimuonMuonPtErrOverPtM200->Fill(ptError(tk0)/tk0->pt(), _madgraphWeight*_kFactor*_puWeight);
	  DimuonMuonPtErrOverPtM200->Fill(ptError(tk1)/tk1->pt(), _madgraphWeight*_kFactor*_puWeight);
	}
	if (mass > 500.) {
	  DimuonMuonPtErrOverPtM500->Fill(ptError(tk0)/tk0->pt(), _madgraphWeight*_kFactor*_puWeight);
	  DimuonMuonPtErrOverPtM500->Fill(ptError(tk1)/tk1->pt(), _madgraphWeight*_kFactor*_puWeight);
	}
      }
    } 

    const pat::Electron* ele0 = toConcretePtr<pat::Electron>(lep0);
    const pat::Electron* ele1 = toConcretePtr<pat::Electron>(lep1);
    if (ele0 && ele1) {
  	_eleMCFac = 1;
	if (fill_gen_info){
		bool e1_pass_trigger = true;
		bool e2_pass_trigger = true;
		bool e1_pass_l1 = true;
		bool e2_pass_l1 = true;
		if     (year_info==2016)e1_pass_trigger=trigEle_2016::passTrig(ele0->et(), ele0->superCluster()->eta()) ;
		else if(year_info==2017){
			e1_pass_trigger=trigEle_2017::passTrig(ele0->et(), ele0->superCluster()->eta());
			e1_pass_l1 = trigEle33l1::passTrig(ele0->et(),ele0->superCluster()->eta(),"Run_all", false);
		}
		else if(year_info==2018)e1_pass_trigger=trigEle_2018::passTrig(ele0->et(), ele0->superCluster()->eta(), "Run_all" , true) ;
		if     (year_info==2016)e2_pass_trigger=trigEle_2016::passTrig(ele1->et(), ele1->superCluster()->eta()) ;
		else if(year_info==2017){
			e2_pass_trigger=trigEle_2017::passTrig(ele1->et(), ele1->superCluster()->eta());
			e2_pass_l1 = trigEle33l1::passTrig(ele1->et(),ele1->superCluster()->eta(),"Run_all", false);
		}
		else if(year_info==2018)e2_pass_trigger=trigEle_2018::passTrig(ele1->et(), ele1->superCluster()->eta(), "Run_all" , true) ;
	
		_eleMCFac = _prefireWeight;
		if (!(e1_pass_trigger && e2_pass_trigger)) _eleMCFac = 0;
		if (!(e1_pass_l1 || e2_pass_l1)) _eleMCFac = 0;
		if (!(ele0->userFloat("genMatch") && ele1->userFloat("genMatch") )) _eleMCFac = 0;
		//std::cout << ele0->userFloat("genMatch") << " " << ele1->userFloat("genMatch") << std::endl;
	}
	if (_eleMCFac != 0){
   		GenMass->Fill(gM); 
		double massScaleUp = 1.;
		double massScaleDown = 1.;
		double ttFac = 1.;
		
		if (_useTTBarWeight){
                 edm::Handle<LHEEventProduct> lheInfoHandle;
                 event.getByToken(LHEEventToken_ , lheInfoHandle);
		if (lheInfoHandle.isValid()) {
			lhef::HEPEUP lheParticleInfo = lheInfoHandle->hepeup();
			std::vector<lhef::HEPEUP::FiveVector> allParticles = lheParticleInfo.PUP;
			std::vector<int> statusCodes = lheParticleInfo.ISTUP;
			for (unsigned int i = 0; i < statusCodes.size(); i++) {
				if (statusCodes[i] == 1) {
					if (abs(lheParticleInfo.IDUP[i]) == 11 || abs(lheParticleInfo.IDUP[i]) == 13 ||  abs(lheParticleInfo.IDUP[i]) == 15) {
						std::cout << sqrt(pow(allParticles[i][0], 2) + pow(allParticles[i][1], 2)) << " "  << abs(lheParticleInfo.IDUP[i]) <<  std::endl;
					}
				}
			}
			}
			std::string ttFileName = std::string(std::getenv("CMSSW_BASE")) + std::string("/src/SUSYBSMAnalysis/Zprime2muAnalysis/data/xgao_ttbar_ratio.root");

			TFile *f_NNPDF = new TFile(ttFileName.c_str(),"read");
			TH1D  *h_tt_BB = (TH1D*)f_NNPDF->Get("gen_M_emu_BB_log");
			TH1D  *h_tt_BE = (TH1D*)f_NNPDF->Get("gen_M_emu_BE_log");
			TAxis *tt_BB_axis = h_tt_BB->GetXaxis();
			TAxis *tt_BE_axis = h_tt_BE->GetXaxis();
			float tt_BB_axis_min = tt_BB_axis->GetXmin();
			float tt_BB_axis_max = tt_BB_axis->GetXmax();
			float tt_BE_axis_min = tt_BE_axis->GetXmin();
			float tt_BE_axis_max = tt_BE_axis->GetXmax();
			//std::cout << gM << std::endl;
			if(gM>500){
				if(fabs(ele0->superCluster()->eta()) < 1.4442 && fabs(ele1->superCluster()->eta()) < 1.4442){
					Int_t bin_tt = tt_BB_axis->FindBin(gM);
					if(gM > tt_BB_axis_max)      bin_tt=tt_BB_axis->FindBin(tt_BB_axis_max-0.01);
					else if(gM < tt_BB_axis_min) bin_tt=tt_BB_axis->FindBin(tt_BB_axis_min+0.01);
					ttFac = h_tt_BB->GetBinContent(bin_tt);// get the ratio of NNPDF3.0 to NNPDF3.1
				}
				else{
					Int_t bin_tt = tt_BE_axis->FindBin(gM);
					if(gM > tt_BE_axis_max)      bin_tt=tt_BE_axis->FindBin(tt_BE_axis_max-0.01);
					else if(gM < tt_BE_axis_min) bin_tt=tt_BE_axis->FindBin(tt_BE_axis_min+0.01);
					ttFac = h_tt_BE->GetBinContent(bin_tt);// get the ratio of NNPDF3.0 to NNPDF3.1
				}
			}	
			//std::cout << ttFac << std::endl;
			f_NNPDF->Close();
		}
		if (fabs(ele0->superCluster()->eta()) < 1.4442 && fabs(ele1->superCluster()->eta()) < 1.4442) {
			_kFactor = _kFactor_bb*ttFac;
			massScaleUp = 1+_scaleUncertEleBB;
			massScaleDown = 1-_scaleUncertEleBB;
		}
		else{
			_kFactor = _kFactor_be*ttFac;
			massScaleUp = 1+_scaleUncertEleBE;
			massScaleDown = 1-_scaleUncertEleBE;
		}
		//std::cout << _madgraphWeight << " " << _kFactor << " " << _eleMCFac << " " << _puWeight << std::endl;	
		DielectronMass->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
		DielectronMassVsCS->Fill(dil.mass(),cos_cs, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
		DielectronMassScaleUp->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
		DielectronMassScaleDown->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
		DielectronMassPrefireUp->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
		DielectronMassPrefireDown->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
		DielectronMassPUScaleUp->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
		DielectronMassPUScaleDown->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
		if (cos_cs >= 0){
			 DielectronMass_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			 DielectronMassScaleUp_CSPos->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			 DielectronMassScaleDown_CSPos->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			 DielectronMassPrefireUp_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
			 DielectronMassPrefireDown_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
			 DielectronMassPUScaleUp_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
			 DielectronMassPUScaleDown_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
		}
		else{
			 DielectronMass_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			 DielectronMassScaleUp_CSNeg->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			 DielectronMassScaleDown_CSNeg->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			 DielectronMassPrefireUp_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
			 DielectronMassPrefireDown_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
			 DielectronMassPUScaleUp_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
			 DielectronMassPUScaleDown_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
		}

		if (fabs(ele0->superCluster()->eta()) < 1.4442 && fabs(ele1->superCluster()->eta()) < 1.4442) {
			//std::cout << "filling BB " << _madgraphWeight << " " << _kFactor << " " <<  _eleMCFac << " "  << _puWeight << " " << gM <<  std::endl;
			DielectronMass_bb->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMass_gen_bb->Fill(gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronResponse_bb->Fill(dil.mass(),gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassVsCS_bb->Fill(dil.mass(),cos_cs, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleUp_bb->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleDown_bb->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassPrefireUp_bb->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
			DielectronMassPrefireDown_bb->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
			DielectronMassPUScaleUp_bb->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
			DielectronMassPUScaleDown_bb->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			DielectronMass_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassVsCS_bbbe->Fill(dil.mass(),cos_cs, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleUp_bbbe->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleDown_bbbe->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassPrefireUp_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
			DielectronMassPrefireDown_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
			DielectronMassPUScaleUp_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
			DielectronMassPUScaleDown_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			DielectronResponseMassScaleUp_bb->Fill(dil.mass()*massScaleUp,gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronResponseMassScaleDown_bb->Fill(dil.mass()*massScaleDown,gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			if (cos_cs >= 0){
				 DielectronMass_bb_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleUp_bb_CSPos->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleDown_bb_CSPos->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassPrefireUp_bb_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				 DielectronMassPrefireDown_bb_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				 DielectronMassPUScaleUp_bb_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				 DielectronMassPUScaleDown_bb_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);

				 DielectronMass_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleUp_bbbe_CSPos->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleDown_bbbe_CSPos->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassPrefireUp_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				 DielectronMassPrefireDown_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				 DielectronMassPUScaleDown_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				 DielectronMassPUScaleUp_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			}
			else {
				DielectronMass_bb_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleUp_bb_CSNeg->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleDown_bb_CSNeg->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassPrefireUp_bb_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				DielectronMassPrefireDown_bb_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				DielectronMassPUScaleUp_bb_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				DielectronMassPUScaleDown_bb_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);

				DielectronMass_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleUp_bbbe_CSNeg->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleDown_bbbe_CSNeg->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassPrefireUp_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				DielectronMassPrefireDown_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				DielectronMassPUScaleUp_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				DielectronMassPUScaleDown_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			}

		}
		else if ((fabs(ele0->superCluster()->eta()) < 1.4442 && fabs(ele1->superCluster()->eta()) > 1.566) ||(fabs(ele0->superCluster()->eta()) > 1.566 && fabs(ele1->superCluster()->eta()) < 1.4442)) {
			//std::cout << "filling BE " << _madgraphWeight << " " << _kFactor << " " <<  _eleMCFac << " "  << _puWeight << " " << gM <<  std::endl;
			DielectronMass_be->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMass_gen_be->Fill(gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
                        DielectronResponse_be->Fill(dil.mass(),gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassVsCS_be->Fill(dil.mass(),cos_cs, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleUp_be->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleDown_be->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassPrefireUp_be->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
			DielectronMassPrefireDown_be->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
			DielectronMassPUScaleUp_be->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
			DielectronMassPUScaleDown_be->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);

			DielectronMass_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassVsCS_bbbe->Fill(dil.mass(),cos_cs, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleUp_bbbe->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleDown_bbbe->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassPrefireUp_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
			DielectronMassPrefireDown_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
			DielectronMassPUScaleUp_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
			DielectronMassPUScaleDown_bbbe->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			DielectronResponseMassScaleUp_be->Fill(dil.mass()*massScaleUp,gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
                        DielectronResponseMassScaleDown_be->Fill(dil.mass()*massScaleDown,gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			if (cos_cs >= 0) {
				DielectronMass_be_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleUp_be_CSPos->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleDown_be_CSPos->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassPrefireUp_be_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				DielectronMassPrefireDown_be_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				DielectronMassPUScaleUp_be_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				DielectronMassPUScaleDown_be_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);

				DielectronMass_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleUp_bbbe_CSPos->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleDown_bbbe_CSPos->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassPrefireUp_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				DielectronMassPrefireDown_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				DielectronMassPUScaleUp_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				DielectronMassPUScaleDown_bbbe_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			}
			else{
				 DielectronMass_be_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleUp_be_CSNeg->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleDown_be_CSNeg->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
	 			 DielectronMassPrefireUp_be_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				 DielectronMassPrefireDown_be_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				 DielectronMassPUScaleUp_be_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				 DielectronMassPUScaleDown_be_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);

				 DielectronMass_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleUp_bbbe_CSNeg->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleDown_bbbe_CSNeg->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassPrefireUp_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				 DielectronMassPrefireDown_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				 DielectronMassPUScaleUp_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				 DielectronMassPUScaleDown_bbbe_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			}

		}
		else if (fabs(ele0->superCluster()->eta()) > 1.566 && fabs(ele1->superCluster()->eta()) > 1.566) {
			DielectronMass_ee->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMass_gen_ee->Fill(gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
                        DielectronResponse_ee->Fill(dil.mass(),gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassVsCS_ee->Fill(dil.mass(),cos_cs, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleUp_ee->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassScaleDown_ee->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			DielectronMassPrefireUp_ee->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
			DielectronMassPrefireDown_ee->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
			DielectronMassPUScaleUp_ee->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
			DielectronMassPUScaleDown_ee->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			DielectronResponseMassScaleUp_ee->Fill(dil.mass()*massScaleUp,gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
                        DielectronResponseMassScaleDown_ee->Fill(dil.mass()*massScaleDown,gM, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
			if (cos_cs >= 0) {
				DielectronMass_ee_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleUp_ee_CSPos->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassScaleDown_ee_CSPos->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				DielectronMassPrefireUp_ee_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				DielectronMassPrefireDown_ee_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				DielectronMassPUScaleUp_ee_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				DielectronMassPUScaleDown_ee_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			}
			else
			{ 
				 DielectronMass_ee_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleUp_ee_CSNeg->Fill(dil.mass()*massScaleUp, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassScaleDown_ee_CSNeg->Fill(dil.mass()*massScaleDown, _madgraphWeight*_kFactor*_eleMCFac*_puWeight);
				 DielectronMassPrefireUp_ee_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightUp*_puWeight);
				 DielectronMassPrefireDown_ee_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_prefireWeightDown*_puWeight);
				 DielectronMassPUScaleUp_ee_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleUp);
				 DielectronMassPUScaleDown_ee_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_eleMCFac*_puWeight_scaleDown);
			}

		}

	}

    }

  }

  DileptonDaughterIds->Fill(dil.daughter(0)->pdgId(), dil.daughter(1)->pdgId(), _madgraphWeight*_kFactor*_puWeight);

  DileptonDaughterDeltaR->Fill(reco::deltaR(*dil.daughter(0), *dil.daughter(1)), _madgraphWeight*_kFactor*_puWeight);
  DileptonDaughterDeltaPhi->Fill(reco::deltaPhi(dil.daughter(0)->phi(), dil.daughter(1)->phi()), _madgraphWeight*_kFactor*_puWeight);

  if (dil.hasUserFloat("vertexM") && dil.hasUserFloat("vertexMError")) {

   GenMass->Fill(gM); 
    float vertex_mass = dil.userFloat("vertexM");
    float vertex_mass_err = dil.userFloat("vertexMError");
 
    double SF1 = 1.0;
    double SF2 = 1.0;
    if (fabs(dil.daughter(0)->eta()) <= 1.6 && dil.daughter(0)->p() > 100)
    	SF1 = (0.994 - 4.08e-6 * dil.daughter(0)->p())/(0.994 - 4.08e-6 * 100);
    else if (fabs(dil.daughter(0)->eta()) > 1.6 && dil.daughter(0)->p() > 200)
    	SF1 = ((0.9784 - 4.73e-5 * dil.daughter(0)->eta())/(0.9908 - 1.26e-5 * dil.daughter(0)->p())) / ((0.9784 - 4.73e-5 * 200)/(0.9908 - 1.26e-5 * 200)) ;
    if (fabs(dil.daughter(1)->eta()) <= 1.6 && dil.daughter(1)->p() > 100)
    	SF2 = (0.994 - 4.08e-6 * dil.daughter(1)->p())/(0.994 - 4.08e-6 * 100);
    else if (fabs(dil.daughter(1)->eta()) > 1.6 && dil.daughter(1)->p() > 200)
	SF2 = ((0.9784 - 4.73e-5 * dil.daughter(1)->p())/(0.9908 - 1.26e-5 * dil.daughter(1)->p())) / ((0.9784 - 4.73e-5 * 200)/(0.9908 - 1.26e-5 * 200) ) ;

    double recoWeight = SF1*SF2;

    if (year_info == 2017 || year_info == 2018){
	if (dil.daughter(0)->eta()<-1.2 || dil.daughter(1)->eta()<-1.2 || dil.daughter(0)->eta()>1.2 || dil.daughter(1)->eta()>1.2) recoWeight = getRecoWeight(vertex_mass,false,year_info);
	else recoWeight = getRecoWeight(vertex_mass,true,year_info);
    }
    
    

     //std::cout<<" filling mass "<<vertex_mass<<std::endl;
    if (fill_gen_info) vertex_mass = getSmearedMass(dil,gM,year_info,false);
    float smearedMass = getSmearedMass(dil,gM,year_info,true);

    //if (year_info == 2017 || year_info == 2018){
    _scaleUncertBB = ScaleUncert(vertex_mass, true, year_info);
    _scaleUncertBE = ScaleUncert(vertex_mass, false, year_info);
    //}
    //std::cout <<  _madgraphWeight << " " << _kFactor << " " << _puWeight << std::endl;
    DimuonMassVertexConstrained->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight);
    DimuonMassVertexConstrainedVsCS->Fill(vertex_mass,cos_cs, _madgraphWeight*_kFactor*_puWeight);
    DimuonMassVertexConstrainedSmear->Fill(smearedMass, _madgraphWeight*_kFactor*_puWeight);
    DimuonMassVertexConstrainedMuonID->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight*recoWeight);
    DimuonMassVertexConstrainedScaleUp->Fill(vertex_mass*(1+_scaleUncertBB), _madgraphWeight*_kFactor*_puWeight);
    DimuonMassVertexConstrainedScaleDown->Fill(vertex_mass*(1-_scaleUncertBB), _madgraphWeight*_kFactor*_puWeight);
    if (cos_cs > -998.){
     if (cos_cs >= 0) DimuonMassVertexConstrained_CSPos->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight);
     else DimuonMassVertexConstrained_CSNeg->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight);
    }
    DimuonMassVtxConstrainedLog->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight);
    DimuonMassConstrainedVsUn->Fill(dil.mass(), vertex_mass, _madgraphWeight*_kFactor*_puWeight);
    DimuonMassVertexConstrainedError->Fill(vertex_mass, vertex_mass_err, _madgraphWeight*_kFactor*_puWeight);
    DimuonMassVertexConstrainedWeight->Fill(vertex_mass,_prescaleWeight*_madgraphWeight*_kFactor*_puWeight);
    DimuonMassVtxConstrainedLogWeight->Fill(vertex_mass,_prescaleWeight*_madgraphWeight*_kFactor*_puWeight);

  
    // plot per categories
  if (dil.daughter(0)->eta()<=1.2 && dil.daughter(1)->eta()<=1.2 && dil.daughter(0)->eta()>=-1.2 && dil.daughter(1)->eta()>=-1.2){
        DimuonMassVertexConstrained_bb->Fill(vertex_mass,_madgraphWeight*_kFactor_bb*_puWeight);
        DimuonMassVertexConstrainedVsCS_bb->Fill(vertex_mass,cos_cs,_madgraphWeight*_kFactor_bb*_puWeight);
        DimuonMassVtxConstrainedLog_bb->Fill(vertex_mass, _madgraphWeight*_kFactor_bb*_puWeight);
        DimuonMassVertexConstrainedSmear_bb->Fill(smearedMass, _madgraphWeight*_kFactor_bb*_puWeight);
    	DimuonMassVertexConstrainedMuonID_bb->Fill(vertex_mass, _madgraphWeight*_kFactor_bb*_puWeight*recoWeight);
        DimuonMassVertexConstrainedScaleUp_bb->Fill(vertex_mass*(1+_scaleUncertBB), _madgraphWeight*_kFactor_bb*_puWeight);
        DimuonMassVertexConstrainedScaleDown_bb->Fill(vertex_mass*(1-_scaleUncertBB), _madgraphWeight*_kFactor_bb*_puWeight);
        DileptonMass_bb->Fill(dil.mass(), _madgraphWeight*_kFactor_bb*_puWeight);
	DimuonMassVertexConstrained_gen_bb->Fill(gM,_madgraphWeight*_kFactor_bb*_puWeight);
	DimuonMassVtxConstrainedLog_gen_bb->Fill(gM, _madgraphWeight*_kFactor_bb*_puWeight);
	DimuonResponse_bb->Fill(vertex_mass,gM,_madgraphWeight*_kFactor_bb*_puWeight);
	DimuonResponseSmear_bb->Fill(smearedMass,gM,_madgraphWeight*_kFactor_bb*_puWeight);
	DimuonResponseMassScaleUp_bb->Fill(vertex_mass*(1+_scaleUncertBB),gM,_madgraphWeight*_kFactor_bb*_puWeight);
	DimuonResponseMassScaleDown_bb->Fill(vertex_mass*(1-_scaleUncertBB),gM,_madgraphWeight*_kFactor_bb*_puWeight);
	DimuonResponseMassScaleUpSmear_bb->Fill(smearedMass*(1+_scaleUncertBB),gM,_madgraphWeight*_kFactor_bb*_puWeight);
        DimuonResponseMassScaleDownSmear_bb->Fill(smearedMass*(1-_scaleUncertBB),gM,_madgraphWeight*_kFactor_bb*_puWeight);
        if (cos_cs > -998.){
    		if (cos_cs >= 0){
			 DimuonMassVertexConstrained_bb_CSPos->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight);
			 DimuonMassVertexConstrainedSmear_bb_CSPos->Fill(smearedMass, _madgraphWeight*_kFactor*_puWeight);
    			 DimuonMassVertexConstrainedMuonID_bb_CSPos->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight*recoWeight);
    			 DimuonMassVertexConstrainedScaleUp_bb_CSPos->Fill(vertex_mass*(1+_scaleUncertBB), _madgraphWeight*_kFactor*_puWeight);
    			 DimuonMassVertexConstrainedScaleDown_bb_CSPos->Fill(vertex_mass*(1-_scaleUncertBB), _madgraphWeight*_kFactor*_puWeight);
		}
     		else{
			 DimuonMassVertexConstrained_bb_CSNeg->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight);
		         DimuonMassVertexConstrainedSmear_bb_CSNeg->Fill(smearedMass, _madgraphWeight*_kFactor*_puWeight);
    			 DimuonMassVertexConstrainedMuonID_bb_CSNeg->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight*recoWeight);
    			 DimuonMassVertexConstrainedScaleUp_bb_CSNeg->Fill(vertex_mass*(1+_scaleUncertBB), _madgraphWeight*_kFactor*_puWeight);
    			 DimuonMassVertexConstrainedScaleDown_bb_CSNeg->Fill(vertex_mass*(1-_scaleUncertBB), _madgraphWeight*_kFactor*_puWeight);
		}
   		if (cos_cs >= 0) DileptonMass_bb_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_puWeight);
     		else DileptonMass_bb_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_puWeight);
   	 }



}

 if (dil.daughter(0)->eta()<-1.2 || dil.daughter(1)->eta()<-1.2 || dil.daughter(0)->eta()>1.2 || dil.daughter(1)->eta()>1.2){
        DimuonMassVertexConstrained_be->Fill(vertex_mass,_madgraphWeight*_kFactor_be*_puWeight);
        DimuonMassVertexConstrainedVsCS_be->Fill(vertex_mass,cos_cs,_madgraphWeight*_kFactor_be*_puWeight);
        DimuonMassVtxConstrainedLog_be->Fill(vertex_mass, _madgraphWeight*_kFactor_be*_puWeight);
        DimuonMassVertexConstrainedSmear_be->Fill(smearedMass, _madgraphWeight*_kFactor_be*_puWeight);
    	DimuonMassVertexConstrainedMuonID_be->Fill(vertex_mass, _madgraphWeight*_kFactor_be*_puWeight*recoWeight);
        DimuonMassVertexConstrainedScaleUp_be->Fill(vertex_mass*(1+_scaleUncertBE), _madgraphWeight*_kFactor_be*_puWeight);
        DimuonMassVertexConstrainedScaleDown_be->Fill(vertex_mass*(1-_scaleUncertBE), _madgraphWeight*_kFactor_be*_puWeight);
	DimuonMassVertexConstrained_gen_be->Fill(gM,_madgraphWeight*_kFactor_be*_puWeight);
        DimuonMassVtxConstrainedLog_gen_be->Fill(gM, _madgraphWeight*_kFactor_be*_puWeight);
        DimuonResponse_be->Fill(vertex_mass,gM,_madgraphWeight*_kFactor_be*_puWeight);
        DimuonResponseSmear_be->Fill(smearedMass,gM,_madgraphWeight*_kFactor_be*_puWeight);
        DimuonResponseMassScaleUp_be->Fill(vertex_mass*(1+_scaleUncertBE),gM,_madgraphWeight*_kFactor_be*_puWeight);
        DimuonResponseMassScaleDown_be->Fill(vertex_mass*(1-_scaleUncertBE),gM,_madgraphWeight*_kFactor_be*_puWeight);
        DimuonResponseMassScaleUpSmear_be->Fill(smearedMass*(1+_scaleUncertBE),gM,_madgraphWeight*_kFactor_be*_puWeight);
        DimuonResponseMassScaleDownSmear_be->Fill(smearedMass*(1-_scaleUncertBE),gM,_madgraphWeight*_kFactor_be*_puWeight);
       DileptonMass_be->Fill(dil.mass(), _madgraphWeight*_kFactor_be*_puWeight);
        if (cos_cs > -998.){
   		if (cos_cs >= 0){
			DimuonMassVertexConstrained_be_CSPos->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight);
		        DimuonMassVertexConstrainedSmear_be_CSPos->Fill(smearedMass, _madgraphWeight*_kFactor*_puWeight);
    			DimuonMassVertexConstrainedMuonID_be_CSPos->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight*recoWeight);
                        DimuonMassVertexConstrainedScaleUp_be_CSPos->Fill(vertex_mass*(1+_scaleUncertBE), _madgraphWeight*_kFactor*_puWeight);
   			DimuonMassVertexConstrainedScaleDown_be_CSPos->Fill(vertex_mass*(1-_scaleUncertBE), _madgraphWeight*_kFactor*_puWeight);
		}
     		else{
			 DimuonMassVertexConstrained_be_CSNeg->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight);
		         DimuonMassVertexConstrainedSmear_be_CSNeg->Fill(smearedMass, _madgraphWeight*_kFactor*_puWeight);
        		 DimuonMassVertexConstrainedScaleUp_be_CSNeg->Fill(vertex_mass*(1+_scaleUncertBE), _madgraphWeight*_kFactor*_puWeight);
    			 DimuonMassVertexConstrainedMuonID_be_CSNeg->Fill(vertex_mass, _madgraphWeight*_kFactor*_puWeight*recoWeight);
    			 DimuonMassVertexConstrainedScaleDown_be_CSNeg->Fill(vertex_mass*(1-_scaleUncertBE), _madgraphWeight*_kFactor*_puWeight);
		}
    		if (cos_cs >= 0) DileptonMass_be_CSPos->Fill(dil.mass(), _madgraphWeight*_kFactor*_puWeight);
     		else DileptonMass_be_CSNeg->Fill(dil.mass(), _madgraphWeight*_kFactor*_puWeight);
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
    if (fill_gen_info) hardInteraction->Fill(event);

    if (_usePrescaleWeight) {
        edm::Handle<int> totalPrescale;
        event.getByLabel(edm::InputTag("getPrescales","TotalPrescale","Zprime2muAnalysis"), totalPrescale);
        _prescaleWeight = *totalPrescale;
    }
    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid()) _LRWeight = lrWeightProducer.calculateWeight(event,hardInteraction,gen_ev_info->alphaQED()); 
    WeightLR->Fill(_LRWeight);
    if (_useMadgraphWeight) {
        eventWeight = 1.;
	_madgraphWeight = 1.; 
	if (gen_ev_info.isValid()){
        	eventWeight = gen_ev_info->weight();
        	_madgraphWeight = ( eventWeight > 0 ) ? 1 : -1;
	}

        WeightMadGraph->Fill(_madgraphWeight);
    }
    
  _madgraphWeight = _madgraphWeight*_LRWeight;
  if (use_bs_and_pv)
    getBSandPV(event);

  edm::Handle<edm::View<reco::Candidate> > leptons;


  if (fill_gen_info && doElectrons && (year_info == 2016 || year_info == 2017)){
  	edm::Handle< double > theprefweight;
	event.getByToken(prefweight_token, theprefweight ) ;
	_prefireWeight  =(*theprefweight);

	edm::Handle< double > theprefweightup;
	event.getByToken(prefweightup_token, theprefweightup ) ;
	_prefireWeightUp =(*theprefweightup);

	edm::Handle< double > theprefweightdown;
	event.getByToken(prefweightdown_token, theprefweightdown ) ;
	_prefireWeightDown =(*theprefweightdown);
  }
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
  if (fill_gen_info) {
    	hardInteraction->Fill(event);
 	if(hardInteraction->IsValidForRes()) gM = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).mass();
  }
  _kFactor = 1.;
  _kFactor_bb = 1.;
  _kFactor_be = 1.;

  if (!dileptons.isValid())
    edm::LogWarning("DileptonHandleInvalid") << "tried to get " << dilepton_src << " and failed!";
  else {
    if (leptonsFromDileptons){
    	if (fill_gen_info) {
	        if (_usekFactor){
			if(hardInteraction->IsValidForRes()){
    				gM = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).mass();
				double NNPDFFac = 1.;
				double NNPDFFac_bb = 1.;
				double NNPDFFac_be = 1.;
				double gen_lead_pt = 0.;
				if (year_info == 2017){
					if (gM < 120){
						gen_lead_pt = (hardInteraction->lepPlusNoIB->pt()*(hardInteraction->lepPlusNoIB->pt()>hardInteraction->lepMinusNoIB->pt()) + hardInteraction->lepMinusNoIB->pt()*(hardInteraction->lepMinusNoIB->pt()>hardInteraction->lepPlusNoIB->pt()));
						NNPDFFac =((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.8245728118) + (-0.0537728412909)*pow(gen_lead_pt,1) + (0.000731365981935)*pow(gen_lead_pt,2) + (7.16669312495e-06)*pow(gen_lead_pt,3) + (-1.99723894101e-07)*pow(gen_lead_pt,4) + (1.0112316789e-09)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.01849023288));
						NNPDFFac_bb = ((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.91383074609) + (-0.0596201865777)*pow(gen_lead_pt,1) + (0.000811074027001)*pow(gen_lead_pt,2) + (7.90677720686e-06)*pow(gen_lead_pt,3) + (-2.21489848717e-07)*pow(gen_lead_pt,4) + (1.12700571973e-09)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.00484010198));                                                                                            
						NNPDFFac_be = ((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.71913319508) + (-0.0481243962238)*pow(gen_lead_pt,1) + (0.000666286154366)*pow(gen_lead_pt,2) + (6.45776405133e-06)*pow(gen_lead_pt,3) + (-1.82202504311e-07)*pow(gen_lead_pt,4) + (9.24567381899e-10)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.02790393101));                                                                                                                                                          

					}
					else{
						NNPDFFac = ((0.918129) + (6.92702e-05)*pow(gM,1) + (1.62175e-08)*pow(gM,2) + (-2.47833e-11)*pow(gM,3) + (8.75707e-15)*pow(gM,4) + (-7.53019e-19)*pow(gM,5));
						NNPDFFac_bb = ((0.914053) + (7.91618e-05)*pow(gM,1) + (2.19722e-08)*pow(gM,2) + (-3.49212e-11)*pow(gM,3) + (1.22504e-14)*pow(gM,4) + (-1.07347e-18)*pow(gM,5));
						NNPDFFac_be = ((0.933214) + (3.76813e-05)*pow(gM,1) + (1.95612e-08)*pow(gM,2) + (-1.2688e-11)*pow(gM,3) + (3.69867e-15)*pow(gM,4) + (-2.62212e-19)*pow(gM,5)); 
					}

				}
				if (year_info == 2018){
					if (gM < 120){
                                                gen_lead_pt = (hardInteraction->lepPlusNoIB->pt()*(hardInteraction->lepPlusNoIB->pt()>hardInteraction->lepMinusNoIB->pt()) + hardInteraction->lepMinusNoIB->pt()*(hardInteraction->lepMinusNoIB->pt()>hardInteraction->lepPlusNoIB->pt()));
						NNPDFFac = ((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.69147781688) + (-0.0473286496053)*pow(gen_lead_pt,1) + (0.000661599919558)*pow(gen_lead_pt,2) + (6.33324308996e-06)*pow(gen_lead_pt,3) + (-1.80459280586e-07)*pow(gen_lead_pt,4) + (9.19632449685e-10)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.02344217328));
						NNPDFFac_bb = ((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.65477513925) + (-0.0472097707001)*pow(gen_lead_pt,1) + (0.000681831627146)*pow(gen_lead_pt,2) + (6.15645344304e-06)*pow(gen_lead_pt,3) + (-1.82810037593e-07)*pow(gen_lead_pt,4) + (9.43667804224e-10)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.01489199674));
						NNPDFFac_be =((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.60977951604) + (-0.0426122819079)*pow(gen_lead_pt,1) + (0.000599273084801)*pow(gen_lead_pt,2) + (5.88395881526e-06)*pow(gen_lead_pt,3) + (-1.66414436738e-07)*pow(gen_lead_pt,4) + (8.4690800397e-10)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.02846360871));
					}
					else{
						NNPDFFac = ((0.919027) + (5.98337e-05)*pow(gM,1) + (2.56077e-08)*pow(gM,2) + (-2.82876e-11)*pow(gM,3) + (9.2782e-15)*pow(gM,4) + (-7.77529e-19)*pow(gM,5));                                  
						NNPDFFac_bb = ((0.911563) + (0.000113313)*pow(gM,1) + (-2.35833e-08)*pow(gM,2) + (-1.44584e-11)*pow(gM,3) + (8.41748e-15)*pow(gM,4) + (-8.16574e-19)*pow(gM,5));
						NNPDFFac_be = ((0.934502) + (2.21259e-05)*pow(gM,1) + (4.14656e-08)*pow(gM,2) + (-2.26011e-11)*pow(gM,3) + (5.58804e-15)*pow(gM,4) + (-3.92687e-19)*pow(gM,5));
					}
				}		
 			 	if (doElectrons){
					if (year_info == 2017 or year_info == 2018){
						if (gM < 120){
							gen_lead_pt = (hardInteraction->lepPlusNoIB->et()*(hardInteraction->lepPlusNoIB->et()>hardInteraction->lepMinusNoIB->et()) + hardInteraction->lepMinusNoIB->et()*(hardInteraction->lepMinusNoIB->et()>hardInteraction->lepPlusNoIB->et()));
							NNPDFFac = 1.;
							NNPDFFac_bb = (gen_lead_pt<150) ? 3.596-0.2076 *gen_lead_pt+0.005795*pow(gen_lead_pt,2)-7.421e-05*pow(gen_lead_pt,3)+4.447e-07*pow(gen_lead_pt,4)-1.008e-9 *pow(gen_lead_pt,5) : 0.969125;
							NNPDFFac_be = (gen_lead_pt<150) ? 2.066-0.09495*gen_lead_pt+0.002664*pow(gen_lead_pt,2)-3.242e-05*pow(gen_lead_pt,3)+1.755e-07*pow(gen_lead_pt,4)-3.424e-10*pow(gen_lead_pt,5) : 1.191875;
						}
						else{
							NNPDFFac = 1.;	
							NNPDFFac_bb=(gM<5000) ? 0.8934+0.0002193 *gM-1.961e-7*pow(gM,2)+8.704e-11*pow(gM,3)-1.551e-14*pow(gM,4)+1.112e-18*pow(gM,5) : 1.74865;
							NNPDFFac_be=(gM<5000) ? 0.8989+0.000182  *gM-1.839e-7*pow(gM,2)+1.026e-10*pow(gM,3)-2.361e-14*pow(gM,4)+1.927e-18*pow(gM,5) : 1.302025;
						}

					}
					else{
						NNPDFFac=1.;
						NNPDFFac_bb=1.;
						NNPDFFac_be=1.;
					}

					_kFactor = (gM > 120) ? (1.0678 - 0.000120666 * gM + 3.22646e-08 * pow(gM,2) - 3.94886e-12 * pow(gM,3))*NNPDFFac : 1.*NNPDFFac;
					_kFactor_bb = (gM > 120) ? (1.0678 - 0.000120666 * gM + 3.22646e-08 * pow(gM,2) - 3.94886e-12 * pow(gM,3))*NNPDFFac_bb : 1*NNPDFFac_bb;
					_kFactor_be = (gM > 120) ? (1.0678 - 0.000120666 * gM + 3.22646e-08 * pow(gM,2) - 3.94886e-12 * pow(gM,3))*NNPDFFac_be: 1*NNPDFFac_be;
				}
				else{
					if(gM < 150){
						_kFactor = 1*NNPDFFac;
						_kFactor_bb = 1*NNPDFFac_bb;
						_kFactor_be = 1*NNPDFFac_be;
					}
					if(gM > 150){
			 			   	_kFactor = (1.053 - 0.0001552 * gM + 5.661e-08 * pow(gM,2) - 8.382e-12 * pow(gM,3))*NNPDFFac;
			 		   		_kFactor_bb = (1.032 - 0.000138 * gM + 4.827e-08 * pow(gM,2) - 7.321e-12 * pow(gM,3))*NNPDFFac_bb;
			 		   		_kFactor_be = (1.064 - 0.0001674 * gM + 6.599e-08 * pow(gM,2) - 9.657e-12 * pow(gM,3))*NNPDFFac_be;
			 		}
				}
			 	
 	   		} // hardInter
    			else std::cout<<"problems"<<std::endl;
    		} //k-Factor

	        else if (_useTTBarWeight && !(doElectrons)){
    			hardInteraction->Fill(event);
			if(hardInteraction->IsValidForRes()){
 	   			gM = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).mass();
				double NNPDFFac = 1.;
				double NNPDFFac_bb = 1.;
				double NNPDFFac_be = 1.;
				if (year_info == 2017){
					NNPDFFac = ((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005));
					NNPDFFac_bb = ((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005));
					NNPDFFac_be =((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005)); 
				}
				if (year_info == 2018){
					NNPDFFac = ((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005));
					NNPDFFac_bb =((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005));
					NNPDFFac_be =((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005));
				}
				_kFactor = 1*NNPDFFac;
				_kFactor_bb = 1*NNPDFFac_bb;
				_kFactor_be = 1*NNPDFFac_be;
 	   		} // hardInter
    			//else std::cout<<"problems"<<std::endl;
    		} //gen_info
    	} //kFactor
      fillLeptonHistosFromDileptons(*dileptons);
    }
    fillDileptonHistos(*dileptons, event, gM);
  }
}

DEFINE_FWK_MODULE(Zprime2muHistosFromPAT);
