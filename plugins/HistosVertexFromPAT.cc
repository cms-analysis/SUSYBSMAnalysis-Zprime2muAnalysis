// JMTBAD make this fill lepton histos from all leptons + those that
// make it into dileptons always, not just depending on flag.

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"


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

#include "DataFormats/TrackReco/interface/TrackBase.h"


class Zprime2muHistosVertexFromPAT : public edm::EDAnalyzer {
 public:
  explicit Zprime2muHistosVertexFromPAT(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void getBSandPV(const edm::Event&);
  void fillBasicLeptonHistos(const reco::CandidateBaseRef&);
  void fillOfflineMuonHistos(const pat::Muon*);
  void fillOfflineElectronHistos(const pat::Electron*);
  void fillLeptonHistos(const reco::CandidateBaseRef&);
  void fillLeptonHistos(const edm::View<reco::Candidate>&);
  void fillLeptonHistosFromDileptons(const pat::CompositeCandidateCollection&);
  void fillDileptonHistos(const pat::CompositeCandidate&);
  void fillDileptonHistos(const pat::CompositeCandidateCollection&);

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
  TH1F* DileptonMass;
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
  TH1F* DimuonMassVertexConstrainedWeight;
  TH1F* DimuonMassVtxConstrainedLog;
  TH1F* DimuonMassVtxConstrainedLogWeight;
  TH2F* DimuonMassConstrainedVsUn;
  TH2F* DimuonMassVertexConstrainedError;
    //special
    TH1F* DimuonMassVtx_chi2;
    TH1F* DimuonMassVtx_prob;
    TH2F* DimuonMassVertexProbMass;
    TH2F* DimuonMassVertexProbPt;
    TH2F* DimuonMassVertexProbEta;
    TH2F* DimuonMassVertexProbDilEta;
    TH2F* DimuonMassVertexProbAbsEta;
    TH2F* DimuonMassVertexProbAbsDilEta;
    TH1F* DimuonMuonPtProb;
    TH2F* DimuonMuonPtProbVertexProbPt;
    /////////
    TH1F* DimuonMassVtx_LeptonEta;
    TH1F* DimuonMassVtx_LeptonEta_MassRange;
    
    TH1F* DimuonMassVtx_prob_BarrelBarrel;
    TH1F* DimuonMassVtx_prob_EndcapEndcap;
    TH1F* DimuonMassVtx_prob_EndcapBarrel;
    
    TH1F* DimuonMuonPtProb_Barrel;
    TH1F* DimuonMuonPtError_Barrel;
    TH1F* DimuonMuonPtErrOverPt_Barrel;

    TH1F* DimuonMuonPtProb_Endcap;
    TH1F* DimuonMuonPtError_Endcap;
    TH1F* DimuonMuonPtErrOverPt_Endcap;
    
    TH1F* DimuonMassVtx_prob_BarrelBarrel_MassRange;
    TH1F* DimuonMassVtx_prob_EndcapEndcap_MassRange;
    TH1F* DimuonMassVtx_prob_EndcapBarrel_MassRange;
    
    TH1F* DimuonMuonPtProb_Barrel_MassRange;
    TH1F* DimuonMuonPtError_Barrel_MassRange;
    TH1F* DimuonMuonPtErrOverPt_Barrel_MassRange;
    
    TH1F* DimuonMuonPtProb_Endcap_MassRange;
    TH1F* DimuonMuonPtError_Endcap_MassRange;
    TH1F* DimuonMuonPtErrOverPt_Endcap_MassRange;
    
    TH1F* TrackD0BS_Barrel;
    TH1F* TrackD0PV_Barrel;
    TH1F* TrackD0BS_Endcap;
    TH1F* TrackD0PV_Endcap;
    TH1F* TrackD0BS_Barrel_MassRange;
    TH1F* TrackD0PV_Barrel_MassRange;
    TH1F* TrackD0BS_Endcap_MassRange;
    TH1F* TrackD0PV_Endcap_MassRange;
    
    TH1F* TrackIteration;
    TH1F* TrackIteration_Barrel;
    TH1F* TrackIteration_Endcap;
    TH1F* TrackIteration_Barrel_MassRange;
    TH1F* TrackIteration_Endcap_MassRange;
    TH2F* TrackIterationVSPtProb;
    TH2F* TrackIterationVSPtProb_Barrel;
    TH2F* TrackIterationVSPtProb_Endcap;
    TH2F* TrackIterationVSPtProb_Barrel_MassRange;
    TH2F* TrackIterationVSPtProb_Endcap_MassRange;
    
    TH2F* TrackIterationVSDimuonMassVtx_prob_BarrelBarrel;
    TH2F* TrackIterationVSDimuonMassVtx_prob_EndcapEndcap;
    TH2F* TrackIterationVSDimuonMassVtx_prob_EndcapBarrel;
    TH2F* TrackIterationVSDimuonMassVtx_prob1314_BarrelBarrel;
    TH2F* TrackIterationVSDimuonMassVtx_prob1314_EndcapEndcap;
    TH2F* TrackIterationVSDimuonMassVtx_prob1314_EndcapBarrel;
    
    TH2F* TrackIterationVSIsoSumPt;
    TH2F* TrackIterationVSRelIsoSumPt;
    TH2F* TrackIterationVSIsoNTracks;
    TH2F* PtProbVSIsoSumPt;
    TH2F* PtProbVSRelIsoSumPt;
    TH2F* PtProbVSIsoNTracks;
    TH2F* TrackIterationVSLeptonPt;
    TH2F* TrackIterationVSLeptonTrackerPt;
    TH2F* PtProbVSLeptonPt;
    TH2F* PtProbVSLeptonTrackerPt;
    TH2F* PtTrackerProbVSLeptonTrackerPt;


    
};

Zprime2muHistosVertexFromPAT::Zprime2muHistosVertexFromPAT(const edm::ParameterSet& cfg)
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
    _prescaleWeight(1)
{
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
  
  // Dilepton invariant mass.
  DileptonMass            = fs->make<TH1F>("DileptonMass",            titlePrefix + "dil. mass", 2000, 0, 20000);
  DileptonMassWeight      = fs->make<TH1F>("DileptonMassWeight",      titlePrefix + "dil. mass", 2000, 0, 20000);
  DileptonWithPhotonsMass = fs->make<TH1F>("DileptonWithPhotonsMass", titlePrefix + "res. mass", 2000, 0, 20000);
  
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
  DimuonMassVertexConstrained = fs->make<TH1F>("DimuonMassVertexConstrained", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000);
  DimuonMassVertexConstrainedWeight = fs->make<TH1F>("DimuonMassVertexConstrainedWeight", titlePrefix + "dimu. vertex-constrained mass", 2000, 0, 20000);
  // Mass plot in bins of log(mass)
  const int    NMBINS = 100;
  const double MMIN = 60., MMAX = 2100.;
  double logMbins[NMBINS+1];
  for (int ibin = 0; ibin <= NMBINS; ibin++)
    logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
  DimuonMassVtxConstrainedLog = fs->make<TH1F>("DimuonMassVtxConstrainedLog", titlePrefix + "dimu vtx-constrained mass in log bins", NMBINS, logMbins);
  DimuonMassVtxConstrainedLogWeight = fs->make<TH1F>("DimuonMassVtxConstrainedLogWeight", titlePrefix + "dimu vtx-constrained mass in log bins", NMBINS, logMbins);
  DimuonMassConstrainedVsUn = fs->make<TH2F>("DimuonMassConstrainedVsUn", titlePrefix + "dimu. vertex-constrained vs. non-constrained mass", 200, 0, 3000, 200, 0, 3000);
  DimuonMassVertexConstrainedError = fs->make<TH2F>("DimuonMassVertexConstrainedError", titlePrefix + "dimu. vertex-constrained mass error vs. mass", 100, 0, 3000, 100, 0, 400);
    //special
    DimuonMassVtx_chi2 = fs->make<TH1F>("DimuonMassVtx_chi2", titlePrefix + "dimu. vertex #chi^{2}/dof", 300, 0, 30);
    DimuonMassVtx_prob = fs->make<TH1F>("DimuonMassVtx_prob", titlePrefix + "dimu. vertex probability", 1000, 0., 1.);
    DimuonMassVertexProbMass = fs->make<TH2F>("DimuonMassVertexProbMass", titlePrefix + "dimu. vertex prob vs. mass", 2000, 0, 20000, 1000, 0., 1);
    DimuonMassVertexProbPt = fs->make<TH2F>("DimuonMassVertexProbPt", titlePrefix + "dimu. vertex prob vs. pt", 500, 0, 5000, 1000, 0., 1);
    DimuonMassVertexProbEta = fs->make<TH2F>("DimuonMassVertexProbEta", titlePrefix + "dimu. vertex prob vs. muon eta", 100, -6, 6, 1000, 0., 1);
    DimuonMassVertexProbDilEta = fs->make<TH2F>("DimuonMassVertexProbDilEta", titlePrefix + "dimu. vertex prob vs. dilepton eta", 100, -6, 6, 1000, 0., 1);
    DimuonMassVertexProbAbsEta = fs->make<TH2F>("DimuonMassVertexProbAbsEta", titlePrefix + "dimu. vertex prob vs. muon |eta|", 100, -6, 6, 1000, 0., 1);
    DimuonMassVertexProbAbsDilEta = fs->make<TH2F>("DimuonMassVertexProbAbsDilEta", titlePrefix + "dimu. vertex prob vs. dilepton |eta|", 100, -6, 6, 1000, 0., 1);
    DimuonMuonPtProb = fs->make<TH1F>("DimuonMuonPtProb", titlePrefix + "muon track prob", 1000, 0., 1.);
    DimuonMuonPtProbVertexProbPt = fs->make<TH2F>("DimuonMuonPtProbVertexProbPt", titlePrefix + "dimu. vertex prob vs. muon track prob", 1000, 0., 1., 1000, 0., 1);
    
    /////////
    DimuonMassVtx_LeptonEta = fs->make<TH1F>("DimuonMassVtx_LeptonEta", titlePrefix + "#eta", 100, -5, 5);
    DimuonMassVtx_LeptonEta_MassRange = fs->make<TH1F>("DimuonMassVtx_LeptonEta_MassRange", titlePrefix + "#eta 800<M<1200 GeV", 100, -5, 5);
    
    DimuonMassVtx_prob_BarrelBarrel= fs->make<TH1F>("DimuonMassVtx_prob_BarrelBarrel", titlePrefix + "dimu. vertex probability both muons |eta|<0.8", 1000, 0., 1.);
    DimuonMassVtx_prob_EndcapEndcap= fs->make<TH1F>("DimuonMassVtx_prob_EndcapEndcap", titlePrefix + "dimu. vertex probability both muons |eta|>1.", 1000, 0., 1.);
    DimuonMassVtx_prob_EndcapBarrel= fs->make<TH1F>("DimuonMassVtx_prob_EndcapBarrel", titlePrefix + "dimu. vertex probability one muon |eta|<0.8 && one muon |eta|>1.", 1000, 0., 1.);
    
    DimuonMuonPtProb_Barrel= fs->make<TH1F>("DimuonMuonPtProb_Barrel", titlePrefix + "muon track probability |eta|<0.8", 1000, 0., 1.);
    DimuonMuonPtError_Barrel= fs->make<TH1F>("DimuonMuonPtError_Barrel",        titlePrefix + "dil. #sigma_{pT}^{1} |eta|<0.8", 100, 0, 100);
    DimuonMuonPtErrOverPt_Barrel= fs->make<TH1F>("DimuonMuonPtErrOverPt_Barrel",     titlePrefix + "muon #sigma_{pT}/pT |eta|<0.8", 200, 0., 1.);
   
    DimuonMuonPtProb_Endcap= fs->make<TH1F>("DimuonMuonPtProb_Endcap", titlePrefix + "muon track probability |eta|>1", 1000, 0., 1.);
    DimuonMuonPtError_Endcap= fs->make<TH1F>("DimuonMuonPtError_Endcap",        titlePrefix + "dil. #sigma_{pT}^{1} |eta|>1", 100, 0, 100);
    DimuonMuonPtErrOverPt_Endcap= fs->make<TH1F>("DimuonMuonPtErrOverPt_Endcap",     titlePrefix + "muon #sigma_{pT}/pT |eta|>1", 200, 0., 1.);
    
    DimuonMassVtx_prob_BarrelBarrel_MassRange= fs->make<TH1F>("DimuonMassVtx_prob_BarrelBarrel_MassRange", titlePrefix + "dimu. vertex probability both muons |eta|<0.8 800<M<1200 GeV", 1000, 0., 1.);
    DimuonMassVtx_prob_EndcapEndcap_MassRange= fs->make<TH1F>("DimuonMassVtx_prob_EndcapEndcap_MassRange", titlePrefix + "dimu. vertex probability both muons |eta|>1. 800<M<1200 GeV", 1000, 0., 1.);
    DimuonMassVtx_prob_EndcapBarrel_MassRange= fs->make<TH1F>("DimuonMassVtx_prob_EndcapBarrel_MassRange", titlePrefix + "dimu. vertex probability one muon |eta|<0.8 && one muon |eta|>1. 800<M<1200 GeV", 1000, 0., 1.);
    
    DimuonMuonPtProb_Barrel_MassRange= fs->make<TH1F>("DimuonMuonPtProb_Barrel_MassRange", titlePrefix + "muon track probability |eta|<0.8. 800<M<1200 GeV", 1000, 0., 1.);
    DimuonMuonPtError_Barrel_MassRange= fs->make<TH1F>("DimuonMuonPtError_Barrel_MassRange",        titlePrefix + "dil. #sigma_{pT}^{1} |eta|<0.8 800<M<1200 GeV", 100, 0, 100);
    DimuonMuonPtErrOverPt_Barrel_MassRange= fs->make<TH1F>("DimuonMuonPtErrOverPt_Barrel_MassRange",     titlePrefix + "muon #sigma_{pT}/pT |eta|<0.8 800<M<1200 GeV", 200, 0., 1.);
  
    DimuonMuonPtProb_Endcap_MassRange= fs->make<TH1F>("DimuonMuonPtProb_Endcap_MassRange", titlePrefix + "muon track probability |eta|>1. 800<M<1200 GeV", 1000, 0., 1.);
    DimuonMuonPtError_Endcap_MassRange= fs->make<TH1F>("DimuonMuonPtError_Endcap_MassRange",        titlePrefix + "dil. #sigma_{pT}^{1} |eta|>1 800<M<1200 GeV", 100, 0, 100);
    DimuonMuonPtErrOverPt_Endcap_MassRange= fs->make<TH1F>("DimuonMuonPtErrOverPt_Endcap_MassRange",     titlePrefix + "muon #sigma_{pT}/pT |eta|>1 800<M<1200 GeV", 200, 0., 1.);
    
    TrackD0BS_Barrel = fs->make<TH1F>("TrackD0BS_Barrel", titlePrefix + "|d0 wrt BS| Barrel", 100, 0, 0.2);
    TrackD0PV_Barrel = fs->make<TH1F>("TrackD0PV_Barrel", titlePrefix + "|d0 wrt PV| Barrel", 100, 0, 0.2);
    TrackD0BS_Endcap = fs->make<TH1F>("TrackD0BS_Endcap", titlePrefix + "|d0 wrt BS| Endcap", 100, 0, 0.2);
    TrackD0PV_Endcap = fs->make<TH1F>("TrackD0PV_Endcap", titlePrefix + "|d0 wrt PV| Endcap", 100, 0, 0.2);
    TrackD0BS_Barrel_MassRange = fs->make<TH1F>("TrackD0BS_Barrel_MassRange", titlePrefix + "|d0 wrt BS| Barrel 800<M<1200", 100, 0, 0.2);
    TrackD0PV_Barrel_MassRange = fs->make<TH1F>("TrackD0PV_Barrel_MassRange", titlePrefix + "|d0 wrt PV| Barrel 800<M<1200", 100, 0, 0.2);
    TrackD0BS_Endcap_MassRange = fs->make<TH1F>("TrackD0BS_Endcap_MassRange", titlePrefix + "|d0 wrt BS| Endcap 800<M<1200", 100, 0, 0.2);
    TrackD0PV_Endcap_MassRange = fs->make<TH1F>("TrackD0PV_Endcap_MassRange", titlePrefix + "|d0 wrt PV| Endcap 800<M<1200", 100, 0, 0.2);
    
    TrackIteration = fs->make<TH1F>("TrackIteration", titlePrefix + "inner track algo()", 16, -0.5, 15.5 );
    TrackIteration_Barrel = fs->make<TH1F>("TrackIteration_Barrel", titlePrefix + "inner track algo() barrel", 16, -0.5, 15.5 );
    TrackIteration_Endcap = fs->make<TH1F>("TrackIteration_Endcap", titlePrefix + "inner track algo() endcap", 16, -0.5, 15.5 );
    TrackIteration_Barrel_MassRange = fs->make<TH1F>("TrackIteration_Barrel_MassRange", titlePrefix + "inner track algo() barrel 800<M<1200", 16, -0.5, 15.5 );
    TrackIteration_Endcap_MassRange = fs->make<TH1F>("TrackIteration_Endcap_MassRange", titlePrefix + "inner track algo() endcap 800<M<1200", 16, -0.5, 15.5 );
    TrackIterationVSPtProb = fs->make<TH2F>("TrackIterationVSPtProb", titlePrefix + "inner track algo() vs muon track probability", 16, -0.5, 15.5,  1000, 0., 1.);
    TrackIterationVSPtProb_Barrel = fs->make<TH2F>("TrackIterationVSPtProb_Barrel", titlePrefix + "inner track algo() vs muon track probability barrel", 16, -0.5, 15.5, 1000, 0., 1.);
    TrackIterationVSPtProb_Endcap = fs->make<TH2F>("TrackIterationVSPtProb_Endcap", titlePrefix + "inner track algo() vs muon track probability  endcap", 16, -0.5, 15.5, 1000, 0., 1.);
    TrackIterationVSPtProb_Barrel_MassRange = fs->make<TH2F>("TrackIterationVSPtProb_Barrel_MassRange", titlePrefix + "inner track algo() vs muon track probability  barrel 800<M<1200", 16, -0.5, 15.5,  1000, 0., 1.);
    TrackIterationVSPtProb_Endcap_MassRange = fs->make<TH2F>("TrackIterationVSPtProb_Endcap_MassRange", titlePrefix + "inner track algo() vs muon track probability endcap 800<M<1200", 16, -0.5, 15.5,  1000, 0., 1.);
    
    TrackIterationVSDimuonMassVtx_prob_BarrelBarrel= fs->make<TH2F>("TrackIterationVSDimuonMassVtx_prob_BarrelBarrel ", titlePrefix + "inner track algo() vs dimu. vertex probability both muons |eta|<0.8", 16, -0.5, 15.5, 1000, 0., 1.);
    TrackIterationVSDimuonMassVtx_prob_EndcapEndcap = fs->make<TH2F>("TrackIterationVSDimuonMassVtx_prob_EndcapEndcap ", titlePrefix + "inner track algo() vs dimu. vertex probability both muons |eta|>1.", 16, -0.5, 15.5, 1000, 0., 1.);
    TrackIterationVSDimuonMassVtx_prob_EndcapBarrel = fs->make<TH2F>("TrackIterationVSDimuonMassVtx_prob_EndcapBarrel ", titlePrefix + "inner track algo() vs dimu. vertex probability one muon |eta|<0.8 && one muon |eta|>1.", 16, -0.5, 15.5, 1000, 0., 1.);
    TrackIterationVSDimuonMassVtx_prob1314_BarrelBarrel= fs->make<TH2F>("TrackIterationVSDimuonMassVtx_prob1314_BarrelBarrel ", titlePrefix + "inner track algo() =13,14 vs dimu. vertex probability both muons |eta|<0.8", 16, -0.5, 15.5, 1000, 0., 1.);
    TrackIterationVSDimuonMassVtx_prob1314_EndcapEndcap = fs->make<TH2F>("TrackIterationVSDimuonMassVtx_prob1314_EndcapEndcap ", titlePrefix + "inner track algo() =13,14 vs dimu. vertex probability both muons |eta|>1.", 16, -0.5, 15.5, 1000, 0., 1.);
    TrackIterationVSDimuonMassVtx_prob1314_EndcapBarrel = fs->make<TH2F>("TrackIterationVSDimuonMassVtx_prob1314_EndcapBarrel ", titlePrefix + "inner track algo() =13,14 vs dimu. vertex probability one muon |eta|<0.8 && one muon |eta|>1.", 16, -0.5, 15.5, 1000, 0., 1.);
    
    ///////special lepton
    TrackIterationVSIsoSumPt = fs->make<TH2F>("TrackIterationVSIsoSumPt", titlePrefix + "inner track algo() vs Iso. (#Delta R < 0.3) #Sigma pT", 16, -0.5, 15.5, 50, 0, 50);
    TrackIterationVSRelIsoSumPt = fs->make<TH2F>("TrackIterationVSRelIsoSumPt", titlePrefix + "inner track algo() vs Iso. (#Delta R < 0.3) #Sigma pT / tk. pT", 16, -0.5, 15.5, 50, 0, 1);
    TrackIterationVSIsoNTracks = fs->make<TH2F>("TrackIterationVSIsoNTracks", titlePrefix + "inner track algo() vs Iso. (#Delta R < 0.3) nTracks", 16, -0.5, 15.5, 10, 0, 10);

    PtProbVSIsoSumPt = fs->make<TH2F>("PtProbVSIsoSumPt", titlePrefix + "muon track probability vs Iso. (#Delta R < 0.3) #Sigma pT", 1000, 0., 1., 50, 0, 50);
    PtProbVSRelIsoSumPt = fs->make<TH2F>("PtProbVSRelIsoSumPt", titlePrefix + "muon track probability vs Iso. (#Delta R < 0.3) #Sigma pT / tk. pT", 1000, 0., 1., 50, 0, 1);
    PtProbVSIsoNTracks = fs->make<TH2F>("PtProbVSIsoNTracks", titlePrefix + "muon track probability vs Iso. (#Delta R < 0.3) nTracks", 1000, 0., 1., 10, 0, 10);

    TrackIterationVSLeptonPt = fs->make<TH2F>("TrackIterationVSLeptonPt", titlePrefix + "inner track algo() vs pT", 16, -0.5, 15.5, 5000, 0, 5000);
    TrackIterationVSLeptonTrackerPt = fs->make<TH2F>("TrackIterationVSLeptonTrackerPt", titlePrefix + "inner track algo() vs tracker trk pT", 16, -0.5, 15.5, 5000, 0, 5000);
    PtProbVSLeptonPt = fs->make<TH2F>("PtProbVSLeptonPt", titlePrefix + "muon track probability vs pT", 1000, 0., 1., 5000, 0, 5000);
    PtProbVSLeptonTrackerPt = fs->make<TH2F>("PtProbVSLeptonTrackerPt", titlePrefix + "muon track probability vs tracker trk pT", 1000, 0., 1., 5000, 0, 5000);
    PtTrackerProbVSLeptonTrackerPt = fs->make<TH2F>("PtTrackerProbVSLeptonTrackerPt", titlePrefix + "muon tracker track probability vs tracker trk pT", 1000, 0., 1., 5000, 0, 5000);

}

void Zprime2muHistosVertexFromPAT::getBSandPV(const edm::Event& event) {
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
  NVertices->Fill(vertex_count);
}

void Zprime2muHistosVertexFromPAT::fillBasicLeptonHistos(const reco::CandidateBaseRef& lep) {
  LeptonEta->Fill(lep->eta());
  LeptonRap->Fill(lep->rapidity());
  LeptonPhi->Fill(lep->phi());

  LeptonPt->Fill(lep->pt());
  LeptonPz->Fill(fabs(lep->pz()));
  LeptonP ->Fill(lep->p());

  LeptonPtVsEta->Fill(lep->eta(), lep->pt());
  LeptonPVsEta ->Fill(lep->eta(), lep->p());
}

void Zprime2muHistosVertexFromPAT::fillOfflineMuonHistos(const pat::Muon* mu) {
  const reco::MuonIsolation& iso = mu->isolationR03();
  IsoSumPt   ->Fill(iso.sumPt);
  RelIsoSumPt->Fill(iso.sumPt / mu->innerTrack()->pt());
  IsoEcal    ->Fill(iso.emEt);
  IsoHcal    ->Fill(iso.hadEt + iso.hoEt);
  CombIso    ->Fill( iso.sumPt + iso.emEt + iso.hadEt + iso.hoEt);
  RelCombIso ->Fill((iso.sumPt + iso.emEt + iso.hadEt + iso.hoEt) / mu->innerTrack()->pt());
  IsoNTracks ->Fill(iso.nTracks);
  IsoNJets   ->Fill(iso.nJets);

  CombIsoNoECAL   ->Fill( iso.sumPt + iso.hadEt + iso.hoEt);
  RelCombIsoNoECAL->Fill((iso.sumPt + iso.hadEt + iso.hoEt) / mu->innerTrack()->pt());
    
    ///////////special
    const reco::TrackRef tk = patmuon::getPickedTrack(*mu);
    if (tk.isAvailable()) {
    TrackIterationVSIsoSumPt->Fill(mu->innerTrack()->algo(),iso.sumPt);
    TrackIterationVSRelIsoSumPt->Fill(mu->innerTrack()->algo(),iso.sumPt / mu->innerTrack()->pt());
    TrackIterationVSIsoNTracks ->Fill(mu->innerTrack()->algo(),iso.nTracks);

    PtProbVSIsoSumPt->Fill(TMath::Prob(tk->chi2(), tk->ndof()),iso.sumPt);
    PtProbVSRelIsoSumPt->Fill(TMath::Prob(tk->chi2(), tk->ndof()),iso.sumPt / mu->innerTrack()->pt());
    PtProbVSIsoNTracks ->Fill(TMath::Prob(tk->chi2(), tk->ndof()),iso.nTracks);
    
    TrackIterationVSLeptonPt->Fill(mu->innerTrack()->algo(),tk->pt());
    TrackIterationVSLeptonTrackerPt->Fill(mu->innerTrack()->algo(),mu->innerTrack()->pt());
    PtProbVSLeptonPt->Fill(TMath::Prob(tk->chi2(), tk->ndof()),tk->pt());
    PtProbVSLeptonTrackerPt->Fill(TMath::Prob(tk->chi2(), tk->ndof()),mu->innerTrack()->pt());
    PtTrackerProbVSLeptonTrackerPt->Fill(TMath::Prob(mu->innerTrack()->chi2(), mu->innerTrack()->ndof()),mu->innerTrack()->pt());
    }
    //////////////////

  const reco::TrackRef track = patmuon::getPickedTrack(*mu);
  if (track.isAvailable()) {
    Chi2dof->Fill(track->normalizedChi2());

    if (beamspot != 0) {
      TrackD0BS->Fill(fabs(track->dxy(beamspot->position())));
      TrackDZBS->Fill(fabs(track->dz (beamspot->position())));
    }

    if (vertex != 0) {
      TrackD0PV->Fill(fabs(track->dxy(vertex->position())));
      TrackDZPV->Fill(fabs(track->dz (vertex->position())));
    }

    const reco::HitPattern& hp = track->hitPattern();
    NPxHits->Fill(hp.numberOfValidPixelHits());
    NStHits->Fill(hp.numberOfValidStripHits());
    NTkHits->Fill(hp.numberOfValidTrackerHits());
    NMuHits->Fill(hp.numberOfValidMuonHits());

    NHits->Fill(hp.numberOfValidHits());
    NInvalidHits->Fill(hp.numberOfHits(reco::HitPattern::TRACK_HITS) - hp.numberOfValidHits());
    //NInvalidHits->Fill(hp.numberOfHits() - hp.numberOfValidHits());
    
    NPxLayers->Fill(hp.pixelLayersWithMeasurement());
    NStLayers->Fill(hp.stripLayersWithMeasurement());
    NTkLayers->Fill(hp.trackerLayersWithMeasurement());
  }
}

void Zprime2muHistosVertexFromPAT::fillOfflineElectronHistos(const pat::Electron* lep) {
  // Can add electron quantities here.
}

void Zprime2muHistosVertexFromPAT::fillLeptonHistos(const reco::CandidateBaseRef& lep) {
  fillBasicLeptonHistos(lep);
  
  const pat::Muon* muon = toConcretePtr<pat::Muon>(lep);
  if (muon) fillOfflineMuonHistos(muon);
  
  const pat::Electron* electron = toConcretePtr<pat::Electron>(lep);
  if (electron) fillOfflineElectronHistos(electron);
}

void Zprime2muHistosVertexFromPAT::fillLeptonHistos(const edm::View<reco::Candidate>& leptons) {
  NLeptons->Fill(leptons.size());

 // JMTBAD this should use leptonsPassingCuts or whatever
  int total_q = 0;
  for (size_t i = 0; i < leptons.size(); ++i) {
    total_q += leptons[i].charge();
    fillLeptonHistos(leptons.refAt(i));
  }

  LeptonSigns->Fill(leptons.size(), total_q);
}

void Zprime2muHistosVertexFromPAT::fillLeptonHistosFromDileptons(const pat::CompositeCandidateCollection& dileptons) {
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
  NLeptons->Fill(nleptons);
  LeptonSigns->Fill(nleptons, total_q);
}

void Zprime2muHistosVertexFromPAT::fillDileptonHistos(const pat::CompositeCandidate& dil) {
  if (dbg_tree) {
    dbg_t.mass = dil.mass();
    dbg_t.id = dil.daughter(0)->pdgId() + dil.daughter(1)->pdgId();
    dbg_tree->Fill();
  }
  DileptonEta->Fill(dil.eta());
  DileptonRap->Fill(dil.rapidity());
  DileptonPhi->Fill(dil.phi());

  DileptonPt->Fill(dil.pt());
  DileptonPz->Fill(fabs(dil.pz()));
  DileptonP ->Fill(dil.p());

  DileptonPtVsEta->Fill(dil.eta(), dil.pt());
  DileptonPVsEta ->Fill(dil.eta(), dil.p());

  DileptonMass->Fill(dil.mass());
  DileptonMassWeight->Fill(dil.mass(),_prescaleWeight);
  DileptonWithPhotonsMass->Fill(resonanceP4(dil).mass());

  const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
  const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);

  if (lep0.isNonnull() && lep1.isNonnull()) {
    DileptonDeltaPt->Fill(fabs(lep0->pt()) - fabs(lep1->pt()));
    DileptonDeltaP ->Fill(fabs(lep0->p())  - fabs(lep1->p()));

    const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
    const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
    if (mu0 && mu1) {
      const reco::Track* tk0 = patmuon::getPickedTrack(*mu0).get();
      const reco::Track* tk1 = patmuon::getPickedTrack(*mu1).get();
      if (tk0 && tk1) {
	DimuonMuonPtErrors->Fill(ptError(tk0), ptError(tk1));
	DimuonMuonPtErrOverPt->Fill(ptError(tk0)/tk0->pt());
	DimuonMuonPtErrOverPt->Fill(ptError(tk1)/tk1->pt());
	float mass = -999.;
	// Use mass calculated with the vertex constraint when available
	if (dil.hasUserFloat("vertexM"))
	  mass = dil.userFloat("vertexM");
	else
	  mass = dil.mass();
	if (mass > 200.) {
	  DimuonMuonPtErrOverPtM200->Fill(ptError(tk0)/tk0->pt());
	  DimuonMuonPtErrOverPtM200->Fill(ptError(tk1)/tk1->pt());
	}
	if (mass > 500.) {
	  DimuonMuonPtErrOverPtM500->Fill(ptError(tk0)/tk0->pt());
	  DimuonMuonPtErrOverPtM500->Fill(ptError(tk1)/tk1->pt());
	}
      }
    }
  }

  DileptonDaughterIds->Fill(dil.daughter(0)->pdgId(), dil.daughter(1)->pdgId());

  DileptonDaughterDeltaR->Fill(reco::deltaR(*dil.daughter(0), *dil.daughter(1)));
  DileptonDaughterDeltaPhi->Fill(reco::deltaPhi(dil.daughter(0)->phi(), dil.daughter(1)->phi())); 

  if (dil.hasUserFloat("vertexM") && dil.hasUserFloat("vertexMError")) {
    float vertex_mass = dil.userFloat("vertexM");
    float vertex_mass_err = dil.userFloat("vertexMError");
      std::cout<<"***************** filling dil mass "<<dil.mass()<<" vtx constrain "<<vertex_mass<<std::endl;
    DimuonMassVertexConstrained->Fill(vertex_mass);
    DimuonMassVtxConstrainedLog->Fill(vertex_mass);
    DimuonMassConstrainedVsUn->Fill(dil.mass(), vertex_mass);
    DimuonMassVertexConstrainedError->Fill(vertex_mass, vertex_mass_err);
    DimuonMassVertexConstrainedWeight->Fill(vertex_mass,_prescaleWeight);
    DimuonMassVtxConstrainedLogWeight->Fill(vertex_mass,_prescaleWeight);
      
    // special
    float vertex_chi2 = dil.userFloat("vertex_chi2");
      DimuonMassVtx_chi2->Fill(vertex_chi2);
    if (vertex_chi2 > 0 ) {
        float vertex_ndof = dil.userFloat("vertex_ndof");
        float vertex_chi2_noNormalized = vertex_chi2*vertex_ndof;
        float vertex_prob = TMath::Prob(vertex_chi2_noNormalized, vertex_ndof);
        DimuonMassVtx_prob->Fill(vertex_prob);
        DimuonMassVertexProbMass->Fill(dil.mass(), vertex_prob);
        DimuonMassVertexProbPt->Fill(dil.daughter(0)->pt(), vertex_prob);
        DimuonMassVertexProbPt->Fill(dil.daughter(1)->pt(), vertex_prob);
        DimuonMassVertexProbEta->Fill(dil.daughter(0)->eta(), vertex_prob);
        DimuonMassVertexProbEta->Fill(dil.daughter(1)->eta(), vertex_prob);
        DimuonMassVertexProbDilEta->Fill(dil.eta(), vertex_prob);
        DimuonMassVertexProbAbsEta->Fill(fabs(dil.daughter(0)->eta()), vertex_prob);
        DimuonMassVertexProbAbsEta->Fill(fabs(dil.daughter(1)->eta()), vertex_prob);
        DimuonMassVertexProbAbsDilEta->Fill(fabs(dil.eta()), vertex_prob);
        if (lep0.isNonnull() && lep1.isNonnull()) {
            const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
            const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
            if (mu0 && mu1) {
                const reco::TrackRef inner_tk0 = mu0->innerTrack();
                const reco::TrackRef inner_tk1 = mu1->innerTrack();
                const reco::Track* tk0 = patmuon::getPickedTrack(*mu0).get();
                const reco::Track* tk1 = patmuon::getPickedTrack(*mu1).get();
                const reco::TrackRef track0 = patmuon::getPickedTrack(*mu0);
                const reco::TrackRef track1 = patmuon::getPickedTrack(*mu1);
                if (tk0 && tk1) {
                //if(dil.mass()<=200 && dil.mass()>=100){//solo x dy120
                    //if(inner_tk0 && inner_tk1){
                        std::cout<<"ALGO NAME "<<inner_tk0->algo()<<std::endl;
                        std::cout<<inner_tk0->algoName()<<std::endl;//}

                    DimuonMuonPtProb->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()));
                    DimuonMuonPtProb->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()));
                    DimuonMuonPtProbVertexProbPt->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()), vertex_prob);
                    DimuonMuonPtProbVertexProbPt->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()), vertex_prob);
                    TrackIteration->Fill(inner_tk0->algo());
                    TrackIteration->Fill(inner_tk1->algo());
                    TrackIterationVSPtProb->Fill(inner_tk0->algo(),TMath::Prob(tk0->chi2(), tk0->ndof()));
                    TrackIterationVSPtProb->Fill(inner_tk1->algo(),TMath::Prob(tk1->chi2(), tk1->ndof()));
                    
                    DimuonMassVtx_LeptonEta->Fill(tk0->eta());
                    DimuonMassVtx_LeptonEta->Fill(tk1->eta());
                    if(dil.mass()<=1200 && dil.mass()>=800){
                        DimuonMassVtx_LeptonEta_MassRange->Fill(tk0->eta());
                        DimuonMassVtx_LeptonEta_MassRange->Fill(tk1->eta());
                    }
                    //BARREL
                    if(fabs(tk0->eta())<0.8){
                        DimuonMuonPtProb_Barrel->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()));
                        DimuonMuonPtError_Barrel->Fill(ptError(tk0));
                        DimuonMuonPtErrOverPt_Barrel->Fill(ptError(tk0)/tk0->pt());
                        TrackIteration_Barrel->Fill(inner_tk0->algo());
                        TrackIterationVSPtProb_Barrel->Fill(inner_tk0->algo(),TMath::Prob(tk0->chi2(), tk0->ndof()));
                        if (beamspot != 0) {TrackD0BS_Barrel->Fill(fabs(track0->dxy(beamspot->position())));}
                        if (vertex != 0) {TrackD0PV_Barrel->Fill(fabs(track0->dxy(vertex->position())));}
                    }
                    if(fabs(tk1->eta())<0.8){
                        DimuonMuonPtProb_Barrel->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()));
                        DimuonMuonPtError_Barrel->Fill(ptError(tk1));
                        DimuonMuonPtErrOverPt_Barrel->Fill(ptError(tk1)/tk1->pt());
                        TrackIteration_Barrel->Fill(inner_tk1->algo());
                        TrackIterationVSPtProb_Barrel->Fill(inner_tk1->algo(),TMath::Prob(tk1->chi2(), tk1->ndof()));
                        if (beamspot != 0) {TrackD0BS_Barrel->Fill(fabs(track1->dxy(beamspot->position())));}
                        if (vertex != 0) {TrackD0PV_Barrel->Fill(fabs(track1->dxy(vertex->position())));}
                    }
                    if(fabs(tk1->eta())<0.8 && fabs(tk0->eta())<0.8){
                    DimuonMassVtx_prob_BarrelBarrel->Fill(vertex_prob);
                        if(inner_tk0->algo()>12)TrackIterationVSDimuonMassVtx_prob1314_BarrelBarrel->Fill(inner_tk0->algo(),vertex_prob);
                        if(inner_tk1->algo()>12)TrackIterationVSDimuonMassVtx_prob1314_BarrelBarrel->Fill(inner_tk1->algo(),vertex_prob);
                        if((inner_tk0->algo()<=12) && (inner_tk1->algo()<=12)){
                            TrackIterationVSDimuonMassVtx_prob_BarrelBarrel->Fill(inner_tk0->algo(),vertex_prob);
                            TrackIterationVSDimuonMassVtx_prob_BarrelBarrel->Fill(inner_tk1->algo(),vertex_prob);
                        }
                 
                    }
                    if(dil.mass()<1200 && dil.mass()>800){
                        if(fabs(tk0->eta())<0.8){
                            DimuonMuonPtProb_Barrel_MassRange->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()));
                            DimuonMuonPtError_Barrel_MassRange->Fill(ptError(tk0));
                            DimuonMuonPtErrOverPt_Barrel_MassRange->Fill(ptError(tk0)/tk0->pt());
                            TrackIteration_Barrel_MassRange->Fill(inner_tk0->algo());
                            TrackIterationVSPtProb_Barrel_MassRange->Fill(inner_tk0->algo(),TMath::Prob(tk0->chi2(), tk0->ndof()));
                            if (beamspot != 0) {TrackD0BS_Barrel_MassRange->Fill(fabs(track0->dxy(beamspot->position())));}
                            if (vertex != 0) {TrackD0PV_Barrel_MassRange->Fill(fabs(track0->dxy(vertex->position())));}
                        }
                        if(fabs(tk1->eta())<0.8){
                            DimuonMuonPtProb_Barrel_MassRange->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()));
                            DimuonMuonPtError_Barrel_MassRange->Fill(ptError(tk1));
                            DimuonMuonPtErrOverPt_Barrel_MassRange->Fill(ptError(tk1)/tk1->pt());
                            TrackIteration_Barrel_MassRange->Fill(inner_tk1->algo());
                            TrackIterationVSPtProb_Barrel_MassRange->Fill(inner_tk1->algo(),TMath::Prob(tk1->chi2(), tk1->ndof()));
                            if (beamspot != 0) {TrackD0BS_Barrel_MassRange->Fill(fabs(track1->dxy(beamspot->position())));}
                            if (vertex != 0) {TrackD0PV_Barrel_MassRange->Fill(fabs(track1->dxy(vertex->position())));}
                        }
                        if(fabs(tk1->eta())<0.8 && fabs(tk0->eta())<0.8){
                            DimuonMassVtx_prob_BarrelBarrel_MassRange->Fill(vertex_prob);
                        }
                    }
                    //ENDCAP
                    if(fabs(tk0->eta())>1.){
                        DimuonMuonPtProb_Endcap->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()));
                        DimuonMuonPtError_Endcap->Fill(ptError(tk0));
                        DimuonMuonPtErrOverPt_Endcap->Fill(ptError(tk0)/tk0->pt());
                        TrackIteration_Endcap->Fill(inner_tk0->algo());
                        TrackIterationVSPtProb_Endcap->Fill(inner_tk0->algo(),TMath::Prob(tk0->chi2(), tk0->ndof()));
                        if (beamspot != 0) {TrackD0BS_Endcap->Fill(fabs(track0->dxy(beamspot->position())));}
                        if (vertex != 0) {TrackD0PV_Endcap->Fill(fabs(track0->dxy(vertex->position())));}
                    }
                    if(fabs(tk1->eta())>1.){
                        DimuonMuonPtProb_Endcap->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()));
                        DimuonMuonPtError_Endcap->Fill(ptError(tk1));
                        DimuonMuonPtErrOverPt_Endcap->Fill(ptError(tk1)/tk1->pt());
                        TrackIteration_Endcap->Fill(inner_tk1->algo());
                        TrackIterationVSPtProb_Endcap->Fill(inner_tk1->algo(),TMath::Prob(tk1->chi2(), tk1->ndof()));
                        if (beamspot != 0) {TrackD0BS_Endcap->Fill(fabs(track1->dxy(beamspot->position())));}
                        if (vertex != 0) {TrackD0PV_Endcap->Fill(fabs(track1->dxy(vertex->position())));}
                    }
                    if(fabs(tk1->eta())>1. && fabs(tk0->eta())>1.){
                        DimuonMassVtx_prob_EndcapEndcap->Fill(vertex_prob);
                        if(inner_tk0->algo()>12)TrackIterationVSDimuonMassVtx_prob1314_EndcapEndcap->Fill(inner_tk0->algo(),vertex_prob);
                        if(inner_tk1->algo()>12)TrackIterationVSDimuonMassVtx_prob1314_EndcapEndcap->Fill(inner_tk1->algo(),vertex_prob);
                        if((inner_tk0->algo()<=12) && (inner_tk1->algo()<=12)){
                            TrackIterationVSDimuonMassVtx_prob_EndcapEndcap->Fill(inner_tk0->algo(),vertex_prob);
                            TrackIterationVSDimuonMassVtx_prob_EndcapEndcap->Fill(inner_tk1->algo(),vertex_prob);
                        }
                    }
                    //ENDCAPBARREL
                    if((fabs(tk1->eta())>1. && fabs(tk0->eta())<0.8) || (fabs(tk0->eta())>1. && fabs(tk1->eta())<0.8)){
                        DimuonMassVtx_prob_EndcapBarrel->Fill(vertex_prob);
                        if(inner_tk0->algo()>12)TrackIterationVSDimuonMassVtx_prob1314_EndcapBarrel->Fill(inner_tk0->algo(),vertex_prob);
                        if(inner_tk1->algo()>12)TrackIterationVSDimuonMassVtx_prob1314_EndcapBarrel->Fill(inner_tk1->algo(),vertex_prob);
                        if((inner_tk0->algo()<=12) && (inner_tk1->algo()<=12)){
                            TrackIterationVSDimuonMassVtx_prob_EndcapBarrel->Fill(inner_tk0->algo(),vertex_prob);
                            TrackIterationVSDimuonMassVtx_prob_EndcapBarrel->Fill(inner_tk1->algo(),vertex_prob);
                        }
                    }
                    if(dil.mass()<=1200 && dil.mass()>=800){
                        if(fabs(tk0->eta())>1.){
                            DimuonMuonPtProb_Endcap_MassRange->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()));
                            DimuonMuonPtError_Endcap_MassRange->Fill(ptError(tk0));
                            DimuonMuonPtErrOverPt_Endcap_MassRange->Fill(ptError(tk0)/tk0->pt());
                            TrackIteration_Endcap_MassRange->Fill(inner_tk0->algo());
                            TrackIterationVSPtProb_Endcap_MassRange->Fill(inner_tk0->algo(),TMath::Prob(tk0->chi2(), tk0->ndof()));
                            if (beamspot != 0) {TrackD0BS_Endcap_MassRange->Fill(fabs(track0->dxy(beamspot->position())));}
                            if (vertex != 0) {TrackD0PV_Endcap_MassRange->Fill(fabs(track0->dxy(vertex->position())));}
                        }
                        if(fabs(tk1->eta())>1.){
                            DimuonMuonPtProb_Endcap_MassRange->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()));
                            DimuonMuonPtError_Endcap_MassRange->Fill(ptError(tk1));
                            DimuonMuonPtErrOverPt_Endcap_MassRange->Fill(ptError(tk1)/tk1->pt());
                            TrackIteration_Endcap_MassRange->Fill(inner_tk1->algo());
                            TrackIterationVSPtProb_Endcap_MassRange->Fill(inner_tk1->algo(),TMath::Prob(tk1->chi2(), tk1->ndof()));
                            if (beamspot != 0) {TrackD0BS_Endcap_MassRange->Fill(fabs(track1->dxy(beamspot->position())));}
                            if (vertex != 0) {TrackD0PV_Endcap_MassRange->Fill(fabs(track1->dxy(vertex->position())));}
                        }
                        if(fabs(tk1->eta())>1. && fabs(tk0->eta())>1.){
                            DimuonMassVtx_prob_EndcapEndcap_MassRange->Fill(vertex_prob);
                        }
                        //ENDCAPBARREL
                        if((fabs(tk1->eta())>1. && fabs(tk0->eta())<0.8) || (fabs(tk0->eta())>1. && fabs(tk1->eta())<0.8)){
                            DimuonMassVtx_prob_EndcapBarrel_MassRange->Fill(vertex_prob);
                        }
                    }
                    
                    
                }//}//solo per dy120
            }
        
        }
        
        
    }

  }
}

void Zprime2muHistosVertexFromPAT::fillDileptonHistos(const pat::CompositeCandidateCollection& dileptons) {
  NDileptons->Fill(dileptons.size());

  pat::CompositeCandidateCollection::const_iterator dil = dileptons.begin(), dile = dileptons.end();
  for ( ; dil != dile; ++dil)
    fillDileptonHistos(*dil);
}

void Zprime2muHistosVertexFromPAT::analyze(const edm::Event& event, const edm::EventSetup& setup) {
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
  
  if (!dileptons.isValid())
    edm::LogWarning("DileptonHandleInvalid") << "tried to get " << dilepton_src << " and failed!";
  else {
    if (leptonsFromDileptons)
      fillLeptonHistosFromDileptons(*dileptons);
    
    fillDileptonHistos(*dileptons);
  }
}

DEFINE_FWK_MODULE(Zprime2muHistosVertexFromPAT);
