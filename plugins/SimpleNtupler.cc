#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

class SimpleNtupler : public edm::EDAnalyzer {
public:
  explicit SimpleNtupler(const edm::ParameterSet&);
  ~SimpleNtupler() { delete hardInteraction; }
  void analyze(const edm::Event&, const edm::EventSetup&);

private:
  struct tree_t {
   
    unsigned run;
    unsigned lumi;
    unsigned event;
    float genWeight;
    float beamspot_x;
    float beamspot_x_err;
    float beamspot_y;
    float beamspot_y_err;
    float beamspot_z;
    float beamspot_z_err;
    int nvertices;
    int dil_chosen;
    float dil_mass;
    float dil_pt;
    float dil_rap;
    float dil_eta;
    float dil_phi;
    float dil_dR;
    float dil_dPhi;
    float cos_angle;
    float vertex_chi2;
    float cos_cs;
    float phi_cs;
    float vertex_m;
    float vertex_m_err;
    float vertex_x;
    float vertex_x_err;
    float vertex_y;
    float vertex_y_err;
    float vertex_z;
    float vertex_z_err;
    int lep_id[2];
    float lep_p[2];
    float lep_pt[2];
    float lep_pt_err[2];
    float lep_px[2];
    float lep_py[2];
    float lep_pz[2];
    float lep_E[2];
    float lep_eta[2];
    float lep_phi[2];
    float lep_qOverPt[2];
    float lep_tk_p[2];
    float lep_tk_pt[2];
    float lep_tk_pt_err[2];
    float lep_tk_px[2];
    float lep_tk_py[2];
    float lep_tk_pz[2];
    float lep_tk_eta[2];
    float lep_tk_phi[2];
    float lep_tk_dz[2];
    float lep_tk_vz[2];
    float lep_tk_chi2[2];
    float lep_tk_ndf[2];
    float lep_tk_qOverPt[2];
    float lep_glb_p[2];
    float lep_glb_pt[2];
    float lep_glb_pt_err[2];
    float lep_glb_px[2];
    float lep_glb_py[2];
    float lep_glb_pz[2];
    float lep_glb_eta[2];
    float lep_glb_phi[2];
    float lep_glb_chi2[2];
    float lep_glb_ndf[2];
    float lep_glb_qOverPt[2];
    float lep_tpfms_p[2];
    float lep_tpfms_pt[2];
    float lep_tpfms_pt_err[2];
    float lep_tpfms_px[2];
    float lep_tpfms_py[2];
    float lep_tpfms_pz[2];
    float lep_tpfms_eta[2];
    float lep_tpfms_phi[2];
    float lep_tpfms_chi2[2];
    float lep_tpfms_ndf[2];
    float lep_tpfms_qOverPt[2];
    float lep_picky_p[2];
    float lep_picky_pt[2];
    float lep_picky_pt_err[2];
    float lep_picky_px[2];
    float lep_picky_py[2];
    float lep_picky_pz[2];
    float lep_picky_eta[2];
    float lep_picky_phi[2];
    float lep_picky_chi2[2];
    float lep_picky_ndf[2];
    float lep_picky_qOverPt[2];
    float lep_cocktail_p[2];
    float lep_cocktail_pt[2];
    float lep_cocktail_pt_err[2];
    float lep_cocktail_px[2];
    float lep_cocktail_py[2];
    float lep_cocktail_pz[2];
    float lep_cocktail_eta[2];
    float lep_cocktail_phi[2];
    float lep_cocktail_chi2[2];
    float lep_cocktail_ndf[2];
    float lep_cocktail_qOverPt[2];
    short lep_cocktail_choice[2];
    float lep_tuneP_p[2];
    float lep_tuneP_pt[2];
    float lep_tuneP_pt_err[2];
    float lep_tuneP_px[2];
    float lep_tuneP_py[2];
    float lep_tuneP_pz[2];
    float lep_tuneP_eta[2];
    float lep_tuneP_phi[2];
    float lep_tuneP_dz[2];
    float lep_tuneP_vz[2];
    float lep_tuneP_chi2[2];
    float lep_tuneP_ndf[2];
    float lep_tuneP_qOverPt[2];
    float lep_triggerMatchPt[2];
    float lep_triggerMatchEta[2];
    float lep_chi2dof[2];
    float lep_dB[2];
    float lep_sumPt[2];
    float lep_emEt[2];
    float lep_hadEt[2];
    float lep_hoEt[2];
    float lep_pfIso[2];
    float lep_pfIsoDB[2];
    int lep_timeNdof[2];
    float lep_timeInOut[2];
    float lep_timeOutIn[2];
    float lep_timeInOutErr[2];
    float lep_timeOutInErr[2];
    int lep_heep_id[2];
    float lep_min_muon_dR[2];
    short lep_tk_numberOfValidTrackerHits[2];
    short lep_tk_numberOfValidTrackerLayers[2];
    short lep_tk_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidTrackerHits[2]; 
    short lep_glb_numberOfValidTrackerLayers[2]; 
    short lep_glb_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidMuonHits[2];
    short lep_glb_numberOfValidMuonDTHits[2];
    short lep_glb_numberOfValidMuonCSCHits[2];
    short lep_glb_numberOfValidMuonRPCHits[2];
    short lep_glb_muonStationsWithValidHits[2];
    short lep_glb_dtStationsWithValidHits[2];
    short lep_glb_cscStationsWithValidHits[2];
    short lep_glb_rpcStationsWithValidHits[2];
    short lep_glb_innermostMuonStationWithValidHits[2];
    short lep_glb_outermostMuonStationWithValidHits[2];
    short lep_numberOfMatches[2];
    short lep_numberOfMatchedStations[2];
    unsigned int lep_stationMask[2];
    int lep_numberOfChambers[2];
    int lep_numberOfChambersNoRPC[2];
    unsigned int lep_stationGapMaskDistance[2];
    unsigned int lep_stationGapMaskPull[2];
    bool lep_isGlobalMuon[2];
    bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool GoodVtx;
    bool METFilter;
    int nJets;
    float met_pt;
    float met_phi;
    float jet_pt[4];
    float jet_eta[4];
    float jet_phi[4];
    float gen_res_mass;
    float gen_res_pt;
    float gen_res_rap;
    float gen_res_eta;
    float gen_res_phi;
    float gen_dil_mass;
    float gen_dil_pt;
    float gen_dil_rap;
    float gen_dil_eta;
    float gen_dil_phi;
    float gen_dil_dR;
    float gen_dil_dPhi;
    float gen_lep_p[2];
    float gen_lep_pt[2];
    float gen_lep_px[2];
    float gen_lep_py[2];
    float gen_lep_pz[2];
    float gen_lep_E[2];
    float gen_lep_eta[2];
    float gen_lep_phi[2];
    float gen_lep_qOverPt[2];
    float gen_lep_noib_p[2];
    float gen_lep_noib_pt[2];
    float gen_lep_noib_px[2];
    float gen_lep_noib_py[2];
    float gen_lep_noib_pz[2];
    float gen_lep_noib_E[2];
    float gen_lep_noib_eta[2];
    float gen_lep_noib_phi[2];
    float gen_lep_noib_qOverPt[2];
  };

  tree_t t;
  TTree* tree;

  const edm::InputTag dimu_src;
  const edm::InputTag beamspot_src;
  const edm::InputTag met_src;
  const edm::InputTag jet_src;
  const edm::InputTag vertices_src;
  const bool fill_gen_info;
  const edm::InputTag TriggerResults_src;
  const edm::InputTag genEventInfo_;
  HardInteraction* hardInteraction;
  //const edm::InputTag TriggerResults_src; 
};

TString replace_all(const TString& a, const TString& b, const TString& c) {
  TString ret = a;
  ret.ReplaceAll(b, c);
  return ret;
}

SimpleNtupler::SimpleNtupler(const edm::ParameterSet& cfg)
  : dimu_src(cfg.getParameter<edm::InputTag>("dimu_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    met_src(cfg.getParameter<edm::InputTag>("met_src")),
    jet_src(cfg.getParameter<edm::InputTag>("jet_src")),
    vertices_src(cfg.getParameter<edm::InputTag>("vertices_src")),
    fill_gen_info(cfg.existsAs<edm::ParameterSet>("hardInteraction")),
    TriggerResults_src(cfg.getParameter<edm::InputTag>("TriggerResults_src")),
    genEventInfo_(cfg.getUntrackedParameter<edm::InputTag>("genEventInfo")),
    hardInteraction(fill_gen_info ? new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")) : 0)
   
{
 
  consumes<pat::CompositeCandidateCollection>(dimu_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<std::vector<pat::MET>>(met_src);
  consumes<std::vector<pat::Jet>>(jet_src);
  consumes<reco::VertexCollection>(vertices_src);
  consumes<edm::TriggerResults>(TriggerResults_src);
  consumes<GenEventInfoProduct>(genEventInfo_);
  if (fill_gen_info) consumes<reco::GenParticleCollection>(hardInteraction->src);
  
 
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("run", &t.run, "run/i");
  tree->Branch("lumi", &t.lumi, "lumi/i");
  tree->Branch("event", &t.event, "event/i");
  tree->Branch("beamspot_x", &t.beamspot_x, "beamspot_x/F");
  tree->Branch("beamspot_x_err", &t.beamspot_x_err, "beamspot_x_err/F");
  tree->Branch("beamspot_y", &t.beamspot_y, "beamspot_y/F");
  tree->Branch("beamspot_y_err", &t.beamspot_y_err, "beamspot_y_err/F");
  tree->Branch("beamspot_z", &t.beamspot_z, "beamspot_z/F");
  tree->Branch("beamspot_z_err", &t.beamspot_z_err, "beamspot_z_err/F");
  tree->Branch("nvertices", &t.nvertices, "nvertices/I");
  tree->Branch("dil_chosen", &t.dil_chosen, "dil_chosen/I");
  tree->Branch("dil_mass", &t.dil_mass, "dil_mass/F");
  tree->Branch("dil_pt", &t.dil_pt, "dil_pt/F");
  tree->Branch("dil_rap", &t.dil_rap, "dil_rap/F");
  tree->Branch("dil_eta", &t.dil_eta, "dil_eta/F");
  tree->Branch("dil_phi", &t.dil_phi, "dil_phi/F");
  tree->Branch("dil_dR", &t.dil_dR, "dil_dR/F");
  tree->Branch("dil_dPhi", &t.dil_dPhi, "dil_dPhi/F");
  tree->Branch("cos_angle", &t.cos_angle, "cos_angle/F");
  tree->Branch("vertex_chi2", &t.vertex_chi2, "vertex_chi2/F");
  tree->Branch("cos_cs", &t.cos_cs, "cos_cs/F");
  tree->Branch("phi_cs", &t.phi_cs, "phi_cs/F");
  tree->Branch("vertex_m", &t.vertex_m, "vertex_m/F");
  tree->Branch("vertex_m_err", &t.vertex_m_err, "vertex_m_err/F");
  tree->Branch("vertex_x", &t.vertex_x, "vertex_x/F");
  tree->Branch("vertex_x_err", &t.vertex_x_err, "vertex_x_err/F");
  tree->Branch("vertex_y", &t.vertex_y, "vertex_y/F");
  tree->Branch("vertex_y_err", &t.vertex_y_err, "vertex_y_err/F");
  tree->Branch("vertex_z", &t.vertex_z, "vertex_z/F");
  tree->Branch("vertex_z_err", &t.vertex_z_err, "vertex_z_err/F");
  tree->Branch("lep_id", t.lep_id, "lep_id[2]/I");
  tree->Branch("lep_heep_id", t.lep_heep_id, "lep_heep_id[2]/I");
  tree->Branch("lep_p", t.lep_p, "lep_p[2]/F");
  tree->Branch("lep_pt", t.lep_pt, "lep_pt[2]/F");
  tree->Branch("lep_pt_err", t.lep_pt_err, "lep_pt_err[2]/F");
  tree->Branch("lep_px", t.lep_px, "lep_px[2]/F");
  tree->Branch("lep_py", t.lep_py, "lep_py[2]/F");
  tree->Branch("lep_pz", t.lep_pz, "lep_pz[2]/F");
  tree->Branch("lep_E", t.lep_E, "lep_E[2]/F");
  tree->Branch("lep_eta", t.lep_eta, "lep_eta[2]/F");
  tree->Branch("lep_phi", t.lep_phi, "lep_phi[2]/F");
  tree->Branch("lep_qOverPt", t.lep_qOverPt, "lep_qOverPt[2]/F");
  tree->Branch("lep_tk_p", t.lep_tk_p, "lep_tk_p[2]/F");
  tree->Branch("lep_tk_pt", t.lep_tk_pt, "lep_tk_pt[2]/F");
  tree->Branch("lep_tk_pt_err", t.lep_tk_pt_err, "lep_tk_pt_err[2]/F");
  tree->Branch("lep_tk_px", t.lep_tk_px, "lep_tk_px[2]/F");
  tree->Branch("lep_tk_py", t.lep_tk_py, "lep_tk_py[2]/F");
  tree->Branch("lep_tk_pz", t.lep_tk_pz, "lep_tk_pz[2]/F");
  tree->Branch("lep_tk_eta", t.lep_tk_eta, "lep_tk_eta[2]/F");
  tree->Branch("lep_tk_phi", t.lep_tk_phi, "lep_tk_phi[2]/F");
  tree->Branch("lep_tk_dz", t.lep_tk_dz, "lep_tk_dz[2]/F");
  tree->Branch("lep_tk_vz", t.lep_tk_vz, "lep_tk_vz[2]/F");
  tree->Branch("lep_tk_chi2", t.lep_tk_chi2, "lep_tk_chi2[2]/F");
  tree->Branch("lep_tk_ndf", t.lep_tk_ndf, "lep_tk_ndf[2]/F");
  tree->Branch("lep_tk_qOverPt", t.lep_tk_qOverPt, "lep_tk_qOverPt[2]/F");
  tree->Branch("lep_glb_p", t.lep_glb_p, "lep_glb_p[2]/F");
  tree->Branch("lep_glb_pt", t.lep_glb_pt, "lep_glb_pt[2]/F");
  tree->Branch("lep_glb_pt_err", t.lep_glb_pt_err, "lep_glb_pt_err[2]/F");
  tree->Branch("lep_glb_px", t.lep_glb_px, "lep_glb_px[2]/F");
  tree->Branch("lep_glb_py", t.lep_glb_py, "lep_glb_py[2]/F");
  tree->Branch("lep_glb_pz", t.lep_glb_pz, "lep_glb_pz[2]/F");
  tree->Branch("lep_glb_eta", t.lep_glb_eta, "lep_glb_eta[2]/F");
  tree->Branch("lep_glb_phi", t.lep_glb_phi, "lep_glb_phi[2]/F");
  tree->Branch("lep_glb_chi2", t.lep_glb_chi2, "lep_glb_chi2[2]/F");
  tree->Branch("lep_glb_ndf", t.lep_glb_ndf, "lep_glb_ndf[2]/F");
  tree->Branch("lep_glb_qOverPt", t.lep_glb_qOverPt, "lep_glb_qOverPt[2]/F");
  tree->Branch("lep_tpfms_p", t.lep_tpfms_p, "lep_tpfms_p[2]/F");
  tree->Branch("lep_tpfms_pt", t.lep_tpfms_pt, "lep_tpfms_pt[2]/F");
  tree->Branch("lep_tpfms_pt_err", t.lep_tpfms_pt_err, "lep_tpfms_pt_err[2]/F");
  tree->Branch("lep_tpfms_px", t.lep_tpfms_px, "lep_tpfms_px[2]/F");
  tree->Branch("lep_tpfms_py", t.lep_tpfms_py, "lep_tpfms_py[2]/F");
  tree->Branch("lep_tpfms_pz", t.lep_tpfms_pz, "lep_tpfms_pz[2]/F");
  tree->Branch("lep_tpfms_eta", t.lep_tpfms_eta, "lep_tpfms_eta[2]/F");
  tree->Branch("lep_tpfms_phi", t.lep_tpfms_phi, "lep_tpfms_phi[2]/F");
  tree->Branch("lep_tpfms_chi2", t.lep_tpfms_chi2, "lep_tpfms_chi2[2]/F");
  tree->Branch("lep_tpfms_ndf", t.lep_tpfms_ndf, "lep_tpfms_ndf[2]/F");
  tree->Branch("lep_tpfms_qOverPt", t.lep_tpfms_qOverPt, "lep_tpfms_qOverPt[2]/F");
  tree->Branch("lep_picky_p", t.lep_picky_p, "lep_picky_p[2]/F");
  tree->Branch("lep_picky_pt", t.lep_picky_pt, "lep_picky_pt[2]/F");
  tree->Branch("lep_picky_pt_err", t.lep_picky_pt_err, "lep_picky_pt_err[2]/F");
  tree->Branch("lep_picky_px", t.lep_picky_px, "lep_picky_px[2]/F");
  tree->Branch("lep_picky_py", t.lep_picky_py, "lep_picky_py[2]/F");
  tree->Branch("lep_picky_pz", t.lep_picky_pz, "lep_picky_pz[2]/F");
  tree->Branch("lep_picky_eta", t.lep_picky_eta, "lep_picky_eta[2]/F");
  tree->Branch("lep_picky_phi", t.lep_picky_phi, "lep_picky_phi[2]/F");
  tree->Branch("lep_picky_chi2", t.lep_picky_chi2, "lep_picky_chi2[2]/F");
  tree->Branch("lep_picky_ndf", t.lep_picky_ndf, "lep_picky_ndf[2]/F");
  tree->Branch("lep_picky_qOverPt", t.lep_picky_qOverPt, "lep_picky_qOverPt[2]/F");
  tree->Branch("lep_cocktail_p", t.lep_cocktail_p, "lep_cocktail_p[2]/F");
  tree->Branch("lep_cocktail_pt", t.lep_cocktail_pt, "lep_cocktail_pt[2]/F");
  tree->Branch("lep_cocktail_pt_err", t.lep_cocktail_pt_err, "lep_cocktail_pt_err[2]/F");
  tree->Branch("lep_cocktail_px", t.lep_cocktail_px, "lep_cocktail_px[2]/F");
  tree->Branch("lep_cocktail_py", t.lep_cocktail_py, "lep_cocktail_py[2]/F");
  tree->Branch("lep_cocktail_pz", t.lep_cocktail_pz, "lep_cocktail_pz[2]/F");
  tree->Branch("lep_cocktail_eta", t.lep_cocktail_eta, "lep_cocktail_eta[2]/F");
  tree->Branch("lep_cocktail_phi", t.lep_cocktail_phi, "lep_cocktail_phi[2]/F");
  tree->Branch("lep_cocktail_chi2", t.lep_cocktail_chi2, "lep_cocktail_chi2[2]/F");
  tree->Branch("lep_cocktail_ndf", t.lep_cocktail_ndf, "lep_cocktail_ndf[2]/F");
  tree->Branch("lep_cocktail_qOverPt", t.lep_cocktail_qOverPt, "lep_cocktail_qOverPt[2]/F");
  tree->Branch("lep_cocktail_choice", t.lep_cocktail_choice, "lep_cocktail_choice[2]/S");
  tree->Branch("lep_tuneP_p", t.lep_tuneP_p, "lep_tuneP_p[2]/F");
  tree->Branch("lep_tuneP_pt", t.lep_tuneP_pt, "lep_tuneP_pt[2]/F");
  tree->Branch("lep_tuneP_pt_err", t.lep_tuneP_pt_err, "lep_tuneP_pt_err[2]/F");
  tree->Branch("lep_tuneP_px", t.lep_tuneP_px, "lep_tuneP_px[2]/F");
  tree->Branch("lep_tuneP_py", t.lep_tuneP_py, "lep_tuneP_py[2]/F");
  tree->Branch("lep_tuneP_pz", t.lep_tuneP_pz, "lep_tuneP_pz[2]/F");
  tree->Branch("lep_tuneP_eta", t.lep_tuneP_eta, "lep_tuneP_eta[2]/F");
  tree->Branch("lep_tuneP_phi", t.lep_tuneP_phi, "lep_tuneP_phi[2]/F");
  tree->Branch("lep_tuneP_chi2", t.lep_tuneP_chi2, "lep_tuneP_chi2[2]/F");
  tree->Branch("lep_tuneP_ndf", t.lep_tuneP_ndf, "lep_tuneP_ndf[2]/F");
  tree->Branch("lep_tuneP_qOverPt", t.lep_tuneP_qOverPt, "lep_tuneP_qOverPt[2]/F");
  tree->Branch("lep_triggerMatchPt", t.lep_triggerMatchPt, "lep_triggerMatchPt[2]/F");
  tree->Branch("lep_triggerMatchEta", t.lep_triggerMatchEta, "lep_triggerMatchEta[2]/F");
  tree->Branch("lep_chi2dof", t.lep_chi2dof, "lep_chi2dof[2]/F");
  tree->Branch("lep_dB", t.lep_dB, "lep_dB[2]/F");
  tree->Branch("lep_sumPt", t.lep_sumPt, "lep_sumPt[2]/F");
  tree->Branch("lep_emEt", t.lep_emEt, "lep_emEt[2]/F");
  tree->Branch("lep_hadEt", t.lep_hadEt, "lep_hadEt[2]/F");
  tree->Branch("lep_hoEt", t.lep_hoEt, "lep_hoEt[2]/F");
  tree->Branch("lep_pfIso", t.lep_pfIso, "lep_pfIso[2]/F");
  tree->Branch("lep_pfIsoDB", t.lep_pfIsoDB, "lep_pfIsoDB[2]/F");
  tree->Branch("lep_timeNdof", t.lep_timeNdof, "lep_timeNdof[2]/I");
  tree->Branch("lep_timeInOut", t.lep_timeInOut, "lep_timeInOut[2]/F");
  tree->Branch("lep_timeOutIn", t.lep_timeOutIn, "lep_timeOutIn[2]/F");
  tree->Branch("lep_timeInOutErr", t.lep_timeInOutErr, "lep_timeInOutErr[2]/F");
  tree->Branch("lep_timeOutInErr", t.lep_timeOutInErr, "lep_timeOutInErr[2]/F");
  tree->Branch("lep_min_muon_dR", t.lep_min_muon_dR, "lep_min_muon_dR[2]/F");
  tree->Branch("lep_tk_numberOfValidTrackerHits", t.lep_tk_numberOfValidTrackerHits, "lep_tk_numberOfValidTrackerHits[2]/S");
  tree->Branch("lep_tk_numberOfValidTrackerLayers", t.lep_tk_numberOfValidTrackerLayers, "lep_tk_numberOfValidTrackerLayers[2]/S");
  tree->Branch("lep_tk_numberOfValidPixelHits", t.lep_tk_numberOfValidPixelHits, "lep_tk_numberOfValidPixelHits[2]/S");
  tree->Branch("lep_glb_numberOfValidTrackerHits", t.lep_glb_numberOfValidTrackerHits, "lep_glb_numberOfValidTrackerHits[2]/S");
  tree->Branch("lep_glb_numberOfValidTrackerLayers", t.lep_glb_numberOfValidTrackerLayers, "lep_glb_numberOfValidTrackerLayers[2]/S");
  tree->Branch("lep_glb_numberOfValidPixelHits", t.lep_glb_numberOfValidPixelHits, "lep_glb_numberOfValidPixelHits[2]/S");
  tree->Branch("lep_glb_numberOfValidMuonHits", t.lep_glb_numberOfValidMuonHits, "lep_glb_numberOfValidMuonHits[2]/S");
  tree->Branch("lep_glb_numberOfValidMuonDTHits", t.lep_glb_numberOfValidMuonDTHits, "lep_glb_numberOfValidMuonDTHits[2]/S");
  tree->Branch("lep_glb_numberOfValidMuonCSCHits", t.lep_glb_numberOfValidMuonCSCHits, "lep_glb_numberOfValidMuonCSCHits[2]/S");
  tree->Branch("lep_glb_numberOfValidMuonRPCHits", t.lep_glb_numberOfValidMuonRPCHits, "lep_glb_numberOfValidMuonRPCHits[2]/S");
  tree->Branch("lep_glb_muonStationsWithValidHits", t.lep_glb_muonStationsWithValidHits, "lep_glb_muonStationsWithValidHits[2]/S");
  tree->Branch("lep_glb_dtStationsWithValidHits", t.lep_glb_dtStationsWithValidHits, "lep_glb_dtStationsWithValidHits[2]/S");
  tree->Branch("lep_glb_cscStationsWithValidHits", t.lep_glb_cscStationsWithValidHits, "lep_glb_cscStationsWithValidHits[2]/S");
  tree->Branch("lep_glb_rpcStationsWithValidHits", t.lep_glb_rpcStationsWithValidHits, "lep_glb_rpcStationsWithValidHits[2]/S");
  tree->Branch("lep_glb_innermostMuonStationWithValidHits", t.lep_glb_innermostMuonStationWithValidHits, "lep_glb_innermostMuonStationWithValidHits[2]/S");
  tree->Branch("lep_glb_outermostMuonStationWithValidHits", t.lep_glb_outermostMuonStationWithValidHits, "lep_glb_outermostMuonStationWithValidHits[2]/S");
  tree->Branch("lep_numberOfMatches", t.lep_numberOfMatches, "lep_numberOfMatches[2]/S");
  tree->Branch("lep_numberOfMatchedStations", t.lep_numberOfMatchedStations, "lep_numberOfMatchedStations[2]/S");
  tree->Branch("lep_stationMask", t.lep_stationMask, "lep_stationMask[2]/I");
  tree->Branch("lep_numberOfChambers", t.lep_numberOfChambers, "lep_numberOfChambers[2]/I");
  tree->Branch("lep_numberOfChambersNoRPC", t.lep_numberOfChambersNoRPC, "lep_numberOfChambersNoRPC[2]/I");
  tree->Branch("lep_stationGapMaskDistance", t.lep_stationGapMaskDistance, "lep_stationGapMaskDistance[2]/I");
  tree->Branch("lep_stationGapMaskPull", t.lep_stationGapMaskPull, "lep_stationGapMaskPull[2]/I");
  tree->Branch("lep_isGlobalMuon", t.lep_isGlobalMuon, "lep_isGlobalMuon[2]/O");
  tree->Branch("lep_isTrackerMuon", t.lep_isTrackerMuon, "lep_isTrackerMuon[2]/O");
  tree->Branch("GoodDataRan", &t.GoodDataRan, "GoodDataRan/O");
  tree->Branch("GoodVtx", &t.GoodVtx, "GoodVtx/O");
  tree->Branch("METFilter", &t.METFilter, "METFilter/O");
  tree->Branch("met_pt", &t.met_pt, "met_pt/F");
  tree->Branch("met_phi", &t.met_phi, "met_phi/F");
  tree->Branch("nJets", &t.nJets, "nJets/I");
  tree->Branch("jet_pt", t.jet_pt, "jet_pt[4]/F");
  tree->Branch("jet_eta", t.jet_eta, "jet_eta[4]/F");
  tree->Branch("jet_phi", t.jet_phi, "jet_phi[4]/F");
  if (fill_gen_info) {
    tree->Branch("genWeight", &t.genWeight, "genWeight/F");
    tree->Branch("gen_res_mass", &t.gen_res_mass, "gen_res_mass/F");
    tree->Branch("gen_res_pt", &t.gen_res_pt, "gen_res_pt/F");
    tree->Branch("gen_res_rap", &t.gen_res_rap, "gen_res_rap/F");
    tree->Branch("gen_res_eta", &t.gen_res_eta, "gen_res_eta/F");
    tree->Branch("gen_res_phi", &t.gen_res_phi, "gen_res_phi/F");
    tree->Branch("gen_dil_mass", &t.gen_dil_mass, "gen_dil_mass/F");
    tree->Branch("gen_dil_pt", &t.gen_dil_pt, "gen_dil_pt/F");
    tree->Branch("gen_dil_rap", &t.gen_dil_rap, "gen_dil_rap/F");
    tree->Branch("gen_dil_eta", &t.gen_dil_eta, "gen_dil_eta/F");
    tree->Branch("gen_dil_phi", &t.gen_dil_phi, "gen_dil_phi/F");
    tree->Branch("gen_dil_dR", &t.gen_dil_dR, "gen_dil_dR/F");
    tree->Branch("gen_dil_dPhi", &t.gen_dil_dPhi, "gen_dil_dPhi/F");
    tree->Branch("gen_lep_p", t.gen_lep_p, "gen_lep_p[2]/F");
    tree->Branch("gen_lep_pt", t.gen_lep_pt, "gen_lep_pt[2]/F");
    tree->Branch("gen_lep_px", t.gen_lep_px, "gen_lep_px[2]/F");
    tree->Branch("gen_lep_py", t.gen_lep_py, "gen_lep_py[2]/F");
    tree->Branch("gen_lep_pz", t.gen_lep_pz, "gen_lep_pz[2]/F");
    tree->Branch("gen_lep_E", t.gen_lep_E, "gen_lep_E[2]/F");
    tree->Branch("gen_lep_eta", t.gen_lep_eta, "gen_lep_eta[2]/F");
    tree->Branch("gen_lep_phi", t.gen_lep_phi, "gen_lep_phi[2]/F");
    tree->Branch("gen_lep_qOverPt", t.gen_lep_qOverPt, "gen_lep_qOverPt[2]/F");
    tree->Branch("gen_lep_noib_pt", t.gen_lep_noib_pt, "gen_lep_noib_pt[2]/F");
    tree->Branch("gen_lep_noib_px", t.gen_lep_noib_px, "gen_lep_noib_px[2]/F");
    tree->Branch("gen_lep_noib_py", t.gen_lep_noib_py, "gen_lep_noib_py[2]/F");
    tree->Branch("gen_lep_noib_pz", t.gen_lep_noib_pz, "gen_lep_noib_pz[2]/F");
    tree->Branch("gen_lep_noib_E", t.gen_lep_noib_E, "gen_lep_noib_E[2]/F");
    tree->Branch("gen_lep_noib_eta", t.gen_lep_noib_eta, "gen_lep_noib_eta[2]/F");
    tree->Branch("gen_lep_noib_phi", t.gen_lep_noib_phi, "gen_lep_noib_phi[2]/F");
    tree->Branch("gen_lep_noib_qOverPt", t.gen_lep_noib_qOverPt, "gen_lep_noib_qOverPt[2]/F");
  }

  tree->SetAlias("OppSign",  "lep_id[0]*lep_id[1] < 0");
  tree->SetAlias("SameSign", "lep_id[0]*lep_id[1] > 0");
  tree->SetAlias("Dimu",     "abs(lep_id[0]*lep_id[1]) == 169");
  tree->SetAlias("Emu",      "abs(lep_id[0]*lep_id[1]) == 143");

#define offlineMinPt "53"
#define triggerMatchMinPt "50"
#define triggerMatchMaxEta "2.1"

//  tree->SetAlias("trigger_match_0", "lep_triggerMatchPt[0] > " triggerMatchMinPt " && abs(lep_triggerMatchEta[0]) < " triggerMatchMaxEta);
//  tree->SetAlias("trigger_match_1", "lep_triggerMatchPt[1] > " triggerMatchMinPt " && abs(lep_triggerMatchEta[1]) < " triggerMatchMaxEta);
  tree->SetAlias("trigger_match_0", "lep_triggerMatchPt[0] > " triggerMatchMinPt );
  tree->SetAlias("trigger_match_1", "lep_triggerMatchPt[1] > " triggerMatchMinPt );
  tree->SetAlias("triggerMatched", "trigger_match_0 || trigger_match_1");

  // tree->SetAlias("GoodData", "GoodDataRan && HLTPhysicsDeclared && NoScraping && GoodVtx");
  tree->SetAlias("GoodData", "GoodDataRan && HLTPhysicsDeclared && GoodVtx");

  tree->SetAlias("extraDimuonCuts", "cos_angle > -0.9998 && vertex_chi2 < 20");

  TString loose_2010 =
    "lep_isGlobalMuon[X] && "						\
    "lep_pt[X] > " offlineMinPt " && "					\
    "lep_tk_numberOfValidTrackerHits[X] >= 10 && "			\
    "lep_sumPt[X] / lep_tk_pt[X] < 0.1";

  TString tight_2010 =
    "abs(lep_dB[X]) < 0.2 && "						\
    "lep_chi2dof[X] < 10 && "						\
    "lep_tk_numberOfValidPixelHits[X] >= 1 && "				\
    "lep_glb_muonStationsWithValidHits[X] >= 2 && "			\
    "lep_isTrackerMuon[X] && "						\
    "lep_triggerMatchPt[X] > " triggerMatchMinPt;

  TString vbtf =
    "lep_isGlobalMuon[X] && "						\
    "lep_isTrackerMuon[X] && "						\
    "lep_tk_pt[X] > " offlineMinPt " && "				\
    "abs(lep_tk_eta[X]) < 2.1 && "					\
    "abs(lep_dB[X]) < 0.2 && "						\
    "lep_sumPt[X] < 3 && "						\
    "lep_glb_numberOfValidTrackerHits[X] > 10 && "			\
    "lep_glb_numberOfValidPixelHits[X] > 0 && "				\
    "lep_glb_numberOfValidMuonHits[X] > 0 && "				\
    "lep_numberOfMatches[X] > 1";

  TString loose_no_iso =
    "lep_isGlobalMuon[X] && "						\
    "lep_isTrackerMuon[X] && "						\
    "lep_pt[X] > " offlineMinPt " && "					\
    "abs(lep_dB[X]) < 0.2 && "						\
    "lep_glb_numberOfValidTrackerLayers[X] > 5 && "			\
    "lep_glb_numberOfValidPixelHits[X] >= 1 && "			\
    "lep_glb_numberOfValidMuonHits[X] > 0 && "				\
    "lep_numberOfMatchedStations[X] > 1 && "                            \
    "lep_pt_err[X] / lep_pt[X] < 0.3";

  TString tight_2015 =
    "lep_isGlobalMuon[X] && "						\
    "lep_isTrackerMuon[X] && "						\
    "lep_tuneP_pt[X] > " offlineMinPt " && "				\
    "abs(lep_dB[X]) < 0.2 && "						\
    "lep_glb_numberOfValidTrackerLayers[X] > 5 && "			\
    "lep_glb_numberOfValidPixelHits[X] >= 1 && "			\
    "lep_glb_numberOfValidMuonHits[X] > 0 && "				\
    "lep_numberOfMatchedStations[X] > 1 && "                            \
    "lep_tuneP_pt_err[X] / lep_tuneP_pt[X] < 0.3 && " 			\
    "lep_sumPt[X] / lep_tk_pt[X] < 0.1";

  TString loose_no_iso_pt_ptErr = 
    "lep_isGlobalMuon[X] && "						\
    "lep_isTrackerMuon[X] && " 						\
    "lep_glb_numberOfValidTrackerLayers[X] > 5 && "			\
    "lep_glb_numberOfValidPixelHits[X] >= 1 && "			\
    "lep_glb_numberOfValidMuonHits[X] > 0 && " 				\
    "lep_numberOfMatchedStations[X] > 1 && "				\
    "abs(lep_dB[X]) < 0.2";

  TString loose_2012 = loose_no_iso + " && lep_sumPt[X] / lep_tk_pt[X] < 0.1";

  TString loose_new(loose_2012);
  loose_new.ReplaceAll("lep_glb_numberOfValidTrackerLayers[X] > 5", 
		       "lep_glb_numberOfValidTrackerLayers[X] > 8");
  loose_new.ReplaceAll(" && lep_pt_err[X] / lep_pt[X] < 0.3", "");

  TString loose_2011eps(loose_new);
  loose_2011eps.ReplaceAll("lep_glb_numberOfValidTrackerLayers", 
			   "lep_glb_numberOfValidTrackerHits");

  tree->SetAlias("loose_2010_0",    replace_all(loose_2010,    "[X]", "[0]"));
  tree->SetAlias("loose_2010_1",    replace_all(loose_2010,    "[X]", "[1]"));
  tree->SetAlias("tight_2010_0",    replace_all(tight_2010,    "[X]", "[0]"));
  tree->SetAlias("tight_2010_1",    replace_all(tight_2010,    "[X]", "[1]"));
  tree->SetAlias("vbtf_0",          replace_all(vbtf,          "[X]", "[0]"));
  tree->SetAlias("vbtf_1",          replace_all(vbtf,          "[X]", "[1]"));
  tree->SetAlias("loose_2011eps_0", replace_all(loose_2011eps, "[X]", "[0]"));
  tree->SetAlias("loose_2011eps_1", replace_all(loose_2011eps, "[X]", "[1]"));
  tree->SetAlias("loose_new_0",     replace_all(loose_new,     "[X]", "[0]"));
  tree->SetAlias("loose_new_1",     replace_all(loose_new,     "[X]", "[1]"));
  tree->SetAlias("loose_no_iso_0",  replace_all(loose_no_iso,  "[X]", "[0]"));
  tree->SetAlias("loose_no_iso_1",  replace_all(loose_no_iso,  "[X]", "[1]"));
  tree->SetAlias("loose_2012_0",    replace_all(loose_2012,    "[X]", "[0]"));
  tree->SetAlias("loose_2012_1",    replace_all(loose_2012,    "[X]", "[1]"));
  tree->SetAlias("tight_2015_0",    replace_all(tight_2015,    "[X]", "[0]"));
  tree->SetAlias("tight_2015_1",    replace_all(tight_2015,    "[X]", "[1]"));
  tree->SetAlias("loose_no_iso_pt_ptErr_0", replace_all(loose_no_iso_pt_ptErr, "[X]", "[0]"));
  tree->SetAlias("loose_no_iso_pt_ptErr_1", replace_all(loose_no_iso_pt_ptErr, "[X]", "[1]"));

  tree->SetAlias("OurSel2010",
		 "loose_2010_0 && loose_2010_1 && "			\
		 "(tight_2010_0 || tight_2010_1) && "			\
		 "OppSign && "						\
		 "extraDimuonCuts && "					\
		 "GoodData");

  tree->SetAlias("VBTFSel",
		 "vbtf_0 && vbtf_1 && "					\
		 "triggerMatched && "					\
		 "OppSign");

  tree->SetAlias("OurSel2011EPS",
		 "loose_2011eps_0 && loose_2011eps_1 && "		\
		 "triggerMatched && "					\
		 "OppSign && "						\
		 "extraDimuonCuts && "					\
		 "GoodData");

  tree->SetAlias("OurSelNewNoSign",
		 "loose_new_0 && loose_new_1 && "			\
		 "triggerMatched && "					\
		 "extraDimuonCuts && "					\
		 "GoodData");

  tree->SetAlias("OurSelNew",   "OurSelNewNoSign && OppSign");
  tree->SetAlias("OurSelNewSS", "OurSelNewNoSign && SameSign");

  tree->SetAlias("OurSel2012NoSign",
		 "loose_2012_0 && loose_2012_1 && "			\
		 "triggerMatched && "					\
		 "extraDimuonCuts && "					\
		 "GoodData");

  tree->SetAlias("OurSel2012",   "OurSel2012NoSign && OppSign");
  tree->SetAlias("OurSel2012SS", "OurSel2012NoSign && SameSign");

  tree->SetAlias("OurSel2015",
                 "OppSign && "						\
                 "tight_2015_0 && tight_2015_1 && "			\
                 "triggerMatched && " 					\
                 "extraDimuonCuts && " 					\
                 "GoodData");

  // For e-mu dileptons, below we always put the muon in [0] and the
  // electron in [1], so don't have to check the other combination.
  tree->SetAlias("EmuSelNoSign",
		 "abs(lep_id[1]) == 11 && "				\
		 "lep_heep_id[1] == 0 && "				\
		 "loose_2012_0 && "					\
		 "trigger_match_0 && "					\
		 "GoodData");

  tree->SetAlias("EmuSel", "EmuSelNoSign && OppSign");
}

template <typename T>
float userFloat(const T& patobj, const char* name, float def=-999.) {
  return patobj.hasUserFloat(name) ? patobj.userFloat(name) : def;
}

template <typename T>
int userInt(const T& patobj, const char* name, int def=-999) {
  return patobj.hasUserInt(name) ? patobj.userInt(name) : def;
}

void SimpleNtupler::analyze(const edm::Event& event, const edm::EventSetup&) {
  memset(&t, 0, sizeof(tree_t));

  //
  // Event Information
  //
  t.run = event.id().run();
  t.lumi = event.luminosityBlock();
  t.event = event.id().event();

  // Get Trigger information




  edm::Handle<edm::TriggerResults> respat;
  //event.getByLabel(edm::InputTag("TriggerResults", "", "PAT"), respat);
  event.getByLabel(TriggerResults_src, respat);
  
  const edm::TriggerNames& namespat = event.triggerNames(*respat);
  
  if (namespat.triggerIndex("goodDataPrimaryVertexFilter") < respat->size()) {
    t.GoodDataRan = 1;
    t.GoodVtx = respat->accept(namespat.triggerIndex("goodDataPrimaryVertexFilter"));
    t.METFilter = respat->accept(namespat.triggerIndex("goodDataMETFilter"));
  }

  // Get Beamspot information
  edm::Handle<reco::BeamSpot> bs;
  event.getByLabel(beamspot_src, bs);
  t.beamspot_x     = bs->x0();
  t.beamspot_x_err = bs->x0Error();
  t.beamspot_y     = bs->y0();
  t.beamspot_y_err = bs->y0Error();
  t.beamspot_z     = bs->z0();
  t.beamspot_z_err = bs->z0Error();

  // Get Vertex information
  edm::Handle<reco::VertexCollection> pvs;
  event.getByLabel(vertices_src, pvs);
  t.nvertices = 0;
  BOOST_FOREACH(const reco::Vertex& vtx, *pvs)
    if (vtx.ndof() > 4 && fabs(vtx.z()) <= 24 && fabs(vtx.position().rho()) <= 2)
      t.nvertices += 1;
  

  if (fill_gen_info) {

    // This only works for DY/Z'/RSG events, and really just for PYTHIA!
    hardInteraction->Fill(event);

   int EventWeight = 1.;
   edm::Handle<GenEventInfoProduct> gen_ev_info;
   event.getByLabel(genEventInfo_, gen_ev_info);
   EventWeight = gen_ev_info->weight();
   t.genWeight = ( EventWeight > 0 ) ? 1 : -1;


    //
    // Store Generator Level information
    //
    if(hardInteraction->IsValid()){

        t.gen_res_mass = hardInteraction->resonance->mass();
        t.gen_res_pt   = hardInteraction->resonance->pt();
        t.gen_res_rap  = hardInteraction->resonance->rapidity();
        t.gen_res_eta  = hardInteraction->resonance->eta();
        t.gen_res_phi  = hardInteraction->resonance->phi();

        t.gen_dil_mass = hardInteraction->dilepton().mass();
        t.gen_dil_pt   = hardInteraction->dilepton().pt();
        t.gen_dil_rap  = hardInteraction->dilepton().Rapidity();
        t.gen_dil_eta  = hardInteraction->dilepton().eta();
        t.gen_dil_phi  = hardInteraction->dilepton().phi();
        t.gen_dil_dR   = deltaR(*hardInteraction->lepMinus, *hardInteraction->lepPlus);
        t.gen_dil_dPhi = deltaPhi(*hardInteraction->lepMinus, *hardInteraction->lepPlus);

        t.gen_lep_p[0]  = hardInteraction->lepMinus->p();
        t.gen_lep_pt[0]  = hardInteraction->lepMinus->pt();
        t.gen_lep_px[0]  = hardInteraction->lepMinus->px();
        t.gen_lep_py[0]  = hardInteraction->lepMinus->py();
        t.gen_lep_pz[0]  = hardInteraction->lepMinus->pz();
        t.gen_lep_E[0]  = hardInteraction->lepMinus->energy();
        t.gen_lep_eta[0] = hardInteraction->lepMinus->eta();
        t.gen_lep_phi[0] = hardInteraction->lepMinus->phi();
        t.gen_lep_qOverPt[0] = hardInteraction->lepMinus->charge() / hardInteraction->lepMinus->pt();

        t.gen_lep_p[1]  = hardInteraction->lepPlus->p();
        t.gen_lep_pt[1]  = hardInteraction->lepPlus->pt();
        t.gen_lep_px[1]  = hardInteraction->lepMinus->px();
        t.gen_lep_py[1]  = hardInteraction->lepMinus->py();
        t.gen_lep_pz[1]  = hardInteraction->lepMinus->pz();
        t.gen_lep_E[1]  = hardInteraction->lepMinus->energy();
        t.gen_lep_eta[1] = hardInteraction->lepPlus->eta();
        t.gen_lep_phi[1] = hardInteraction->lepPlus->phi();
        t.gen_lep_qOverPt[1] = hardInteraction->lepPlus->charge() / hardInteraction->lepPlus->pt();

        /*
        t.gen_lep_noib_pt[0]  = hardInteraction->lepMinusNoIB->pt();
        t.gen_lep_noib_px[0]  = hardInteraction->lepMinusNoIB->px();
        t.gen_lep_noib_py[0]  = hardInteraction->lepMinusNoIB->py();
        t.gen_lep_noib_pz[0]  = hardInteraction->lepMinusNoIB->pz();
        t.gen_lep_noib_e[0]  = hardInteraction->lepMinusNoIB->energy();
        t.gen_lep_noib_eta[0] = hardInteraction->lepMinusNoIB->eta();
        t.gen_lep_noib_phi[0] = hardInteraction->lepMinusNoIB->phi();

        t.gen_lep_noib_pt[1]  = hardInteraction->lepPlusNoIB->pt();
        t.gen_lep_noib_px[1]  = hardInteraction->lepMinusNoIB->px();
        t.gen_lep_noib_py[1]  = hardInteraction->lepMinusNoIB->py();
        t.gen_lep_noib_pz[1]  = hardInteraction->lepMinusNoIB->pz();
        t.gen_lep_noib_e[1]  = hardInteraction->lepMinusNoIB->energy();
        t.gen_lep_noib_eta[1] = hardInteraction->lepPlusNoIB->eta();
        t.gen_lep_noib_phi[1] = hardInteraction->lepPlusNoIB->phi();
        */
    } // end if hardInteraction->IsValid()

  } // end if fill_gen_info

  //
  // Get dilepton collection
  //
  edm::Handle<pat::CompositeCandidateCollection> dils;
  event.getByLabel(dimu_src, dils);

  //
  // Loop over dil candidates in dils
  //
  int chosen = 0;
  BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils) {

    
    // The dils come pre-sorted so that the first in the list is the one to use
    t.dil_chosen = chosen;
    chosen++;
    t.dil_mass = dil.mass();
    t.dil_pt = dil.pt();
    t.dil_rap = dil.rapidity();
    t.dil_eta = dil.eta();
    t.dil_phi = dil.phi();
    t.dil_dR = deltaR(*dil.daughter(0), *dil.daughter(1));
    t.dil_dPhi = deltaPhi(*dil.daughter(0), *dil.daughter(1));

    // Only deal with dileptons composed of e,mu for now.
    assert(dil.numberOfDaughters() == 2);
    assert(abs(dil.daughter(0)->pdgId()) == 11 || abs(dil.daughter(0)->pdgId()) == 13);
    assert(abs(dil.daughter(1)->pdgId()) == 11 || abs(dil.daughter(1)->pdgId()) == 13);

    // set opp_sign and diff_flavor
    const bool opp_sign = dil.daughter(0)->charge() + dil.daughter(1)->charge() == 0;
    const bool diff_flavor = abs(dil.daughter(0)->pdgId()) != abs(dil.daughter(1)->pdgId());
    //const bool dimuon = abs(dil.daughter(0)->pdgId()) == 13 && abs(dil.daughter(1)->pdgId()) == 13;
    
    //
    // Loop over dil.daughters
    //
    for (size_t i = 0; i < 2; ++i) {

      // For e-mu dileptons, put the muon first. For opposite-sign
      // dileptons, always put the negative lepton first. Otherwise
      // don't mess with the order.
      size_t w = i;
      if (diff_flavor)
	w = abs(dil.daughter(i)->pdgId()) == 13 ? 0 : 1;
      else if (opp_sign)
	w = dil.daughter(i)->charge() < 0 ? 0 : 1;

      // Set lepton information
      t.lep_id[w] = dil.daughter(i)->pdgId();
      t.lep_eta[w] = dil.daughter(i)->eta();
      t.lep_phi[w] = dil.daughter(i)->phi();

      //
      // Non-muon information
      //
      if (abs(t.lep_id[w]) != 13) {
	t.lep_pt_err[w] = -999;
	t.lep_tk_p[w] = -999;
	t.lep_tk_pt[w] = -999;
	t.lep_tk_pt_err[w] = -999;
	t.lep_tk_px[w] = -999;
	t.lep_tk_py[w] = -999;
	t.lep_tk_pz[w] = -999;
	t.lep_tk_eta[w] = -999;
	t.lep_tk_phi[w] = -999;
	t.lep_tk_vz[w] = -999;
	t.lep_tk_dz[w] = -999;
	t.lep_tk_chi2[w] = -999;
	t.lep_tk_ndf[w] = -999;
	t.lep_tk_qOverPt[w] = -999;
	t.lep_glb_p[w] = -999;
	t.lep_glb_pt[w] = -999;
	t.lep_glb_pt_err[w] = -999;
	t.lep_glb_px[w] = -999;
	t.lep_glb_py[w] = -999;
	t.lep_glb_pz[w] = -999;
	t.lep_glb_eta[w] = -999;
	t.lep_glb_phi[w] = -999;
	t.lep_glb_chi2[w] = -999;
	t.lep_glb_ndf[w] = -999;
	t.lep_glb_qOverPt[w] = -999;
	t.lep_tpfms_p[w] = -999;
	t.lep_tpfms_pt[w] = -999;
	t.lep_tpfms_pt_err[w] = -999;
	t.lep_tpfms_px[w] = -999;
	t.lep_tpfms_py[w] = -999;
	t.lep_tpfms_pz[w] = -999;
	t.lep_tpfms_eta[w] = -999;
	t.lep_tpfms_phi[w] = -999;
	t.lep_tpfms_chi2[w] = -999;
	t.lep_tpfms_ndf[w] = -999;
	t.lep_tpfms_qOverPt[w] = -999;
	t.lep_picky_p[w] = -999;
	t.lep_picky_pt[w] = -999;
	t.lep_picky_pt_err[w] = -999;
	t.lep_picky_px[w] = -999;
	t.lep_picky_py[w] = -999;
	t.lep_picky_pz[w] = -999;
	t.lep_picky_eta[w] = -999;
	t.lep_picky_phi[w] = -999;
	t.lep_picky_chi2[w] = -999;
	t.lep_picky_ndf[w] = -999;
	t.lep_picky_qOverPt[w] = -999;
	t.lep_cocktail_p[w] = -999;
	t.lep_cocktail_pt[w] = -999;
	t.lep_cocktail_pt_err[w] = -999;
	t.lep_cocktail_px[w] = -999;
	t.lep_cocktail_py[w] = -999;
	t.lep_cocktail_pz[w] = -999;
	t.lep_cocktail_eta[w] = -999;
	t.lep_cocktail_phi[w] = -999;
	t.lep_cocktail_chi2[w] = -999;
	t.lep_cocktail_ndf[w] = -999;
	t.lep_cocktail_qOverPt[w] = -999;
	t.lep_tuneP_p[w] = -999;
	t.lep_tuneP_pt[w] = -999;
	t.lep_tuneP_pt_err[w] = -999;
	t.lep_tuneP_px[w] = -999;
	t.lep_tuneP_py[w] = -999;
	t.lep_tuneP_pz[w] = -999;
	t.lep_tuneP_eta[w] = -999;
	t.lep_tuneP_phi[w] = -999;
	t.lep_tuneP_chi2[w] = -999;
	t.lep_tuneP_ndf[w] = -999;
	t.lep_tuneP_qOverPt[w] = -999;
	t.lep_triggerMatchPt[w] = -999;
	t.lep_triggerMatchEta[w] = -999;
	t.lep_chi2dof[w] = -999;
	t.lep_dB[w] = -999;
	t.lep_sumPt[w] = -999;
	t.lep_emEt[w] = -999;
	t.lep_hadEt[w] = -999;
	t.lep_hoEt[w] = -999;
	t.lep_pfIso[w] = -999;
	t.lep_pfIsoDB[w] = -999;
	t.lep_timeNdof[w] = -999;
	t.lep_timeInOut[w] = -999;
	t.lep_timeOutIn[w] = -999;
	t.lep_timeInOutErr[w] = -999;
	t.lep_timeOutInErr[w] = -999;
	t.lep_tk_numberOfValidTrackerHits[w] = -999; 
	t.lep_tk_numberOfValidTrackerLayers[w] = -999; 
	t.lep_tk_numberOfValidPixelHits[w] = -999;
	t.lep_glb_numberOfValidTrackerHits[w] = -999; 
	t.lep_glb_numberOfValidTrackerLayers[w] = -999; 
	t.lep_glb_numberOfValidPixelHits[w] = -999;
	t.lep_glb_muonStationsWithValidHits[w] = -999;
	t.lep_glb_numberOfValidMuonHits[w] = -999;
	t.lep_glb_numberOfValidMuonDTHits[w] = -999;
	t.lep_glb_numberOfValidMuonCSCHits[w] = -999;
	t.lep_glb_numberOfValidMuonRPCHits[w] = -999;
	t.lep_glb_muonStationsWithValidHits[w] = -999;
	t.lep_glb_dtStationsWithValidHits[w] = -999;
	t.lep_glb_cscStationsWithValidHits[w] = -999;
	t.lep_glb_rpcStationsWithValidHits[w] = -999;
        t.lep_glb_innermostMuonStationWithValidHits[w] = -999;
        t.lep_glb_outermostMuonStationWithValidHits[w] = -999;
	t.lep_numberOfMatches[w] = -999;
	t.lep_numberOfMatchedStations[w] = -999;
        t.lep_stationMask[w] = 999;
	t.lep_isGlobalMuon[w] = false;
	t.lep_isTrackerMuon[w] = false;

        // 
        // Electron Information
        //
	if (abs(t.lep_id[w]) == 11) {
	  const pat::Electron* el = toConcretePtr<pat::Electron>(dileptonDaughter(dil, i));
	  assert(el);

	  t.lep_p[w] = dil.daughter(i)->p();
	  t.lep_pt[w] = dil.daughter(i)->pt();
	  t.lep_px[w] = dil.daughter(i)->px();
	  t.lep_py[w] = dil.daughter(i)->py();
	  t.lep_pz[w] = dil.daughter(i)->pz();
          t.lep_E[w] = dil.daughter(i)->energy();
	  t.lep_heep_id[w] = userInt(*el, "HEEPId", 999);
	  t.lep_min_muon_dR[w] = userFloat(*el, "min_muon_dR", 999);

	} // end if electron 

      } // end if !muon

      //
      // Muon Information
      //
      else { // else of if (abs(t.lep_id[w]) != 13) 
	t.lep_heep_id[w] = 999;
	t.lep_min_muon_dR[w] = 999;

        // 
        // Muon is always from dilepton object
        //
	const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
	assert(mu);

        //
        // Default Muon info (tuneP)
        //
	const reco::Track* tk = patmuon::getPickedTrack(*mu).get();
	assert (tk);
	t.lep_p[w]     = tk->p();
	t.lep_pt[w]     = tk->pt();
	t.lep_px[w]     = tk->px();
	t.lep_py[w]     = tk->py();
	t.lep_pz[w]     = tk->pz();
	t.lep_qOverPt[w] = tk->charge() / tk->pt();
	t.lep_pt_err[w] = ptError(tk);
        
        //
        // Tracker Track Muon Information
        //
	t.lep_tk_p[w] = mu->innerTrack()->p();
	t.lep_tk_pt[w] = mu->innerTrack()->pt();
	t.lep_tk_pt_err[w] = mu->innerTrack()->ptError();
	t.lep_tk_px[w] = mu->innerTrack()->px();
	t.lep_tk_py[w] = mu->innerTrack()->py();
	t.lep_tk_pz[w] = mu->innerTrack()->pz();
	t.lep_tk_eta[w] = mu->innerTrack()->eta();
	t.lep_tk_phi[w] = mu->innerTrack()->phi();
	t.lep_tk_vz[w] = mu->innerTrack()->vz();
	t.lep_tk_dz[w] = mu->innerTrack()->dz();
	t.lep_tk_chi2[w] = mu->innerTrack()->chi2();
	t.lep_tk_ndf[w] = mu->innerTrack()->ndof();
	t.lep_tk_qOverPt[w] = (mu->charge())/(mu->innerTrack()->pt());

        // 
        // Global Muon Information
        //
	t.lep_glb_p[w] = mu->globalTrack()->p();
	t.lep_glb_pt[w] = mu->globalTrack()->pt();
	t.lep_glb_pt_err[w] = mu->globalTrack()->ptError();
	t.lep_glb_px[w] = mu->innerTrack()->px();
	t.lep_glb_py[w] = mu->innerTrack()->py();
	t.lep_glb_pz[w] = mu->innerTrack()->pz();
	t.lep_glb_eta[w] = mu->globalTrack()->eta();
	t.lep_glb_phi[w] = mu->globalTrack()->phi();
	t.lep_glb_chi2[w] = mu->globalTrack()->chi2();
	t.lep_glb_ndf[w] = mu->globalTrack()->ndof();
	t.lep_glb_qOverPt[w] = (mu->charge())/(mu->globalTrack()->pt());

        //
        // Tracker Plus First Muon Station Muon Information
        //
	if (mu->tpfmsMuon().isNull()) {
	  t.lep_tpfms_p[w] = -999;
	  t.lep_tpfms_pt[w] = -999;
	  t.lep_tpfms_pt_err[w] = -999;
	  t.lep_tpfms_px[w] = -999;
	  t.lep_tpfms_py[w] = -999;
	  t.lep_tpfms_pz[w] = -999;
	  t.lep_tpfms_eta[w] = -999;
	  t.lep_tpfms_phi[w] = -999;
	  t.lep_tpfms_chi2[w] = -999;
	  t.lep_tpfms_ndf[w] = -999;
	  t.lep_tpfms_qOverPt[w] = -999;
	}
	else {
	  t.lep_tpfms_p[w] = mu->tpfmsMuon()->p();
	  t.lep_tpfms_pt[w] = mu->tpfmsMuon()->pt();
	  t.lep_tpfms_pt_err[w] = mu->tpfmsMuon()->ptError();
	  t.lep_tpfms_px[w] = mu->tpfmsMuon()->px();
	  t.lep_tpfms_py[w] = mu->tpfmsMuon()->py();
	  t.lep_tpfms_pz[w] = mu->tpfmsMuon()->pz();
	  t.lep_tpfms_eta[w] = mu->tpfmsMuon()->eta();
	  t.lep_tpfms_phi[w] = mu->tpfmsMuon()->phi();
	  t.lep_tpfms_chi2[w] = mu->tpfmsMuon()->chi2();
	  t.lep_tpfms_ndf[w] = mu->tpfmsMuon()->ndof();
	  t.lep_tpfms_qOverPt[w] = (mu->charge())/(mu->tpfmsMuon()->pt());
	}

        //
        // Picky Muon Information
        //
	if (mu->pickyMuon().isNull()) {
	  t.lep_picky_p[w] = -999;
	  t.lep_picky_pt[w] = -999;
	  t.lep_picky_pt_err[w] = -999;
	  t.lep_picky_px[w] = -999;
	  t.lep_picky_py[w] = -999;
	  t.lep_picky_pz[w] = -999;
	  t.lep_picky_eta[w] = -999;
	  t.lep_picky_phi[w] = -999;
	  t.lep_picky_chi2[w] = -999;
	  t.lep_picky_ndf[w] = -999;
	  t.lep_picky_qOverPt[w] = -999;
	}
	else {
	  t.lep_picky_p[w] = mu->pickyMuon()->p();
	  t.lep_picky_pt[w] = mu->pickyMuon()->pt();
	  t.lep_picky_pt_err[w] = mu->pickyMuon()->ptError();
	  t.lep_picky_px[w] = mu->pickyMuon()->px();
	  t.lep_picky_py[w] = mu->pickyMuon()->py();
	  t.lep_picky_pz[w] = mu->pickyMuon()->pz();
	  t.lep_picky_eta[w] = mu->pickyMuon()->eta();
	  t.lep_picky_phi[w] = mu->pickyMuon()->phi();
	  t.lep_picky_chi2[w] = mu->pickyMuon()->chi2();
	  t.lep_picky_ndf[w] = mu->pickyMuon()->ndof();
	  t.lep_picky_qOverPt[w] = (mu->charge())/(mu->pickyMuon()->pt());
	}

        //
        // Cocktail Muon Information
        //
	reco::TrackRef cocktail = muon::tevOptimized(*mu, 200, 17, 40, 0.25).first;
	if (cocktail.isNull()) {
	  t.lep_cocktail_p[w] = -999;
	  t.lep_cocktail_pt[w] = -999;
	  t.lep_cocktail_pt_err[w] = -999;
	  t.lep_cocktail_px[w] = -999;
	  t.lep_cocktail_py[w] = -999;
	  t.lep_cocktail_pz[w] = -999;
	  t.lep_cocktail_eta[w] = -999;
	  t.lep_cocktail_phi[w] = -999;
	  t.lep_cocktail_chi2[w] = -999;
	  t.lep_cocktail_ndf[w] = -999;
	  t.lep_cocktail_qOverPt[w] = -999;
	  t.lep_cocktail_choice[w] = -999;
	}
	else {
	  t.lep_cocktail_p[w] = cocktail->p();
	  t.lep_cocktail_pt[w] = cocktail->pt();
	  t.lep_cocktail_pt_err[w] = cocktail->ptError();
	  t.lep_cocktail_px[w] = cocktail->px();
	  t.lep_cocktail_py[w] = cocktail->py();
	  t.lep_cocktail_pz[w] = cocktail->pz();
	  t.lep_cocktail_eta[w] = cocktail->eta();
	  t.lep_cocktail_phi[w] = cocktail->phi();
	  t.lep_cocktail_chi2[w] = cocktail->chi2();
	  t.lep_cocktail_ndf[w] = cocktail->ndof();
	  t.lep_cocktail_qOverPt[w] = (mu->charge())/(cocktail->pt());
	  t.lep_cocktail_choice[w] = short(patmuon::whichTrack(*mu, cocktail));
        }
	if (mu->tunePMuonBestTrack().isNull()) {
	  t.lep_tuneP_p[w] = -999;
	  t.lep_tuneP_pt[w] = -999;
	  t.lep_tuneP_pt_err[w] = -999;
	  t.lep_tuneP_px[w] = -999;
	  t.lep_tuneP_py[w] = -999;
	  t.lep_tuneP_pz[w] = -999;
	  t.lep_tuneP_eta[w] = -999;
	  t.lep_tuneP_phi[w] = -999;
	  t.lep_tuneP_chi2[w] = -999;
	  t.lep_tuneP_ndf[w] = -999;
	  t.lep_tuneP_qOverPt[w] = -999;
	}
	else {
	  t.lep_tuneP_p[w] = mu->tunePMuonBestTrack()->p();
	  t.lep_tuneP_pt[w] = mu->tunePMuonBestTrack()->pt();
	  t.lep_tuneP_pt_err[w] = mu->tunePMuonBestTrack()->ptError();
	  t.lep_tuneP_px[w] = mu->tunePMuonBestTrack()->px();
	  t.lep_tuneP_py[w] = mu->tunePMuonBestTrack()->py();
	  t.lep_tuneP_pz[w] = mu->tunePMuonBestTrack()->pz();
	  t.lep_tuneP_eta[w] = mu->tunePMuonBestTrack()->eta();
	  t.lep_tuneP_phi[w] = mu->tunePMuonBestTrack()->phi();
	  t.lep_tuneP_chi2[w] = mu->tunePMuonBestTrack()->chi2();
	  t.lep_tuneP_ndf[w] = mu->tunePMuonBestTrack()->ndof();
	  t.lep_tuneP_qOverPt[w] = (mu->charge())/(mu->tunePMuonBestTrack()->pt());
	}

        //
        // TuneP Track information
        //
	if (mu->tunePMuonBestTrack().isNull()) {
	  t.lep_tuneP_p[w] = -999;
	  t.lep_tuneP_pt[w] = -999;
	  t.lep_tuneP_pt_err[w] = -999;
	  t.lep_tuneP_px[w] = -999;
	  t.lep_tuneP_py[w] = -999;
	  t.lep_tuneP_pz[w] = -999;
	  t.lep_tuneP_eta[w] = -999;
	  t.lep_tuneP_phi[w] = -999;
	  t.lep_tuneP_chi2[w] = -999;
	  t.lep_tuneP_ndf[w] = -999;
	  t.lep_tuneP_qOverPt[w] = -999;
	}
	else {
	  t.lep_tuneP_p[w] = mu->tunePMuonBestTrack()->p();
	  t.lep_tuneP_pt[w] = mu->tunePMuonBestTrack()->pt();
	  t.lep_tuneP_pt_err[w] = mu->tunePMuonBestTrack()->ptError();
	  t.lep_tuneP_px[w] = mu->tunePMuonBestTrack()->px();
	  t.lep_tuneP_py[w] = mu->tunePMuonBestTrack()->py();
	  t.lep_tuneP_pz[w] = mu->tunePMuonBestTrack()->pz();
	  t.lep_tuneP_eta[w] = mu->tunePMuonBestTrack()->eta();
	  t.lep_tuneP_phi[w] = mu->tunePMuonBestTrack()->phi();
	  t.lep_tuneP_chi2[w] = mu->tunePMuonBestTrack()->chi2();
	  t.lep_tuneP_ndf[w] = mu->tunePMuonBestTrack()->ndof();
	  t.lep_tuneP_qOverPt[w] = (mu->charge())/(mu->tunePMuonBestTrack()->pt());
	}

        //
        // Trigger Match Information
        //
	t.lep_triggerMatchPt[w]  = userFloat(*mu, "TriggerMatchPt",  -999);
	t.lep_triggerMatchEta[w] = userFloat(*mu, "TriggerMatchEta", -999);

        //
        // Misc. event quantities
        //
	t.lep_chi2dof[w] = mu->globalTrack()->normalizedChi2();
	t.lep_dB[w] = mu->dB();
	t.lep_sumPt[w] = mu->isolationR03().sumPt;
	t.lep_emEt[w] = mu->isolationR03().emEt;
	t.lep_hadEt[w] = mu->isolationR03().hadEt;
	t.lep_hoEt[w] = mu->isolationR03().hoEt;
        t.lep_pfIso[w] = mu->pfIsolationR04().sumChargedHadronPt + mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt;
        t.lep_pfIsoDB[w] = mu->pfIsolationR04().sumChargedHadronPt + std::max(mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt,0.0);
	if (mu->isTimeValid()) {
	  t.lep_timeNdof[w] = mu->time().nDof;
	  t.lep_timeInOut[w] = mu->time().timeAtIpInOut;
	  t.lep_timeOutIn[w] = mu->time().timeAtIpOutIn;
	  t.lep_timeInOutErr[w] = mu->time().timeAtIpInOutErr;
	  t.lep_timeOutInErr[w] = mu->time().timeAtIpOutInErr;
	}
	else {
	  t.lep_timeNdof[w] = -999;
	  t.lep_timeInOut[w] = -999;
	  t.lep_timeOutIn[w] = -999;
	  t.lep_timeInOutErr[w] = -999;
	  t.lep_timeOutInErr[w] = -999;
	}

        //
        // Tracker track
        // 
	t.lep_isTrackerMuon[w] = mu->isTrackerMuon();
        // Tracker Track tracker information
	t.lep_tk_numberOfValidTrackerHits[w] = mu->innerTrack()->hitPattern().numberOfValidTrackerHits();
	t.lep_tk_numberOfValidTrackerLayers[w] = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
	t.lep_tk_numberOfValidPixelHits[w] = mu->innerTrack()->hitPattern().numberOfValidPixelHits();

        // 
        // Global track
        //
	t.lep_isGlobalMuon[w] = mu->isGlobalMuon();
        // Global track tracker information
	t.lep_glb_numberOfValidTrackerHits[w] = mu->globalTrack()->hitPattern().numberOfValidTrackerHits();
	t.lep_glb_numberOfValidTrackerLayers[w] = mu->globalTrack()->hitPattern().trackerLayersWithMeasurement();
	t.lep_glb_numberOfValidPixelHits[w] = mu->globalTrack()->hitPattern().numberOfValidPixelHits();
        // Valid Muon (all), DT, CSC, RPC hits
	t.lep_glb_numberOfValidMuonHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
	t.lep_glb_numberOfValidMuonDTHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonDTHits();
	t.lep_glb_numberOfValidMuonCSCHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonCSCHits();
	t.lep_glb_numberOfValidMuonRPCHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonRPCHits();
        // Valid Muon, DT, CSC, RPC, innermost, outermost Station Hits
	t.lep_glb_muonStationsWithValidHits[w] = mu->globalTrack()->hitPattern().muonStationsWithValidHits();
	t.lep_glb_dtStationsWithValidHits[w] = mu->globalTrack()->hitPattern().dtStationsWithValidHits();
	t.lep_glb_cscStationsWithValidHits[w] = mu->globalTrack()->hitPattern().cscStationsWithValidHits();
	t.lep_glb_rpcStationsWithValidHits[w] = mu->globalTrack()->hitPattern().rpcStationsWithValidHits();
        t.lep_glb_innermostMuonStationWithValidHits[w] = mu->globalTrack()->hitPattern().innermostMuonStationWithValidHits();
        t.lep_glb_outermostMuonStationWithValidHits[w] = mu->globalTrack()->hitPattern().outermostMuonStationWithValidHits();
        // number of chambers with matched segments
	t.lep_numberOfMatches[w] = mu->numberOfMatches();
        // number of stations with matched segments
	t.lep_numberOfMatchedStations[w] = mu->numberOfMatchedStations();
        // get bit map of stations with matched segments
        // bits 0-1-2-3 = DT stations 1-2-3-4
        // bits 4-5-6-7 = CSC stations 1-2-3-4
        t.lep_stationMask[w] = mu->stationMask();
        // number of chambers
        t.lep_numberOfChambers[w] = mu->numberOfChambers();
        // number of chambers not including RPC matches
        t.lep_numberOfChambersNoRPC[w] = mu->numberOfChambersNoRPC();
        // distanceCut = 10cm by default (distance in cm)
        t.lep_stationGapMaskDistance[w] = mu->stationGapMaskDistance();
        // sigmaCut = 3 by default (in # sigmas)
        t.lep_stationGapMaskPull[w] = mu->stationGapMaskPull();

      } // end else of if (abs(t.lep_id[w]) != 13) 

    } // end for (int i = 0; i < 2; i++) // Loop over dilepton leptons

    // 
    // more event quantites
    //
    t.cos_angle    = userFloat(dil, "cos_angle", 999);
    t.vertex_chi2  = userFloat(dil, "vertex_chi2");
    t.vertex_m     = userFloat(dil, "vertexM");
    t.vertex_m_err = userFloat(dil, "vertexMError");
    t.vertex_x     = userFloat(dil, "vertexX");
    t.vertex_x_err = userFloat(dil, "vertexXError");
    t.vertex_y     = userFloat(dil, "vertexY");
    t.vertex_y_err = userFloat(dil, "vertexYError");
    t.vertex_z     = userFloat(dil, "vertexZ");
    t.vertex_z_err = userFloat(dil, "vertexZError");

    if (opp_sign) {
      const reco::CandidateBaseRef mum = dileptonDaughterByCharge(dil, -1);
      const reco::CandidateBaseRef mup = dileptonDaughterByCharge(dil, +1);

      t.cos_cs = calcCosThetaCSAnal(mum->pz(), mum->energy(), mup->pz(), mup->energy(), dil.pt(), dil.pz(), dil.mass());
      t.phi_cs = calcPhiCSAnal(mum->px(), mum->py(), mup->px(), mup->py(), dil.pt(), dil.eta(), dil.phi(), dil.mass(), true);
    } // end if opp_sign
    else {
      t.cos_cs = -999;
      t.phi_cs = -999;
    } // end if !opp_sign


  edm::Handle< std::vector< pat::MET > > mets;
  event.getByLabel(met_src, mets);
  t.met_pt = mets->front().pt();
  t.met_phi = mets->front().phi();

  edm::Handle< std::vector< pat::Jet > > jets;
  event.getByLabel(jet_src, jets);


  int nJets = 0;

  for (std::vector<pat::Jet>::const_iterator itJet = jets->begin(); itJet != jets->end(); itJet++) {
        if (fabs(itJet->eta()) < 2.4 && itJet->pt() > 30 && itJet->neutralHadronEnergyFraction() < 0.99 && itJet->neutralEmEnergyFraction() < 0.99 && itJet->chargedHadronEnergyFraction() > 0 && itJet->muonEnergyFraction() < 0.8 && itJet->chargedEmEnergyFraction() < 0.99 && (itJet->chargedMultiplicity()+itJet->neutralMultiplicity()) > 1  && itJet->chargedMultiplicity() > 0 && deltaR((*itJet),dil.daughter(0)->p4()) > 0.4 && deltaR((*itJet),dil.daughter(1)->p4())){
                if (nJets < 4){
                        t.jet_pt[nJets] = itJet->pt();
                        t.jet_eta[nJets] = itJet->eta();
                        t.jet_phi[nJets] = itJet->phi();
                }
                nJets++;
        }
  }
    t.nJets = nJets;
    if (nJets < 4){
    	if (nJets < 3){
    		if (nJets < 2){
			if (nJets < 1){
				t.jet_pt[0] = -999.;
				t.jet_eta[0] = -999.;
				t.jet_phi[0] = -999.;
			}
                        t.jet_pt[1] = -999.;
                        t.jet_eta[1] = -999.;
			t.jet_phi[1] = -999;
		}
		t.jet_pt[2] = -999.;
                t.jet_eta[2] = -999.;
                t.jet_phi[2] = -999.;
	  }
	  t.jet_pt[3] = -999.;                                                                                                                                                                                   t.jet_eta[3] = -999.;
          t.jet_phi[3] = -999.;
    }




    //
    // Fill tree
    //
    tree->Fill();

  } // end BOOST_FOREACH(dil, *dils)
    
  // 
  // Put additional event info here
  // MET, Jets, etc.
  //

} // end SimpleNtupler::analyze

DEFINE_FWK_MODULE(SimpleNtupler);
