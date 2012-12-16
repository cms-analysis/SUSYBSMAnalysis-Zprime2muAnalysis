#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
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
    float beamspot_x;
    float beamspot_x_err;
    float beamspot_y;
    float beamspot_y_err;
    float beamspot_z;
    float beamspot_z_err;
    int nvertices;
    float dil_mass;
    float dil_pt;
    float dil_rap;
    float dil_eta;
    float dil_phi;
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
    float lep_pt[2];
    float lep_pt_err[2];
    float lep_eta[2];
    float lep_phi[2];
    float lep_tk_pt[2];
    float lep_tk_pt_err[2];
    float lep_tk_eta[2];
    float lep_tk_phi[2];
    float lep_tk_dz[2];
    float lep_tk_vz[2];
    float lep_tk_chi2[2];
    float lep_tk_ndf[2];
    float lep_glb_pt[2];
    float lep_glb_pt_err[2];
    float lep_glb_eta[2];
    float lep_glb_phi[2];
    float lep_glb_chi2[2];
    float lep_glb_ndf[2];
    float lep_tpfms_pt[2];
    float lep_tpfms_pt_err[2];
    float lep_tpfms_eta[2];
    float lep_tpfms_phi[2];
    float lep_tpfms_chi2[2];
    float lep_tpfms_ndf[2];
    float lep_picky_pt[2];
    float lep_picky_pt_err[2];
    float lep_picky_eta[2];
    float lep_picky_phi[2];
    float lep_picky_chi2[2];
    float lep_picky_ndf[2];
    float lep_cocktail_pt[2];
    float lep_cocktail_pt_err[2];
    float lep_cocktail_eta[2];
    float lep_cocktail_phi[2];
    float lep_cocktail_chi2[2];
    float lep_cocktail_ndf[2];
    short lep_cocktail_choice[2];
    float lep_triggerMatchPt[2];
    float lep_triggerMatchEta[2];
    float lep_chi2dof[2];
    float lep_dB[2];
    float lep_sumPt[2];
    float lep_emEt[2];
    float lep_hadEt[2];
    float lep_hoEt[2];
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
    short lep_glb_muonStationsWithValidHits[2];
    short lep_numberOfMatches[2];
    short lep_numberOfMatchedStations[2];
    bool lep_isGlobalMuon[2];
    bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool HLTPhysicsDeclared;
    bool GoodVtx;
    bool NoScraping;
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
    float gen_lep_pt[2];
    float gen_lep_eta[2];
    float gen_lep_phi[2];
    float gen_lep_noib_pt[2];
    float gen_lep_noib_eta[2];
    float gen_lep_noib_phi[2];
  };

  tree_t t;
  TTree* tree;

  const edm::InputTag dimu_src;
  const edm::InputTag beamspot_src;
  const edm::InputTag vertices_src;
  const bool fill_gen_info;
  HardInteraction* hardInteraction;
};

TString replace_all(const TString& a, const TString& b, const TString& c) {
  TString ret = a;
  ret.ReplaceAll(b, c);
  return ret;
}

SimpleNtupler::SimpleNtupler(const edm::ParameterSet& cfg)
  : dimu_src(cfg.getParameter<edm::InputTag>("dimu_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertices_src(cfg.getParameter<edm::InputTag>("vertices_src")),
    fill_gen_info(cfg.existsAs<edm::ParameterSet>("hardInteraction")),
    hardInteraction(fill_gen_info ? new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")) : 0)
{
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
  tree->Branch("dil_mass", &t.dil_mass, "dil_mass/F");
  tree->Branch("dil_pt", &t.dil_pt, "dil_pt/F");
  tree->Branch("dil_rap", &t.dil_rap, "dil_rap/F");
  tree->Branch("dil_eta", &t.dil_eta, "dil_eta/F");
  tree->Branch("dil_phi", &t.dil_phi, "dil_phi/F");
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
  tree->Branch("lep_pt", t.lep_pt, "lep_pt[2]/F");
  tree->Branch("lep_pt_err", t.lep_pt_err, "lep_pt_err[2]/F");
  tree->Branch("lep_eta", t.lep_eta, "lep_eta[2]/F");
  tree->Branch("lep_phi", t.lep_phi, "lep_phi[2]/F");
  tree->Branch("lep_tk_pt", t.lep_tk_pt, "lep_tk_pt[2]/F");
  tree->Branch("lep_tk_pt_err", t.lep_tk_pt_err, "lep_tk_pt_err[2]/F");
  tree->Branch("lep_tk_eta", t.lep_tk_eta, "lep_tk_eta[2]/F");
  tree->Branch("lep_tk_phi", t.lep_tk_phi, "lep_tk_phi[2]/F");
  tree->Branch("lep_tk_dz", t.lep_tk_dz, "lep_tk_dz[2]/F");
  tree->Branch("lep_tk_vz", t.lep_tk_vz, "lep_tk_vz[2]/F");
  tree->Branch("lep_tk_chi2", t.lep_tk_chi2, "lep_tk_chi2[2]/F");
  tree->Branch("lep_tk_ndf", t.lep_tk_ndf, "lep_tk_ndf[2]/F");
  tree->Branch("lep_glb_pt", t.lep_glb_pt, "lep_glb_pt[2]/F");
  tree->Branch("lep_glb_pt_err", t.lep_glb_pt_err, "lep_glb_pt_err[2]/F");
  tree->Branch("lep_glb_eta", t.lep_glb_eta, "lep_glb_eta[2]/F");
  tree->Branch("lep_glb_phi", t.lep_glb_phi, "lep_glb_phi[2]/F");
  tree->Branch("lep_glb_chi2", t.lep_glb_chi2, "lep_glb_chi2[2]/F");
  tree->Branch("lep_glb_ndf", t.lep_glb_ndf, "lep_glb_ndf[2]/F");
  tree->Branch("lep_tpfms_pt", t.lep_tpfms_pt, "lep_tpfms_pt[2]/F");
  tree->Branch("lep_tpfms_pt_err", t.lep_tpfms_pt_err, "lep_tpfms_pt_err[2]/F");
  tree->Branch("lep_tpfms_eta", t.lep_tpfms_eta, "lep_tpfms_eta[2]/F");
  tree->Branch("lep_tpfms_phi", t.lep_tpfms_phi, "lep_tpfms_phi[2]/F");
  tree->Branch("lep_tpfms_chi2", t.lep_tpfms_chi2, "lep_tpfms_chi2[2]/F");
  tree->Branch("lep_tpfms_ndf", t.lep_tpfms_ndf, "lep_tpfms_ndf[2]/F");
  tree->Branch("lep_picky_pt", t.lep_picky_pt, "lep_picky_pt[2]/F");
  tree->Branch("lep_picky_pt_err", t.lep_picky_pt_err, "lep_picky_pt_err[2]/F");
  tree->Branch("lep_picky_eta", t.lep_picky_eta, "lep_picky_eta[2]/F");
  tree->Branch("lep_picky_phi", t.lep_picky_phi, "lep_picky_phi[2]/F");
  tree->Branch("lep_picky_chi2", t.lep_picky_chi2, "lep_picky_chi2[2]/F");
  tree->Branch("lep_picky_ndf", t.lep_picky_ndf, "lep_picky_ndf[2]/F");
  tree->Branch("lep_cocktail_pt", t.lep_cocktail_pt, "lep_cocktail_pt[2]/F");
  tree->Branch("lep_cocktail_pt_err", t.lep_cocktail_pt_err, "lep_cocktail_pt_err[2]/F");
  tree->Branch("lep_cocktail_eta", t.lep_cocktail_eta, "lep_cocktail_eta[2]/F");
  tree->Branch("lep_cocktail_phi", t.lep_cocktail_phi, "lep_cocktail_phi[2]/F");
  tree->Branch("lep_cocktail_chi2", t.lep_cocktail_chi2, "lep_cocktail_chi2[2]/F");
  tree->Branch("lep_cocktail_ndf", t.lep_cocktail_ndf, "lep_cocktail_ndf[2]/F");
  tree->Branch("lep_cocktail_choice", t.lep_cocktail_choice, "lep_cocktail_choice[2]/S");
  tree->Branch("lep_triggerMatchPt", t.lep_triggerMatchPt, "lep_triggerMatchPt[2]/F");
  tree->Branch("lep_triggerMatchEta", t.lep_triggerMatchEta, "lep_triggerMatchEta[2]/F");
  tree->Branch("lep_chi2dof", t.lep_chi2dof, "lep_chi2dof[2]/F");
  tree->Branch("lep_dB", t.lep_dB, "lep_dB[2]/F");
  tree->Branch("lep_sumPt", t.lep_sumPt, "lep_sumPt[2]/F");
  tree->Branch("lep_emEt", t.lep_emEt, "lep_emEt[2]/F");
  tree->Branch("lep_hadEt", t.lep_hadEt, "lep_hadEt[2]/F");
  tree->Branch("lep_hoEt", t.lep_hoEt, "lep_hoEt[2]/F");
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
  tree->Branch("lep_glb_muonStationsWithValidHits", t.lep_glb_muonStationsWithValidHits, "lep_glb_muonStationsWithValidHits[2]/S");
  tree->Branch("lep_numberOfMatches", t.lep_numberOfMatches, "lep_numberOfMatches[2]/S");
  tree->Branch("lep_numberOfMatchedStations", t.lep_numberOfMatchedStations, "lep_numberOfMatchedStations[2]/S");
  tree->Branch("lep_isGlobalMuon", t.lep_isGlobalMuon, "lep_isGlobalMuon[2]/O");
  tree->Branch("lep_isTrackerMuon", t.lep_isTrackerMuon, "lep_isTrackerMuon[2]/O");
  tree->Branch("GoodDataRan", &t.GoodDataRan, "GoodDataRan/O");
  tree->Branch("HLTPhysicsDeclared", &t.HLTPhysicsDeclared, "HLTPhysicsDeclared/O");
  tree->Branch("GoodVtx", &t.GoodVtx, "GoodVtx/O");
  tree->Branch("NoScraping", &t.NoScraping, "NoScraping/O");
  if (fill_gen_info) {
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
    tree->Branch("gen_lep_pt", t.gen_lep_pt, "gen_lep_pt[2]/F");
    tree->Branch("gen_lep_eta", t.gen_lep_eta, "gen_lep_eta[2]/F");
    tree->Branch("gen_lep_phi", t.gen_lep_phi, "gen_lep_phi[2]/F");
    tree->Branch("gen_lep_noib_pt", t.gen_lep_noib_pt, "gen_lep_noib_pt[2]/F");
    tree->Branch("gen_lep_noib_eta", t.gen_lep_noib_eta, "gen_lep_noib_eta[2]/F");
    tree->Branch("gen_lep_noib_phi", t.gen_lep_noib_phi, "gen_lep_noib_phi[2]/F");
  }

  tree->SetAlias("OppSign",  "lep_id[0]*lep_id[1] < 0");
  tree->SetAlias("SameSign", "lep_id[0]*lep_id[1] > 0");
  tree->SetAlias("Dimu",     "abs(lep_id[0]*lep_id[1]) == 169");
  tree->SetAlias("Emu",      "abs(lep_id[0]*lep_id[1]) == 143");

#define offlineMinPt "45"
#define triggerMatchMinPt "40"
#define triggerMatchMaxEta "2.1"

  tree->SetAlias("trigger_match_0", "lep_triggerMatchPt[0] > " triggerMatchMinPt " && abs(lep_triggerMatchEta[0]) < " triggerMatchMaxEta);
  tree->SetAlias("trigger_match_1", "lep_triggerMatchPt[1] > " triggerMatchMinPt " && abs(lep_triggerMatchEta[1]) < " triggerMatchMaxEta);
  tree->SetAlias("triggerMatched", "trigger_match_0 || trigger_match_1");

  tree->SetAlias("GoodData", "GoodDataRan && HLTPhysicsDeclared && NoScraping && GoodVtx");

  tree->SetAlias("extraDimuonCuts", "cos_angle > -0.9998 && vertex_chi2 < 10");

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

  t.run = event.id().run();
  t.lumi = event.luminosityBlock();
  t.event = event.id().event();

  edm::Handle<edm::TriggerResults> respat;
  event.getByLabel(edm::InputTag("TriggerResults", "", "PAT"), respat);
  const edm::TriggerNames& namespat = event.triggerNames(*respat);
  if (namespat.triggerIndex("goodDataHLTPhysicsDeclared") < respat->size()) {
    t.GoodDataRan = 1;
    t.HLTPhysicsDeclared = respat->accept(namespat.triggerIndex("goodDataHLTPhysicsDeclared"));
    t.GoodVtx = respat->accept(namespat.triggerIndex("goodDataPrimaryVertexFilter"));
    t.NoScraping = respat->accept(namespat.triggerIndex("goodDataNoScraping"));
  }

  edm::Handle<reco::BeamSpot> bs;
  event.getByLabel(beamspot_src, bs);
  t.beamspot_x     = bs->x0();
  t.beamspot_x_err = bs->x0Error();
  t.beamspot_y     = bs->y0();
  t.beamspot_y_err = bs->y0Error();
  t.beamspot_z     = bs->z0();
  t.beamspot_z_err = bs->z0Error();

  edm::Handle<reco::VertexCollection> pvs;
  event.getByLabel(vertices_src, pvs);
  t.nvertices = 0;
  BOOST_FOREACH(const reco::Vertex& vtx, *pvs)
    if (vtx.ndof() > 4 && fabs(vtx.z()) <= 24 && fabs(vtx.position().rho()) <= 2)
      t.nvertices += 1;
  
  if (fill_gen_info) {
    // This only works for DY/Z'/RSG events, and really just for PYTHIA!
    hardInteraction->Fill(event);

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

    t.gen_lep_pt[0]  = hardInteraction->lepMinus->pt();
    t.gen_lep_eta[0] = hardInteraction->lepMinus->eta();
    t.gen_lep_phi[0] = hardInteraction->lepMinus->phi();

    t.gen_lep_pt[1]  = hardInteraction->lepPlus->pt();
    t.gen_lep_eta[1] = hardInteraction->lepPlus->eta();
    t.gen_lep_phi[1] = hardInteraction->lepPlus->phi();

    t.gen_lep_noib_pt[0]  = hardInteraction->lepMinusNoIB->pt();
    t.gen_lep_noib_eta[0] = hardInteraction->lepMinusNoIB->eta();
    t.gen_lep_noib_phi[0] = hardInteraction->lepMinusNoIB->phi();

    t.gen_lep_noib_pt[1]  = hardInteraction->lepPlusNoIB->pt();
    t.gen_lep_noib_eta[1] = hardInteraction->lepPlusNoIB->eta();
    t.gen_lep_noib_phi[1] = hardInteraction->lepPlusNoIB->phi();
  }

  edm::Handle<pat::CompositeCandidateCollection> dils;
  event.getByLabel(dimu_src, dils);

  BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils) {
    t.dil_mass = dil.mass();
    t.dil_pt = dil.pt();
    t.dil_rap = dil.rapidity();
    t.dil_eta = dil.eta();
    t.dil_phi = dil.phi();

    // Only deal with dileptons composed of e,mu for now.
    assert(dil.numberOfDaughters() == 2);
    assert(abs(dil.daughter(0)->pdgId()) == 11 || abs(dil.daughter(0)->pdgId()) == 13);
    assert(abs(dil.daughter(1)->pdgId()) == 11 || abs(dil.daughter(1)->pdgId()) == 13);

    const bool opp_sign = dil.daughter(0)->charge() + dil.daughter(1)->charge() == 0;
    const bool diff_flavor = abs(dil.daughter(0)->pdgId()) != abs(dil.daughter(1)->pdgId());
    //const bool dimuon = abs(dil.daughter(0)->pdgId()) == 13 && abs(dil.daughter(1)->pdgId()) == 13;

    for (size_t i = 0; i < 2; ++i) {
      // For e-mu dileptons, put the muon first. For opposite-sign
      // dileptons, always put the negative lepton first. Otherwise
      // don't mess with the order.
      size_t w = i;
      if (diff_flavor)
	w = abs(dil.daughter(i)->pdgId()) == 13 ? 0 : 1;
      else if (opp_sign)
	w = dil.daughter(i)->charge() < 0 ? 0 : 1;

      t.lep_id[w] = dil.daughter(i)->pdgId();
      t.lep_eta[w] = dil.daughter(i)->eta();
      t.lep_phi[w] = dil.daughter(i)->phi();

      if (abs(t.lep_id[w]) != 13) {
	t.lep_pt_err[w] = -999;
	t.lep_tk_pt[w] = -999;
	t.lep_tk_pt_err[w] = -999;
	t.lep_tk_eta[w] = -999;
	t.lep_tk_phi[w] = -999;
	t.lep_tk_vz[w] = -999;
	t.lep_tk_dz[w] = -999;
	t.lep_tk_chi2[w] = -999;
	t.lep_tk_ndf[w] = -999;
	t.lep_glb_pt[w] = -999;
	t.lep_glb_pt_err[w] = -999;
	t.lep_glb_eta[w] = -999;
	t.lep_glb_phi[w] = -999;
	t.lep_glb_chi2[w] = -999;
	t.lep_glb_ndf[w] = -999;
	t.lep_tpfms_pt[w] = -999;
	t.lep_tpfms_pt_err[w] = -999;
	t.lep_tpfms_eta[w] = -999;
	t.lep_tpfms_phi[w] = -999;
	t.lep_tpfms_chi2[w] = -999;
	t.lep_tpfms_ndf[w] = -999;
	t.lep_picky_pt[w] = -999;
	t.lep_picky_pt_err[w] = -999;
	t.lep_picky_eta[w] = -999;
	t.lep_picky_phi[w] = -999;
	t.lep_picky_chi2[w] = -999;
	t.lep_picky_ndf[w] = -999;
	t.lep_cocktail_pt[w] = -999;
	t.lep_cocktail_pt_err[w] = -999;
	t.lep_cocktail_eta[w] = -999;
	t.lep_cocktail_phi[w] = -999;
	t.lep_cocktail_chi2[w] = -999;
	t.lep_cocktail_ndf[w] = -999;
	t.lep_triggerMatchPt[w] = -999;
	t.lep_triggerMatchEta[w] = -999;
	t.lep_chi2dof[w] = -999;
	t.lep_dB[w] = -999;
	t.lep_sumPt[w] = -999;
	t.lep_emEt[w] = -999;
	t.lep_hadEt[w] = -999;
	t.lep_hoEt[w] = -999;
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
	t.lep_glb_numberOfValidMuonHits[w] = -999;
	t.lep_glb_muonStationsWithValidHits[w] = -999;
	t.lep_numberOfMatches[w] = -999;
	t.lep_numberOfMatchedStations[w] = -999;
	t.lep_isGlobalMuon[w] = false;
	t.lep_isTrackerMuon[w] = false;

	if (abs(t.lep_id[w]) == 11) {
	  const pat::Electron* el = toConcretePtr<pat::Electron>(dileptonDaughter(dil, i));
	  assert(el);

	  t.lep_pt[w] = dil.daughter(i)->pt();
	  t.lep_heep_id[w] = userInt(*el, "HEEPId", 999);
	  t.lep_min_muon_dR[w] = userFloat(*el, "min_muon_dR", 999);
	}
      }
      else {
	t.lep_heep_id[w] = 999;
	t.lep_min_muon_dR[w] = 999;

	const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
	assert(mu);

	const reco::Track* tk = patmuon::getPickedTrack(*mu).get();
	assert (tk);

	t.lep_pt[w]     = tk->pt();
	t.lep_pt_err[w] = ptError(tk);
	t.lep_tk_pt[w] = mu->innerTrack()->pt();
	t.lep_tk_pt_err[w] = mu->innerTrack()->ptError();
	t.lep_tk_eta[w] = mu->innerTrack()->eta();
	t.lep_tk_phi[w] = mu->innerTrack()->phi();
	t.lep_tk_vz[w] = mu->innerTrack()->vz();
	t.lep_tk_dz[w] = mu->innerTrack()->dz();
	t.lep_tk_chi2[w] = mu->innerTrack()->chi2();
	t.lep_tk_ndf[w] = mu->innerTrack()->ndof();
	t.lep_glb_pt[w] = mu->globalTrack()->pt();
	t.lep_glb_pt_err[w] = mu->globalTrack()->ptError();
	t.lep_glb_eta[w] = mu->globalTrack()->eta();
	t.lep_glb_phi[w] = mu->globalTrack()->phi();
	t.lep_glb_chi2[w] = mu->globalTrack()->chi2();
	t.lep_glb_ndf[w] = mu->globalTrack()->ndof();
	if (mu->tpfmsMuon().isNull()) {
	  t.lep_tpfms_pt[w] = -999;
	  t.lep_tpfms_pt_err[w] = -999;
	  t.lep_tpfms_eta[w] = -999;
	  t.lep_tpfms_phi[w] = -999;
	  t.lep_tpfms_chi2[w] = -999;
	  t.lep_tpfms_ndf[w] = -999;
	}
	else {
	  t.lep_tpfms_pt[w] = mu->tpfmsMuon()->pt();
	  t.lep_tpfms_pt_err[w] = mu->tpfmsMuon()->ptError();
	  t.lep_tpfms_eta[w] = mu->tpfmsMuon()->eta();
	  t.lep_tpfms_phi[w] = mu->tpfmsMuon()->phi();
	  t.lep_tpfms_chi2[w] = mu->tpfmsMuon()->chi2();
	  t.lep_tpfms_ndf[w] = mu->tpfmsMuon()->ndof();
	}
	if (mu->pickyMuon().isNull()) {
	  t.lep_picky_pt[w] = -999;
	  t.lep_picky_pt_err[w] = -999;
	  t.lep_picky_eta[w] = -999;
	  t.lep_picky_phi[w] = -999;
	  t.lep_picky_chi2[w] = -999;
	  t.lep_picky_ndf[w] = -999;
	}
	else {
	  t.lep_picky_pt[w] = mu->pickyMuon()->pt();
	  t.lep_picky_pt_err[w] = mu->pickyMuon()->ptError();
	  t.lep_picky_eta[w] = mu->pickyMuon()->eta();
	  t.lep_picky_phi[w] = mu->pickyMuon()->phi();
	  t.lep_picky_chi2[w] = mu->pickyMuon()->chi2();
	  t.lep_picky_ndf[w] = mu->pickyMuon()->ndof();
	}

	reco::TrackRef cocktail = muon::tevOptimized(*mu, 200, 17, 40, 0.25).first;
	if (cocktail.isNull()) {
	  t.lep_cocktail_pt[w] = -999;
	  t.lep_cocktail_pt_err[w] = -999;
	  t.lep_cocktail_eta[w] = -999;
	  t.lep_cocktail_phi[w] = -999;
	  t.lep_cocktail_chi2[w] = -999;
	  t.lep_cocktail_ndf[w] = -999;
	  t.lep_cocktail_choice[w] = -999;
	}
	else {
	  t.lep_cocktail_pt[w] = cocktail->pt();
	  t.lep_cocktail_pt_err[w] = cocktail->ptError();
	  t.lep_cocktail_eta[w] = cocktail->eta();
	  t.lep_cocktail_phi[w] = cocktail->phi();
	  t.lep_cocktail_chi2[w] = cocktail->chi2();
	  t.lep_cocktail_ndf[w] = cocktail->ndof();
	  t.lep_cocktail_choice[w] = short(patmuon::whichTrack(*mu, cocktail));
	}

	t.lep_triggerMatchPt[w]  = userFloat(*mu, "TriggerMatchPt",  -999);
	t.lep_triggerMatchEta[w] = userFloat(*mu, "TriggerMatchEta", -999);

	t.lep_chi2dof[w] = mu->globalTrack()->normalizedChi2();
	t.lep_dB[w] = mu->dB();
	t.lep_sumPt[w] = mu->isolationR03().sumPt;
	t.lep_emEt[w] = mu->isolationR03().emEt;
	t.lep_hadEt[w] = mu->isolationR03().hadEt;
	t.lep_hoEt[w] = mu->isolationR03().hoEt;
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
	    
	t.lep_tk_numberOfValidTrackerHits[w] = mu->innerTrack()->hitPattern().numberOfValidTrackerHits();
	t.lep_tk_numberOfValidTrackerLayers[w] = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
	t.lep_tk_numberOfValidPixelHits[w] = mu->innerTrack()->hitPattern().numberOfValidPixelHits();
	t.lep_glb_numberOfValidTrackerHits[w] = mu->globalTrack()->hitPattern().numberOfValidTrackerHits();
	t.lep_glb_numberOfValidTrackerLayers[w] = mu->globalTrack()->hitPattern().trackerLayersWithMeasurement();
	t.lep_glb_numberOfValidPixelHits[w] = mu->globalTrack()->hitPattern().numberOfValidPixelHits();
	t.lep_glb_numberOfValidMuonHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
	t.lep_glb_muonStationsWithValidHits[w] = mu->globalTrack()->hitPattern().muonStationsWithValidHits();
	t.lep_numberOfMatches[w] = mu->numberOfMatches();
	t.lep_numberOfMatchedStations[w] = mu->numberOfMatchedStations();
	t.lep_isGlobalMuon[w] = mu->isGlobalMuon();
	t.lep_isTrackerMuon[w] = mu->isTrackerMuon();
      }
    }

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
    }
    else {
      t.cos_cs = -999;
      t.phi_cs = -999;
    }

    tree->Fill();
  }
}

DEFINE_FWK_MODULE(SimpleNtupler);
