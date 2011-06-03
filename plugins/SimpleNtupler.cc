#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

class SimpleNtupler : public edm::EDAnalyzer {
public:
  explicit SimpleNtupler(const edm::ParameterSet&);
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
    float lep_glb_pt[2];
    float lep_glb_pt_err[2];
    float lep_glb_eta[2];
    float lep_glb_phi[2];
    float lep_tpfms_pt[2];
    float lep_tpfms_pt_err[2];
    float lep_tpfms_eta[2];
    float lep_tpfms_phi[2];
    float lep_picky_pt[2];
    float lep_picky_pt_err[2];
    float lep_picky_eta[2];
    float lep_picky_phi[2];
    float lep_triggerMatchPt[2];
    float lep_chi2dof[2];
    float lep_dB[2];
    float lep_sumPt[2];
    float lep_emEt[2];
    float lep_hadEt[2];
    float lep_hoEt[2];
    short lep_tk_numberOfValidTrackerHits[2]; 
    short lep_tk_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidTrackerHits[2]; 
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
    bool firstOppDimu;
  };

  tree_t t;
  TTree* tree;

  edm::InputTag hlt_src;
  edm::InputTag dimu_src;
  edm::InputTag beamspot_src;
  edm::InputTag vertices_src;
};

SimpleNtupler::SimpleNtupler(const edm::ParameterSet& cfg)
  : hlt_src(cfg.getParameter<edm::InputTag>("hlt_src")),
    dimu_src(cfg.getParameter<edm::InputTag>("dimu_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertices_src(cfg.getParameter<edm::InputTag>("vertices_src"))
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
  tree->Branch("lep_pt", t.lep_pt, "lep_pt[2]/F");
  tree->Branch("lep_eta", t.lep_eta, "lep_eta[2]/F");
  tree->Branch("lep_phi", t.lep_phi, "lep_phi[2]/F");
  tree->Branch("lep_tk_pt", t.lep_tk_pt, "lep_tk_pt[2]/F");
  tree->Branch("lep_tk_pt_err", t.lep_tk_pt_err, "lep_tk_pt_err[2]/F");
  tree->Branch("lep_tk_eta", t.lep_tk_eta, "lep_tk_eta[2]/F");
  tree->Branch("lep_tk_phi", t.lep_tk_phi, "lep_tk_phi[2]/F");
  tree->Branch("lep_glb_pt", t.lep_glb_pt, "lep_glb_pt[2]/F");
  tree->Branch("lep_glb_pt_err", t.lep_glb_pt_err, "lep_glb_pt_err[2]/F");
  tree->Branch("lep_glb_eta", t.lep_glb_eta, "lep_glb_eta[2]/F");
  tree->Branch("lep_glb_phi", t.lep_glb_phi, "lep_glb_phi[2]/F");
  tree->Branch("lep_tpfms_pt", t.lep_tpfms_pt, "lep_tpfms_pt[2]/F");
  tree->Branch("lep_tpfms_pt_err", t.lep_tpfms_pt_err, "lep_tpfms_pt_err[2]/F");
  tree->Branch("lep_tpfms_eta", t.lep_tpfms_eta, "lep_tpfms_eta[2]/F");
  tree->Branch("lep_tpfms_phi", t.lep_tpfms_phi, "lep_tpfms_phi[2]/F");
  tree->Branch("lep_picky_pt", t.lep_picky_pt, "lep_picky_pt[2]/F");
  tree->Branch("lep_picky_pt_err", t.lep_picky_pt_err, "lep_picky_pt_err[2]/F");
  tree->Branch("lep_picky_eta", t.lep_picky_eta, "lep_picky_eta[2]/F");
  tree->Branch("lep_picky_phi", t.lep_picky_phi, "lep_picky_phi[2]/F");
  tree->Branch("lep_triggerMatchPt", t.lep_triggerMatchPt, "lep_triggerMatchPt[2]/F");
  tree->Branch("lep_chi2dof", t.lep_chi2dof, "lep_chi2dof[2]/F");
  tree->Branch("lep_dB", t.lep_dB, "lep_dB[2]/F");
  tree->Branch("lep_sumPt", t.lep_sumPt, "lep_sumPt[2]/F");
  tree->Branch("lep_emEt", t.lep_emEt, "lep_emEt[2]/F");
  tree->Branch("lep_hadEt", t.lep_hadEt, "lep_hadEt[2]/F");
  tree->Branch("lep_hoEt", t.lep_hoEt, "lep_hoEt[2]/F");
  tree->Branch("lep_tk_numberOfValidTrackerHits", t.lep_tk_numberOfValidTrackerHits, "lep_tk_numberOfValidTrackerHits[2]/S");
  tree->Branch("lep_tk_numberOfValidPixelHits", t.lep_tk_numberOfValidPixelHits, "lep_tk_numberOfValidPixelHits[2]/S");
  tree->Branch("lep_glb_numberOfValidTrackerHits", t.lep_glb_numberOfValidTrackerHits, "lep_glb_numberOfValidTrackerHits[2]/S");
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
  tree->Branch("firstOppDimu", &t.firstOppDimu, "firstOppDimu/O");

  tree->SetAlias("OurSel",
		 "("							\
		 "lep_isGlobalMuon[0] && "				\
		 "lep_pt[0] > 35 && "					\
		 "lep_tk_numberOfValidTrackerHits[0] >= 10 && "		\
		 "lep_sumPt[0] / lep_tk_pt[0] < 0.1"			\
		 ") && ("						\
		 "lep_isGlobalMuon[1] && "				\
		 "lep_pt[1] > 35 && "					\
		 "lep_tk_numberOfValidTrackerHits[1] >= 10 && "		\
		 "lep_sumPt[1] / lep_tk_pt[1] < 0.1"			\
		 ") && ( ("						\
		 "abs(lep_dB[0]) < 0.2 && "				\
		 "lep_chi2dof[0] < 10 && "				\
		 "lep_tk_numberOfValidPixelHits[0] >= 1 && "		\
		 "lep_glb_muonStationsWithValidHits[0] >= 2 && "	\
		 "lep_isTrackerMuon[0] && "				\
		 "lep_triggerMatchPt[0] > 30"				\
		 ") || ("						\
		 "abs(lep_dB[1]) < 0.2 && "				\
		 "lep_chi2dof[1] < 10 && "				\
		 "lep_tk_numberOfValidPixelHits[1] >= 1 && "		\
		 "lep_glb_muonStationsWithValidHits[1] >= 2 && "	\
		 "lep_isTrackerMuon[1] && "				\
		 "lep_triggerMatchPt[1] > 30"				\
		 ") ) && "						\
		 "lep_id[0] + lep_id[1] == 0 && "			\
		 "cos_angle > -0.9998 && "				\
		 "vertex_chi2 < 10 && "					\
		 "GoodDataRan && "					\
		 "HLTPhysicsDeclared && "				\
		 "NoScraping && "					\
		 "GoodVtx && "						\
		 "firstOppDimu");

  tree->SetAlias("VBTFSel",
		 "lep_isGlobalMuon[0] && "				\
		 "lep_isTrackerMuon[0] && "				\
		 "lep_tk_pt[0] > 35 && "				\
		 "abs(lep_tk_eta[0]) < 2.1 && "				\
		 "abs(lep_dB[0]) < 0.2 && "				\
		 "lep_sumPt[0] < 3 && "					\
		 "lep_glb_numberOfValidTrackerHits[0] > 10 && "		\
		 "lep_glb_numberOfValidPixelHits[0] > 0 && "		\
		 "lep_glb_numberOfValidMuonHits[0] > 0 && "		\
		 "lep_numberOfMatches[0] > 1 && "			\
		 "lep_isGlobalMuon[1] && "				\
		 "lep_isTrackerMuon[1] && "				\
		 "lep_tk_pt[1] > 35 && "				\
		 "abs(lep_tk_eta[1]) < 2.1 && "				\
		 "abs(lep_dB[1]) < 0.2 && "				\
		 "lep_sumPt[1] < 3 && "					\
		 "lep_glb_numberOfValidTrackerHits[1] > 10 && "		\
		 "lep_glb_numberOfValidPixelHits[1] > 0 && "		\
		 "lep_glb_numberOfValidMuonHits[1] > 0 && "		\
		 "lep_numberOfMatches[1] > 1 && "			\
		 "(lep_triggerMatchPt[0] > 30 || lep_triggerMatchPt[1] > 30) && " \
		 "lep_id[0] + lep_id[1] == 0");

  tree->SetAlias("OurNewSel",
		 "lep_isGlobalMuon[0] && "				\
		 "lep_isTrackerMuon[0] && "				\
		 "lep_pt[0] > 35 && "					\
		 "abs(lep_dB[0]) < 0.2 && "				\
		 "lep_chi2dof[0] < 10 && "				\
		 "lep_sumPt[0] / lep_tk_pt[0] < 0.1 && "		\
		 "lep_glb_numberOfValidTrackerHits[0] > 10 && "		\
		 "lep_glb_numberOfValidPixelHits[0] >= 1 && "		\
		 "lep_glb_numberOfValidMuonHits[0] > 0 && "		\
		 "lep_numberOfMatchedStations[0] > 1 && "		\
		 "lep_isGlobalMuon[1] && "				\
		 "lep_isTrackerMuon[1] && "				\
		 "lep_pt[1] > 35 && "					\
		 "abs(lep_dB[1]) < 0.2 && "				\
		 "lep_chi2dof[1] < 10 && "				\
		 "lep_sumPt[1] / lep_tk_pt[1] < 0.1 && "		\
		 "lep_glb_numberOfValidTrackerHits[1] > 10 && "		\
		 "lep_glb_numberOfValidPixelHits[1] >= 1 && "		\
		 "lep_glb_numberOfValidMuonHits[1] > 0 && "		\
		 "lep_numberOfMatchedStations[1] > 1 && "		\
		 "(lep_triggerMatchPt[0] > 30 || lep_triggerMatchPt[1] > 30) && " \
		 "lep_id[0] + lep_id[1] == 0 && "			\
		 "cos_angle > -0.9998 && "				\
		 "vertex_chi2 < 10 && "					\
		 "GoodDataRan && "					\
		 "HLTPhysicsDeclared && "				\
		 "NoScraping && "					\
		 "GoodVtx && "						\
		 "firstOppDimu");
}

float userFloat(const pat::CompositeCandidate& dil, const char* name, float def=-999.) {
  if (dil.hasUserFloat(name))
    return dil.userFloat(name);
  else
    return def;
}

int userInt(const pat::CompositeCandidate& dil, const char* name, int def=-999) {
  if (dil.hasUserInt(name))
    return dil.userInt(name);
  else
    return def;
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

  edm::Handle<pat::CompositeCandidateCollection> dils;
  event.getByLabel(dimu_src, dils);

  bool seen_first_oppsign_dimu = false;

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
    const bool dimuon = abs(dil.daughter(0)->pdgId()) == 13 && abs(dil.daughter(1)->pdgId()) == 13;

    if (opp_sign && dimuon) {
      if (!seen_first_oppsign_dimu) {
	t.firstOppDimu = true;
	seen_first_oppsign_dimu = true;
      }
      else
	t.firstOppDimu = false;
    }
      
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
      t.lep_pt[w] = dil.daughter(i)->pt();
      t.lep_eta[w] = dil.daughter(i)->eta();
      t.lep_phi[w] = dil.daughter(i)->phi();

      if (abs(t.lep_id[w]) != 13) {
	t.lep_pt_err[w] = -999;
	t.lep_tk_pt[w] = -999;
	t.lep_tk_pt_err[w] = -999;
	t.lep_tk_eta[w] = -999;
	t.lep_tk_phi[w] = -999;
	t.lep_glb_pt[w] = -999;
	t.lep_glb_pt_err[w] = -999;
	t.lep_glb_eta[w] = -999;
	t.lep_glb_phi[w] = -999;
	t.lep_tpfms_pt[w] = -999;
	t.lep_tpfms_pt_err[w] = -999;
	t.lep_tpfms_eta[w] = -999;
	t.lep_tpfms_phi[w] = -999;
	t.lep_picky_pt[w] = -999;
	t.lep_picky_pt_err[w] = -999;
	t.lep_picky_eta[w] = -999;
	t.lep_picky_phi[w] = -999;
	t.lep_triggerMatchPt[w] = -999;
	t.lep_chi2dof[w] = -999;
	t.lep_dB[w] = -999;
	t.lep_sumPt[w] = -999;
	t.lep_emEt[w] = -999;
	t.lep_hadEt[w] = -999;
	t.lep_hoEt[w] = -999;
	t.lep_tk_numberOfValidTrackerHits[w] = -999; 
	t.lep_tk_numberOfValidPixelHits[w] = -999;
	t.lep_glb_numberOfValidTrackerHits[w] = -999; 
	t.lep_glb_numberOfValidPixelHits[w] = -999;
	t.lep_glb_numberOfValidMuonHits[w] = -999;
	t.lep_glb_muonStationsWithValidHits[w] = -999;
	t.lep_numberOfMatches[w] = -999;
	t.lep_numberOfMatchedStations[w] = -999;
	t.lep_isGlobalMuon[w] = false;
	t.lep_isTrackerMuon[w] = false;
      }
      else {
	const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
	assert(mu);

	t.lep_tk_pt[w] = mu->innerTrack()->pt();
	t.lep_tk_pt_err[w] = mu->innerTrack()->ptError();
	t.lep_tk_eta[w] = mu->innerTrack()->eta();
	t.lep_tk_phi[w] = mu->innerTrack()->phi();
	t.lep_glb_pt[w] = mu->globalTrack()->pt();
	t.lep_glb_pt_err[w] = mu->globalTrack()->ptError();
	t.lep_glb_eta[w] = mu->globalTrack()->eta();
	t.lep_glb_phi[w] = mu->globalTrack()->phi();
	if (mu->tpfmsMuon().isNull()) {
	  t.lep_tpfms_pt[w] = -999;
	  t.lep_tpfms_pt_err[w] = -999;
	  t.lep_tpfms_eta[w] = -999;
	  t.lep_tpfms_phi[w] = -999;
	}
	else {
	  t.lep_tpfms_pt[w] = mu->tpfmsMuon()->pt();
	  t.lep_tpfms_pt_err[w] = mu->tpfmsMuon()->ptError();
	  t.lep_tpfms_eta[w] = mu->tpfmsMuon()->eta();
	  t.lep_tpfms_phi[w] = mu->tpfmsMuon()->phi();
	}
	if (mu->pickyMuon().isNull()) {
	  t.lep_picky_pt[w] = -999;
	  t.lep_picky_pt_err[w] = -999;
	  t.lep_picky_eta[w] = -999;
	  t.lep_picky_phi[w] = -999;
	}
	else {
	  t.lep_picky_pt[w] = mu->pickyMuon()->pt();
	  t.lep_picky_pt_err[w] = mu->pickyMuon()->ptError();
	  t.lep_picky_eta[w] = mu->pickyMuon()->eta();
	  t.lep_picky_phi[w] = mu->pickyMuon()->phi();
	}

	static const size_t n_single_mu_path_names = 8;
	static const char* single_mu_path_names[n_single_mu_path_names] = {"HLT_Mu30_v3", "HLT_Mu30_v2", "HLT_Mu30_v1", "HLT_Mu24_v2", "HLT_Mu24_v1", "HLT_Mu15_v2", "HLT_Mu15_v1", "HLT_Mu9"};
	t.lep_triggerMatchPt[w] = -999;
	for (size_t j = 0; j < n_single_mu_path_names; ++j) {
	  if (!mu->triggerObjectMatchesByPath(single_mu_path_names[j]).empty()) { 
	    t.lep_triggerMatchPt[w] = mu->triggerObjectMatchesByPath(single_mu_path_names[j]).at(0).pt();
	    break;
	  }
	}

	t.lep_chi2dof[w] = mu->globalTrack()->normalizedChi2();
	t.lep_dB[w] = mu->dB();
	t.lep_sumPt[w] = mu->isolationR03().sumPt;
	t.lep_emEt[w] = mu->isolationR03().emEt;
	t.lep_hadEt[w] = mu->isolationR03().hadEt;
	t.lep_hoEt[w] = mu->isolationR03().hoEt;
	t.lep_tk_numberOfValidTrackerHits[w] = mu->innerTrack()->hitPattern().numberOfValidTrackerHits();
	t.lep_tk_numberOfValidPixelHits[w] = mu->innerTrack()->hitPattern().numberOfValidPixelHits();
	t.lep_glb_numberOfValidTrackerHits[w] = mu->globalTrack()->hitPattern().numberOfValidTrackerHits();
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
