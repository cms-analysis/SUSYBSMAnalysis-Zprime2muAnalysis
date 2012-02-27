#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

class HardInteractionNtuple : public edm::EDAnalyzer {
 public:
  explicit HardInteractionNtuple(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  struct tree_t {
    unsigned run;
    unsigned lumi;
    unsigned event;
    float res_mass;
    float res_pt;
    float res_rap;
    float res_eta;
    float res_phi;
    float dil_mass;
    float dil_pt;
    float dil_rap;
    float dil_eta;
    float dil_phi;
    float lep_pt[2];
    float lep_eta[2];
    float lep_phi[2];
    float lep_noib_pt[2];
    float lep_noib_eta[2];
    float lep_noib_phi[2];
    float reco_lep_pt[2];
    float reco_lep_pt_err[2];
    float reco_lep_eta[2];
    float reco_lep_phi[2];
    short reco_lep_q[2];
    bool reco_lep_passes[2];
  };

  tree_t t;
  TTree* tree;

  HardInteraction hi;
  const bool fill_reco;
  const edm::InputTag muon_src, picky_src;
};

HardInteractionNtuple::HardInteractionNtuple(const edm::ParameterSet& cfg)
  : hi(cfg.getParameter<edm::ParameterSet>("hardInteraction")),
    fill_reco(cfg.existsAs<edm::InputTag>("muon_src")),
    muon_src (fill_reco ? cfg.getParameter<edm::InputTag>("muon_src")  : edm::InputTag()),
    picky_src(fill_reco ? cfg.getParameter<edm::InputTag>("picky_src") : edm::InputTag())
{
  
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  TString branch_spec = "run/i:lumi:event:res_mass/F:res_pt:res_rap:res_eta:res_phi:dil_mass:dil_pt:dil_rap:dil_eta:dil_phi:lep_pt[2]:lep_eta[2]:lep_phi[2]:lep_noib_pt[2]:lep_noib_eta[2]:lep_noib_phi[2]";
  if (fill_reco)
    branch_spec += ":reco_lep_pt[2]:reco_lep_pt_err[2]:reco_lep_eta[2]:reco_lep_phi[2]:reco_lep_q[2]/S:reco_lep_passes[2]/O";
  tree->Branch("tt", &t, branch_spec);
}

bool passes(const reco::Muon& mu) {
  return
    mu.isGlobalMuon() &&
    mu.isTrackerMuon() &&
    mu.isolationR03().sumPt / mu.innerTrack()->pt() < 0.1 &&
    mu.globalTrack()->hitPattern().trackerLayersWithMeasurement() > 8 &&
    mu.globalTrack()->hitPattern().numberOfValidPixelHits() >= 1 &&
    mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
    mu.numberOfMatchedStations() > 1;
}

void HardInteractionNtuple::analyze(const edm::Event& event, const edm::EventSetup&) {
  memset(&t, 0, sizeof(tree_t));

  t.run = event.id().run();
  t.lumi = event.luminosityBlock();
  t.event = event.id().event();

  hi.Fill(event);

  t.res_mass = hi.resonance->mass();
  t.res_pt = hi.resonance->pt();
  t.res_rap = hi.resonance->rapidity();
  t.res_eta = hi.resonance->eta();
  t.res_phi = hi.resonance->phi();

  t.dil_mass = hi.dilepton().mass();
  t.dil_pt = hi.dilepton().pt();
  t.dil_rap = hi.dilepton().Rapidity();
  t.dil_eta = hi.dilepton().eta();
  t.dil_phi = hi.dilepton().phi();

  t.lep_pt[0] = hi.lepMinus->pt();
  t.lep_eta[0] = hi.lepMinus->eta();
  t.lep_phi[0] = hi.lepMinus->phi();

  t.lep_pt[1] = hi.lepPlus->pt();
  t.lep_eta[1] = hi.lepPlus->eta();
  t.lep_phi[1] = hi.lepPlus->phi();

  t.lep_noib_pt[0] = hi.lepMinusNoIB->pt();
  t.lep_noib_eta[0] = hi.lepMinusNoIB->eta();
  t.lep_noib_phi[0] = hi.lepMinusNoIB->phi();

  t.lep_noib_pt[1] = hi.lepPlusNoIB->pt();
  t.lep_noib_eta[1] = hi.lepPlusNoIB->eta();
  t.lep_noib_phi[1] = hi.lepPlusNoIB->phi();

  if (fill_reco) {
    edm::Handle<reco::MuonCollection> muons;
    event.getByLabel(muon_src, muons);

    edm::Handle<reco::TrackToTrackMap> picky;
    event.getByLabel(picky_src, picky);

    double min_dR[2] = {1e99, 1e99};

    t.reco_lep_q[0] = t.reco_lep_q[1] = 0; // Indicates no reco muon found (e.g. out of acceptance).

    // For each muon of the generated pair, try to find a reconstructed
    // muon within deltaR < 0.2 that is the best match to it. Don't try
    // to arbitrate.
    for (int i = 0; i < 2; ++i) {
      BOOST_FOREACH(const reco::Muon& mu, *muons) {
	if (!mu.isGlobalMuon())
	  continue;

	reco::TrackRef pk;
	reco::TrackToTrackMap::const_iterator it = picky->find(mu.globalTrack());
	if (it != picky->end()) pk = it->val;
	if (pk.isNull())
	  continue;

	double dR = reco::deltaR(pk->eta(), pk->phi(), t.lep_eta[i], t.lep_phi[i]);
	if (dR > 0.2 || dR > min_dR[i])
	  continue;
      
	min_dR[i] = dR;
      
	t.reco_lep_q[i] = pk->charge();
	t.reco_lep_pt[i] = pk->pt();
	t.reco_lep_pt_err[i] = pk->ptError();
	t.reco_lep_eta[i] = pk->eta();
	t.reco_lep_phi[i] = pk->phi();
	t.reco_lep_passes[i] = passes(mu);
      }
    }
  }

  tree->Fill();
}

DEFINE_FWK_MODULE(HardInteractionNtuple);
