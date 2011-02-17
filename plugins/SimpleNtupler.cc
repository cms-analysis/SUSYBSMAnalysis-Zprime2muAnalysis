#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
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
    float dil_mass;
    float dil_pt;
    float dil_rap;
    float dil_eta;
    float dil_phi;
    int lep_id[2];
    float lep_pt[2];
    float lep_eta[2];
    float lep_phi[2];
    float lep_tk_pt[2];
    float lep_tk_eta[2];
    float lep_tk_phi[2];
    float lep_glb_pt[2];
    float lep_glb_eta[2];
    float lep_glb_phi[2];
    float lep_triggerMatchPt[2];
    float lep_chi2dof[2];
    float lep_dB[2];
    float lep_sumPt[2];
    float lep_emEt[2];
    float lep_hadEt[2];
    float lep_hoEt[2];
    short lep_numberOfValidTrackerHits[2]; 
    short lep_numberOfValidPixelHits[2];
    short lep_numberOfValidMuonHits[2];
    short lep_muonStationsWithValidHits[2];
    bool lep_isGlobalMuon[2];
    bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool HLTPhysicsDeclared;
    bool GoodVtx;
    bool NoScraping;
    bool HLT_Single;
    bool HLT_Double;
  };

  tree_t t;
  TTree* tree;

  edm::InputTag hlt_src;
  edm::InputTag dimu_src;
};

SimpleNtupler::SimpleNtupler(const edm::ParameterSet& cfg)
  : hlt_src(cfg.getParameter<edm::InputTag>("hlt_src")),
    dimu_src(cfg.getParameter<edm::InputTag>("dimu_src"))
{
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("tt", &t, "run/i:lumi:event:dil_mass/F:dil_pt:dil_rap:dil_eta:dil_phi:lep_id[2]/I:lep_pt[2]/F:lep_eta[2]:lep_phi[2]:lep_tk_pt[2]:lep_tk_eta[2]:lep_tk_phi[2]:lep_glb_pt[2]:lep_glb_eta[2]:lep_glb_phi[2]:lep_triggerMatchPt[2]:lep_chi2dof[2]:lep_dB[2]:lep_sumPt[2]:lep_emEt[2]:lep_hadEt[2]:lep_hoEt[2]:lep_numberOfValidTrackerHits[2]/S:lep_numberOfValidPixelHits[2]:lep_numberOfValidMuonHits[2]:lep_muonStationsWithValidHits[2]:lep_isGlobalMuon[2]/O:lep_isTrackerMuon[2]:GoodDataRan:HLTPhysicsDeclared:GoodVtx:NoScraping:HLT_Single:HLT_Double");
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

  edm::Handle<edm::TriggerResults> reshlt;
  event.getByLabel(hlt_src, reshlt);
  const edm::TriggerNames& nameshlt = event.triggerNames(*reshlt);

  const unsigned r = event.id().run();
  if (r <= 147119) {
    t.HLT_Single = reshlt->accept(nameshlt.triggerIndex("HLT_Mu9")); // changing this to pt > 15 is taken care of by the VBTF selection
    t.HLT_Double = reshlt->accept(nameshlt.triggerIndex("HLT_DoubleMu3"));
  }
  else {
    t.HLT_Single = reshlt->accept(nameshlt.triggerIndex("HLT_Mu15_v1"));
    t.HLT_Double = reshlt->accept(nameshlt.triggerIndex("HLT_DoubleMu3_v2"));
  }

  edm::Handle<pat::CompositeCandidateCollection> dils;
  event.getByLabel(dimu_src, dils);

  BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils) {
    t.dil_mass = dil.mass();
    t.dil_pt = dil.pt();
    t.dil_rap = dil.rapidity();
    t.dil_eta = dil.eta();
    t.dil_phi = dil.phi();

    for (size_t i = 0; i < 2; ++i) {
      t.lep_id[i] = dil.daughter(i)->pdgId();
      t.lep_pt[i] = dil.daughter(i)->pt();
      t.lep_eta[i] = dil.daughter(i)->eta();
      t.lep_phi[i] = dil.daughter(i)->phi();

      if (abs(t.lep_id[i]) != 13) {
	t.lep_tk_pt[i] = -999;
	t.lep_tk_eta[i] = -999;
	t.lep_tk_phi[i] = -999;
	t.lep_glb_pt[i] = -999;
	t.lep_glb_eta[i] = -999;
	t.lep_glb_phi[i] = -999;
	t.lep_triggerMatchPt[i] = -999;
	t.lep_chi2dof[i] = -999;
	t.lep_dB[i] = -999;
	t.lep_sumPt[i] = -999;
	t.lep_emEt[i] = -999;
	t.lep_hadEt[i] = -999;
	t.lep_hoEt[i] = -999;
	t.lep_numberOfValidTrackerHits[i] = -999; 
	t.lep_numberOfValidPixelHits[i] = -999;
	t.lep_numberOfValidMuonHits[i] = -999;
	t.lep_muonStationsWithValidHits[i] = -999;
      }
      else {
	const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
	assert(mu);

	t.lep_tk_pt[i] = mu->innerTrack()->pt();
	t.lep_tk_eta[i] = mu->innerTrack()->eta();
	t.lep_tk_phi[i] = mu->innerTrack()->phi();
	t.lep_glb_pt[i] = mu->globalTrack()->pt();
	t.lep_glb_eta[i] = mu->globalTrack()->eta();
	t.lep_glb_phi[i] = mu->globalTrack()->phi();
	if (!mu->triggerObjectMatchesByPath("HLT_Mu15_v1").empty())
	  t.lep_triggerMatchPt[i] = mu->triggerObjectMatchesByPath("HLT_Mu15_v1").at(0).pt();
	else if (!mu->triggerObjectMatchesByPath("HLT_Mu9").empty())
	  t.lep_triggerMatchPt[i] = mu->triggerObjectMatchesByPath("HLT_Mu9").at(0).pt();
	else
	  t.lep_triggerMatchPt[i] = -999;
	t.lep_chi2dof[i] = mu->globalTrack()->normalizedChi2();
	t.lep_dB[i] = mu->dB();
	t.lep_sumPt[i] = mu->isolationR03().sumPt;
	t.lep_emEt[i] = mu->isolationR03().emEt;
	t.lep_hadEt[i] = mu->isolationR03().hadEt;
	t.lep_hoEt[i] = mu->isolationR03().hoEt;
	t.lep_numberOfValidTrackerHits[i] = mu->innerTrack()->hitPattern().numberOfValidTrackerHits();
	t.lep_numberOfValidPixelHits[i] = mu->innerTrack()->hitPattern().numberOfValidPixelHits();
	t.lep_numberOfValidMuonHits[i] = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
	t.lep_muonStationsWithValidHits[i] = mu->globalTrack()->hitPattern().muonStationsWithValidHits();
      }
    }

    tree->Fill();
  }
}

DEFINE_FWK_MODULE(SimpleNtupler);
