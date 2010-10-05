#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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
  //printf("tree size %i\n", sizeof(tree_t));
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("tt", &t, "run/i:lumi:event:dil_mass/F:dil_pt:dil_rap:dil_eta:dil_phi:lep_id[2]/I:lep_pt[2]/F:lep_eta[2]:lep_phi[2]:GoodDataRan/O:HLTPhysicsDeclared:GoodVtx:NoScraping:HLT_Single:HLT_Double");
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

  t.HLT_Single = reshlt->accept(nameshlt.triggerIndex("HLT_Mu9"));
  t.HLT_Double = reshlt->accept(nameshlt.triggerIndex("HLT_DoubleMu3"));

  edm::Handle<pat::CompositeCandidateCollection> dils;
  event.getByLabel(dimu_src, dils);

  BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils) {
    t.dil_mass = dil.mass();
    t.dil_pt = dil.pt();
    t.dil_rap = dil.rapidity();
    t.dil_eta = dil.eta();
    t.dil_phi = dil.phi();
    t.lep_id[0] = dil.daughter(0)->pdgId();
    t.lep_pt[0] = dil.daughter(0)->pt();
    t.lep_eta[0] = dil.daughter(0)->eta();
    t.lep_phi[0] = dil.daughter(0)->phi();
    t.lep_id[1] = dil.daughter(1)->pdgId();
    t.lep_pt[1] = dil.daughter(1)->pt();
    t.lep_eta[1] = dil.daughter(1)->eta();
    t.lep_phi[1] = dil.daughter(1)->phi();
    tree->Fill();
  }
}

DEFINE_FWK_MODULE(SimpleNtupler);
