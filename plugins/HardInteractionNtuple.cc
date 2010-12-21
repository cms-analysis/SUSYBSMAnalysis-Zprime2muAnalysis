#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
  };

  tree_t t;
  TTree* tree;

  HardInteraction hi;
};

HardInteractionNtuple::HardInteractionNtuple(const edm::ParameterSet& cfg)
  : hi(cfg.getParameter<edm::ParameterSet>("hardInteraction"))
{
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("tt", &t, "run/i:lumi:event:res_mass/F:res_pt:res_rap:res_eta:res_phi:dil_mass:dil_pt:dil_rap:dil_eta:dil_phi:lep_pt[2]:lep_eta[2]:lep_phi[2]:lep_noib_pt[2]:lep_noib_eta[2]:lep_noib_phi[2]");
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

  tree->Fill();
}

DEFINE_FWK_MODULE(HardInteractionNtuple);
