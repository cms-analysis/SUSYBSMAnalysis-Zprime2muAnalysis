#include "TH1F.h"
#include "TMath.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

class GenPlots : public edm::EDAnalyzer {
 public:
  explicit GenPlots(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  HardInteraction hardInteraction;
  const bool both_in_acc;

  TH1F* res_mass;
  TH1F* dil_mass;
  TH1F* res_pt;
  TH1F* dil_pt;
  TH1F* res_rap;
  TH1F* dil_rap;
  TH1F* res_eta;
  TH1F* dil_eta;
  TH1F* res_phi;
  TH1F* dil_phi;
};

GenPlots::GenPlots(const edm::ParameterSet& cfg)
  : hardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")),
    both_in_acc(cfg.getParameter<bool>("both_in_acc"))
{
  edm::Service<TFileService> fs;
  res_mass = fs->make<TH1F>("res_mass", "", 2000, 0, 2000);
  dil_mass = fs->make<TH1F>("dil_mass", "", 2000, 0, 2000);
  res_pt = fs->make<TH1F>("res_pt", "", 2000, 0, 2000);
  dil_pt = fs->make<TH1F>("dil_pt", "", 2000, 0, 2000);
  res_rap = fs->make<TH1F>("res_rap", "", 100, -5, 5);
  dil_rap = fs->make<TH1F>("dil_rap", "", 100, -5, 5);
  res_eta = fs->make<TH1F>("res_eta", "", 100, -5, 5);
  dil_eta = fs->make<TH1F>("dil_eta", "", 100, -5, 5);
  res_phi = fs->make<TH1F>("res_phi", "", 100, -TMath::Pi(), TMath::Pi());
  dil_phi = fs->make<TH1F>("dil_phi", "", 100, -TMath::Pi(), TMath::Pi());
}

void GenPlots::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  hardInteraction.Fill(event);

  if (both_in_acc && !(fabs(hardInteraction.lepMinus->eta()) < 2.4 && fabs(hardInteraction.lepPlus->eta()) < 2.4))
    return;
  
  res_mass->Fill(hardInteraction.resonance->mass());
  res_pt->Fill(hardInteraction.resonance->pt());
  res_eta->Fill(hardInteraction.resonance->eta());
  res_phi->Fill(hardInteraction.resonance->phi());
  res_rap->Fill(hardInteraction.resonance->rapidity());

  reco::Particle::LorentzVector dil = hardInteraction.dilepton();
  dil_mass->Fill(dil.mass());
  dil_pt->Fill(dil.pt());
  dil_eta->Fill(dil.eta());
  dil_phi->Fill(dil.phi());
  dil_rap->Fill(dil.Rapidity());
}

DEFINE_FWK_MODULE(GenPlots);
