#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"

class ResolutionAtZ : public edm::EDAnalyzer {
 public:
  explicit ResolutionAtZ(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void fillDileptonMassResolution(const reco::CompositeCandidate&);
  void fillDileptonHistos(const pat::CompositeCandidateCollection&);

  edm::InputTag lepton_src;
  edm::InputTag dilepton_src;
  const bool leptonsFromDileptons;
  const bool doQoverP;
  
  TH1F* DileptonMass;
  TH1F* DileptonMass_BB;
  TH1F* DileptonMass_BE;
  TH2F* DileptonMass_2d_vsPt;
  TH2F* DileptonMass_2d_vsPt_BB;
  TH2F* DileptonMass_2d_vsPt_BE;

};

ResolutionAtZ::ResolutionAtZ(const edm::ParameterSet& cfg)
  : lepton_src(cfg.getParameter<edm::InputTag>("lepton_src")),
    dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    leptonsFromDileptons(cfg.getParameter<bool>("leptonsFromDileptons")),
    doQoverP(cfg.getParameter<bool>("doQoverP"))
{

   consumes<edm::View<reco::Candidate>>(lepton_src);
   consumes<pat::CompositeCandidateCollection>(dilepton_src);
   

   std::string title_prefix = cfg.getUntrackedParameter<std::string>("titlePrefix", "");
   if (title_prefix.size() && title_prefix[title_prefix.size()-1] != ' ')
     title_prefix += " ";
   const TString titlePrefix(title_prefix.c_str());
   
   edm::Service<TFileService> fs;
   //   double ptbins[10]  = { 20, 30, 50, 70, 100, 150, 200, 250, 300, 600};   
   DileptonMass    = fs->make<TH1F>("DileptonMass"   , titlePrefix + "dil. mass", 2000, 0, 2000);
   DileptonMass_BB = fs->make<TH1F>("DileptonMass_BB", titlePrefix + "dil. mass", 2000, 0, 2000);
   DileptonMass_BE = fs->make<TH1F>("DileptonMass_BE", titlePrefix + "dil. mass", 2000, 0, 2000);
   
   DileptonMass_2d_vsPt   = fs->make<TH2F>("DileptonMass_2d_vsPt", titlePrefix + " dil. mass vs pt"    , 200, 50., 150., 500., 0., 2000.);
   DileptonMass_2d_vsPt_BB = fs->make<TH2F>("DileptonMass_2d_vsPt_BB", titlePrefix + " dil. mass vs pt", 200, 50., 150., 500., 0., 2000.);
   DileptonMass_2d_vsPt_BE = fs->make<TH2F>("DileptonMass_2d_vsPt_BE", titlePrefix + " dil. mass vs pt", 200, 50., 150., 500., 0., 2000.);
}

void ResolutionAtZ::fillDileptonMassResolution(const reco::CompositeCandidate& dil) {
  const double mass         = dil.mass();    

  DileptonMass        ->Fill(mass);
  DileptonMass_2d_vsPt->Fill(mass, dil.daughter(0)->pt());
  DileptonMass_2d_vsPt->Fill(mass, dil.daughter(1)->pt());
  
  if (abs(dil.daughter(0)->eta()) < 1.2 && abs(dil.daughter(1)->eta()) < 1.2) { 
    DileptonMass_BB->Fill(mass);
    DileptonMass_2d_vsPt_BB->Fill(mass, dil.daughter(0)->pt());
    DileptonMass_2d_vsPt_BB->Fill(mass, dil.daughter(1)->pt());
  }
  else { 
    DileptonMass_BE->Fill(mass);
    DileptonMass_2d_vsPt_BE->Fill(mass, dil.daughter(0)->pt());
    DileptonMass_2d_vsPt_BE->Fill(mass, dil.daughter(1)->pt());
  }
}


void ResolutionAtZ::fillDileptonHistos(const pat::CompositeCandidateCollection& dileptons) {
  for (pat::CompositeCandidateCollection::const_iterator dil = dileptons.begin(), dile = dileptons.end(); dil != dile; ++dil)
    fillDileptonMassResolution(*dil);
}

void ResolutionAtZ::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  
  edm::Handle<edm::View<reco::Candidate> > leptons;
  event.getByLabel(lepton_src, leptons);

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);

  if (!dileptons.isValid())
    edm::LogWarning("DileptonHandleInvalid") << "tried to get " << dilepton_src << " and failed!";
  else {
    fillDileptonHistos(*dileptons);
  }
}

DEFINE_FWK_MODULE(ResolutionAtZ);
