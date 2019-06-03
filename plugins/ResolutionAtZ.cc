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
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

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
  TH2F* DileptonMass_2d_vsPtPlus;
  TH2F* DileptonMass_2d_vsPtPlus_BB;
  TH2F* DileptonMass_2d_vsPtPlus_BE;
  TH2F* DileptonMass_2d_vsPtMinus;
  TH2F* DileptonMass_2d_vsPtMinus_BB;
  TH2F* DileptonMass_2d_vsPtMinus_BE;

  TH2F* DileptonMass_2d_vsP;
  TH2F* DileptonMass_2d_vsP_BB;
  TH2F* DileptonMass_2d_vsP_BE;
  TH2F* DileptonMass_2d_vsPPlus;
  TH2F* DileptonMass_2d_vsPPlus_BB;
  TH2F* DileptonMass_2d_vsPPlus_BE;
  TH2F* DileptonMass_2d_vsPMinus;
  TH2F* DileptonMass_2d_vsPMinus_BB;
  TH2F* DileptonMass_2d_vsPMinus_BE;

  TH2F* DileptonMass_2d_vsPt_BB_neweta;
  TH2F* DileptonMass_2d_vsPt_BE_neweta;

  TH2F* DileptonMass_2d_vsPt_BB_neweta_negPhi;
  TH2F* DileptonMass_2d_vsPt_BE_neweta_negPhi;
  TH2F* DileptonMass_2d_vsPt_BB_neweta_midPhi;
  TH2F* DileptonMass_2d_vsPt_BE_neweta_midPhi;
  TH2F* DileptonMass_2d_vsPt_BB_neweta_posPhi;
  TH2F* DileptonMass_2d_vsPt_BE_neweta_posPhi;

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

   DileptonMass_2d_vsPt         = fs->make<TH2F>("DileptonMass_2d_vsPt",         titlePrefix + " dil. mass vs pt"    , 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPt_BB      = fs->make<TH2F>("DileptonMass_2d_vsPt_BB",      titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPt_BE      = fs->make<TH2F>("DileptonMass_2d_vsPt_BE",      titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPtPlus     = fs->make<TH2F>("DileptonMass_2d_vsPtPlus",     titlePrefix + " dil. mass vs pt"    , 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPtPlus_BB  = fs->make<TH2F>("DileptonMass_2d_vsPtPlus_BB",  titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPtPlus_BE  = fs->make<TH2F>("DileptonMass_2d_vsPtPlus_BE",  titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPtMinus    = fs->make<TH2F>("DileptonMass_2d_vsPtMinus",    titlePrefix + " dil. mass vs pt"    , 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPtMinus_BB = fs->make<TH2F>("DileptonMass_2d_vsPtMinus_BB", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPtMinus_BE = fs->make<TH2F>("DileptonMass_2d_vsPtMinus_BE", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);

   DileptonMass_2d_vsP         = fs->make<TH2F>("DileptonMass_2d_vsP",         titlePrefix + " dil. mass vs p"    , 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsP_BB      = fs->make<TH2F>("DileptonMass_2d_vsP_BB",      titlePrefix + " dil. mass vs p", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsP_BE      = fs->make<TH2F>("DileptonMass_2d_vsP_BE",      titlePrefix + " dil. mass vs p", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPPlus     = fs->make<TH2F>("DileptonMass_2d_vsPPlus",     titlePrefix + " dil. mass vs p"    , 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPPlus_BB  = fs->make<TH2F>("DileptonMass_2d_vsPPlus_BB",  titlePrefix + " dil. mass vs p", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPPlus_BE  = fs->make<TH2F>("DileptonMass_2d_vsPPlus_BE",  titlePrefix + " dil. mass vs p", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPMinus    = fs->make<TH2F>("DileptonMass_2d_vsPMinus",    titlePrefix + " dil. mass vs p"    , 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPMinus_BB = fs->make<TH2F>("DileptonMass_2d_vsPMinus_BB", titlePrefix + " dil. mass vs p", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPMinus_BE = fs->make<TH2F>("DileptonMass_2d_vsPMinus_BE", titlePrefix + " dil. mass vs p", 200, 50., 150., 2000., 0., 2000.);

   DileptonMass_2d_vsPt_BB_neweta = fs->make<TH2F>("DileptonMass_2d_vsPt_BB_neweta", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPt_BE_neweta = fs->make<TH2F>("DileptonMass_2d_vsPt_BE_neweta", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);

   DileptonMass_2d_vsPt_BB_neweta_negPhi = fs->make<TH2F>("DileptonMass_2d_vsPt_BB_neweta_midPhi", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPt_BE_neweta_negPhi = fs->make<TH2F>("DileptonMass_2d_vsPt_BE_neweta_midPhi", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPt_BB_neweta_midPhi = fs->make<TH2F>("DileptonMass_2d_vsPt_BB_neweta_midPhi", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPt_BE_neweta_midPhi = fs->make<TH2F>("DileptonMass_2d_vsPt_BE_neweta_midPhi", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPt_BB_neweta_midPhi = fs->make<TH2F>("DileptonMass_2d_vsPt_BB_neweta_midPhi", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);
   DileptonMass_2d_vsPt_BE_neweta_midPhi = fs->make<TH2F>("DileptonMass_2d_vsPt_BE_neweta_midPhi", titlePrefix + " dil. mass vs pt", 200, 50., 150., 2000., 0., 2000.);

}

void ResolutionAtZ::fillDileptonMassResolution(const reco::CompositeCandidate& dil) {
  const double mass         = dil.mass();    

  DileptonMass        ->Fill(mass);
  DileptonMass_2d_vsPt->Fill(mass, dil.daughter(0)->pt());
  DileptonMass_2d_vsPt->Fill(mass, dil.daughter(1)->pt());
  DileptonMass_2d_vsP->Fill(mass, dil.daughter(0)->p());
  DileptonMass_2d_vsP->Fill(mass, dil.daughter(1)->p());
  if (dil.daughter(0)->charge() > 0){
    	DileptonMass_2d_vsPtPlus->Fill(mass, dil.daughter(0)->pt());
    	DileptonMass_2d_vsPPlus->Fill(mass, dil.daughter(0)->p());
  }
  if (dil.daughter(1)->charge() > 0){
    	DileptonMass_2d_vsPtPlus->Fill(mass, dil.daughter(1)->pt());
    	DileptonMass_2d_vsPPlus->Fill(mass, dil.daughter(1)->p());
  }
  if (dil.daughter(0)->charge() < 0){
    	DileptonMass_2d_vsPtMinus->Fill(mass, dil.daughter(0)->pt());
    	DileptonMass_2d_vsPMinus->Fill(mass, dil.daughter(0)->p());
  }
  if (dil.daughter(1)->charge() < 0){
   	DileptonMass_2d_vsPtMinus->Fill(mass, dil.daughter(1)->pt());
   	DileptonMass_2d_vsPMinus->Fill(mass, dil.daughter(1)->p());
  }

  if (fabs(dil.daughter(0)->eta()) < 1.2 && fabs(dil.daughter(1)->eta()) < 1.2) { 
    DileptonMass_BB->Fill(mass);
    DileptonMass_2d_vsPt_BB->Fill(mass, dil.daughter(0)->pt());
    DileptonMass_2d_vsPt_BB->Fill(mass, dil.daughter(1)->pt());
    DileptonMass_2d_vsP_BB->Fill(mass, dil.daughter(0)->p());
    DileptonMass_2d_vsP_BB->Fill(mass, dil.daughter(1)->p());
    if (dil.daughter(0)->charge() > 0){
    	DileptonMass_2d_vsPtPlus_BB->Fill(mass, dil.daughter(0)->pt());
    	DileptonMass_2d_vsPPlus_BB->Fill(mass, dil.daughter(0)->p());
    }
    if (dil.daughter(1)->charge() > 0){
    	DileptonMass_2d_vsPtPlus_BB->Fill(mass, dil.daughter(1)->pt());
    	DileptonMass_2d_vsPPlus_BB->Fill(mass, dil.daughter(1)->p());
    }
    if (dil.daughter(0)->charge() < 0){
    	DileptonMass_2d_vsPtMinus_BB->Fill(mass, dil.daughter(0)->pt());
    	DileptonMass_2d_vsPMinus_BB->Fill(mass, dil.daughter(0)->p());
    }
    if (dil.daughter(1)->charge() < 0){
    	DileptonMass_2d_vsPtMinus_BB->Fill(mass, dil.daughter(1)->pt());
    	DileptonMass_2d_vsPMinus_BB->Fill(mass, dil.daughter(1)->p());
    }
  }
  else { 
    DileptonMass_BE->Fill(mass);
    DileptonMass_2d_vsPt_BE->Fill(mass, dil.daughter(0)->pt());
    DileptonMass_2d_vsPt_BE->Fill(mass, dil.daughter(1)->pt());
    DileptonMass_2d_vsP_BE->Fill(mass, dil.daughter(0)->p());
    DileptonMass_2d_vsP_BE->Fill(mass, dil.daughter(1)->p());
    if (dil.daughter(0)->charge() > 0){
    	DileptonMass_2d_vsPtPlus_BE->Fill(mass, dil.daughter(0)->pt());
    	DileptonMass_2d_vsPPlus_BE->Fill(mass, dil.daughter(0)->p());
    }
    if (dil.daughter(1)->charge() > 0){
    	DileptonMass_2d_vsPtPlus_BE->Fill(mass, dil.daughter(1)->pt());
    	DileptonMass_2d_vsPPlus_BE->Fill(mass, dil.daughter(1)->p());
    }
    if (dil.daughter(0)->charge() < 0){
    	DileptonMass_2d_vsPtMinus_BE->Fill(mass, dil.daughter(0)->pt());
    	DileptonMass_2d_vsPMinus_BE->Fill(mass, dil.daughter(0)->p());
    }
    if (dil.daughter(1)->charge() < 0){
    	DileptonMass_2d_vsPtMinus_BE->Fill(mass, dil.daughter(1)->pt());
    	DileptonMass_2d_vsPMinus_BE->Fill(mass, dil.daughter(1)->p());
    }

  }

  if (fabs(dil.daughter(0)->eta()) < 1.6 && fabs(dil.daughter(1)->eta()) < 1.6) {
    DileptonMass_2d_vsPt_BB_neweta->Fill(mass, dil.daughter(0)->pt());
    DileptonMass_2d_vsPt_BB_neweta->Fill(mass, dil.daughter(1)->pt());
  } else {
    DileptonMass_2d_vsPt_BE_neweta->Fill(mass, dil.daughter(0)->pt());
    DileptonMass_2d_vsPt_BE_neweta->Fill(mass, dil.daughter(1)->pt());

    if (dil.daughter(0)->phi()* (180.0/3.141592653589793238463) < -60.) DileptonMass_2d_vsPt_BE_neweta_negPhi->Fill(mass, dil.daughter(0)->pt());
    if (dil.daughter(1)->phi()* (180.0/3.141592653589793238463) < -60.) DileptonMass_2d_vsPt_BE_neweta_negPhi->Fill(mass, dil.daughter(1)->pt());
    if (dil.daughter(0)->phi()* (180.0/3.141592653589793238463) >= -60. && dil.daughter(0)->phi()* (180.0/3.141592653589793238463) < 60.) DileptonMass_2d_vsPt_BE_neweta_midPhi->Fill(mass, dil.daughter(0)->pt());
    if (dil.daughter(1)->phi()* (180.0/3.141592653589793238463) >= -60. && dil.daughter(1)->phi()* (180.0/3.141592653589793238463) < 60.) DileptonMass_2d_vsPt_BE_neweta_midPhi->Fill(mass, dil.daughter(1)->pt());
    if (dil.daughter(0)->phi()* (180.0/3.141592653589793238463) >= 60.) DileptonMass_2d_vsPt_BE_neweta_posPhi->Fill(mass, dil.daughter(0)->pt());
    if (dil.daughter(1)->phi()* (180.0/3.141592653589793238463) >= 60.) DileptonMass_2d_vsPt_BE_neweta_posPhi->Fill(mass, dil.daughter(1)->pt());
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
