#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "math.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"

class ResolutionUsingMC : public edm::EDAnalyzer {
 public:
  explicit ResolutionUsingMC(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void fillLeptonResolution(const reco::GenParticle*, const reco::CandidateBaseRef&);
  void fillLeptonExtraMomentumResolution(const reco::GenParticle*, const reco::CandidateBaseRef&);
  void fillChargeResolution(const reco::GenParticle*, const reco::CandidateBaseRef&);
  void fillDileptonMassResolution(const reco::CompositeCandidate&);
  void fillLeptonHistos(const reco::CandidateBaseRef&);
  void fillLeptonHistos(const edm::View<reco::Candidate>&);
  void fillLeptonHistosFromDileptons(const pat::CompositeCandidateCollection&);
  void fillDileptonHistos(const pat::CompositeCandidateCollection&);

  HardInteraction hardInteraction;
  edm::InputTag lepton_src;
  edm::InputTag dilepton_src;
  const bool leptonsFromDileptons;
  const bool doQoverP;

  TH1F* weights;

  TH1F* LeptonEtaDiff;
  TH1F* LeptonPhiDiff;
  TH1F* LeptonPtDiff;
  TH2F* LeptonPtScatter;
  TH1F* LeptonPtRes;
  TH1F* LeptonPRes;
  TH1F* LeptonInvPtRes;
  TH1F* LeptonInvPRes;
  
  TH2F* LeptonPtVSEta;
  TH2F* LeptonPtVSTheta;

  TH2F* LeptonInvPResVPGen_2d;
  TH2F* LeptonInvPResVPGen_2d_B;
  TH2F* LeptonInvPResVPGen_2d_BO;
  TH2F* LeptonInvPResVPGen_2d_E;
  TH2F* LeptonInvPResVPGen_2d_EO;
  TH2F* LeptonInvPResVPGen_2d_O;

  TH2F* LeptonInvPtResVPtWeighted_2d;
  TH2F* LeptonInvPtResVPtWeighted_2d_B;
  TH2F* LeptonInvPtResVPtWeighted_2d_BO;
  TH2F* LeptonInvPtResVPtWeighted_2d_E;
  TH2F* LeptonInvPtResVPtWeighted_2d_EO;
  TH2F* LeptonInvPtResVPtWeighted_2d_O;
   
  TH2F* LeptonInvPtResVPtGen_2d;
  TH2F* LeptonInvPtResVPtGen_2d_B;
  TH2F* LeptonInvPtResVPtGen_2d_BO;
  TH2F* LeptonInvPtResVPtGen_2d_E;
  TH2F* LeptonInvPtResVPtGen_2d_EO;
  TH2F* LeptonInvPtResVPtGen_2d_E2;
  TH2F* LeptonInvPtResVPtGen_2d_O;
  TProfile* LeptonInvPtResVPtGen;
  TProfile* LeptonInvPResVPGen;
  TH1F* LeptonInvPtPull;
  TH1F* LeptonInvPPull;

  TH1F* LeptonInvPtResBy[W_L_MAX];
  TH1F* LeptonInvPResBy[W_L_MAX];
  TH1F* LeptonInvPtPullBy[W_L_MAX];
  TH1F* LeptonInvPPullBy[W_L_MAX];

  TH1F* ChargeDiff;
  TH1F* ChargeRightVInvPt;
  TH1F* ChargeWrongVInvPt;



  TH1F* LeptonPtLogB;
  TH1F* LeptonPtLogO;
  TH1F* LeptonPtLogE;


  TH1F* DileptonMassReco;
  TH1F* DileptonMassGen;
  TH1F* DileptonMassRes;
  TH1F* DileptonResMassRes;
  TH1F* ResonanceMassRes;
  TH1F* DileptonMassInvRes;

  TH2F* DileptonMassResVMass_2d;
  TH2F* DileptonMassResVMass_2d_BBTight;
  TH2F* DileptonMassResVMass_2d_BB;
  TH2F* DileptonMassResVMass_2d_BE;
  TH2F* DileptonMassResVMass_2d_EE;
  TH2F* DileptonMassResVMass_2d_OO;
  TH2F* DileptonMassResVMass_2d_vsPtsum;
  TH2F* DileptonMassResVMass_2d_vsLeadPt;
  TH2F* DileptonMassResVMass_2d_vsPt_B;
  TH2F* DileptonMassResVMass_2d_vsPt_E;

  TProfile* DileptonMassResVMass;
  TProfile* DileptonResMassResVMass;
  TProfile* ResonanceMassResVMass;
  TProfile* DileptonInvMassResVMass;

  
//  TH1F* DileptonMassResBy[W_D_MAX];
//  TH1F* DileptonResMassResBy[W_D_MAX];
//  TH1F* ResonanceMassResBy[W_D_MAX];
//
//  TH1F* DileptonMassResVsMassBy[M_MAX][W_D_MAX];
};

ResolutionUsingMC::ResolutionUsingMC(const edm::ParameterSet& cfg)
  : hardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")),
    lepton_src(cfg.getParameter<edm::InputTag>("lepton_src")),
    dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    leptonsFromDileptons(cfg.getParameter<bool>("leptonsFromDileptons")),
    doQoverP(cfg.getParameter<bool>("doQoverP"))
{

   consumes<reco::GenParticleCollection>(hardInteraction.src);
   consumes<edm::View<reco::Candidate>>(lepton_src);
   consumes<pat::CompositeCandidateCollection>(dilepton_src);
   


  std::string title_prefix = cfg.getUntrackedParameter<std::string>("titlePrefix", "");
  if (title_prefix.size() && title_prefix[title_prefix.size()-1] != ' ')
    title_prefix += " ";
  const TString titlePrefix(title_prefix.c_str());
  
  edm::Service<TFileService> fs;
  
  double massmax = 5000;
  int nbinsmass=250;
  
  LeptonEtaDiff = fs->make<TH1F>("LeptonEtaDiff", titlePrefix + "#eta - gen #eta", 100, -0.002, 0.002);
  LeptonPhiDiff = fs->make<TH1F>("LeptonPhiDiff", titlePrefix + "#phi - gen #phi", 100, -0.002, 0.002);
  LeptonPtDiff  = fs->make<TH1F>("LeptonPtDiff",  titlePrefix + "pT - gen pT", 100, -20, 20);
  LeptonPtScatter = fs->make<TH2F>("LeptonPtScatter", titlePrefix + "pT vs. gen pT", 100, 0, 2000, 100, 0, 2000);

  LeptonPtRes = fs->make<TH1F>("LeptonPtRes", titlePrefix + "(pT - gen pT)/(gen pT)", 100, -0.5, 0.5);
  LeptonPRes  = fs->make<TH1F>("LeptonPRes",  titlePrefix + "(p - gen p)/(gen p)",    100, -0.5, 0.5);

  LeptonInvPtRes = fs->make<TH1F>("LeptonInvPtRes", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT)", 100, -0.5, 0.5);
  LeptonInvPRes  = fs->make<TH1F>("LeptonInvPRes",  titlePrefix + "(1/p - 1/gen p)/(1/gen p)",    100, -0.5, 0.5);

  LeptonPtVSEta = fs->make<TH2F>("LeptonPtVSEta", titlePrefix + "pt vs. eta",  300,0, 3000,100, -5, 5);
  LeptonPtVSTheta = fs->make<TH2F>("LeptonPtVSTheta", titlePrefix + "pt vs. theta",  300,0, 3000,100, -TMath::Pi(), TMath::Pi());

  LeptonInvPResVPGen_2d    = fs->make<TH2F>("LeptonInvPResVPGen_2d",    titlePrefix + "(1/p - 1/gen p)/(1/gen p) vs. gen p",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPResVPGen_2d_B  = fs->make<TH2F>("LeptonInvPResVPGen_2d_B",  titlePrefix + "(1/p - 1/gen p)/(1/gen p) vs. gen p",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPResVPGen_2d_BO = fs->make<TH2F>("LeptonInvPResVPGen_2d_BO", titlePrefix + "(1/p - 1/gen p)/(1/gen p) vs. gen p",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPResVPGen_2d_E  = fs->make<TH2F>("LeptonInvPResVPGen_2d_E",  titlePrefix + "(1/p - 1/gen p)/(1/gen p) vs. gen p",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPResVPGen_2d_EO = fs->make<TH2F>("LeptonInvPResVPGen_2d_EO", titlePrefix + "(1/p - 1/gen p)/(1/gen p) vs. gen p",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPResVPGen_2d_O = fs->make<TH2F>("LeptonInvPResVPGen_2d_O", titlePrefix + "(1/p - 1/gen p)/(1/gen p) vs. gen p",  300,0, 3000,1500, -1, 1.5);

  LeptonInvPtResVPtWeighted_2d    = fs->make<TH2F>("LeptonInvPtResVPtWeighted_2d",    titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtWeighted_2d_B  = fs->make<TH2F>("LeptonInvPtResVPtWeighted_2d_B",  titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtWeighted_2d_BO = fs->make<TH2F>("LeptonInvPtResVPtWeighted_2d_BO", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtWeighted_2d_E  = fs->make<TH2F>("LeptonInvPtResVPtWeighted_2d_E",  titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtWeighted_2d_EO = fs->make<TH2F>("LeptonInvPtResVPtWeighted_2d_EO", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  300,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtWeighted_2d_O = fs->make<TH2F>("LeptonInvPtResVPtWeighted_2d_O", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  300,0, 3000,1500, -1, 1.5);
    
  LeptonInvPtResVPtGen_2d    = fs->make<TH2F>("LeptonInvPtResVPtGen_2d",    titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  600,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtGen_2d_B  = fs->make<TH2F>("LeptonInvPtResVPtGen_2d_B",  titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  600,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtGen_2d_BO = fs->make<TH2F>("LeptonInvPtResVPtGen_2d_BO", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  600,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtGen_2d_E2  = fs->make<TH2F>("LeptonInvPtResVPtGen_2d_E2",  titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  600,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtGen_2d_E  = fs->make<TH2F>("LeptonInvPtResVPtGen_2d_E",  titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  600,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtGen_2d_EO = fs->make<TH2F>("LeptonInvPtResVPtGen_2d_EO", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  600,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtGen_2d_O = fs->make<TH2F>("LeptonInvPtResVPtGen_2d_O", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT",  600,0, 3000,1500, -1, 1.5);
  LeptonInvPtResVPtGen       = fs->make<TProfile>("LeptonInvPtResVPtGen",   titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT", 50, 0, 1000, -0.5, 0.5);
  LeptonInvPResVPGen         = fs->make<TProfile>("LeptonInvPResVPGen",     titlePrefix + "(1/p - 1/gen p)/(1/gen p) vs. gen p",     50, 0, 1000, -0.5, 0.5);
  
  LeptonInvPtPull = fs->make<TH1F>("LeptonInvPtPull", titlePrefix + "(1/pT - 1/gen pT)/#sigma_{1/pT}", 100, -10, 10);
  LeptonInvPPull  = fs->make<TH1F>("LeptonInvPPull",  titlePrefix + "(1/p - 1/gen p)/#sigma_{1/p}",    100, -10, 10);

  static const TString lepton_where_names[W_L_MAX] = {"Barrel", "Overlap", "Endcap", "Outside"};
  for (size_t i = 0; i < W_L_MAX; ++i) {
    LeptonInvPtResBy [i] = fs->make<TH1F>("LeptonInvPtRes"  + lepton_where_names[i], titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT), "    + lepton_where_names[i], 100, -0.5, 0.5);
    LeptonInvPResBy  [i] = fs->make<TH1F>("LeptonInvPRes"   + lepton_where_names[i], titlePrefix + "(1/p - 1/gen p)/(1/gen p), "       + lepton_where_names[i], 100, -0.5, 0.5);
    LeptonInvPtPullBy[i] = fs->make<TH1F>("LeptonInvPtPull" + lepton_where_names[i], titlePrefix + "(1/pT - 1/gen pT)/#sigma_{1/pT}, " + lepton_where_names[i], 100, -10, 10);
    LeptonInvPPullBy [i] = fs->make<TH1F>("LeptonInvPPull"  + lepton_where_names[i], titlePrefix + "(1/p - 1/gen p)/#sigma_{1/p}, "    + lepton_where_names[i], 100, -10, 10);
  }

  ChargeDiff = fs->make<TH1F>("ChargeDiff", titlePrefix + "q - gen q", 7, -3.5, 3.5);  
  ChargeRightVInvPt = fs->make<TH1F>("ChargeRightVInvPt", titlePrefix + "right q vs. 1/(gen pT)", 50, 0, 0.01);
  ChargeWrongVInvPt = fs->make<TH1F>("ChargeWrongVInvPt", titlePrefix + "wrong q vs. 1/(gen pT)", 50, 0, 0.01);
  
  ChargeRightVInvPt->Sumw2();
  ChargeWrongVInvPt->Sumw2();
    
  LeptonPtLogB   = fs->make<TH1F>("LeptonPtLogB",   titlePrefix + "lepton pT", 90, 1, 4);
  LeptonPtLogO   = fs->make<TH1F>("LeptonPtLogO",   titlePrefix + "lepton pT", 90, 1, 4);
  LeptonPtLogE   = fs->make<TH1F>("LeptonPtLogE",   titlePrefix + "lepton pT", 90, 1, 4);



  DileptonMassReco   = fs->make<TH1F>("DileptonMassReco",   titlePrefix + "dil. mass", 20000, 0, 20000);
  DileptonMassGen    = fs->make<TH1F>("DileptonMassGen",    titlePrefix + "dil. mass gen", 20000, 0, 20000);
  DileptonMassRes    = fs->make<TH1F>("DileptonMassRes",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 100, -0.5, 0.5);
  DileptonResMassRes = fs->make<TH1F>("DileptonResMassRes", titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass)", 100, -0.5, 0.5);
  ResonanceMassRes   = fs->make<TH1F>("ResonanceMassRes",   titlePrefix + "(res. mass - gen res. mass)/(gen res. mass)", 100, -0.5, 0.5);
  DileptonMassInvRes = fs->make<TH1F>("DileptonMassInvRes", titlePrefix + "(1/dil. mass - 1/gen dil. mass)/(1/gen dil. mass)", 100, -0.5, 0.5);

  DileptonMassResVMass_2d          = fs->make<TH2F>("DileptonMassResVMass_2d",          titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_BB       = fs->make<TH2F>("DileptonMassResVMass_2d_BB",       titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_BBTight  = fs->make<TH2F>("DileptonMassResVMass_2d_BBTight",  titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_BE       = fs->make<TH2F>("DileptonMassResVMass_2d_BE",       titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_EE       = fs->make<TH2F>("DileptonMassResVMass_2d_EE",       titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_OO       = fs->make<TH2F>("DileptonMassResVMass_2d_OO",       titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_vsPtsum  = fs->make<TH2F>("DileptonMassResVMass_2d_vsPtsum",  titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_vsLeadPt = fs->make<TH2F>("DileptonMassResVMass_2d_vsLeadPt", titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_vsPt_B   = fs->make<TH2F>("DileptonMassResVMass_2d_vsPt_B",   titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);
  DileptonMassResVMass_2d_vsPt_E   = fs->make<TH2F>("DileptonMassResVMass_2d_vsPt_E",   titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 120, 0., 6000., 1000, -1., 1.0);

  DileptonMassResVMass    = fs->make<TProfile>("DileptonMassResVMass",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", nbinsmass,0, massmax, -1, 1);
  DileptonResMassResVMass = fs->make<TProfile>("DileptonResMassResVMass", titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass)", nbinsmass,0, massmax, -1, 1);
  ResonanceMassResVMass   = fs->make<TProfile>("ResonanceMassResVMass",   titlePrefix + "(res. mass - gen res. mass)/(gen res. mass)", nbinsmass,0, massmax, -1, 1);
  DileptonInvMassResVMass = fs->make<TProfile>("DileptonInvMassResVMass", titlePrefix + "(./dil. mass - 1/gen dil. mass)/(1/gen dil. mass)", nbinsmass,0, massmax, -1, 1);

  weights = fs->make<TH1F>("weights","weights",90,1,4);
weights->SetBinContent(1,0.0);                                                                                                                                                                                                                
weights->SetBinContent(2,0.0);                                                                                                                                                                                                                
weights->SetBinContent(3,0.0);                                                                                                                                                                                                                
weights->SetBinContent(4,0.0);                                                                                                                                                                                                                
weights->SetBinContent(5,0.0);                                                                                                                                                                                                                
weights->SetBinContent(6,0.0);                                                                                                                                                                                                                
weights->SetBinContent(7,0.0);                                                                                                                                                                                                                
weights->SetBinContent(8,0.0);                                                                                                                                                                                                                
weights->SetBinContent(9,0.0);                                                                                                                                                                                                                
weights->SetBinContent(10,0.0);                                                                                                                                                                                                               
weights->SetBinContent(11,0.0);                                                                                                                                                                                                               
weights->SetBinContent(12,0.0);                                                                                                                                                                                                               
weights->SetBinContent(13,0.0);                                                                                                                                                                                                               
weights->SetBinContent(14,0.0);                                                                                                                                                                                                               
weights->SetBinContent(15,0.0);                                                                                                                                                                                                               
weights->SetBinContent(16,0.0);                                                                                                                                                                                                               
weights->SetBinContent(17,0.0);                                                                                                                                                                                                               
weights->SetBinContent(18,0.0);                                                                                                                                                                                                               
weights->SetBinContent(19,0.0);                                                                                                                                                                                                               
weights->SetBinContent(20,0.0);                                                                                                                                                                                                               
weights->SetBinContent(21,0.0);                                                                                                                                                                                                               
weights->SetBinContent(22,1.064515471458435);                                                                                                                                                                                                 
weights->SetBinContent(23,0.3278546929359436);                                                                                                                                                                                                
weights->SetBinContent(24,0.344069242477417);                                                                                                                                                                                                 
weights->SetBinContent(25,0.419121652841568);                                                                                                                                                                                                 
weights->SetBinContent(26,0.3980500400066376);                                                                                                                                                                                                
weights->SetBinContent(27,0.5091821551322937);                                                                                                                                                                                                
weights->SetBinContent(28,0.555397093296051);                                                                                                                                                                                                 
weights->SetBinContent(29,0.7060061693191528);                                                                                                                                                                                                
weights->SetBinContent(30,0.6712177395820618);                                                                                                                                                                                                
weights->SetBinContent(31,0.7814756035804749);                                                                                                                                                                                                
weights->SetBinContent(32,0.8268739581108093);                                                                                                                                                                                                
weights->SetBinContent(33,1.0383925437927246);                                                                                                                                                                                                
weights->SetBinContent(34,1.2656561136245728);                                                                                                                                                                                                
weights->SetBinContent(35,1.217139482498169);                                                                                                                                                                                                 
weights->SetBinContent(36,1.3628764152526855);                                                                                                                                                                                                
weights->SetBinContent(37,1.9920384883880615);                                                                                                                                                                                                
weights->SetBinContent(38,1.828977108001709);                                                                                                                                                                                                 
weights->SetBinContent(39,1.9653478860855103);
weights->SetBinContent(40,2.5994772911071777);
weights->SetBinContent(41,3.3766162395477295);
weights->SetBinContent(42,2.8907573223114014);
weights->SetBinContent(43,4.098905563354492);
weights->SetBinContent(44,4.869252681732178);
weights->SetBinContent(45,4.353845596313477);
weights->SetBinContent(46,4.309980869293213);
weights->SetBinContent(47,3.7432472705841064);
weights->SetBinContent(48,12.955918312072754);
weights->SetBinContent(49,7.140853404998779);
weights->SetBinContent(50,5.972690105438232);
weights->SetBinContent(51,14.11471939086914);
weights->SetBinContent(52,16.666122436523438);
weights->SetBinContent(53,17.403629302978516);
weights->SetBinContent(54,36.072357177734375);
weights->SetBinContent(55,36.21830749511719);
weights->SetBinContent(56,13.891311645507812);
weights->SetBinContent(57,75.06332397460938);
weights->SetBinContent(58,63.87181091308594);
weights->SetBinContent(59,84.57400512695312);
weights->SetBinContent(60,107.95161437988281);
weights->SetBinContent(61,212.1295928955078);
weights->SetBinContent(62,251.29901123046875);
weights->SetBinContent(63,203.99586486816406);
weights->SetBinContent(64,193.7982635498047);
weights->SetBinContent(65,215.29000854492188);
weights->SetBinContent(66,364.6211242675781);
weights->SetBinContent(67,570.5145874023438);
weights->SetBinContent(68,373.2234191894531);
weights->SetBinContent(69,2058.141845703125);
weights->SetBinContent(70,2607.290771484375);
weights->SetBinContent(71,5086.4833984375);
weights->SetBinContent(72,5189.6787109375);
weights->SetBinContent(73,8548.02734375);
weights->SetBinContent(74,0.0);
weights->SetBinContent(75,0.0);
weights->SetBinContent(76,0.0);
weights->SetBinContent(77,0.0);
weights->SetBinContent(78,0.0);
weights->SetBinContent(79,0.0);
weights->SetBinContent(80,0.0);
weights->SetBinContent(81,0.0);
weights->SetBinContent(82,0.0);
weights->SetBinContent(83,0.0);
weights->SetBinContent(84,0.0);
weights->SetBinContent(85,0.0);
weights->SetBinContent(86,0.0);
weights->SetBinContent(87,0.0);
weights->SetBinContent(88,0.0);
weights->SetBinContent(89,0.0);


//  static const TString dilepton_where_names[W_D_MAX] = {"BB", "BE"};
//  static const TString dilepton_mass_names[M_MAX] = {"0to500", "500to1000", "1000to1500", "1500to2000", "2000to2500", "2500to3000", "3000to3500", "3500to4000", "4000to4500", "4500to5000","5000toInf"};
//  
//  for (size_t i = 0; i < W_D_MAX; ++i) {
//    DileptonMassResBy   [i] = fs->make<TH1F>("DileptonMassRes"    + dilepton_where_names[i], titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
//    DileptonResMassResBy[i] = fs->make<TH1F>("DileptonResMassRes" + dilepton_where_names[i], titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
//    ResonanceMassResBy  [i] = fs->make<TH1F>("ResonanceMassRes"   + dilepton_where_names[i], titlePrefix + "(res. mass - gen res. mass)/(gen res. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
//    
//    for (size_t j = 0; j < M_MAX; ++j) { 
//      DileptonMassResVsMassBy[j][i]    = fs->make<TH1F>("DileptonMassResVMass_" + dilepton_mass_names[j] + dilepton_where_names[i],    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 100, -0.5, 0.5);
//    }
//  }  
}

void ResolutionUsingMC::fillLeptonResolution(const reco::GenParticle* gen_lep, const reco::CandidateBaseRef& lep) {
  // Angular diffs.
  LeptonEtaDiff->Fill(lep->eta() - gen_lep->eta());
  LeptonPhiDiff->Fill(lep->phi() - gen_lep->phi());

  const double gen_pt  = gen_lep->pt();
  const double gen_p   = gen_lep->p();

  // Momentum diffs/resolutions.
  LeptonPtDiff->Fill(lep->pt() - gen_pt);
  LeptonPtScatter->Fill(gen_pt, lep->pt());

  LeptonPtRes->Fill((lep->pt() - gen_pt)/gen_pt);
  LeptonPRes ->Fill((lep->p()  - gen_p) /gen_p);

  const double inv_pt_diff = 1/lep->pt() - 1/gen_pt;
  const double inv_p_diff  = 1/lep->p()  - 1/gen_p;

  // Inverse momentum resolutions.
  LeptonInvPtRes->Fill(inv_pt_diff/(1/gen_pt));
  LeptonInvPRes ->Fill(inv_p_diff /(1/gen_p));

  

}

void ResolutionUsingMC::fillLeptonExtraMomentumResolution(const reco::GenParticle* gen_lep, const reco::CandidateBaseRef& lep) {
  // More histograms (pulls, profiles) for offline reconstructed
  // leptons.
  double inv_gen_pt, inv_gen_p, inv_pt, inv_p;
  
  if (doQoverP) {
    const int gen_q = gen_lep->charge();
    inv_gen_pt  = gen_q/gen_lep->pt();
    inv_gen_p   = gen_q/gen_lep->p();
    
    const int q = lep->charge();
    inv_pt = q/lep->pt();
    inv_p  = q/lep->p();
  }
  else {
    inv_gen_pt  = 1/gen_lep->pt();
    inv_gen_p   = 1/gen_lep->p();
    
    inv_pt = 1/lep->pt();
    inv_p  = 1/lep->p();
  }
  if(fabs(lep->eta())<0.9)LeptonPtLogB->Fill(log10(lep->pt()));
  if(fabs(lep->eta())>1.2)LeptonPtLogE->Fill(log10(lep->pt()));
  if(fabs(lep->eta())>0.9 && fabs(gen_lep->eta()) < 1.2)LeptonPtLogO->Fill(log10(lep->pt()));



  const double inv_pt_diff = inv_pt - inv_gen_pt;
  const double inv_p_diff  = inv_p  - inv_gen_p;

  const double inv_pt_res = inv_pt_diff/inv_gen_pt;
  const double inv_p_res  = inv_p_diff /inv_gen_p;
    
  LeptonPtVSEta->Fill(lep->pt(), lep->eta());
  LeptonPtVSTheta->Fill(lep->pt(), lep->theta());
  
  // Inverse momentum resolutions as a function of generated momenta.


  LeptonInvPtResVPtWeighted_2d->Fill(lep->pt(), inv_pt_res,weights->GetBinContent(weights->FindBin(log10(lep->pt()))));
  if(fabs(lep->eta())<0.9)LeptonInvPtResVPtWeighted_2d_B->Fill(lep->pt(), inv_pt_res,weights->GetBinContent(weights->FindBin(log10(lep->pt()))));
  if(fabs(lep->eta())<1.2)LeptonInvPtResVPtWeighted_2d_BO->Fill(lep->pt(), inv_pt_res,weights->GetBinContent(weights->FindBin(log10(lep->pt()))));
  if(fabs(lep->eta())>1.2)LeptonInvPtResVPtWeighted_2d_E->Fill(lep->pt(), inv_pt_res,weights->GetBinContent(weights->FindBin(log10(lep->pt()))));
  if(fabs(lep->eta())>0.9)LeptonInvPtResVPtWeighted_2d_EO->Fill(lep->pt(), inv_pt_res,weights->GetBinContent(weights->FindBin(log10(lep->pt()))));
  if(fabs(lep->eta())>0.9 && fabs(lep->eta()) < 1.2)LeptonInvPtResVPtWeighted_2d_O->Fill(lep->pt(), inv_pt_res,weights->GetBinContent(weights->FindBin(log10(lep->pt()))));


  LeptonInvPtResVPtGen_2d->Fill(gen_lep->pt(), inv_pt_res);
  if(fabs(gen_lep->eta())<0.9)LeptonInvPtResVPtGen_2d_B->Fill(gen_lep->pt(), inv_pt_res);
  if(fabs(gen_lep->eta())<1.2)LeptonInvPtResVPtGen_2d_BO->Fill(gen_lep->pt(), inv_pt_res);
  if(fabs(gen_lep->eta())>1.2)LeptonInvPtResVPtGen_2d_E->Fill(gen_lep->pt(), inv_pt_res);
  if(fabs(gen_lep->eta())>0.9)LeptonInvPtResVPtGen_2d_EO->Fill(gen_lep->pt(), inv_pt_res);
  if(fabs(gen_lep->eta())>0.9 && fabs(gen_lep->eta()) < 1.2)LeptonInvPtResVPtGen_2d_O->Fill(gen_lep->pt(), inv_pt_res);
  if(fabs(gen_lep->eta())>1.2 && fabs(gen_lep->eta()) < 1.6)LeptonInvPtResVPtGen_2d_E2->Fill(gen_lep->pt(), inv_pt_res);

  LeptonInvPResVPGen_2d->Fill(gen_lep->p(), inv_p_res);
  if(fabs(gen_lep->eta())<0.9)LeptonInvPResVPGen_2d_B->Fill(gen_lep->p(), inv_p_res);
  if(fabs(gen_lep->eta())<1.2)LeptonInvPResVPGen_2d_BO->Fill(gen_lep->p(), inv_p_res);
  if(fabs(gen_lep->eta())>1.2)LeptonInvPResVPGen_2d_E->Fill(gen_lep->p(), inv_p_res);
  if(fabs(gen_lep->eta())>0.9)LeptonInvPResVPGen_2d_EO->Fill(gen_lep->p(), inv_p_res);
  if(fabs(gen_lep->eta())>0.9 && fabs(gen_lep->eta()) < 1.2)LeptonInvPResVPGen_2d_O->Fill(gen_lep->p(), inv_p_res);
 
  LeptonInvPtResVPtGen->Fill(gen_lep->pt(), inv_pt_res*inv_pt_res);
  LeptonInvPResVPGen  ->Fill(gen_lep->p(),  inv_p_res*inv_p_res);
  
  // Try to get the reconstructed momentum errors for pulls.
  double inv_pt_error, inv_p_error;
  inv_pt_error = inv_p_error = -999;
  bool errorOK = false;

  const reco::Track* tk = 0;
  const pat::Muon* mu = toConcretePtr<pat::Muon>(lep);
  if (mu) tk = patmuon::getPickedTrack(*mu).get();
    
  if (tk) {
    errorOK = true;
    inv_pt_error = invPtError(&*tk);
    inv_p_error  = invPError(&*tk);
  }

  // Inverse momentum pulls.
  if (errorOK) {
    LeptonInvPtPull->Fill(inv_pt_diff/inv_pt_error);
    LeptonInvPPull ->Fill(inv_p_diff /inv_p_error);
  }
  
  // The above inverse momentum resolutions and pulls, except
  // separately for barrel, overlap, and endcap.
  size_t w = whereIsLepton(lep);
  LeptonInvPtResBy[w]->Fill(inv_pt_res);
  LeptonInvPResBy [w]->Fill(inv_p_res);
  if (errorOK) {
    LeptonInvPtPullBy[w]->Fill(inv_pt_diff/inv_pt_error);
    LeptonInvPPullBy [w]->Fill(inv_p_diff /inv_p_error);
  }
}

void ResolutionUsingMC::fillChargeResolution(const reco::GenParticle* gen_lep, const reco::CandidateBaseRef& lep) {
  const int delta_q = lep->charge() - gen_lep->charge();
  ChargeDiff->Fill(delta_q);
  if (delta_q == 0)
    ChargeRightVInvPt->Fill(1/gen_lep->pt());
  else //if (delta_q == 2), also don't throw exception?
    ChargeWrongVInvPt->Fill(1/gen_lep->pt());
}

void ResolutionUsingMC::fillDileptonMassResolution(const reco::CompositeCandidate& dil) {
  if (!hardInteraction.IsValidForRes())
    return;
  
  const double mass         = dil.mass();    
  const double gen_mass     = (hardInteraction.lepPlusNoIB->p4() + hardInteraction.lepMinusNoIB->p4()).mass();

  const double res_mass     = resonanceP4(dil).mass();
  const double gen_res_mass = hardInteraction.resonance->mass();

  const double rdil    = mass    /gen_mass     - 1;
  const double rdilres = mass    /gen_res_mass - 1;
  const double rres    = res_mass/gen_res_mass - 1;
  double invmres = (1./mass-1./gen_mass)/(1./gen_mass);
 
  DileptonMassReco  ->Fill(mass);
  DileptonMassGen   ->Fill(gen_mass);
  DileptonMassRes   ->Fill(rdil);
  DileptonResMassRes->Fill(rdilres);
  ResonanceMassRes  ->Fill(rres);
  DileptonMassInvRes->Fill(invmres);

  DileptonMassResVMass_2d   ->Fill(gen_mass,     rdil);
  if (dil.daughter(0)->eta()<=0.9 && dil.daughter(1)->eta()<=0.9 && dil.daughter(0)->eta()>=-0.9 && dil.daughter(1)->eta()>=-0.9) 
    //    if   (abs(dil.daughter(0)->eta())<1.2 && abs(dil.daughter(1)->eta())<1.2 )   
    DileptonMassResVMass_2d_BBTight ->Fill(gen_mass,     rdil);
  if (dil.daughter(0)->eta()<=1.2 && dil.daughter(1)->eta()<=1.2 && dil.daughter(0)->eta()>=-1.2 && dil.daughter(1)->eta()>=-1.2) 
    //    if   (abs(dil.daughter(0)->eta())<1.2 && abs(dil.daughter(1)->eta())<1.2 )   
    DileptonMassResVMass_2d_BB ->Fill(gen_mass,     rdil);
  else                                                                         
    DileptonMassResVMass_2d_BE ->Fill(gen_mass,     rdil);
  if (fabs(dil.daughter(0)->eta())>=1.2 && fabs(dil.daughter(1)->eta())>=1.2) 
    //    if   (abs(dil.daughter(0)->eta())<1.2 && abs(dil.daughter(1)->eta())<1.2 )   
    DileptonMassResVMass_2d_EE ->Fill(gen_mass,     rdil);
  if (fabs(dil.daughter(0)->eta()) >=0.9 && fabs(dil.daughter(1)->eta())>=0.9 && fabs(dil.daughter(0)->eta())<=1.2 && fabs(dil.daughter(1)->eta())<=1.2) 
    //    if   (abs(dil.daughter(0)->eta())<1.2 && abs(dil.daughter(1)->eta())<1.2 )   
    DileptonMassResVMass_2d_OO ->Fill(gen_mass,     rdil);
 
  DileptonMassResVMass_2d_vsPtsum   ->Fill(hardInteraction.lepPlusNoIB->pt()+hardInteraction.lepMinusNoIB->pt(),     rdil);
  if (hardInteraction.lepPlusNoIB->pt() > hardInteraction.lepMinusNoIB->pt())  DileptonMassResVMass_2d_vsLeadPt   ->Fill(hardInteraction.lepPlusNoIB->pt(),     rdil);
  else                                                                 DileptonMassResVMass_2d_vsLeadPt   ->Fill(hardInteraction.lepMinusNoIB->pt(),     rdil);
  
  if   (abs(hardInteraction.lepPlusNoIB->eta()) < 1.2) DileptonMassResVMass_2d_vsPt_B ->Fill(hardInteraction.lepPlusNoIB->pt(),     rdil);
  else                                             DileptonMassResVMass_2d_vsPt_E ->Fill(hardInteraction.lepPlusNoIB->pt(),     rdil);
  if   (abs(hardInteraction.lepMinusNoIB->eta()) < 1.2) DileptonMassResVMass_2d_vsPt_B ->Fill(hardInteraction.lepMinusNoIB->pt(),     rdil);
  else                                              DileptonMassResVMass_2d_vsPt_E ->Fill(hardInteraction.lepMinusNoIB->pt(),     rdil);
  
  DileptonMassResVMass   ->Fill(gen_mass,     rdil*rdil);
  DileptonResMassResVMass->Fill(gen_res_mass, rdilres*rdilres);
  ResonanceMassResVMass  ->Fill(gen_res_mass, rres*rres);
  DileptonInvMassResVMass   ->Fill(gen_mass,  invmres*invmres);   
}

const reco::GenParticle* getGenParticle(const reco::CandidateBaseRef& lep) {
  // Unfortunately, cannot cast to PATObject since it is templated to
  // the Muon or Electron classes anyway, so just try each of the two.
  const pat::Muon* mu = toConcretePtr<pat::Muon>(lep);
  if (mu) return mu->genParticle();
  const pat::Electron* el = toConcretePtr<pat::Electron>(lep);
  if (el) return el->genParticle();
  return 0;
}

void ResolutionUsingMC::fillLeptonHistos(const reco::CandidateBaseRef& lep) {
  const reco::GenParticle* gen_lep = getGenParticle(lep);
  if (gen_lep != 0) {
    fillLeptonResolution(gen_lep, lep);
    fillLeptonExtraMomentumResolution(gen_lep, lep);
    fillChargeResolution(gen_lep, lep);
  }
}

void ResolutionUsingMC::fillLeptonHistos(const edm::View<reco::Candidate>& leptons) {
  // Fill lepton histos from all leptons.
  for (size_t i = 0, n = leptons.size(); i < n; ++i)
    fillLeptonHistos(leptons.refAt(i));
}

void ResolutionUsingMC::fillLeptonHistosFromDileptons(const pat::CompositeCandidateCollection& dileptons) {
  for (pat::CompositeCandidateCollection::const_iterator dil = dileptons.begin(), dile = dileptons.end(); dil != dile; ++dil)
    for (size_t i = 0; i < dil->numberOfDaughters(); ++i)
      fillLeptonHistos(dil->daughter(i)->masterClone());
}

void ResolutionUsingMC::fillDileptonHistos(const pat::CompositeCandidateCollection& dileptons) {
  for (pat::CompositeCandidateCollection::const_iterator dil = dileptons.begin(), dile = dileptons.end(); dil != dile; ++dil)
    fillDileptonMassResolution(*dil);
}

void ResolutionUsingMC::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  hardInteraction.Fill(event);
  
  edm::Handle<edm::View<reco::Candidate> > leptons;
  event.getByLabel(lepton_src, leptons);

  if (!leptons.isValid())
    edm::LogWarning("LeptonHandleInvalid") << "tried to get " << lepton_src << " with edm::Handle<edm::View<reco::Candidate> > and failed!";
  else {
    if (!leptonsFromDileptons)
      fillLeptonHistos(*leptons);
  }

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);
  
  if (!dileptons.isValid())
    edm::LogWarning("DileptonHandleInvalid") << "tried to get " << dilepton_src << " and failed!";
  else {
    if (leptonsFromDileptons)
      fillLeptonHistosFromDileptons(*dileptons);
    
    fillDileptonHistos(*dileptons);
  }
}

DEFINE_FWK_MODULE(ResolutionUsingMC);
