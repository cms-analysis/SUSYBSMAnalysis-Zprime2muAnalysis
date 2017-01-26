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
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"

class ResolutionVertexUsingMC : public edm::EDAnalyzer {
 public:
  explicit ResolutionVertexUsingMC(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void fillLeptonResolution(const reco::GenParticle*, const reco::CandidateBaseRef&);
  void fillLeptonExtraMomentumResolution(const reco::GenParticle*, const reco::CandidateBaseRef&);
  void fillChargeResolution(const reco::GenParticle*, const reco::CandidateBaseRef&);
  void fillDileptonMassResolution(const pat::CompositeCandidate&);
  void fillLeptonHistos(const reco::CandidateBaseRef&);
  void fillLeptonHistos(const edm::View<reco::Candidate>&);
  void fillLeptonHistosFromDileptons(const pat::CompositeCandidateCollection&);
  void fillDileptonHistos(const pat::CompositeCandidateCollection&);

  HardInteraction hardInteraction;
  edm::InputTag lepton_src;
  edm::InputTag dilepton_src;
  const bool leptonsFromDileptons;
  const bool doQoverP;

  TH1F* LeptonEtaDiff;
  TH1F* LeptonPhiDiff;
  TH1F* LeptonPtDiff;
  TH2F* LeptonPtScatter;
  TH1F* LeptonPtRes;
  TH1F* LeptonPRes;
  TH1F* LeptonInvPtRes;
  TH1F* LeptonInvPRes;
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

  TH1F* DileptonMassReco;
  TH1F* DileptonMassGen;
  TH1F* DileptonMassRes;
  TH1F* DileptonResMassRes;
  TH1F* ResonanceMassRes;
  TH1F* DileptonMassInvRes;

  TProfile* DileptonMassResVMass;
  TProfile* DileptonResMassResVMass;
  TProfile* ResonanceMassResVMass;
  TProfile* DileptonInvMassResVMass;

  TH1F* DileptonMassResBy[W_D_MAX];
  TH1F* DileptonResMassResBy[W_D_MAX];
  TH1F* ResonanceMassResBy[W_D_MAX];
  
  /////////////////
  TH1F* DileptonMassRes_BarrelBarrel;
  TH1F* DileptonMassRes_EndcapBarrel;
  TH1F* DileptonMassRes_EndcapEndcap;
    TH1F* LeptonPt_BarrelBarrel;
    TH1F* LeptonPt_EndcapBarrel;
    TH1F* LeptonPt_EndcapEndcap;
    TH1F* LeptonPtRes_BarrelBarrel;
    TH1F* LeptonPtRes_EndcapBarrel;
    TH1F* LeptonPtRes_EndcapEndcap;
    TH1F* LeptonInvPtRes_BarrelBarrel;
    TH1F* LeptonInvPtRes_EndcapBarrel;
    TH1F* LeptonInvPtRes_EndcapEndcap;
    TH1F* LeptonInvPtPull_BarrelBarrel;
    TH1F* LeptonInvPtPull_EndcapBarrel;
    TH1F* LeptonInvPtPull_EndcapEndcap;
    TH1F* LeptonPtPull_BarrelBarrel;
    TH1F* LeptonPtPull_EndcapBarrel;
    TH1F* LeptonPtPull_EndcapEndcap;
    TH2F* DileptonMassResVSvtxChi2;
    TH2F* DileptonMassVtxResVSvtxChi2;
    TH2F* LeptonPtResVSvtxChi2;
    TH2F* LeptonInvPtResVSvtxChi2;
    TH2F* DileptonMassResVSMuonPtProb;
    TH2F* DileptonMassResVSMinMuonPtProb;
    TH2F* DileptonMassVtxResVSMinMuonPtProb;
    TH2F* LeptonPtResVSMuonPtProb;
    TH2F* LeptonInvPtResVSMuonPtProb;
    TH2F* vtxChi2VSMuonPtProb;
    TH2F* vtxChi2VSMinMuonPtProb;
    TH2F* LeptonRecoGenDxyVSvtxChi2;
    TH2F* LeptonRecoGenDzVSvtxChi2;
    TH2F* LeptonGenDxy2VSvtxChi2;
    TH2F* LeptonGenDz2VSvtxChi2;
    TH2F* vtxChi2VSMuonEta;
    TH2F* vtxChi2VSMuonMDist;
    
};

ResolutionVertexUsingMC::ResolutionVertexUsingMC(const edm::ParameterSet& cfg)
  : hardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")),
    lepton_src(cfg.getParameter<edm::InputTag>("lepton_src")),
    dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    leptonsFromDileptons(cfg.getParameter<bool>("leptonsFromDileptons")),
    doQoverP(cfg.getParameter<bool>("doQoverP"))
{
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

  LeptonInvPtResVPtGen = fs->make<TProfile>("LeptonInvPtResVPtGen", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) vs. gen pT", 50, 0, 1000, -0.5, 0.5);
  LeptonInvPResVPGen   = fs->make<TProfile>("LeptonInvPResVPGen",   titlePrefix + "(1/p - 1/gen p)/(1/gen p) vs. gen p",     50, 0, 1000, -0.5, 0.5);
  
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
    
    
  DileptonMassReco            = fs->make<TH1F>("DileptonMassReco",            titlePrefix + "dil. mass", 20000, 0, 20000);
  DileptonMassGen            = fs->make<TH1F>("DileptonMassGen",            titlePrefix + "dil. mass gen", 20000, 0, 20000);
  DileptonMassRes    = fs->make<TH1F>("DileptonMassRes",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 100, -0.5, 0.5);
  DileptonResMassRes = fs->make<TH1F>("DileptonResMassRes", titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass)", 100, -0.5, 0.5);
  ResonanceMassRes   = fs->make<TH1F>("ResonanceMassRes",   titlePrefix + "(res. mass - gen res. mass)/(gen res. mass)", 100, -0.5, 0.5);
  DileptonMassInvRes    = fs->make<TH1F>("DileptonMassInvRes",    titlePrefix + "(1/dil. mass - 1/gen dil. mass)/(1/gen dil. mass)", 100, -0.5, 0.5);
  
  DileptonMassResVMass    = fs->make<TProfile>("DileptonMassResVMass",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", nbinsmass,0, massmax, -1, 1);
  DileptonResMassResVMass = fs->make<TProfile>("DileptonResMassResVMass", titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass)", nbinsmass,0, massmax, -1, 1);
  ResonanceMassResVMass   = fs->make<TProfile>("ResonanceMassResVMass",   titlePrefix + "(res. mass - gen res. mass)/(gen res. mass)", nbinsmass,0, massmax, -1, 1);
  DileptonInvMassResVMass    = fs->make<TProfile>("DileptonInvMassResVMass",    titlePrefix + "(./dil. mass - 1/gen dil. mass)/(1/gen dil. mass)", nbinsmass,0, massmax, -1, 1);

  static const TString dilepton_where_names[W_D_MAX] = {"BB", "BO", "BE", "BU", "OO", "OE", "OU", "EE", "EU", "UU"};
  for (size_t i = 0; i < W_D_MAX; ++i) {
    DileptonMassResBy   [i] = fs->make<TH1F>("DileptonMassRes"    + dilepton_where_names[i], titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
    DileptonResMassResBy[i] = fs->make<TH1F>("DileptonResMassRes" + dilepton_where_names[i], titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
    ResonanceMassResBy  [i] = fs->make<TH1F>("ResonanceMassRes"   + dilepton_where_names[i], titlePrefix + "(res. mass - gen res. mass)/(gen res. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
  }
    ///////////////
    DileptonMassRes_EndcapEndcap    = fs->make<TH1F>("DileptonMassRes_EndcapEndcap",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 100, -0.5, 0.5);
    DileptonMassRes_EndcapBarrel    = fs->make<TH1F>("DileptonMassRes_EndcapBarrel",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 100, -0.5, 0.5);
    DileptonMassRes_BarrelBarrel    = fs->make<TH1F>("DileptonMassRes_BarrelBarrel",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 100, -0.5, 0.5);

   LeptonPt_BarrelBarrel = fs->make<TH1F>("LeptonPt_BarrelBarrel", titlePrefix + "pT BB", 5000, 0, 5000);
   LeptonPt_EndcapBarrel = fs->make<TH1F>("LeptonPt_EndcapBarrel", titlePrefix + "pT EB", 5000, 0, 5000);
   LeptonPt_EndcapEndcap = fs->make<TH1F>("LeptonPt_EndcapEndcap", titlePrefix + "pT EE", 5000, 0, 5000);
    
   LeptonPtRes_BarrelBarrel = fs->make<TH1F>("LeptonPtRes_BarrelBarrel", titlePrefix + "(pT - gen pT)/(gen pT) BB", 100, -0.5, 0.5);
   LeptonPtRes_EndcapBarrel = fs->make<TH1F>("LeptonPtRes_EndcapBarrel", titlePrefix + "(pT - gen pT)/(gen pT) EB", 100, -0.5, 0.5);
   LeptonPtRes_EndcapEndcap = fs->make<TH1F>("LeptonPtRes_EndcapEndcap", titlePrefix + "(pT - gen pT)/(gen pT) EE", 100, -0.5, 0.5);
   
    LeptonInvPtRes_BarrelBarrel = fs->make<TH1F>("LeptonInvPtRes_BarrelBarrel", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) BB", 100, -0.5, 0.5);
    LeptonInvPtRes_EndcapBarrel = fs->make<TH1F>("LeptonInvPtRes_EndcapBarrel", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) EB", 100, -0.5, 0.5);
    LeptonInvPtRes_EndcapEndcap = fs->make<TH1F>("LeptonInvPtRes_EndcapEndcap", titlePrefix + "(1/pT - 1/gen pT)/(1/gen pT) EE", 100, -0.5, 0.5);
    
    LeptonInvPtPull_BarrelBarrel = fs->make<TH1F>("LeptonInvPtPull_BarrelBarrel", titlePrefix + "(1/pT - 1/gen pT)/#sigma_{1/pT} BB", 100, -10, 10);
    LeptonInvPtPull_EndcapBarrel = fs->make<TH1F>("LeptonInvPtPull_EndcapBarrel", titlePrefix + "(1/pT - 1/gen pT)/#sigma_{1/pT} EB", 100, -10, 10);
    LeptonInvPtPull_EndcapEndcap = fs->make<TH1F>("LeptonInvPtPull_EndcapEndcap", titlePrefix + "(1/pT - 1/gen pT)/#sigma_{1/pT} EE", 100, -10, 10);
    
    LeptonPtPull_BarrelBarrel = fs->make<TH1F>("LeptonPtPull_BarrelBarrel", titlePrefix + "(pT - gen pT)/#sigma_{pT} BB", 100, -10, 10);
    LeptonPtPull_EndcapBarrel = fs->make<TH1F>("LeptonPtPull_EndcapBarrel", titlePrefix + "(pT - gen pT)/#sigma_{pT} EB", 100, -10, 10);
    LeptonPtPull_EndcapEndcap = fs->make<TH1F>("LeptonPtPull_EndcapEndcap", titlePrefix + "(pT - gen pT)/#sigma_{pT} EE", 100, -10, 10);
    
    DileptonMassResVSvtxChi2= fs->make<TH2F>("DileptonMassResVSvtxChi2",    titlePrefix + "dimu. vertex #chi^{2}/dof VS (dil. mass - gen dil. mass)/(gen dil. mass)", 300, 0, 30, 100, -0.5, 0.5);
    DileptonMassVtxResVSvtxChi2= fs->make<TH2F>("DileptonMassVtxResVSvtxChi2",    titlePrefix + "dimu. vertex #chi^{2}/dof VS (dil. mass - gen dil. mass)/(gen dil. mass)", 300, 0, 30, 100, -0.5, 0.5);
    LeptonPtResVSvtxChi2= fs->make<TH2F>("LeptonPtResVSvtxChi2",    titlePrefix + "dimu. vertex #chi^{2}/dof VS (pT - gen pT)/(gen pT)", 300, 0, 30, 100, -0.5, 0.5);
    LeptonInvPtResVSvtxChi2= fs->make<TH2F>("LeptonInvPtResVSvtxChi2",    titlePrefix + "dimu. vertex #chi^{2}/dof VS (1/pT - 1/gen pT)/#sigma_{1/pT}", 300, 0, 30, 100, -0.5, 0.5);
    
    DileptonMassResVSMuonPtProb= fs->make<TH2F>("DileptonMassResVSMuonPtProb",    titlePrefix + "muon track prob. VS (dil. mass - gen dil. mass)/(gen dil. mass)", 1000, 0., 1., 100, -0.5, 0.5);
    DileptonMassResVSMinMuonPtProb= fs->make<TH2F>("DileptonMassResVSMinMuonPtProb",    titlePrefix + "min muon track prob. VS (dil. mass - gen dil. mass)/(gen dil. mass)", 1000, 0., 1., 100, -0.5, 0.5);
    DileptonMassVtxResVSMinMuonPtProb= fs->make<TH2F>("DileptonMassVtxResVSMinMuonPtProb",    titlePrefix + "min muon track prob. VS (dil. mass - gen dil. mass)/(gen dil. mass)", 1000, 0., 1., 100, -0.5, 0.5);
    LeptonPtResVSMuonPtProb= fs->make<TH2F>("LeptonPtResVSMuonPtProb",    titlePrefix + "muon track prob. VS (pT - gen pT)/(gen pT)", 1000, 0., 1., 100, -0.5, 0.5);
    LeptonInvPtResVSMuonPtProb= fs->make<TH2F>("LeptonInvPtResVSMuonPtProb",    titlePrefix + "muon track prob. VS (1/pT - 1/gen pT)/#sigma_{1/pT}", 1000, 0., 1., 100, -0.5, 0.5);
    vtxChi2VSMuonPtProb= fs->make<TH2F>("vtxChi2VSMuonPtProb",    titlePrefix + "muon track prob. VS dimu. vertex #chi^{2}/dof", 1000, 0., 1., 300, 0, 30);
    vtxChi2VSMinMuonPtProb= fs->make<TH2F>("vtxChi2VSMinMuonPtProb",    titlePrefix + "min muon track prob. VS dimu. vertex #chi^{2}/dof", 1000, 0., 1., 300, 0, 30);
    vtxChi2VSMuonEta= fs->make<TH2F>("vtxChi2VSMuonEta",    titlePrefix + "max muon track eta VS dimu. vertex #chi^{2}/dof", 100, -5, 5, 600, 0, 60);
    vtxChi2VSMuonMDist= fs->make<TH2F>("vtxChi2VSMuonMDist",    titlePrefix + "eta (muon max distance from gen vertex) VS dimu. vertex #chi^{2}/dof", 100, -5, 5, 600, 0, 60);
    
    LeptonRecoGenDxyVSvtxChi2= fs->make<TH2F>("LeptonRecoGenDxyVSvtxChi2",    titlePrefix + "dimu. vertex #chi^{2}/dof VS |dxy wrt gen vtx|", 300, 0, 30, 1000, 0, 0.2);
    LeptonRecoGenDzVSvtxChi2= fs->make<TH2F>("LeptonRecoGenDzVSvtxChi2",    titlePrefix + "dimu. vertex #chi^{2}/dof VS |dz wrt gen vtx|", 300, 0, 30, 1000, 0, 0.2);
    LeptonGenDxy2VSvtxChi2= fs->make<TH2F>("LeptonGenDxy2VSvtxChi2",    titlePrefix + "dimu. vertex #chi^{2}/dof VS #Delta(gen_vx)^2 + #Delta(gen_vy)^2", 300, 0, 30, 2000, -0.2, 0.2);
    LeptonGenDz2VSvtxChi2= fs->make<TH2F>("LeptonGenDz2VSvtxChi2",    titlePrefix + "dimu. vertex #chi^{2}/dof VS #Delta(gen_vz)^2", 300, 0, 30, 2000, -0.2, 0.2);
}

void ResolutionVertexUsingMC::fillLeptonResolution(const reco::GenParticle* gen_lep, const reco::CandidateBaseRef& lep) {
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

void ResolutionVertexUsingMC::fillLeptonExtraMomentumResolution(const reco::GenParticle* gen_lep, const reco::CandidateBaseRef& lep) {
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

  const double inv_pt_diff = inv_pt - inv_gen_pt;
  const double inv_p_diff  = inv_p  - inv_gen_p;

  const double inv_pt_res = inv_pt_diff/inv_gen_pt;
  const double inv_p_res  = inv_p_diff /inv_gen_p;

  // Inverse momentum resolutions as a function of generated momenta.
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

void ResolutionVertexUsingMC::fillChargeResolution(const reco::GenParticle* gen_lep, const reco::CandidateBaseRef& lep) {
  const int delta_q = lep->charge() - gen_lep->charge();
  ChargeDiff->Fill(delta_q);
  if (delta_q == 0)
    ChargeRightVInvPt->Fill(1/gen_lep->pt());
  else //if (delta_q == 2), also don't throw exception?
    ChargeWrongVInvPt->Fill(1/gen_lep->pt());
}

void ResolutionVertexUsingMC::fillDileptonMassResolution(const pat::CompositeCandidate& dil) {
  if (!hardInteraction.IsValid())
    return;

  const double mass         = dil.mass();
  const double gen_mass     = (hardInteraction.lepPlus->p4() + hardInteraction.lepMinus->p4()).mass();
  const double mass_vtx         =dil.userFloat("vertexM");

  const double res_mass     = resonanceP4(dil).mass();
  const double gen_res_mass = hardInteraction.resonance->mass();

  const double rdil    = mass    /gen_mass     - 1;
  const double rdil_vtx    = mass_vtx    /gen_mass     - 1;
  const double rdilres = mass    /gen_res_mass - 1;
  const double rres    = res_mass/gen_res_mass - 1;
  double invmres = (1./mass-1./gen_mass)/(1./gen_mass);
 
  DileptonMassReco ->Fill(mass);
  DileptonMassGen->Fill(gen_mass);
  DileptonMassRes   ->Fill(rdil);
  DileptonResMassRes->Fill(rdilres);
  ResonanceMassRes  ->Fill(rres);
  DileptonMassInvRes   ->Fill(invmres);

  DileptonMassResVMass   ->Fill(gen_mass,     rdil*rdil);
  DileptonResMassResVMass->Fill(gen_res_mass, rdilres*rdilres);
  ResonanceMassResVMass  ->Fill(gen_res_mass, rres*rres);
  DileptonInvMassResVMass   ->Fill(gen_mass,  invmres*invmres);

  size_t w = whereIsDilepton(dil);
  DileptonMassResBy   [w]->Fill(rdil);
  DileptonResMassResBy[w]->Fill(rdilres);
  ResonanceMassResBy  [w]->Fill(rres);
  
    
  //////////////////////////////////
  //float vertex_chi2 = dil.userFloat("vertex_chi2");
  const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
  const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);
    if (lep0.isNonnull() && lep1.isNonnull()) {
        const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
        const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
        if (mu0 && mu1) {
            const reco::GenParticle* gen_lep0= mu0->genParticle();
            const reco::GenParticle* gen_lep1= mu1->genParticle();
            const reco::Track* tk0 = patmuon::getPickedTrack(*mu0).get();
            const reco::Track* tk1 = patmuon::getPickedTrack(*mu1).get();
            //const reco::Track* gen_tk0 = gen_lep0->bestTrack(); //BOH
            //std::cout<<gen_tk0<<std::endl;
            //if(gen_tk0) std::cout<<"STI "<<tk0->d0()<<tk0->d0Error()<<gen_tk0->d0()<<gen_tk0->d0Error()<<std::endl;
            
            //std::cout<<dil.userFloat("vertexM")<<std::endl;
            if (tk0 && tk1 && (gen_lep0 != 0) && (gen_lep1 != 0) ) {
            //if(gen_mass<=200 && gen_mass>=100){//solo x dy120
                float vertex_chi2 = dil.userFloat("vertex_chi2");
                DileptonMassResVSvtxChi2->Fill(vertex_chi2,rdil);
                DileptonMassVtxResVSvtxChi2->Fill(vertex_chi2,rdil_vtx);
                LeptonPtResVSvtxChi2->Fill(vertex_chi2,(lep0->pt() - gen_lep0->pt())/gen_lep0->pt());
                LeptonPtResVSvtxChi2->Fill(vertex_chi2,(lep1->pt() - gen_lep1->pt())/gen_lep1->pt());
                LeptonInvPtResVSvtxChi2->Fill(vertex_chi2,(1/lep0->pt() - 1/gen_lep0->pt())/(1/gen_lep0->pt()));
                LeptonInvPtResVSvtxChi2->Fill(vertex_chi2,(1/lep1->pt() - 1/gen_lep1->pt())/(1/gen_lep1->pt()));
                
                DileptonMassResVSMuonPtProb->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()),rdil);
                DileptonMassResVSMuonPtProb->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()),rdil);
                    if (TMath::Prob(tk0->chi2(), tk0->ndof()) < TMath::Prob(tk1->chi2(), tk1->ndof())) DileptonMassResVSMinMuonPtProb->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()),rdil);
                    else DileptonMassResVSMinMuonPtProb->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()),rdil);
                
                    if (TMath::Prob(tk0->chi2(), tk0->ndof()) < TMath::Prob(tk1->chi2(), tk1->ndof())) DileptonMassVtxResVSMinMuonPtProb->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()),rdil_vtx);
                    else DileptonMassVtxResVSMinMuonPtProb->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()),rdil_vtx);
                
                LeptonPtResVSMuonPtProb->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()),(lep0->pt() - gen_lep0->pt())/gen_lep0->pt());
                LeptonPtResVSMuonPtProb->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()),(lep1->pt() - gen_lep1->pt())/gen_lep1->pt());
                LeptonInvPtResVSMuonPtProb->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()),(1/lep0->pt() - 1/gen_lep0->pt())/(1/gen_lep0->pt()));
                LeptonInvPtResVSMuonPtProb->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()),(1/lep1->pt() - 1/gen_lep1->pt())/(1/gen_lep1->pt()));
                vtxChi2VSMuonPtProb->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()),vertex_chi2);
                vtxChi2VSMuonPtProb->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()),vertex_chi2);
                    if (TMath::Prob(tk0->chi2(), tk0->ndof()) < TMath::Prob(tk1->chi2(), tk1->ndof())) vtxChi2VSMinMuonPtProb->Fill(TMath::Prob(tk0->chi2(), tk0->ndof()),vertex_chi2);
                    else vtxChi2VSMinMuonPtProb->Fill(TMath::Prob(tk1->chi2(), tk1->ndof()),vertex_chi2);
                
                if(fabs(tk1->eta())>fabs(tk0->eta())) vtxChi2VSMuonEta->Fill(fabs(tk1->eta()),vertex_chi2);
                else vtxChi2VSMuonEta->Fill(fabs(tk0->eta()),vertex_chi2);
                //std::cout<<"gen_lep0 eta pt "<<gen_lep0->eta()<<"  "<<gen_lep0->pt()<<" lep eta pt "<<lep0->eta()<<"  "<<lep0->pt()<<std::endl;
                //std::cout<<"traccia eta pt "<<gen_lep1->eta()<<"  "<<gen_lep1->pt()<<" lep eta pt "<<lep1->eta()<<"  "<<lep1->pt()<<std::endl;
                //std::cout<<"hardInteraction.lepPlus eta pt "<<hardInteraction.lepPlus->eta()<<"  "<<hardInteraction.lepPlus->pt()<<" hardInteraction.lepMinus eta pt "<<hardInteraction.lepMinus->eta()<<"  "<<hardInteraction.lepMinus->pt()<<std::endl;
                //if((hardInteraction.lepMinus->eta() != gen_lep1->eta() && hardInteraction.lepMinus->eta() != gen_lep0->eta()) || (hardInteraction.lepPlus->eta() != gen_lep1->eta() && hardInteraction.lepPlus->eta() != gen_lep0->eta())) std::cout<<" CASPUR "<<std::endl;

                std::cout<<gen_lep0->vx()<<" "<<gen_lep1->vx()<<std::endl;
                LeptonRecoGenDxyVSvtxChi2->Fill(vertex_chi2,fabs(tk0->dxy(gen_lep0->vertex())));
                LeptonRecoGenDxyVSvtxChi2->Fill(vertex_chi2,fabs(tk1->dxy(gen_lep1->vertex())));
                LeptonRecoGenDzVSvtxChi2->Fill(vertex_chi2,fabs(tk0->dz(gen_lep0->vertex())));
                LeptonRecoGenDzVSvtxChi2->Fill(vertex_chi2,fabs(tk1->dz(gen_lep1->vertex())));
                
                LeptonGenDz2VSvtxChi2->Fill(vertex_chi2,(gen_lep0->vz()-gen_lep1->vz())*(gen_lep0->vz()-gen_lep1->vz()));
                LeptonGenDxy2VSvtxChi2->Fill(vertex_chi2,(gen_lep0->vx()-gen_lep1->vx())*(gen_lep0->vx()-gen_lep1->vx()) + (gen_lep0->vy()-gen_lep1->vy())*(gen_lep0->vy()-gen_lep1->vy()));

                if((fabs(tk0->dxy(gen_lep0->vertex()))+fabs(tk0->dz(gen_lep0->vertex()))) < fabs(tk1->dxy(gen_lep1->vertex()))+fabs(tk1->dz(gen_lep1->vertex()))) vtxChi2VSMuonMDist->Fill(fabs(tk1->eta()),vertex_chi2);
                else vtxChi2VSMuonMDist->Fill(fabs(tk0->eta()),vertex_chi2);
                
                //ENDCAPBARREL
                if((fabs(tk1->eta())>1. && fabs(tk0->eta())<0.8) || (fabs(tk0->eta())>1. && fabs(tk1->eta())<0.8)){
                    DileptonMassRes_EndcapBarrel->Fill(rdil);
                    
                    LeptonPt_EndcapBarrel->Fill(lep0->pt());
                    LeptonPt_EndcapBarrel->Fill(lep1->pt());
                    
                    LeptonPtRes_EndcapBarrel->Fill((lep0->pt() - gen_lep0->pt())/gen_lep0->pt());
                    LeptonPtRes_EndcapBarrel->Fill((lep1->pt() - gen_lep1->pt())/gen_lep1->pt());

                    LeptonInvPtRes_EndcapBarrel->Fill((1/lep0->pt() - 1/gen_lep0->pt())/(1/gen_lep0->pt()));
                    LeptonInvPtRes_EndcapBarrel->Fill((1/lep1->pt() - 1/gen_lep1->pt())/(1/gen_lep1->pt()));
                    
                    LeptonInvPtPull_EndcapBarrel->Fill((1/lep0->pt() - 1/gen_lep0->pt())/invPtError(&*tk0));
                    LeptonInvPtPull_EndcapBarrel->Fill((1/lep1->pt() - 1/gen_lep1->pt())/invPtError(&*tk1));
                    
                    LeptonPtPull_EndcapBarrel->Fill((lep0->pt() - gen_lep0->pt())/ptError(&*tk0));
                    LeptonPtPull_EndcapBarrel->Fill((lep1->pt() - gen_lep1->pt())/ptError(&*tk1));

                }
                if(fabs(tk1->eta())<0.8 && fabs(tk0->eta())<0.8){
                    DileptonMassRes_BarrelBarrel->Fill(rdil);
                    
                    LeptonPt_BarrelBarrel->Fill(lep0->pt());
                    LeptonPt_BarrelBarrel->Fill(lep1->pt());
                    
                    LeptonPtRes_BarrelBarrel->Fill((lep0->pt() - gen_lep0->pt())/gen_lep0->pt());
                    LeptonPtRes_BarrelBarrel->Fill((lep1->pt() - gen_lep1->pt())/gen_lep1->pt());
                    
                    LeptonInvPtRes_BarrelBarrel->Fill((1/lep0->pt() - 1/gen_lep0->pt())/(1/gen_lep0->pt()));
                    LeptonInvPtRes_BarrelBarrel->Fill((1/lep1->pt() - 1/gen_lep1->pt())/(1/gen_lep1->pt()));
                    
                    LeptonInvPtPull_BarrelBarrel->Fill((1/lep0->pt() - 1/gen_lep0->pt())/invPtError(&*tk0));
                    LeptonInvPtPull_BarrelBarrel->Fill((1/lep1->pt() - 1/gen_lep1->pt())/invPtError(&*tk1));
                    
                    LeptonPtPull_BarrelBarrel->Fill((lep0->pt() - gen_lep0->pt())/ptError(&*tk0));
                    LeptonPtPull_BarrelBarrel->Fill((lep1->pt() - gen_lep1->pt())/ptError(&*tk1));

                }
                if(fabs(tk1->eta())>1. && fabs(tk0->eta())>1.){
                    DileptonMassRes_EndcapEndcap->Fill(rdil);
                    
                    LeptonPt_EndcapEndcap->Fill(lep0->pt());
                    LeptonPt_EndcapEndcap->Fill(lep1->pt());
                    
                    LeptonPtRes_EndcapEndcap->Fill((lep0->pt() - gen_lep0->pt())/gen_lep0->pt());
                    LeptonPtRes_EndcapEndcap->Fill((lep1->pt() - gen_lep1->pt())/gen_lep1->pt());
                    
                    LeptonInvPtRes_EndcapEndcap->Fill((1/lep0->pt() - 1/gen_lep0->pt())/(1/gen_lep0->pt()));
                    LeptonInvPtRes_EndcapEndcap->Fill((1/lep1->pt() - 1/gen_lep1->pt())/(1/gen_lep1->pt()));
                    
                    LeptonInvPtPull_EndcapEndcap->Fill((1/lep0->pt() - 1/gen_lep0->pt())/invPtError(&*tk0));
                    LeptonInvPtPull_EndcapEndcap->Fill((1/lep1->pt() - 1/gen_lep1->pt())/invPtError(&*tk1));
                    
                    LeptonPtPull_EndcapEndcap->Fill((lep0->pt() - gen_lep0->pt())/ptError(&*tk0));
                    LeptonPtPull_EndcapEndcap->Fill((lep1->pt() - gen_lep1->pt())/ptError(&*tk1));
                    
                    
                }
            }//}//solo per dy120
        }
    }

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

void ResolutionVertexUsingMC::fillLeptonHistos(const reco::CandidateBaseRef& lep) {
  const reco::GenParticle* gen_lep = getGenParticle(lep);
  if (gen_lep != 0) {
    fillLeptonResolution(gen_lep, lep);
    fillLeptonExtraMomentumResolution(gen_lep, lep);
    fillChargeResolution(gen_lep, lep);
  }
}

void ResolutionVertexUsingMC::fillLeptonHistos(const edm::View<reco::Candidate>& leptons) {
  // Fill lepton histos from all leptons.
  for (size_t i = 0, n = leptons.size(); i < n; ++i)
    fillLeptonHistos(leptons.refAt(i));
}

void ResolutionVertexUsingMC::fillLeptonHistosFromDileptons(const pat::CompositeCandidateCollection& dileptons) {
  for (pat::CompositeCandidateCollection::const_iterator dil = dileptons.begin(), dile = dileptons.end(); dil != dile; ++dil)
    for (size_t i = 0; i < dil->numberOfDaughters(); ++i)
      fillLeptonHistos(dil->daughter(i)->masterClone());
}

void ResolutionVertexUsingMC::fillDileptonHistos(const pat::CompositeCandidateCollection& dileptons) {
  for (pat::CompositeCandidateCollection::const_iterator dil = dileptons.begin(), dile = dileptons.end(); dil != dile; ++dil)
    fillDileptonMassResolution(*dil);
}

void ResolutionVertexUsingMC::analyze(const edm::Event& event, const edm::EventSetup& setup) {
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

DEFINE_FWK_MODULE(ResolutionVertexUsingMC);
