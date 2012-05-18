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

  TH1F* DileptonMassRes;
  TH1F* DileptonResMassRes;
  TH1F* ResonanceMassRes;

  TProfile* DileptonMassResVMass;
  TProfile* DileptonResMassResVMass;
  TProfile* ResonanceMassResVMass;

  TH1F* DileptonMassResBy[W_D_MAX];
  TH1F* DileptonResMassResBy[W_D_MAX];
  TH1F* ResonanceMassResBy[W_D_MAX];
};

ResolutionUsingMC::ResolutionUsingMC(const edm::ParameterSet& cfg)
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

  DileptonMassRes    = fs->make<TH1F>("DileptonMassRes",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 100, -0.5, 0.5);
  DileptonResMassRes = fs->make<TH1F>("DileptonResMassRes", titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass)", 100, -0.5, 0.5);
  ResonanceMassRes   = fs->make<TH1F>("ResonanceMassRes",   titlePrefix + "(res. mass - gen res. mass)/(gen res. mass)", 100, -0.5, 0.5);
  
  DileptonMassResVMass    = fs->make<TProfile>("DileptonMassResVMass",    titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass)", 100, nbinsmass, massmax, -0.5, 0.5);
  DileptonResMassResVMass = fs->make<TProfile>("DileptonResMassResVMass", titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass)", 100, nbinsmass, massmax, -0.5, 0.5);
  ResonanceMassResVMass   = fs->make<TProfile>("ResonanceMassResVMass",   titlePrefix + "(res. mass - gen res. mass)/(gen res. mass)", 100, nbinsmass, massmax, -0.5, 0.5);

  static const TString dilepton_where_names[W_D_MAX] = {"BB", "BO", "BE", "BU", "OO", "OE", "OU", "EE", "EU", "UU"};
  for (size_t i = 0; i < W_D_MAX; ++i) {
    DileptonMassResBy   [i] = fs->make<TH1F>("DileptonMassRes"    + dilepton_where_names[i], titlePrefix + "(dil. mass - gen dil. mass)/(gen dil. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
    DileptonResMassResBy[i] = fs->make<TH1F>("DileptonResMassRes" + dilepton_where_names[i], titlePrefix + "(dil. mass - gen res. mass)/(gen res. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
    ResonanceMassResBy  [i] = fs->make<TH1F>("ResonanceMassRes"   + dilepton_where_names[i], titlePrefix + "(res. mass - gen res. mass)/(gen res. mass), " + dilepton_where_names[i], 100, -0.5, 0.5);
  }
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

void ResolutionUsingMC::fillChargeResolution(const reco::GenParticle* gen_lep, const reco::CandidateBaseRef& lep) {
  const int delta_q = lep->charge() - gen_lep->charge();
  ChargeDiff->Fill(delta_q);
  if (delta_q == 0)
    ChargeRightVInvPt->Fill(1/gen_lep->pt());
  else //if (delta_q == 2), also don't throw exception?
    ChargeWrongVInvPt->Fill(1/gen_lep->pt());
}

void ResolutionUsingMC::fillDileptonMassResolution(const reco::CompositeCandidate& dil) {
  if (!hardInteraction.IsValid())
    return;

  const double mass         = dil.mass();
  const double gen_mass     = (hardInteraction.lepPlus->p4() + hardInteraction.lepMinus->p4()).mass();

  const double res_mass     = resonanceP4(dil).mass();
  const double gen_res_mass = hardInteraction.resonance->mass();

  const double rdil    = mass    /gen_mass     - 1;
  const double rdilres = mass    /gen_res_mass - 1;
  const double rres    = res_mass/gen_res_mass - 1;
  
  DileptonMassRes   ->Fill(rdil);
  DileptonResMassRes->Fill(rdilres);
  ResonanceMassRes  ->Fill(rres);

  DileptonMassResVMass   ->Fill(gen_mass,     rdil*rdil);
  DileptonResMassResVMass->Fill(gen_res_mass, rdilres*rdilres);
  ResonanceMassResVMass  ->Fill(gen_res_mass, rres*rres);

  size_t w = whereIsDilepton(dil);
  DileptonMassResBy   [w]->Fill(rdil);
  DileptonResMassResBy[w]->Fill(rdilres);
  ResonanceMassResBy  [w]->Fill(rres);
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
