/**
  \class    Zprime2muResolution
  \brief    Calculates and plots lepton/dilepton resolutions and efficiencies.

  \author   Jordan Tucker, Slava Valuev
  \version  $Id: Zprime2muResolution.cc,v 1.36 2008/12/11 20:53:59 tucker Exp $
*/

#include "TString.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muResolution.h"

using namespace std;

Zprime2muResolution::Zprime2muResolution(const edm::ParameterSet& config) : Zprime2muAnalysis(config) {
  leptonsFromDileptons = config.getParameter<bool>("leptonsFromDileptons");

  // Get the parameters specific to the data sample on which we are running.
  string dataSet = config.getParameter<string>("dataSet");
  edm::ParameterSet dataSetConfig = config.getParameter<edm::ParameterSet>(dataSet);

  peakMass     = dataSetConfig.getParameter<double>("peakMass");
  lowerMassWin = dataSetConfig.getParameter<double>("lowerMassWin");
  upperMassWin = dataSetConfig.getParameter<double>("upperMassWin");
  binSize      = dataSetConfig.getParameter<int>   ("binSize");
  maxTrigMass  = dataSetConfig.getParameter<double>("maxTrigMass");

  bookGenLevelHistos();
  bookLeptonResolutionHistos();
  bookChargeResolutionHistos();
  bookEfficiencyHistos();
  bookDileptonResolutionHistos();
}

void Zprime2muResolution::bookGenLevelHistos() {
  const int nx = 20;
  LeptonOrigin[0] = fs->make<TH1F>("LeptonOrigin0", "Particle Id of Mother of all leptons",                 nx, 0, nx);
  LeptonOrigin[1] = fs->make<TH1F>("LeptonOrigin1", "Particle Id of Mother of opp-sign dilepton daughters", nx, 0, nx);
  char *mother[nx] = {"  ","pi","K ","K0","eta","rho","c ","b ","tau","  ",
		      "Z ","W ","H ","Z'","G*","  ","  ","  ","  ","  "};
  for (int i = 0; i < nx; i++) {
    LeptonOrigin[0]->GetXaxis()->SetBinLabel(i+1, mother[i]);
    LeptonOrigin[1]->GetXaxis()->SetBinLabel(i+1, mother[i]);
  }

  GenMassAllEvents = fs->make<TH1F>("GenMassAllEvents", "Gen mass, all events",           24, lowerMassWin, upperMassWin);
  GenMassInAccept  = fs->make<TH1F>("GenMassInAccept",  "Gen mass, events in acceptance", 24, lowerMassWin, upperMassWin);

  GenMassAllEvents->Sumw2();
  GenMassInAccept->Sumw2();
}

void Zprime2muResolution::bookLeptonResolutionHistos() {
  double scale;
  for (int rec = lL1; rec < MAX_LEVELS; ++rec) {
    TString level(levelName(rec));
    
    scale = rec == lL1 ? 0.1 : (rec == lL2 ? 0.02 : 0.002);
    LeptonEtaDiff[rec] = fs->make<TH1F>(nameHist("LeptonEtaDiff", rec), level + " #eta - gen #eta", 100, -scale, scale);
    LeptonPhiDiff[rec] = fs->make<TH1F>(nameHist("LeptonPhiDiff", rec), level + " #phi - gen #phi", 100, -scale, scale);

    scale = 0.1*peakMass;
    scale = rec == lL1 ? 700 + 2*scale : (rec == lL2 ? 4*scale : scale);
    LeptonPtDiff[rec] = fs->make<TH1F>(nameHist("LeptonPtDiff", rec), level + " pT - gen pT", 100, -scale, scale);

    scale = rec == lL1 ? 0.6*peakMass/140 : (rec == lL2 ? 2 : 0.3);
    LeptonPtRes[rec] = fs->make<TH1F>(nameHist("LeptonPtRes", rec), level + " (pT - gen pT)/(gen pT)", 100, -scale, scale);
    LeptonPRes[rec]  = fs->make<TH1F>(nameHist("LeptonPRes",  rec), level + " (p - gen p)/(gen p)",    100, -scale, scale);

    LeptonInvPtRes[rec] = fs->make<TH1F>(nameHist("LeptonInvPtRes", rec), level + " (1/pT - 1/gen pT)/(1/gen pT)", 100, -scale, scale);
    LeptonInvPRes[rec]  = fs->make<TH1F>(nameHist("LeptonInvPRes",  rec), level + " (1/p - 1/gen p)/(1/gen p)",    100, -scale, scale);

    if (rec >= lGR) {
      LeptonInvPtResVPtGen[rec] = fs->make<TProfile>(nameHist("LeptonInvPtResVPtGen", rec), level + " (1/pT - 1/gen pT)/(1/gen pT) vs. gen pT", 50, 0, peakMass, -scale, scale);
      LeptonInvPResVPGen[rec]   = fs->make<TProfile>(nameHist("LeptonInvPResVPGen",   rec), level + " (1/p - 1/gen p)/(1/gen p) vs. gen p",     50, 0, peakMass, -scale, scale);
  
      LeptonInvPtPull[rec] = fs->make<TH1F>(nameHist("LeptonInvPtPull", rec), level + " (1/pT - 1/gen pT)/#sigma_{1/pT}", 100, -10, 10);
      LeptonInvPPull[rec]  = fs->make<TH1F>(nameHist("LeptonInvPPull",  rec), level + " (1/p - 1/gen p)/#sigma_{1/p}",    100, -10, 10);
  
      LeptonInvPtResBarrel[rec] = fs->make<TH1F>(nameHist("LeptonInvPtResBarrel", rec), level + " (1/pT - 1/gen pT)/(1/gen pT), barrel", 100, -scale, scale);
      LeptonInvPResBarrel[rec]  = fs->make<TH1F>(nameHist("LeptonInvPResBarrel",  rec), level + " (1/p - 1/gen p)/(1/gen p), barrel",    100, -scale, scale);
    
      LeptonInvPtPullBarrel[rec] = fs->make<TH1F>(nameHist("LeptonInvPtPullBarrel", rec), level + " (1/pT - 1/gen pT)/#sigma_{1/pT}, barrel", 100, -10, 10);
      LeptonInvPPullBarrel[rec]  = fs->make<TH1F>(nameHist("LeptonInvPPullBarrel",  rec), level + " (1/p - 1/gen p)/#sigma_{1/p}, barrel",    100, -10, 10);

      LeptonInvPtResEndcap[rec] = fs->make<TH1F>(nameHist("LeptonInvPtResEndcap", rec), level + " (1/pT - 1/gen pT)/(1/gen pT), endcap", 100, -scale, scale);
      LeptonInvPResEndcap[rec]  = fs->make<TH1F>(nameHist("LeptonInvPResEndcap",  rec), level + " (1/p - 1/gen p)/(1/gen p), endcap",    100, -scale, scale);
    
      LeptonInvPtPullEndcap[rec] = fs->make<TH1F>(nameHist("LeptonInvPtPullEndcap", rec), level + " (1/pT - 1/gen pT)/#sigma_{1/pT}, endcap", 100, -10, 10);
      LeptonInvPPullEndcap[rec]  = fs->make<TH1F>(nameHist("LeptonInvPPullEndcap",  rec), level + " (1/p - 1/gen p)/#sigma_{1/p}, endcap",    100, -10, 10);
    }
  }
}

void Zprime2muResolution::bookChargeResolutionHistos() {
  for (int rec = lGR; rec < MAX_LEVELS; ++rec) {
    TString level(levelName(rec));

    ChargeDiff[rec] = fs->make<TH1F>(nameHist("ChargeDiff", rec), level + " q - gen q", 7, -3.5, 3.5);

    ChargeRightVInvPt[rec] = fs->make<TH1F>(nameHist("ChargeRightVInvPt", rec), level + " right q vs. 1/(gen pT)", 50, 0, 0.01);
    ChargeWrongVInvPt[rec] = fs->make<TH1F>(nameHist("ChargeWrongVInvPt", rec), level + " wrong q vs. 1/(gen pT)", 50, 0, 0.01);

    ChargeRightVInvPt[rec]->Sumw2();
    ChargeWrongVInvPt[rec]->Sumw2();
  }
}

void Zprime2muResolution::bookEfficiencyHistos() {
  const TString type[3] = { "all events", "#eta < 2.4", "#eta < 2.1" };
  for (int rec = 0; rec < TRIG_LEVELS; ++rec) {
    TString level(levelName(rec));
    for (int i = 0; i < 3; ++i) {
      TrigEffVsDilMass[rec][i] = fs->make<TH1F>(nameHist("TrigEffVsDilMass", rec, i), "Gen mass, " + level + ", " + type[i], 27, 0., maxTrigMass);
      TrigEffVsDilMass[rec][i]->Sumw2();
    }
  }

  for (int rec = 0; rec < MAX_LEVELS; ++rec) {
    TString level(levelName(rec));
    EffVsEta[rec] = fs->make<TH1F>(nameHist("EffVsEta", rec), "Gen #eta, "  + level + " leptons", 50, -2.5, 2.5);
    EffVsPhi[rec] = fs->make<TH1F>(nameHist("EffVsPhi", rec), "Gen #phi, "  + level + " leptons", 50, -TMath::Pi(), TMath::Pi());
    EffVsPt[rec]  = fs->make<TH1F>(nameHist("EffVsPt",  rec), "Gen pT, "    + level + " leptons", 50, 0, 3500);

    EffVsEta[rec]->Sumw2();
    EffVsPhi[rec]->Sumw2();
    EffVsPt[rec]->Sumw2();
  }

  DilRecEffVsMass[0][0] = fs->make<TH1F>("DilRecEffVsMass00", "Gen mass, all",           27, 0, maxTrigMass);
  DilRecEffVsMass[0][1] = fs->make<TH1F>("DilRecEffVsMass01", "Gen mass, in acceptance", 27, 0, maxTrigMass);
  const TString dil_type[3] = { "2 l", "dil", "2 l, w/ cuts" };
  for (int rec = lGR; rec < MAX_LEVELS; ++rec) {
    TString level(levelName(rec));
    for (int i = 0; i < 3; ++i) {
      DilRecEffVsMass[rec][i] = fs->make<TH1F>(nameHist("DilRecEffVsMass", rec, i), "Gen mass, " + level + ", " + dil_type[i], 27, 0, maxTrigMass);
      DilRecEffVsMass[rec][i]->Sumw2();
    }
  }
}

void Zprime2muResolution::bookDileptonResolutionHistos() {
  for (int rec = lGR; rec < MAX_LEVELS; ++rec) {
    TString level(levelName(rec));

    DileptonMassRes[rec]    = fs->make<TH1F>(nameHist("DileptonMassRes",    rec), level + " (dil. mass - gen dil. mass)/(gen dil. mass)", 100, -0.3, 0.3);
    DileptonResMassRes[rec] = fs->make<TH1F>(nameHist("DileptonResMassRes", rec), level + " (dil. mass - gen res. mass)/(gen res. mass)", 100, -0.3, 0.3);
    ResonanceMassRes[rec]   = fs->make<TH1F>(nameHist("ResonanceMassRes",   rec), level + " (res. mass - gen res. mass)/(gen res. mass)", 100, -0.3, 0.3);

    DileptonMassResVMass[rec]    = fs->make<TProfile>(nameHist("DileptonMassResVMass",    rec), level + " (dil. mass - gen dil. mass)/(gen dil. mass)", 50, lowerMassWin, upperMassWin, -0.3, 0.3);
    DileptonResMassResVMass[rec] = fs->make<TProfile>(nameHist("DileptonResMassResVMass", rec), level + " (dil. mass - gen res. mass)/(gen res. mass)", 50, lowerMassWin, upperMassWin, -0.3, 0.3);
    ResonanceMassResVMass[rec]   = fs->make<TProfile>(nameHist("ResonanceMassResVMass",   rec), level + " (res. mass - gen res. mass)/(gen res. mass)", 50, lowerMassWin, upperMassWin, -0.3, 0.3);
  }
}

int Zprime2muResolution::encodeLeptonOrigin(const int id) const {
  // Group mother id into broad "origin" types defined similarly to
  // MuonSimtrackAnalyser.cc in ORCA.
  enum {undef=0, pi=1, K=2, KL=3, eta=4, rho=5, c=6, b=7, tau=8, Z=10, W=11, H=12, Zprime=13, G=14};

  switch (abs(id)) {
  case 211: return pi;              // pi+/-
  case 321: return K;               // K+/-
  case 311: return KL;              // K0
  case 221: return eta;             // eta
  case 113: return rho;             // rho0
  case 411: case 421: case 431:     // D+/-, D0, D_s+/-
  case 443: case 4122: case 100443: // j/psi, Lambda_c+, psi'
    return c; 
  case 511: case 521: case 531:     // B0, B+/-, B_s0
  case 553: case 5122: case 100553: // Y, Lambda_b0, Y'
    return b;
  case 15: return tau;              // tau
  case 23: return Z;                // Z0
  case 24: return W;                // W+/-
  case 25: return H;                // SM H
  case 32: return Zprime;           // Z'
  case 39: case 5000039: return G;  // G*

  default: return undef;
  }
}

void Zprime2muResolution::fillGenLevelHistos() {
  // PDG id of mothers of all generated leptons.
  for (reco::CandidateBaseRefVector::const_iterator lep = allLeptons[lGN].begin(); lep != allLeptons[lGN].end(); ++lep)
    LeptonOrigin[0]->Fill(encodeLeptonOrigin(motherId(*lep)));

  // PDG id of mothers of all leptons in all dileptons.
  for (reco::CompositeCandidateCollection::const_iterator dil = allDileptons[lGN].begin(); dil != allDileptons[lGN].end(); ++dil)
    for (unsigned i = 0; i < dil->numberOfDaughters(); ++i)
      LeptonOrigin[1]->Fill(encodeLeptonOrigin(motherId(dileptonDaughter(*dil, i))));

  // Acceptance for the dilepton chosen (i.e. the first one in the vector).
  if (allDileptons[lGN].size() > 0) {
    double mass = allDileptons[lGN][0].mass();
    GenMassAllEvents->Fill(mass);
    if (numDaughtersInAcc(allDileptons[lGN][0]) >= 2)
      GenMassInAccept->Fill(mass);
  }
}

void Zprime2muResolution::fillLeptonEfficiencyHistos() {
  // Efficiency to reconstruct muons from Z' decays at various trigger 
  // levels and by various off-line reconstructors.  L1/HLT efficiencies
  // are not included.

  for (reco::CandidateBaseRefVector::const_iterator gen_lep = allLeptons[lGN].begin(); gen_lep != allLeptons[lGN].end(); ++gen_lep) {
    // Only calculate efficiency for muons coming from the resonance.
    if (!HardInteraction::IsResonance(motherId(*gen_lep))) continue;
    
    double gen_eta = (*gen_lep)->eta();
    double gen_phi = (*gen_lep)->phi();
    double gen_pt  = (*gen_lep)->pt();
    int    gen_id  = recLevelHelper.id(*gen_lep);

    // Fill the denominator histos.
    EffVsEta[0]->Fill(gen_eta);
    EffVsPhi[0]->Fill(gen_phi);
    EffVsPt[0]->Fill(gen_pt);

    for (int rec = lL1; rec < MAX_LEVELS; ++rec) {
      //if (!trigDecision.pass(i)) continue;

      // If there is a matched lepton at this rec level, fill the
      // numerator histos.
      bool matched = false;
      for (reco::CandidateBaseRefVector::const_iterator lep = allLeptons[rec].begin(); lep != allLeptons[rec].end(); ++lep) {
	if (recLevelHelper.genMatchId(*lep) == gen_id) {
	  matched = true;
	  break;
	}
      }

      if (matched) {
	EffVsEta[rec]->Fill(gen_eta);
	EffVsPhi[rec]->Fill(gen_phi);
	EffVsPt[rec]->Fill(gen_pt);
      }
    }
  }
}

void Zprime2muResolution::fillTriggerEfficiencyHistos() {
  // Calculations below assume exactly one generated dilepton.
  if (allDileptons[lGN].size() != 1) return;
  const reco::CompositeCandidate& gen_dil = allDileptons[lGN].at(0);

  // Plots are efficiency versus mass in TeV.
  const double gen_mass = gen_dil.mass()/1000;

  // Plots are split into all events, events with both muons in full
  // eta coverage (|eta| < 2.4), and events with at least one muon in
  // the limited muon trigger acceptance (|eta| < 2.1).
  const bool accept[3] = {
    true,
    numDaughtersInAcc(gen_dil, 2.4) >= 2,
    numDaughtersInAcc(gen_dil, 2.4) >= 2 && numDaughtersInAcc(gen_dil, 2.1) >= 1
  };

  for (int rec = lGN; rec <= lL3; ++rec)
    if (trigDecision.pass(rec)) // defined true for rec == lGN
      for (int i = 0; i < 3; ++i)
	if (accept[i])
	  TrigEffVsDilMass[rec][i]->Fill(gen_mass);

  DilRecEffVsMass[lGN][0]->Fill(gen_mass);
  if (trigDecision.pass()) {
    if (accept[1])
      DilRecEffVsMass[lGN][1]->Fill(gen_mass);

    for (int rec = lGR; rec < MAX_LEVELS; ++rec) {
      if (allLeptons[rec].size() > 1)
	DilRecEffVsMass[rec][0]->Fill(gen_mass);
      if (allDileptons[rec].size() > 0)
	DilRecEffVsMass[rec][1]->Fill(gen_mass);
      
      unsigned passCut = 0;
      for (reco::CandidateBaseRefVector::const_iterator lep = allLeptons[rec].begin(); lep != allLeptons[rec].end(); lep++) {
	if (!cutHelper.leptonIsCut(**lep)) passCut++;
	if (passCut > 1) break;
      }

      if (passCut > 1)
	DilRecEffVsMass[rec][2]->Fill(gen_mass);
    }
  }      
}

void Zprime2muResolution::fillDileptonEfficiencyHistos() {
}

void Zprime2muResolution::fillLeptonResolution(const reco::CandidateBaseRef& gen_lep, const reco::CandidateBaseRef& lep, const int rec) {
  // Angular diffs.
  LeptonEtaDiff[rec]->Fill(lep->eta() - gen_lep->eta());
  LeptonPhiDiff[rec]->Fill(lep->phi() - gen_lep->phi());

  const double gen_pt  = gen_lep->pt();
  const double gen_p   = gen_lep->p();

  // Momentum diffs/resolutions.
  LeptonPtDiff[rec]->Fill(lep->pt() - gen_pt);

  LeptonPtRes[rec]->Fill((lep->pt() - gen_pt)/gen_pt);
  LeptonPRes [rec]->Fill((lep->p()  - gen_p) /gen_p);

  const double inv_pt_diff = 1/lep->pt() - 1/gen_pt;
  const double inv_p_diff  = 1/lep->p()  - 1/gen_p;

  // Inverse momentum resolutions.
  LeptonInvPtRes[rec]->Fill(inv_pt_diff/(1/gen_pt));
  LeptonInvPRes [rec]->Fill(inv_p_diff /(1/gen_p));
}

void Zprime2muResolution::fillLeptonExtraMomentumResolution(const reco::CandidateBaseRef& gen_lep, const reco::CandidateBaseRef& lep, const int rec) {
  // More histograms (pulls, profiles) for offline reconstructed
  // leptons.
  const double inv_gen_pt  = 1/gen_lep->pt();
  const double inv_gen_p   = 1/gen_lep->p();

  const double inv_pt_diff = 1/lep->pt() - inv_gen_pt;
  const double inv_p_diff  = 1/lep->p()  - inv_gen_p;

  const double inv_pt_res = inv_pt_diff/inv_gen_pt;
  const double inv_p_res  = inv_p_diff /inv_gen_p;

  // Inverse momentum resolutions as a function of generated momenta.
  LeptonInvPtResVPtGen[rec]->Fill(1/inv_gen_pt, inv_pt_res*inv_pt_res);
  LeptonInvPResVPGen  [rec]->Fill(1/inv_gen_p,  inv_p_res*inv_p_res);
  
  // Try to get the reconstructed momentum errors for pulls.
  double inv_pt_error, inv_p_error;
  inv_pt_error = inv_p_error = -999;
  bool errorOK = false;
  const reco::RecoCandidate* cand = toConcretePtr<reco::RecoCandidate>(lep);
  if (cand) {
    const reco::Track* tk = cand->combinedMuon().get();
    if (tk) {
      errorOK = true;
      inv_pt_error = invPtError(tk);
      inv_p_error  = invPtError(tk);
    }
  }      

  // Inverse momentum pulls.
  if (errorOK) {
    LeptonInvPtPull[rec]->Fill(inv_pt_diff/inv_pt_error);
    LeptonInvPPull [rec]->Fill(inv_p_diff /inv_p_error);
  }
  
  // The above inverse momentum resolutions and pulls, except
  // separately for barrel and endcap.
  if (fabs(gen_lep->eta()) < 1.04) {
    LeptonInvPtResBarrel[rec]->Fill(inv_pt_res);
    LeptonInvPResBarrel [rec]->Fill(inv_p_res);
    
    if (errorOK) {
      LeptonInvPtPullBarrel[rec]->Fill(inv_pt_diff/inv_pt_error);
      LeptonInvPPullBarrel [rec]->Fill(inv_p_diff /inv_p_error);
    }
  }
  else {
    LeptonInvPtResEndcap[rec]->Fill(inv_pt_res);
    LeptonInvPResEndcap [rec]->Fill(inv_p_res);
    
    if (errorOK) {
      LeptonInvPtPullEndcap[rec]->Fill(inv_pt_diff/inv_pt_error);
      LeptonInvPPullEndcap [rec]->Fill(inv_p_diff /inv_p_error);
    }
  }
}

void Zprime2muResolution::fillChargeResolution(const reco::CandidateBaseRef& gen_lep, const reco::CandidateBaseRef& lep, const int rec) {
  const int delta_q = lep->charge() - gen_lep->charge();
  ChargeDiff[rec]->Fill(delta_q);
  if (delta_q == 0)
    ChargeRightVInvPt[rec]->Fill(1/gen_lep->pt());
  else //if (delta_q == 2), also don't throw exception?
    ChargeWrongVInvPt[rec]->Fill(1/gen_lep->pt());
}

void Zprime2muResolution::fillDileptonMassResolution(const reco::CompositeCandidate& gen_dil, const reco::CompositeCandidate& dil, const int rec) {
  const double mass         = dil.mass();
  const double gen_mass     = gen_dil.mass();

  const double res_mass     = resonanceMass(dil);
  const double gen_res_mass = resonanceMass(gen_dil);

  const double rdil    = mass    /gen_mass     - 1;
  const double rdilres = mass    /gen_res_mass - 1;
  const double rres    = res_mass/gen_res_mass - 1;
  
  DileptonMassRes   [rec]->Fill(rdil);
  DileptonResMassRes[rec]->Fill(rdilres);
  ResonanceMassRes  [rec]->Fill(rres);

  DileptonMassResVMass   [rec]->Fill(gen_mass,     rdil*rdil);
  DileptonResMassResVMass[rec]->Fill(gen_res_mass, rdilres*rdilres);
  ResonanceMassResVMass  [rec]->Fill(gen_res_mass, rres*rres);
}

void Zprime2muResolution::fillLeptonHistos(const reco::CandidateBaseRef& lep, const int rec) {
  const reco::CandidateBaseRef& gen_lep = recLevelHelper.matchGenLepton(lep);

  if (gen_lep.isNonnull()) {
    fillLeptonResolution(gen_lep, lep, rec);

    if (rec >= lGR) {
      fillLeptonExtraMomentumResolution(gen_lep, lep, rec);
      fillChargeResolution(gen_lep, lep, rec);
    }
  }
}

void Zprime2muResolution::fillLeptonHistos(const int rec) {
  if (leptonsFromDileptons) {
    // Only fill lepton histos from the leptons that made it into dileptons.
    for (reco::CompositeCandidateCollection::const_iterator dil = allDileptons[rec].begin(); dil != allDileptons[rec].end(); ++dil)
      for (unsigned i = 0; i < dil->numberOfDaughters(); ++i)
	fillLeptonHistos(dileptonDaughter(*dil, i), rec);
  }
  else {
    // Fill lepton histos from all leptons.
    for (reco::CandidateBaseRefVector::const_iterator lep = allLeptons[rec].begin(); lep != allLeptons[rec].end(); ++lep)
      fillLeptonHistos(*lep, rec);
  }
}

void Zprime2muResolution::fillDileptonHistos(const int rec) {
  if (allDileptons[lGN].size() > 0) {
    // Only one dilepton at generator level to choose.
    const reco::CompositeCandidate& gen_dil = allDileptons[lGN].at(0);

    for (reco::CompositeCandidateCollection::const_iterator dil = allDileptons[rec].begin(); dil != allDileptons[rec].end(); ++dil)
      fillDileptonMassResolution(gen_dil, *dil, rec);
  }
}

void Zprime2muResolution::analyze(const edm::Event& event, const edm::EventSetup& eSetup) {
  // Delegate filling our lepton vectors to the parent class.
  Zprime2muAnalysis::analyze(event, eSetup);

  fillGenLevelHistos();

  fillLeptonEfficiencyHistos();
  fillTriggerEfficiencyHistos();
  fillDileptonEfficiencyHistos();

  for (int rec = lL1; rec < MAX_LEVELS; ++rec) {
    fillLeptonHistos(rec);
    if (rec >= lGR)
      fillDileptonHistos(rec);
  }
}

void Zprime2muResolution::makeRMSHist(const TProfile* prof) const {
  // Produce a histogram whose bins and errors are from RMS of the
  // bins of the TProfile passed in.
  int nbins = prof->GetNbinsX();
  TAxis* axis = prof->GetXaxis();
  TString RMS("RMS");
  TH1F* h = fs->make<TH1F>(RMS + prof->GetName(), RMS + TString(" ") + prof->GetTitle(),
			   nbins, axis->GetXmin(), axis->GetXmax());
  for (int ibin = 1; ibin <= nbins; ibin++) {
    double f_bin   = prof->GetBinContent(ibin);
    double ent_bin = prof->GetBinEntries(ibin);

    f_bin = f_bin > 0 ? sqrt(f_bin) : 0;
    double err_bin = ent_bin > 0 ? f_bin/sqrt(2*ent_bin) : 0;

    h->SetBinContent(ibin, f_bin);
    h->SetBinError(ibin, err_bin);
  }
}

void Zprime2muResolution::endJob() {
  Zprime2muAnalysis::endJob();

  for (int rec = lGR; rec < MAX_LEVELS; ++rec) {
    makeRMSHist(LeptonInvPtResVPtGen[rec]);
    makeRMSHist(LeptonInvPResVPGen[rec]);
    makeRMSHist(DileptonMassResVMass[rec]);
    makeRMSHist(DileptonResMassResVMass[rec]);
    makeRMSHist(ResonanceMassResVMass[rec]);
  }
}

DEFINE_FWK_MODULE(Zprime2muResolution);
