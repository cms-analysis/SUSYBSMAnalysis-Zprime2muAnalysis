#include <fstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymmetryHelpers.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Functions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"

class AsymmetryFrameAnalysis : public edm::EDAnalyzer {
public:
  explicit AsymmetryFrameAnalysis(const edm::ParameterSet&);
  
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();
  
private:
  void fillFrameHistos(const pat::CompositeCandidate&, const reco::CandidateBaseRef&, const reco::CandidateBaseRef&);
  void drawFrameHistos();

  const edm::InputTag dilepton_src;
  const bool debug;

  TH1F* cosGJ;
  TH1F* cosGJTag;
  TH1F* cosCS;
  TH1F* cosCSAn;
  TH1F* cosBoost;
  TH1F* cosW;
  TH2F* rap_vs_cosCS;
  TH1F* FMassGJ;
  TH1F* FMassGJTag;
  TH1F* FMassCS;
  TH1F* FMassCSAn;
  TH1F* FMassBoost;
  TH1F* FMassW;
  TH1F* BMassGJ;
  TH1F* BMassGJTag;
  TH1F* BMassCS;
  TH1F* BMassCSAn;
  TH1F* BMassBoost;
  TH1F* BMassW;
  TH1F* AMassGJ;
  TH1F* AMassGJTag;
  TH1F* AMassCS;
  TH1F* AMassCSAn;
  TH1F* AMassBoost;
  TH1F* AMassW;
  TH1F* FRapGJ;
  TH1F* FRapGJTag;
  TH1F* FRapCS;
  TH1F* FRapCSAn;
  TH1F* FRapBoost;
  TH1F* FRapW;
  TH1F* BRapGJ;
  TH1F* BRapGJTag;
  TH1F* BRapCS;
  TH1F* BRapCSAn;
  TH1F* BRapBoost;
  TH1F* BRapW;
  TH1F* ARapGJ;
  TH1F* ARapGJTag;
  TH1F* ARapCS;
  TH1F* ARapCSAn;
  TH1F* ARapBoost;
  TH1F* ARapW;
  TH1F* FPseudGJ;
  TH1F* FPseudCS;
  TH1F* FPseudBoost;
  TH1F* FPseudW;
  TH1F* BPseudGJ;
  TH1F* BPseudCS;
  TH1F* BPseudBoost;
  TH1F* BPseudW;
  TH1F* FMBoostCut[6];
  TH1F* BMBoostCut[6];
  TH1F* AsymMBoostCut[6];
};

AsymmetryFrameAnalysis::AsymmetryFrameAnalysis(const edm::ParameterSet& cfg)
  : dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    debug(cfg.getUntrackedParameter<bool>("debug"))
{
  asymFitManager.setConstants(cfg);

  edm::Service<TFileService> fs;

  const int NBIN = 9;
  const double DMBINS[NBIN] = {200, 400, 500, 750, 1000, 1250, 1500, 2000, 3000};

  cosGJ    = fs->make<TH1F>("cosGJ",    "cos #theta in Gottfried-Jackson Frame",        100, -1, 1);
  cosGJTag = fs->make<TH1F>("cosGJTag", "cos #theta in tagged Gottfried-Jackson Frame", 100, -1, 1);
  cosCS    = fs->make<TH1F>("cosCS",    "cos #theta in Collins-Soper Frame",            100, -1, 1);
  cosCSAn  = fs->make<TH1F>("cosCSAn",  "cos #theta in analytical Collins-Soper Frame", 100, -1, 1);
  cosBoost = fs->make<TH1F>("cosBoost", "cos #theta in Boost Frame",                    100, -1, 1);
  cosW     = fs->make<TH1F>("cosW",     "cos #theta in Wulz Frame",                     100, -1, 1);

  rap_vs_cosCS = fs->make<TH2F>("rap_vs_cosCS", "Rap vs Cos theta CS", 50, -1, 1, 50, -4, 4);

  FMassGJ    = fs->make<TH1F>("FMassGJ",    "F(M), GJ frame",        NBIN-1, DMBINS);
  FMassGJTag = fs->make<TH1F>("FMassGJTag", "F(M), tagged GJ frame", NBIN-1, DMBINS);
  FMassCS    = fs->make<TH1F>("FMassCS",    "F(M), CS frame",        NBIN-1, DMBINS);
  FMassCSAn  = fs->make<TH1F>("FMassCSAn",  "F(M), analyt CS frame", NBIN-1, DMBINS);
  FMassBoost = fs->make<TH1F>("FMassBoost", "F(M), Boost frame",     NBIN-1, DMBINS);
  FMassW     = fs->make<TH1F>("FMassW",     "F(M), Wulz frame",      NBIN-1, DMBINS);

  BMassGJ    = fs->make<TH1F>("BMassGJ",    "B(M), GJ frame",        NBIN-1, DMBINS);
  BMassGJTag = fs->make<TH1F>("BMassGJTag", "B(M), tagged GJ frame", NBIN-1, DMBINS);
  BMassCS    = fs->make<TH1F>("BMassCS",    "B(M), CS frame",        NBIN-1, DMBINS);
  BMassCSAn  = fs->make<TH1F>("BMassCSAn",  "B(M), analyt CS frame", NBIN-1, DMBINS);
  BMassBoost = fs->make<TH1F>("BMassBoost", "B(M), Boost frame",     NBIN-1, DMBINS);
  BMassW     = fs->make<TH1F>("BMassW",     "B(M), Wulz frame",      NBIN-1, DMBINS);

  AMassGJ    = fs->make<TH1F>("AMassGJ",    "A(M), GJ frame",        NBIN-1, DMBINS);
  AMassGJTag = fs->make<TH1F>("AMassGJTag", "A(M), tagged GJ frame", NBIN-1, DMBINS);
  AMassCS    = fs->make<TH1F>("AMassCS",    "A(M), CS frame",        NBIN-1, DMBINS);
  AMassCSAn  = fs->make<TH1F>("AMassCSAn",  "A(M), analyt CS frame", NBIN-1, DMBINS);
  AMassBoost = fs->make<TH1F>("AMassBoost", "A(M), Boost frame",     NBIN-1, DMBINS);
  AMassW     = fs->make<TH1F>("AMassW",     "A(M), Wulz frame",      NBIN-1, DMBINS);

  FRapGJ    = fs->make<TH1F>("FRapGJ",    "F(y), GJ frame",        20, -2.5, 2.5);
  FRapGJTag = fs->make<TH1F>("FRapGJTag", "F(y), tagged GJ frame", 20, -2.5, 2.5);
  FRapCS    = fs->make<TH1F>("FRapCS",    "F(y), CS frame",        20, -2.5, 2.5);
  FRapCSAn  = fs->make<TH1F>("FRapCSAn",  "F(y), analyt CS frame", 20, -2.5, 2.5);
  FRapBoost = fs->make<TH1F>("FRapBoost", "F(y), Boost frame",     20, -2.5, 2.5);
  FRapW     = fs->make<TH1F>("FRapW",     "F(y), Wulz frame",      20, -2.5, 2.5);

  BRapGJ    = fs->make<TH1F>("BRapGJ",    "B(y), GJ frame",        20, -2.5, 2.5);
  BRapGJTag = fs->make<TH1F>("BRapGJTag", "B(y), tagged GJ frame", 20, -2.5, 2.5);
  BRapCS    = fs->make<TH1F>("BRapCS",    "B(y), CS frame",        20, -2.5, 2.5);
  BRapCSAn  = fs->make<TH1F>("BRapCSAn",  "B(y), analyt CS frame", 20, -2.5, 2.5);
  BRapBoost = fs->make<TH1F>("BRapBoost", "B(y), Boost frame",     20, -2.5, 2.5);
  BRapW     = fs->make<TH1F>("BRapW",     "B(y), Wulz frame",      20, -2.5, 2.5);

  ARapGJ    = fs->make<TH1F>("ARapGJ",    "A(y), GJ frame",        20, -2.5, 2.5);
  ARapGJTag = fs->make<TH1F>("ARapGJTag", "A(y), tagged GJ frame", 20, -2.5, 2.5);
  ARapCS    = fs->make<TH1F>("ARapCS",    "A(y), CS frame",        20, -2.5, 2.5);
  ARapCSAn  = fs->make<TH1F>("ARapCSAn",  "A(y), analyt CS frame", 20, -2.5, 2.5);
  ARapBoost = fs->make<TH1F>("ARapBoost", "A(y), Boost frame",     20, -2.5, 2.5);
  ARapW     = fs->make<TH1F>("ARapW",     "A(y), Wulz frame",      20, -2.5, 2.5);

  FPseudGJ    = fs->make<TH1F>("FPseudGJ",    "F(#eta), tagged GJ frame", 50, -6, 6);
  FPseudCS    = fs->make<TH1F>("FPseudCSAn",  "F(#eta), CS frame",        50, -6, 6);
  FPseudBoost = fs->make<TH1F>("FPseudBoost", "F(#eta), Boost frame",     50, -6, 6);
  FPseudW     = fs->make<TH1F>("FPseudW",     "F(#eta), Wulz frame",      50, -6, 6);

  BPseudGJ    = fs->make<TH1F>("BPseudGJ",    "B(#eta), tagged GJ frame", 50, -6, 6);
  BPseudCS    = fs->make<TH1F>("BPseudCSAn",  "B(#eta), CS frame",        50, -6, 6);
  BPseudBoost = fs->make<TH1F>("BPseudBoost", "B(#eta), Boost frame",     50, -6, 6);
  BPseudW     = fs->make<TH1F>("BPseudW",     "B(#eta), Wulz frame",      50, -6, 6);

  FMBoostCut[0] = fs->make<TH1F>("FMBoostCut0", "F(M), boost, |y|<0.4",        NBIN-1, DMBINS);
  FMBoostCut[1] = fs->make<TH1F>("FMBoostCut1", "F(M), boost, 0.4<|y|<0.8",    NBIN-1, DMBINS);
  FMBoostCut[2] = fs->make<TH1F>("FMBoostCut2", "F(M), boost, 0.8<|y|<2.4",    NBIN-1, DMBINS);
  FMBoostCut[3] = fs->make<TH1F>("FMBoostCut3", "F(M), boost, |#eta|<0.4",     NBIN-1, DMBINS);
  FMBoostCut[4] = fs->make<TH1F>("FMBoostCut4", "F(M), boost, 0.4<|#eta|<0.8", NBIN-1, DMBINS);
  FMBoostCut[5] = fs->make<TH1F>("FMBoostCut5", "F(M), boost, 0.8<|#eta|<2.4", NBIN-1, DMBINS);

  BMBoostCut[0] = fs->make<TH1F>("BMBoostCut0", "B(M), boost, |y|<0.4",        NBIN-1, DMBINS);
  BMBoostCut[1] = fs->make<TH1F>("BMBoostCut1", "B(M), boost, 0.4<|y|<0.8",    NBIN-1, DMBINS);
  BMBoostCut[2] = fs->make<TH1F>("BMBoostCut2", "B(M), boost, 0.8<|y|<2.4",    NBIN-1, DMBINS);
  BMBoostCut[3] = fs->make<TH1F>("BMBoostCut3", "B(M), boost, |#eta|<0.4",     NBIN-1, DMBINS);
  BMBoostCut[4] = fs->make<TH1F>("BMBoostCut4", "B(M), boost, 0.4<|#eta|<0.8", NBIN-1, DMBINS);
  BMBoostCut[5] = fs->make<TH1F>("BMBoostCut5", "B(M), boost, 0.8<|#eta|<2.4", NBIN-1, DMBINS);

  AsymMBoostCut[0] = fs->make<TH1F>("AsymMBoostCut0", "A(M), boost, |y|<0.4",        NBIN-1, DMBINS);
  AsymMBoostCut[1] = fs->make<TH1F>("AsymMBoostCut1", "A(M), boost, 0.4<|y|<0.8",    NBIN-1, DMBINS);
  AsymMBoostCut[2] = fs->make<TH1F>("AsymMBoostCut2", "A(M), boost, 0.8<|y|<2.4",    NBIN-1, DMBINS);
  AsymMBoostCut[3] = fs->make<TH1F>("AsymMBoostCut3", "A(M), boost, |#eta|<0.4",     NBIN-1, DMBINS);
  AsymMBoostCut[4] = fs->make<TH1F>("AsymMBoostCut4", "A(M), boost, 0.4<|#eta|<0.8", NBIN-1, DMBINS);
  AsymMBoostCut[5] = fs->make<TH1F>("AsymMBoostCut5", "A(M), boost, 0.8<|#eta|<2.4", NBIN-1, DMBINS);
}

void AsymmetryFrameAnalysis::fillFrameHistos(const pat::CompositeCandidate& dil, const reco::CandidateBaseRef& mum, const reco::CandidateBaseRef& mup) {
  // Calculate cosine theta^* in various frames.

  static const double beam_energy = asymFitManager.beam_energy();
  static const TLorentzVector pp1(0., 0.,  beam_energy, beam_energy);
  static const TLorentzVector pp2(0., 0., -beam_energy, beam_energy);

  if (debug) std::cout << "pp1: " << pp1 << " pp2: " << pp2 << "\n";

  TLorentzVector vmum, vmup, vdil;
  vmum.SetPtEtaPhiM(mum->pt(), mum->eta(), mum->phi(), mum->mass());
  vmup.SetPtEtaPhiM(mup->pt(), mup->eta(), mup->phi(), mup->mass());
  vdil.SetPtEtaPhiM(dil.pt(),  dil.eta(),  dil.phi(),  dil.mass());

  if (debug) std::cout << "Dilepton:" << "\n"
		       << " l+: q: " << mup->charge() << " p4: " << vmup << "\n"
		       << " l-: q: " << mum->charge() << " p4: " << vmum << "\n"
		       << "dil: q: " << dil.charge()  << " p4: " << vdil << "\n";

  // Do this to prevent dividing by zero further on.
  if (vmup.Pz() == 0 || vmum.Pz() == 0) {
    if (debug) std::cout << "Skipping event since pL == 0.\n";
    return;
  }

  const double mass = vdil.M();
  const double rap  = vdil.Rapidity();
  const double eta  = vdil.PseudoRapidity();

  //-----------------------------------------------------------------
  // Compute cosine theta^* in Gottfried-Jackson frame. This frame
  // is defined as the frame with the z^* axis parallel to the
  // beam momentum. This makes the pT of the proton > 0 in the Z'
  // CMS frame. This is a good frame for p-pbar collisions but a
  // bad one for pp collisions because we don't know which of the
  // protons had the quark that interacted.

  // Lorentz boost mu-/+ and beam to dilepton frame.
  const TLorentzVector pmum_star  = LorentzBoost(vdil, vmum);
  const TLorentzVector pmup_star  = LorentzBoost(vdil, vmup);
  const TLorentzVector pb_star_gj = LorentzBoost(vdil, pp1);

  if (debug) std::cout << "l- in dilepton rest frame: " << pmum_star << "\n"
		       << "l+ in dilepton rest frame: " << pmup_star << "\n"
		       << "pb_star_gj: " << pb_star_gj << "\n";

  const double cos_gj = cos_angle(pb_star_gj, pmum_star);
  if (debug) std::cout << "Cos theta^* in Gottfried-Jackson frame: " << cos_gj << "\n";

  cosGJ->Fill(cos_gj);
  if (cos_gj > 0) {
    FMassGJ->Fill(mass);
    FRapGJ->Fill(rap);
  }
  else {
    BMassGJ->Fill(mass);
    BRapGJ->Fill(rap);
  }

  // Now adopt this frame to the pp collisions, by tagging the "beam"
  // proton using the direction of the dilepton system w.r.t the 
  // pp collision axis.
  const bool dil_in_pp1_dir = pp1.Pz() * vdil.Pz() > 0;
  const TLorentzVector pb_star = LorentzBoost(vdil, dil_in_pp1_dir ? pp1 : pp2);
  const double cos_gj_tag = cos_angle(pb_star, pmum_star);
  if (debug) std::cout << "dil in pp1 dir? " << dil_in_pp1_dir << " pb_star: " << pb_star
		       << "\nCos theta^* in tagged Gottfried-Jackson frame: " << cos_gj_tag << "\n";

  cosGJTag->Fill(cos_gj_tag);
  if (cos_gj_tag > 0) {
    FMassGJTag->Fill(mass);
    FRapGJTag->Fill(rap);
    FPseudGJ->Fill(eta);
  }
  else {
    BMassGJTag->Fill(mass);
    BRapGJTag->Fill(rap);
    BPseudGJ->Fill(eta);
  }

  //---------------------------------------------------------------
  // Compute cosine theta^* in Collins-Soper frame. The z^* axis
  // is chosen so that it bisects the angle between the pBeam and
  // -pTarget. This is the best frame in p-pbar collisions since
  // it reduces the uncertainty introduced by the fact that p and
  // p_bar are not parallel in the dilepton rest frame, and the
  // quark directions are not the same as the p and p_bar
  // directions. Choose the pBeam using the sign of the dilepton
  // system. Identify the beam as the p that is in the same
  // direction as the dilepton. Then boost the beam and target
  // into the dilepton CMS frame.
  const TLorentzVector pt_star = LorentzBoost(vdil, dil_in_pp1_dir ? pp2 : pp1);

  // cos_theta_cs = angle between pb_star_hat - pt_star_hat and pmum_star.
  double cos_cs = cos_angle(pb_star.Vect().Unit() - pt_star.Vect().Unit(), pmum_star.Vect());
  if (debug) std::cout << "pt_star: " << pt_star << "\nCos theta^* in Collins-Soper frame: " << cos_cs << "\n";

  cosCS->Fill(cos_cs);
  rap_vs_cosCS->Fill(cos_cs, rap);

  if (cos_cs > 0.) {
    FMassCS->Fill(mass);
    FRapCS->Fill(rap);
    FPseudCS->Fill(eta);
  }
  else {
    BMassCS->Fill(mass);
    BRapCS->Fill(rap);
    BPseudCS->Fill(eta);
  }

  // Now compute the analytic expression. Why do we get exactly the
  // same answer as for the exact formula?
  const double cos_cs_an = calcCosThetaCSAnal(vmum.Pz(), vmum.E(), vmup.Pz(), vmup.E(), vdil.Pt(), vdil.Pz(), mass);
  if (debug) std::cout << "Cos theta^* in Collins-Soper analytic frame: " << cos_cs_an << "\n";

  cosCSAn->Fill(cos_cs_an);
  if (cos_cs_an > 0) {
    FMassCSAn->Fill(mass);
    FRapCSAn->Fill(rap);
  }
  else {
    BMassCSAn->Fill(mass);
    BRapCSAn->Fill(rap);
  }

  //-----------------------------------------------------
  // Now calculate in Baur-boost frame.  This frame takes the
  // quark direction as the boost direction of th, bool res_oke dilepton
  // system.  This approach is described in M. Dittmar,
  // Phys. Rev. D55 (1997) 161; U. Baur et al., hep-ph/9707301.
  const double cos_boost = cos_angle(vdil, pmum_star);
  if (debug) std::cout << "Cos theta^* in Baur-boost frame: " << cos_boost << "\n";

  // Again fill the histos.
  cosBoost->Fill(cos_boost);
  if (cos_boost > 0) {
    FMassBoost->Fill(mass);
    FRapBoost->Fill(rap);
    FPseudBoost->Fill(eta);
  }
  else {
    BMassBoost->Fill(mass);
    BRapBoost->Fill(rap);
    BPseudBoost->Fill(eta);
  }

  // See how this approximation works in various intervals of rapidity.
  if (cos_boost > 0.) {
    if (fabs(rap) < 0.4)
      FMBoostCut[0]->Fill(mass);
    else if (fabs(rap) < 0.8)
      FMBoostCut[1]->Fill(mass);
    else if (fabs(rap) < 2.4)
      FMBoostCut[2]->Fill(mass);

    if (fabs(eta) < 0.4)
      FMBoostCut[3]->Fill(mass);
    else if (fabs(eta) < 0.8)
      FMBoostCut[4]->Fill(mass);
    else if (fabs(eta) < 2.4)
      FMBoostCut[5]->Fill(mass);
  }
  else {
    if (fabs(rap) < 0.4)
      BMBoostCut[0]->Fill(mass);
    else if (fabs(rap) < 0.8)
      BMBoostCut[1]->Fill(mass);
    else if (fabs(rap) < 2.4)
      BMBoostCut[2]->Fill(mass);

    if (fabs(eta) < 0.4)
      BMBoostCut[3]->Fill(mass);
    else if (fabs(eta) < 0.8)
      BMBoostCut[4]->Fill(mass);
    else if (fabs(eta) < 2.4)
      BMBoostCut[5]->Fill(mass);
  }

  //----------------------------------------------------
  // Wulz Frame: This is a variant of Baur frame as described in
  // CMS-TN/93-107 Cosine theta^* is measured between mu+ and Z'
  // for negative rapidity values, and between mu- and Z' for
  // positive rapidity.  The result is an odd function of rapidity
  // and when integrated over the full rapidity interval, gives
  // zero.
  const double cos_wulz = cos_angle(vdil, rap > 0 ? pmum_star : pmup_star); 
  if (debug) std::cout << "Cos theta^* in Wulz frame: " << cos_wulz << "\n";
  cosW->Fill(cos_wulz);
  if (cos_wulz > 0) {
    FMassW->Fill(mass);
    FRapW->Fill(rap);
    FPseudW->Fill(eta);
  }
  else {
    BMassW->Fill(mass);
    BRapW->Fill(rap);
    BPseudW->Fill(eta);  
  }
}

void AsymmetryFrameAnalysis::analyze(const edm::Event& event, const edm::EventSetup&) {
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);

  if (debug) std::cout << "\nRun: " << event.eventAuxiliary().run() << " Event: " << event.eventAuxiliary().event() << "; " << dileptons->size() << " dilepton(s)\n";

  for (size_t i = 0; i < dileptons->size(); ++i) {
    const pat::CompositeCandidate& dil = dileptons->at(i);
    const reco::CandidateBaseRef& mum = dileptonDaughterByCharge(dil, -1);
    const reco::CandidateBaseRef& mup = dileptonDaughterByCharge(dil, +1);
    fillFrameHistos(dil, mum, mup);
  }
}

void AsymmetryFrameAnalysis::endJob() {
  double AFB, eAFB;

  calcAsymmetry(FMassGJ, BMassGJ, AMassGJ, AFB, eAFB);
  std::cout << "Asymmetry in GJ frame:         F = " << FMassGJ->GetEntries() << " B = " << BMassGJ->GetEntries() << " A_FB = " << AFB << " +/- " << eAFB << "\n";
  calcAsymmetry(FMassGJTag, BMassGJTag, AMassGJTag, AFB, eAFB);
  std::cout << "Asymmetry in tagged GJ frame:  F = " << FMassGJTag->GetEntries() << " B = " << BMassGJTag->GetEntries() << " A_FB = " << AFB << " +/- " << eAFB << "\n";
  calcAsymmetry(FMassCS, BMassCS, AMassCS, AFB, eAFB);
  std::cout << "Asymmetry in CS frame:         F = " << FMassCS->GetEntries() << " B = " << BMassCS->GetEntries() << " A_FB = " << AFB << " +/- " << eAFB << "\n";
  calcAsymmetry(FMassCSAn, BMassCSAn, AMassCSAn, AFB, eAFB);
  std::cout << "Asymmetry in analyt CS frame:  F = " << FMassCSAn->GetEntries() << " B = " << BMassCSAn->GetEntries() << " A_FB = " << AFB << " +/- " << eAFB << "\n";
  calcAsymmetry(FMassBoost, BMassBoost, AMassBoost, AFB, eAFB);
  std::cout << "Asymmetry in Boost frame:      F = " << FMassBoost->GetEntries() << " B = " << BMassBoost->GetEntries() << " A_FB = " << AFB << " +/- " << eAFB << "\n";
  calcAsymmetry(FMassW, BMassW, AMassW, AFB, eAFB);
  std::cout << "Asymmetry in Wulz frame:       F = " << FMassW->GetEntries() << " B = " << BMassW->GetEntries() << " A_FB = " << AFB << " +/- " << eAFB << "\n";
  for (int i = 0; i < 6; i++) {
    calcAsymmetry(FMBoostCut[i], BMBoostCut[i], AsymMBoostCut[i], AFB, eAFB);
    std::cout << "Asymmetry in MBCut[" << i << "] frame:   F = " << FMBoostCut[i]->GetEntries() << " B = " << BMBoostCut[i]->GetEntries() << " A_FB = " << AFB << " +/- " << eAFB << "\n";
  }

  calcAsymmetry(FRapGJ,    BRapGJ,    ARapGJ,    AFB, eAFB);
  calcAsymmetry(FRapGJTag, BRapGJTag, ARapGJTag, AFB, eAFB);
  calcAsymmetry(FRapCS,    BRapCS,    ARapCS,    AFB, eAFB);
  calcAsymmetry(FRapCSAn,  BRapCSAn,  ARapCSAn,  AFB, eAFB);
  calcAsymmetry(FRapBoost, BRapBoost, ARapBoost, AFB, eAFB);
  calcAsymmetry(FRapW,     BRapW,     ARapW,     AFB, eAFB);
    
  drawFrameHistos();
}

void AsymmetryFrameAnalysis::drawFrameHistos() {
  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 700);
  TPostScript *ps = new TPostScript("diffFrameAsym.ps", 111);

  const int NUM_PAGES = 8;
  TPad* pad[NUM_PAGES] = {0};
  for (int i_page = 0; i_page < NUM_PAGES; ++i_page)
    pad[i_page] = new TPad("", "", .05, .05, .95, .93);

  int page = 0;
  TPaveLabel *title = 0;
  TText t;
  gStyle->SetOptStat(1110);

  // Gottfried-Jackson frame
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Gottfried-Jackson Frame");
  title->SetFillColor(10);
  title->Draw();
  t.DrawText(.9, .02, TString::Format("- %i -", ++page));
  pad[page]->Draw();
  pad[page]->Divide(2,2);
  pad[page]->cd(1); cosGJ->Draw(); fitCosTheta(std::cout, cosGJ);
  pad[page]->cd(2); cosGJTag->Draw(); fitCosTheta(std::cout, cosGJTag);
  pad[page]->cd(3); AMassGJ->Draw();
  pad[page]->cd(4); AMassGJTag->Draw();
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title->Draw();
  t.DrawText(.9, .02, TString::Format("- %i -", ++page));
  pad[page]->Draw();
  pad[page]->Divide(2,2);
  pad[page]->cd(1); ARapGJ->Draw();
  pad[page]->cd(2); ARapGJTag->Draw();
  pad[page]->cd(3); FPseudGJ->Draw();
  pad[page]->cd(4); BPseudGJ->Draw();
  c1->Update();

  // Collins-Soper frame
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Collins-Soper Frame");
  title->SetFillColor(10);
  title->Draw();
  t.DrawText(.9, .02, TString::Format("- %i -", ++page));
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1); cosCS->Draw(); fitCosTheta(std::cout, cosCS);
  pad[page]->cd(3); AMassCS->Draw();
  pad[page]->cd(4); ARapCS->Draw();
  pad[page]->cd(5); FPseudCS->Draw();
  pad[page]->cd(6); BPseudCS->Draw();
  c1->Update();

  // Collins-Soper frame, special scatter plots and fits
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Collins-Soper Frame, y vs cos");
  title->SetFillColor(10);
  title->Draw();
  t.DrawText(.9, .02, TString::Format("- %i -", ++page));
  pad[page]->Draw();
  pad[page]->Divide(1,2);
  pad[page]->cd(1);
  rap_vs_cosCS->Draw();
  c1->Update();

  // Quark direction is chosen as the boost direction of the dilepton system
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Baur-Boost Frame");
  title->SetFillColor(10);
  title->Draw();
  t.DrawText(.9, .02, TString::Format("- %i -", ++page));
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1); cosBoost->Draw(); fitCosTheta(std::cout, cosBoost);
  pad[page]->cd(3); AMassBoost->Draw();
  pad[page]->cd(4); ARapBoost->Draw();
  pad[page]->cd(5); FPseudBoost->Draw();
  pad[page]->cd(6); BPseudBoost->Draw();
  c1->Update();

  // Drawing the histos for the Asymmetry cuts
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title->Draw();
  t.DrawText(.9, .02, TString::Format("- %i -", ++page));
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  for (int k = 0; k < 3; k++) {
    pad[page]->cd(2*k+1);  AsymMBoostCut[k]->Draw();
    pad[page]->cd(2*k+2);  AsymMBoostCut[k+3]->Draw();
  }
  c1->Update();

  // Wulz frame
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Wulz Frame");
  title->SetFillColor(10);
  title->Draw();
  t.DrawText(.9, .02, TString::Format("- %i -", ++page));
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);  cosW->Draw();  fitCosTheta(std::cout, cosW);
  pad[page]->cd(3);  AMassW->Draw();
  pad[page]->cd(4);  ARapW->Draw();
  pad[page]->cd(5);  FPseudW->Draw();
  pad[page]->cd(6);  BPseudW->Draw();
  c1->Update();

  ps->Close();

  delete title;
  delete c1;
  //  for (int i_page = 0; i_page < NUM_PAGES; ++i_page)
  //    delete pad[i_page]; // why does this cause a segfault?
  delete ps;
}

DEFINE_FWK_MODULE(AsymmetryFrameAnalysis);

#if 0
  TH2F* rec_rap_vs_gen_rap;
  TH1F* cosCSRes;
  TProfile* cosCS_rec_diffsq_vs_gen;

  // Resolution histograms
  cosCSRes = fs->make<TH1F>("cosCSRes", "Rec cos CS - Gen cos CS", 100, -0.005, 0.005);
  cosCS_rec_diffsq_vs_gen = fs->make<TProfile>("cosCS_rec_diffsq_vs_gen", "Rec cos theta CS diffsq vs Gen cos theta CS", 20, -1., 1., 0., 0.0625);
  rec_rap_vs_gen_rap = fs->make<TH2F>("rec_rap_vs_gen_rap","Rec Rap vs Gen Rap", 50, -4., 4., 50, -4., 4.);

  // A few resolution plots
  if (i_rec == 1 && res_okay) {
    const double res = cos_cs[i_rec][i_dil] - cos_cs[0][i_dil];
    cosCSRes->Fill(res);
    cosCS_rec_diffsq_vs_gen->Fill(cos_cs[0][i_dil], res*res);
    rec_rap_vs_gen_rap->Fill(gen_dileptons.at(i_dil).rapidity(), rap);
  }

  // Resolution of cos(theta) CS
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  tit = " Cos theta CS Resolution";
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
  title->SetFillColor(10);
  title->Draw();
  t.DrawText(.9, .02, TString::Format("- %i -", ++page));
  pad[page]->Draw();
  pad[page]->Divide(2,2);
  // pad[page]->cd(1);  rec_rap_vs_gen_rap->Draw();
  pad[page]->cd(1);  cosCSRes->Draw();  cosCSRes->Fit("gaus","Q");
  Stat_t f_bin;
  int nbins = cosCS_rec_diffsq_vs_gen->GetNbinsX();
  TH1F* cosCS3_diffsq_sqrt
    = fs->make<TH1F>("cosCS3_diffsq_sqrt",
		     "Sqrt(Var(L3-Gen cos theta CS)) vs Gen cos theta CS", nbins,
		     cosCS_rec_diffsq_vs_gen->GetXaxis()->GetXmin(),
		     cosCS_rec_diffsq_vs_gen->GetXaxis()->GetXmax());
  for (int ibin = 1; ibin <= nbins; ibin++) {
    f_bin = cosCS_rec_diffsq_vs_gen->GetBinContent(ibin);
    if (f_bin > 0.) {f_bin = sqrt(f_bin);}
    else            {f_bin = 0.;}
    cosCS3_diffsq_sqrt->SetBinContent(ibin, f_bin);
  }
  pad[page]->cd(4);  cosCS3_diffsq_sqrt->Draw();
  c1->Update();
  delete cosCS3_diffsq_sqrt;
#endif
