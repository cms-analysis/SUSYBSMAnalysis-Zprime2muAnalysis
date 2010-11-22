//#define PTDRSTYLE
#ifdef PTDRSTYLE
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/tdrstyle.h"
#else
#include "TH1.h"
#include "TStyle.h"
#include "TROOT.h"
#endif

#include <ostream>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

using namespace std;

////////////////////////////////////////////////////////////////////
// Pretty-printing.
////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream& out, const TLorentzVector& vect) {
  // ROOT makes error, CMSSW makes exception... grr...
  float eta = vect.Pt() > 0 ? vect.Eta() : 10e10;
  out << "cartesian: (" << vect.X() << "," << vect.Y() << "," << vect.Z() << "," << vect.E() << ") pt,eta,phi,m: (" << vect.Pt() << ", " << eta << ", " << vect.Phi() << ", " << vect.M() << ")";
  return out;
}

ostream& operator<<(ostream& out, const reco::Candidate* par) {
  out << "GenParticle: " << " id:" << par->pdgId()
      << " p4: " << par->p4() << " status:" << par->status()
      << " vertex:" << par->vertex();
  return out;
}

////////////////////////////////////////////////////////////////////
// ROOT helpers.
////////////////////////////////////////////////////////////////////

TLorentzVector makeTLV(const reco::Candidate::LorentzVector& lv) {
  TLorentzVector v;
  v.SetPtEtaPhiM(lv.pt(), lv.eta(), lv.phi(), lv.mass());
  return v;
}

char nameHistChar(Long_t i) {
  if (i < 10) return char(48 + i); // '0', '1', ...
  else        return char(55 + i); // 'A', 'B', ...
}

TString nameHist(const char* name, Long_t i, Long_t j, Long_t k, Long_t l, Long_t m) {
  TString s = TString(name) + nameHistChar(i);
  if (j >= 0) s += nameHistChar(j);
  if (k >= 0) s += nameHistChar(k);
  if (l >= 0) s += nameHistChar(l);
  if (m >= 0) s += nameHistChar(m);
  return s;
}

void InitROOT(bool date_histograms) {
#ifdef PTDRSTYLE
  TH1::AddDirectory(false);
  setTDRStyle();
#else
  if (getenv("ZP2MU_NO_DATE_HISTOGRAMS")) // magic for JMT
    gStyle->SetOptDate(0);
  gROOT->SetStyle("Plain");
  gStyle->SetFillColor(0);
  gStyle->SetOptDate();
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetMarkerSize(.1);
  gStyle->SetMarkerStyle(8);
  gStyle->SetGridStyle(3);
  gStyle->SetPaperSize(TStyle::kA4);
  gStyle->SetStatW(0.25); // width of statistics box; default is 0.19
  gStyle->SetStatFormat("6.4g"); // leave default format for now
  gStyle->SetTitleFont(52,"XY"); // italic font for axis
  gStyle->SetLabelFont(52,"XY"); // italic font for axis labels
  gStyle->SetStatFont(52);       // italic font for stat. box
#endif
}

////////////////////////////////////////////////////////////////////
// Lepton/dilepton acceptance
////////////////////////////////////////////////////////////////////

WhereLepton whereIsLepton(const reco::CandidateBaseRef& lepton,
			  const bool isElectron) {
  double abseta = fabs(lepton->eta());
  if (isElectron) {
    // Electron classification should be from the eta of its supercluster.
    const reco::RecoCandidate* el = toConcretePtr<reco::RecoCandidate>(lepton);
    if (el != 0) {
      const reco::SuperClusterRef& sc = el->superCluster();
      if (sc.isNonnull())
	abseta = fabs(sc->eta());
    }
    // JMTBAD are these numbers for electrons useful?
    if      (abseta < 1.4442) return W_BARREL;
    else if (abseta < 1.560)  return W_CRACK;
    else if (abseta < 2.5)    return W_ENDCAP;
    else                      return W_OUTSIDE;
  }
  else {
    if      (abseta < 0.85)  return W_BARREL;
    else if (abseta < 1.2)  return W_OVERLAP;
    else if (abseta < 2.45) return W_ENDCAP;
    else                    return W_OUTSIDE;
  }
}

WhereDilepton whereIsDilepton(const reco::CompositeCandidate& dil,
			      const bool areElectrons) {
  WhereLepton wlep1 = whereIsLepton(dileptonDaughter(dil, 0), areElectrons);
  WhereLepton wlep2 = whereIsLepton(dileptonDaughter(dil, 1), areElectrons);
  // get a unique ordering so testing below is easier
  if (int(wlep1) > int(wlep2)) {
    WhereLepton tmp = wlep1;
    wlep1 = wlep2;
    wlep2 = tmp;
  }

  // For electrons, OVERLAP here means CRACK.

  if      (wlep1 == W_BARREL  && wlep2 == W_BARREL)   return W_BARRELBARREL;
  else if (wlep1 == W_BARREL  && wlep2 == W_OVERLAP)  return W_BARRELOVERLAP;
  else if (wlep1 == W_BARREL  && wlep2 == W_ENDCAP)   return W_BARRELENDCAP;
  else if (wlep1 == W_BARREL  && wlep2 == W_OUTSIDE)  return W_BARRELOUTSIDE;
  else if (wlep1 == W_OVERLAP && wlep2 == W_OVERLAP)  return W_OVERLAPOVERLAP;
  else if (wlep1 == W_OVERLAP && wlep2 == W_ENDCAP)   return W_OVERLAPENDCAP;
  else if (wlep1 == W_OVERLAP && wlep2 == W_OUTSIDE)  return W_OVERLAPOUTSIDE;
  else if (wlep1 == W_ENDCAP  && wlep2 == W_ENDCAP)   return W_ENDCAPENDCAP;
  else if (wlep1 == W_ENDCAP  && wlep2 == W_OUTSIDE)  return W_ENDCAPOUTSIDE;
  else //if (wlep1 == W_OUTSIDE && wlep2 == W_OUTSIDE)
    return W_OUTSIDEOUTSIDE;
}
