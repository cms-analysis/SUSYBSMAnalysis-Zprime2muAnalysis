//#define PTDRSTYLE
#ifdef PTDRSTYLE
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/tdrstyle.h"
#else
#include "TH1.h"
#include "TStyle.h"
#include "TROOT.h"
#endif

#include <ostream>
#include <sstream>
#include <string>

#include "TLorentzVector.h"

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
  out << "(" << vect.X() << "," << vect.Y() << "," << vect.Z() << "," << vect.E() << ")";
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

string nameHist(const char* s, int i, int j, int k, int l, int m,
		bool extended) {
  ostringstream x;
  x << s;
  if (extended) x << "_";
  x << i;
  if (j >= 0) {
    if (extended) x << "_";
    x << j;
  }
  if (k >= 0) {
    if (extended) x << "_";
    x << k;
  }
  if (l >= 0) {
    if (extended) x << "_";
    x << l;
  }
  if (m >= 0) {
    if (extended) x << "_";
    x << m;
  }
  return x.str();
}

void InitROOT() {
#ifdef PTDRSTYLE
  TH1::AddDirectory(false);
  setTDRStyle();
#else
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
    if      (abseta < 0.9) return W_BARREL;
    else if (abseta < 1.2) return W_OVERLAP;
    else if (abseta < 2.4) return W_ENDCAP;
    else                   return W_OUTSIDE;
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
