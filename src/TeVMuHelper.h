#ifndef TEVMUHELPER_H
#define TEVMUHELPER_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

class TeVMuHelper {
public:
  enum CutResult {PASS, PT=0x01, ISO=0x02, CHI2DOF=0x04, D0=0x08, NSIHITS=0x10}; // To go in bitfields.

  // Default values are what was used in the 2007 AN (only applying pT
  // and isolation cuts).
  TeVMuHelper() :
    ptCut(20), isoCut(10), chi2dofCut(5),
    d0Cut(0.25), nSiHitsCut(7),
    cutMask(PT | ISO) {}

  void setCutMask(const unsigned mask) {
    cutMask = mask;
  }

  unsigned leptonIsCut(const reco::Candidate& lepton) const {
    unsigned retval = 0;

    // pT cut.
    if ((cutMask & PT) && lepton.pt() <= ptCut)
      retval |= PT;

    // Try to cast to reco::Muon.
    const reco::Muon* muon = toConcretePtr<reco::Muon>(lepton);
    if (muon != 0) {
      // Sum pT in cone of dR=0.3 cut
      if ((cutMask & ISO) && muon->isIsolationValid() &&
	  muon->isolationR03().sumPt > isoCut)
	retval |= ISO;
      
      const reco::TrackRef& tk = muon->combinedMuon();
      if (tk.isNonnull()) {
	// Cut on chi2/dof.
	if ((cutMask & CHI2DOF) && tk->normalizedChi2() > chi2dofCut)
	  retval |= CHI2DOF;

	// Cut on d0.
	if ((cutMask & D0) && fabs(tk->d0()) > d0Cut) // 2.5 mm
	  retval |= D0;
      }

      const reco::TrackRef& tktk = muon->track();
      // Cut on number of silicon hits for the tracker track.
      if ((cutMask & NSIHITS) && tktk.isNonnull() &&
	  tktk->numberOfValidHits() < nSiHitsCut)
	retval |= NSIHITS;
    }
  
    return retval & cutMask;
  }

private:
  double ptCut;
  double isoCut;
  double chi2dofCut;
  double d0Cut;
  int nSiHitsCut;

  unsigned cutMask;
};

#endif
