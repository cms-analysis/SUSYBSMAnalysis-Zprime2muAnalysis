//-------------------------------------------------
//
//   Class: BaseMuon
//
//   Description: Implements the base class for muon.
//
//   Author List: S. Valuev
//
//--------------------------------------------------

#include <iomanip>

#include "BaseMuon.h"

namespace zp2mu {
  const double BaseMuon::l1_phi_correction = 0.0218;

  BaseMuon::BaseMuon() {
    clear();
  }

  BaseMuon::BaseMuon(const bool valid, const int id,     const int level,
		     const int charge, const double phi, const double eta,
		     const double pt,  const double p) {
    clear();
    fill(valid, id, level, charge, phi, eta, pt, p);
  }

  void BaseMuon::clear() {
    theValid    = false;
    theId       = 0;
    theRecLevel = 0;
    theCharge   = 0;
    thePhi      = 0.;
    theEta      = 0.;
    thePt       = 0.;
    theP        = 0.;
  }

  void BaseMuon::fill(const bool valid, const int id,     const int level,
		      const int charge, const double phi, const double eta,
		      const double pt,  const double p) {
    setValid(valid);
    setId(id);
    setRecLevel(level);
    setCharge(charge);
    setPhi(level, phi);
    setEta(eta);
    setPt(pt);
    setP(p);
  }

  /** sets phi. */
  void BaseMuon::setPhi(const int level, const double theValue) {
    double twoPi = 2.*TMath::Pi();

    // L1 takes phi from the edge of the bin, not the center.
    // This needs to be corrected by 1.25 degrees or 0.0218 rad.
    if (level == 1) thePhi = theValue + l1_phi_correction;
    else            thePhi = theValue;

    if      (thePhi < 0.)    thePhi += twoPi;
    else if (thePhi > twoPi) thePhi -= twoPi;
  }

  /** Overloaded operators. */
  bool BaseMuon::operator > (const BaseMuon& rhs) const {
    return (p() > rhs.p());
  }

  bool BaseMuon::operator < (const BaseMuon& rhs) const {
    return (p() < rhs.p());
  }

  bool BaseMuon::operator == (const BaseMuon& rhs) const {
    return (isValid()   && rhs.isValid()   &&
	    (recLevel() == rhs.recLevel()) &&
	    (id()       == rhs.id()));
  }

  bool BaseMuon::operator != (const BaseMuon& rhs) const {
    return (isValid() == false || rhs.isValid() == false ||
	    (recLevel() != rhs.recLevel()) ||
	    (id()       != rhs.id()));
  }

  std::ostream& operator << (std::ostream& output,  const BaseMuon& rhs) {
    using namespace std;
    if (rhs.isValid()) {
      output << " Id: "     << setw(3) << rhs.id()
	     << " Level: "  << setw(3) << rhs.recLevel()
	     << " Charge: " << setw(3) << rhs.charge()
	     << " Phi: "    << setw(9) << setprecision(5) << rhs.phi()
	     << " Eta: "    << setw(9) << setprecision(5) << rhs.eta()
	     << " Pt: "     << setw(9) << setprecision(5) << rhs.pt()
	     << " P: "      << setw(9) << setprecision(5) << rhs.p() << endl;
    }
    else {
      output << "+++ Not a valid BaseMuon! +++" << endl;
    }
    return output;
  }
}
