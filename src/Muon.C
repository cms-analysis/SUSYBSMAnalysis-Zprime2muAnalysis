//-------------------------------------------------
//
//   Class: Muon
//
//   Description: Implements the class for any muon.
//
//   Author List: S. Valuev
//
//--------------------------------------------------

#include "Muon.h"

#include <iomanip>
#include <string>

namespace zp2mu {
  Muon::Muon() {
    clear();
  }

  void Muon::clear() {
    BaseMuon::clear();

    theMotherLine = 0;
    theMotherId   = 0;
    theGrandmaId  = 0;

    theQuality   = 0;

    theNpixHits  = 0;
    theNsilHits  = 0;
    theNrecHits  = 0;
    theSeedIndex = 0;

    for (int i = 0; i < REC_LEVELS; i++) {
      theClosestId[i]  = -999;
      theSameSeedId[i] = -999;
    }

    theDof             = 0;
    theChi2            = 0.;
    theBackChi2        = 0.;
    theTrackerChi2     = 0.;
    theMuonFitChi2     = 0.;
    theTrackerChi2Diff = 0.;
    theMuonFitChi2Diff = 0.;

    theForwardPt     = 0.;
    theTrackerPt     = 0.;
    theMuonFitPt     = 0.;
    theEInvPt        = 0.;
    theEForwardInvPt = 0.;
    theETrackerInvPt = 0.;
    theEMuonFitInvPt = 0.;
    theEInvP         = 0.;

    for (int i = 0; i < 3; i++) {
      theVertexXYZ[i]  = 0.;
      theTrackerXYZ[i] = 0.;
    }

    thePhoton.SetPxPyPzE(0., 0., 0., 0.);
  }

  /** returns id of the closest muon at rec. level = level. */
  int Muon::closestId(const int level) const {
    if (level >= 0 && level < REC_LEVELS) return theClosestId[level];
    else                                  return -999;
  }

  /** returns id of the same seed muon at rec. level = level. */
  int Muon::sameSeedId(const int level) const {
    if (level >= 0 && level < REC_LEVELS) return theSameSeedId[level];
    else                                  return -999;
  }

  /** returns id of the matched muon at rec. level = level.
      The seed info is not available for generated, L1 and L2 muons, so the
      matched muon is a closest muon.  For muons at other levels, the matched
      muon is the same seed muon if available, and the closest muon otherwise. */
  int Muon::matchedId(const int level) const {
    if (level < 0 || level >= REC_LEVELS) return -999;
    else if (level <= 2)                  return theClosestId[level];
    else {
      if (theSameSeedId[level] != -999)   return theSameSeedId[level];
      else                                return theClosestId[level];
    }
  }

  /** returns X, Y or Z position of the vertex. */
  double Muon::vertexXYZ(const int index) const {
    if (index < 3) return theVertexXYZ[index];
    else           return 0.;
  }

  /** returns X, Y or Z position of the vertex. */
  math::XYZPoint Muon::vertex() const {
    return math::XYZPoint(theVertexXYZ[0],theVertexXYZ[1],theVertexXYZ[2]);
  }

  /** returns X, Y or Z position of the track at the inner tracker surface. */
  double Muon::trackerXYZ(const int index) const {
    if (index < 3) return theTrackerXYZ[index];
    else           return 0.;
  }

  /** returns error on Pt. */
  double Muon::errPt() const {
    return (pt()*pt()*errInvPt());
  }

  /** returns error on Pt from the forward fit. */
  double Muon::errForwardPt() const {
    return (forwardPt()*forwardPt()*errForwardInvPt());
  }

  /** returns error on Pt from the tracker only fit. */
  double Muon::errTrackerPt() const {
    return (trackerPt()*trackerPt()*errTrackerInvPt());
  }

  /** returns error on Pt from the fit in the muon chambers only. */
  double Muon::errMuonFitPt() const {
    return (muonFitPt()*muonFitPt()*errMuonFitInvPt());
  }

  /** returns error on P. */
  double Muon::errP() const {
    return (p()*p()*errInvP());
  }

  void Muon::fill(const bool valid, const int id,     const int level,
		  const int charge, const double phi, const double eta,
		  const double pt,  const double p) {
    BaseMuon::fill(valid, id, level, charge, phi, eta, pt, p);
  }

  /** sets id of the closest muon at rec. level = level. */
  void Muon::setClosestId(const int level, const int closest_id) {
    if (level >= 0 && level < REC_LEVELS) theClosestId[level] = closest_id;
  }

  /** sets id of the same seed muon at rec. level = level. */
  void Muon::setSameSeedId(const int level, const int sameseed_id) {
    if (level >= 0 && level < REC_LEVELS) theSameSeedId[level] = sameseed_id;
  }

  /** sets number of pixel, silicon and muon hits in the track. */
  void Muon::setHits(const int pixel, const int silicon, const int total) {
    theNpixHits = pixel;
    theNsilHits = silicon;
    theNrecHits = total;
  }

  /** sets chi-squared of the fit in the tracker only, and the difference
      between the forward and the backward fits. */
  void Muon::setTrackerChi2(const double trackerchi2,
			    const double trackerchi2diff) {
    theTrackerChi2     = trackerchi2;
    theTrackerChi2Diff = trackerchi2diff;
  }

  /** sets chi-squared of the fit in the muon chambers only, and the difference
      between the forward and the backward fits. */
  void Muon::setMuonFitChi2(const double muonfitchi2,
			    const double muonfitchi2diff) {
    theMuonFitChi2     = muonfitchi2;
    theMuonFitChi2Diff = muonfitchi2diff;
  }

  /** sets array containing X, Y and Z position of the vertex. */
  void Muon::setVertexXYZ(const double xv, const double yv, const double zv) {
    theVertexXYZ[0] = xv;
    theVertexXYZ[1] = yv;
    theVertexXYZ[2] = zv;
  }

  /** sets array containing X, Y and Z position of the track at the inner
      surface of the tracker. */
  void Muon::setTrackerXYZ(const double xt, const double yt, const double zt) {
    theTrackerXYZ[0] = xt;
    theTrackerXYZ[1] = yt;
    theTrackerXYZ[2] = zt;
  }

  /* overloaded assignment; used to sort Muons. */
  Muon& Muon::operator = (const Muon& rhs) {
    if (this != &rhs) {
      theValid    = rhs.theValid;
      theId       = rhs.theId;
      theRecLevel = rhs.theRecLevel;
      theCharge   = rhs.theCharge;
      thePhi      = rhs.thePhi;
      theEta      = rhs.theEta;
      thePt       = rhs.thePt;
      theP        = rhs.theP;

      theMotherLine = rhs.theMotherLine;
      theMotherId   = rhs.theMotherId;
      theGrandmaId  = rhs.theGrandmaId;
      theQuality    = rhs.theQuality;

      theNpixHits  = rhs.theNpixHits;
      theNsilHits  = rhs.theNsilHits;
      theNrecHits  = rhs.theNrecHits;
      theSeedIndex = rhs.theSeedIndex;

      for (int i = 0; i < REC_LEVELS; i++) {
	theClosestId[i]  = rhs.theClosestId[i];
	theSameSeedId[i] = rhs.theSameSeedId[i];
      }

      theDof             = rhs.theDof;
      theChi2            = rhs.theChi2;
      theBackChi2        = rhs.theBackChi2;
      theTrackerChi2     = rhs.theTrackerChi2;
      theMuonFitChi2     = rhs.theMuonFitChi2;
      theTrackerChi2Diff = rhs.theTrackerChi2Diff;
      theMuonFitChi2Diff = rhs.theMuonFitChi2Diff;

      theForwardPt     = rhs.theForwardPt;
      theTrackerPt     = rhs.theTrackerPt;
      theMuonFitPt     = rhs.theMuonFitPt;
      theEInvPt        = rhs.theEInvPt;
      theEForwardInvPt = rhs.theEForwardInvPt;
      theETrackerInvPt = rhs.theETrackerInvPt;
      theEMuonFitInvPt = rhs.theEMuonFitInvPt;
      theEInvP         = rhs.theEInvP;

      for (int i = 0; i < 3; i++) {
	theVertexXYZ[i]  = rhs.theVertexXYZ[i];
	theTrackerXYZ[i] = rhs.theTrackerXYZ[i];
      }

      thePhoton = rhs.thePhoton;
    }
    return *this;
  }

  std::ostream& operator << (std::ostream& output, const Muon& rhs) {
    using namespace std;
    string str_level[] = {"Gn", "L1", "L2", "L3", "GR", "TK", "FS", "PR"};
    string str_lvl = " ";
    if (rhs.isValid()) {
      int level = rhs.recLevel();
      if (level >= 0 && level < REC_LEVELS) str_lvl = str_level[level];
      else if (level == 99)                 str_lvl = "BEST";
      output << str_lvl
	     << " #"          << setw(1) << rhs.id()
	     << "  Q: "       << setw(2) << rhs.charge();
      if (level == 0) {
	output << " Origin: " << setw(4) << rhs.genMotherId()
	       << "/"         << setw(4) << rhs.genMotherLine()
	       << " Phi: "    << setw(7) << setprecision(4) << rhs.phi()
	       << " Eta: "    << setw(7) << setprecision(4) << rhs.eta()
	       << " Pt: "     << setw(7) << setprecision(5) << rhs.pt()
	       << " P: "      << setw(7) << setprecision(5) << rhs.p();
      }
      else if (level == 1) {
	output << " Quality: "<< setw(2) << rhs.l1Quality()
	       << " Phi: "    << setw(7) << setprecision(4) << rhs.phi()
	       << " Eta: "    << setw(7) << setprecision(4) << rhs.eta()
	       << " Pt: "     << setw(7) << setprecision(5) << rhs.pt()
	       << " P: "      << setw(7) << setprecision(5) << rhs.p();
      }
      else if (level == 2) {
	output << " Phi: "      << setw(8) << setprecision(4) << rhs.phi()
	       << "      Eta: " << setw(8) << setprecision(4) << rhs.eta()
	       << endl;
	output << "   Muhits: "    << setw(3) << rhs.nmuHits()
	       << "   Chi2/Ndof: " << setw(8) << setprecision(4) << rhs.chi2()
	       << "/"              << setw(2) << rhs.dof() << endl;
	output << "   P:     " << setw(7) << setprecision(5) << rhs.p()
	       << " +/- " << rhs.errP() << "   "
	       << "   Pt:    " << setw(7) << setprecision(5) << rhs.pt()
	       << " +/- " << rhs.errPt();
      }
      else {
	output << " Phi: "      << setw(8) << setprecision(4) << rhs.phi()
	       << "      Eta: " << setw(8) << setprecision(4) << rhs.eta()
	       << endl;
	output << "   Pixhits: " << setw(3) << rhs.npixHits()
	       << " Silhits: "   << setw(3) << rhs.nsilHits()
	       << " Rechits: "   << setw(3) << rhs.nrecHits()
	       << " Seed: "      << setw(3) << rhs.seedIndex() << endl;
	output << "   Closest  :";
	for (int i = 0; i < REC_LEVELS; i++) {
	  output << setw(5) << rhs.closestId(i)  << "(" << str_level[i] << ")";
	}
	output << endl;
	output << "   Same-seed:";
	for (int i = 0; i < REC_LEVELS; i++) {
	  output << setw(5) << rhs.sameSeedId(i) << "(" << str_level[i] << ")";
	}
	output << endl;
	output << "   P:          " << setw(7) << setprecision(5) << rhs.p()
	       << " +/- "    << rhs.errP() << endl;
	output << "   Pt:         " << setw(7) << rhs.pt()
	       << " +/- "    << rhs.errPt()
	       << "   Chi2/Ndof: "  << setw(11) << rhs.chi2()
	       << "/"        << setw(2) << rhs.dof() << endl;
	output << "   Forward Pt: " << setw(7) << rhs.forwardPt()
	       << " +/- "    << rhs.errForwardPt()
	       << "   Back Chi2: "  << setw(11) << rhs.backChi2() << endl;
	output << "   Tracker Pt: " << setw(7) << rhs.trackerPt()
	       << " +/- "    << rhs.errTrackerPt()
	       << "   Tracker Chi2: "  << setw(8) << rhs.trackerChi2() << endl;
	output << "   MuonFit Pt: " << setw(7) << rhs.muonFitPt()
	       << " +/- "    << rhs.errMuonFitPt()
	       << "   MuonFit Chi2: "  << setw(8) << rhs.muonFitChi2() << endl;
	output << "   Tracker Chi2 diff.: " << setw(7) << rhs.trackerChi2Diff()
	       << "   MuonFit Chi2 diff.: " << setw(7) << rhs.muonFitChi2Diff()
	       << endl;
	output << "   Vertex position: " << setw(11) << rhs.vertexXYZ(0)
	       << " "                    << setw(11) << rhs.vertexXYZ(1)
	       << " "                    << setw(11) << rhs.vertexXYZ(2) << endl;
	output << "   Track. position: " << setw(11) << rhs.trackerXYZ(0)
	       << " "                    << setw(11) << rhs.trackerXYZ(1)
	       << " "                    << setw(11) << rhs.trackerXYZ(2)<< endl;
	output << "   Closest photon: (" << rhs.closestPhoton().Px()
	       << ", "                   << rhs.closestPhoton().Py()
	       << ", "                   << rhs.closestPhoton().Pz()
	       << "; "                   << rhs.closestPhoton().E() << ")";
      }
      output << endl;
    }
    else {
      output << "+++ Not a valid Muon! +++" << endl;
    }
    return output;
  }
}
