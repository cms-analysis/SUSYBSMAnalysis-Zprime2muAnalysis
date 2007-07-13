//-------------------------------------------------
//
//   Class: DiMuon
//
//   Description: Implements the class for any dimuon.
//
//   Author List: S. Valuev
//
//--------------------------------------------------

#include "DiMuon.h"

#include <iomanip>
#include <string>

namespace zp2mu {
  DiMuon::DiMuon() {
    clear();
  }

  void DiMuon::clear() {
    theValid    = false;
    theId       = 0;
    theRecLevel = 0;

    m_muMinus.clear();
    m_muPlus.clear();

    theDilV.Clear();
    theResV.Clear();
  }

  void DiMuon::fill(const bool valid, const int id, const int level) {
    setValid(valid);
    setId(id);
    setRecLevel(level);
  }

  /** overloaded assignment; used to sort DiMuons. */
  DiMuon& DiMuon::operator = (const DiMuon& rhs) {
    if (this != &rhs) {
      theValid    = rhs.theValid;
      theId       = rhs.theId;
      theRecLevel = rhs.theRecLevel;
      m_muMinus   = rhs.m_muMinus;
      m_muPlus    = rhs.m_muPlus;
      theDilV     = rhs.theDilV;
      theResV     = rhs.theResV;
    }
    return *this;
  }

  bool DiMuon::operator > (const DiMuon& rhs) const {
    return (dimuV().M() > rhs.dimuV().M());
  }

  bool DiMuon::operator < (const DiMuon& rhs) const {
    return (dimuV().M() < rhs.dimuV().M());
  }

  bool DiMuon::operator == (const DiMuon& rhs) const {
    return (recLevel()     == rhs.recLevel()     &&
	    muMinus().id() == rhs.muMinus().id() &&
	    muPlus().id()  == rhs.muPlus().id());
  }

  bool DiMuon::operator != (const DiMuon& rhs) const {
    return (recLevel()     != rhs.recLevel()     ||
	    muMinus().id() != rhs.muMinus().id() ||
	    muPlus().id()  != rhs.muPlus().id());
  }

  std::ostream& operator << (std::ostream& output, const DiMuon& rhs) {
    using namespace std;
    string str_level[] = {"Gn", "L1", "L2", "L3", "GR", "TK", "FS", "PR"};
    string str_lvl = " ";
    if (rhs.isValid()) {
      int level = rhs.recLevel();
      if (level >= 0 && level < REC_LEVELS) str_lvl = str_level[level];
      else if (level == 99)                 str_lvl = "BEST";
      TLorentzVector diV = rhs.dimuV();
      output << str_lvl
	     << " #"     << setw(3) << rhs.id()
	     << " Phi: " << setw(7) << setprecision(4) << diV.Phi()
	     << " Eta: " << setw(7) << setprecision(4) << diV.Eta()
	     << " Pt: "  << setw(7) << setprecision(5) << diV.Pt()
	     << " P: "   << setw(7) << setprecision(5) << diV.P()
	     << " M: "   << setw(7) << setprecision(5) << diV.M() << endl;
      output << "  mu+: " << endl << rhs.muPlus();
      output << "  mu-: " << endl << rhs.muMinus();
    }
    else {
      output << "+++ Not a valid DiMuon! +++" << endl;
    }
    return output;
  }
}
