/**

\class DiMuon

Main class for a dimuon, either generated or reconstructed.

\author S. Valuev, 26 August 2003

*/
#ifndef ZP2MUDIMUON_H
#define ZP2MUDIMUON_H

//---------------------
//-- Constants file  --
//---------------------

//--------------------------------------
//-- Collaborating Class Declarations --
//--------------------------------------

#include "Muon.h"
#include "TLorentzVector.h"

//------------------------
//-- Base Class Headers --
//------------------------

//-----------------
//-- C++ Headers --
//-----------------

//---------------------
//-- Class Interface --
//---------------------
namespace zp2mu {
  class DiMuon
  {
  public:
    /** constructor */
    DiMuon();

    /** destructor */
    // ~DiMuon();

    // Access functions:
    /** returns valid flag. */
    inline bool isValid() const        {return theValid;}

    /** returns dimuon id. */
    inline int id() const              {return theId;}

    /** returns reconstruction level of this object. */
    inline int recLevel() const        {return theRecLevel;}

    /** returns mu- for this dimuon. */
    inline const Muon& muMinus() const {return m_muMinus;}

    /** returns mu+ for this dimuon. */
    inline const Muon& muPlus() const  {return m_muPlus;}

    /** returns TLorentzVector for this dimuon. */
    inline const TLorentzVector& dimuV() const {return theDilV;}

    /** returns TLorentzVector for this dimuon, including photons. */
    inline const TLorentzVector& resV()  const {return theResV;}

    // Set functions:
    /** sets everything to zero */
    void clear();

    /** sets valid flag. */
    inline void setValid(const bool theValue)    {theValid = theValue;}

    /** sets muon id. */
    inline void setId(const int theValue)        {theId = theValue;}

    /** sets reconstruction level. */
    inline void setRecLevel(const int theValue)  {theRecLevel = theValue;}

    /** sets mu- for this dimuon. */
    inline void setMuMinus(const Muon& mum)      {m_muMinus = mum;}

    /** sets mu+ for this dimuon. */
    inline void setMuPlus(const Muon& mup)       {m_muPlus = mup;}

    /** sets TLorentzVector for this dimuon. */
    inline void setDimuV(const TLorentzVector v) {theDilV = v;}

    /** sets TLorentzVector for this dimuon, including photons. */
    inline void setResV(const TLorentzVector v)  {theResV = v;}

    /** sets basic class members from separate variables. */
    void fill(const bool valid, const int id, const int level);

    // Overloaded operators:
    DiMuon& operator = (const DiMuon&);
    bool operator >  (const DiMuon&) const;
    bool operator <  (const DiMuon&) const;
    bool operator == (const DiMuon&) const;
    bool operator != (const DiMuon&) const;

    friend std::ostream& operator << (std::ostream&, const DiMuon&);

  private:
    bool   theValid;

    int    theId;
    int    theRecLevel;

    // pointers to mu- and mu+
    Muon   m_muMinus;
    Muon   m_muPlus;

    // dimuon kinematics
    TLorentzVector theDilV; // reconstructed from two leptons
    TLorentzVector theResV; // includes photons
  };
}

#endif // ZP2MUDIMUON_H
