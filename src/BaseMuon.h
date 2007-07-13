/**

\class BaseMuon

Abstract bass class for a generic muon.

\author S. Valuev, 20 August 2003

*/
#ifndef ZP2MUBASEMUON_H
#define ZP2MUBASEMUON_H

//---------------------
//-- Constants file  --
//---------------------

//--------------------------------------
//-- Collaborating Class Declarations --
//--------------------------------------

#include "TMath.h" // for Pi()

//------------------------
//-- Base Class Headers --
//------------------------

//-----------------
//-- C++ Headers --
//-----------------

#include <iosfwd>

//---------------------
//-- Class Interface --
//---------------------
namespace zp2mu {
  class BaseMuon
  {
  public:

    /** constructors */
    BaseMuon();

    BaseMuon(const bool valid, const int id,     const int level,
	     const int charge, const double phi, const double eta,
	     const double pt,  const double p);

    /** destructor */
    // ~BaseMuon();

    // Access functions:
    /** returns valid flag. */
    inline bool isValid()   const {return theValid;}

    /** returns muon id. */
    inline int id()         const {return theId;}

    /** returns reconstruction level of this object. */
    inline int recLevel()   const {return theRecLevel;}

    /** returns muon charge. */
    inline int charge() const {return theCharge;}

    /** returns muon phi value. */
    inline double phi() const {return thePhi;}

    /** returns muon eta value. */
    inline double eta() const {return theEta;}

    /**return muon px value. */
    inline double px() const {return thePt*cos(thePhi);}

    /**return muon py value. */
    inline double py() const {return thePt*sin(thePhi);}

    /** returns muon longitudinal momentum. */
    inline double pz()  const {
      double abs_val = sqrt((theP*theP) - (thePt*thePt));
      if (theEta < 0.) abs_val *= -1.;
      return abs_val;
    }

    /** returns muon pT value. */
    inline double pt()  const {return thePt;}

    /** returns muon momentum value. */
    inline double p()   const {return theP;}

    /** returns muon energy. */
    inline double energy() const {
      const double mumass = 0.10566; //GeV/c^2
      return sqrt((theP*theP) + (mumass*mumass));
    }

    // Set functions:
    /** sets everything to zero */
    virtual void clear();

    /** sets valid flag. */
    inline void setValid(const bool theValue)   {theValid = theValue;}

    /** sets muon id. */
    inline void setId(const int theValue)       {theId = theValue;}

    /** sets reconstruction level. */
    inline void setRecLevel(const int theValue) {theRecLevel = theValue;}

    /** sets charge. */
    inline void setCharge(const int theValue)   {theCharge = theValue;}

    /** sets phi. */
    void setPhi(const int level, const double theValue);

    /** sets eta. */
    inline void setEta(const double theValue)   {theEta = theValue;}

    /** sets pT. */
    inline void setPt(const double theValue)    {thePt = theValue;}

    /** sets p. */
    inline void setP(const double theValue)     {theP = theValue;}

    /** sets class members from separate variables. */
    virtual void fill(const bool valid, const int id,     const int level,
		      const int charge, const double phi, const double eta,
		      const double pt,  const double p);

    // Overloaded operators:
    bool operator >  (const BaseMuon&) const;
    bool operator <  (const BaseMuon&) const;
    bool operator == (const BaseMuon&) const;
    bool operator != (const BaseMuon&) const;

    friend std::ostream& operator << (std::ostream&, const BaseMuon&);

  protected:
    static const double l1_phi_correction;

    bool   theValid;

    int    theId;
    int    theRecLevel;
    int    theCharge;

    double thePhi;
    double theEta;
    double thePt;
    double theP;
  };
}

#endif // ZP2MUBASEMUON_H
