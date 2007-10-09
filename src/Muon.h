/**

\class Muon

Main class for a muon, either generated or reconstructed.

\author S. Valuev, 25 August 2003

*/

#ifndef ZP2MUMUON_H
#define ZP2MUMUON_H

//--------------------------------------
//-- Collaborating Class Declarations --
//--------------------------------------

//------------------------
//-- Base Class Headers --
//------------------------

#include "BaseMuon.h"

//-----------------
//-- C++ Headers --
//-----------------

#include "TLorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

namespace zp2mu {
  const int REC_LEVELS = 8;

  //---------------------
  //-- Class Interface --
  //---------------------
  class Muon: public BaseMuon
  {
  public:
    /** constructor */
    Muon();

    /** destructor */
    virtual ~Muon(){ }

    // Access functions:
    /** returns line number of muon mother in Pythia list.
	Generated muon only. */
    inline int genMotherLine() const      {return theMotherLine;}

    /** returns pythia id of muon mother (32 for Z', 5000039 for G*).
	Generated muon only. */
    inline int genMotherId() const        {return theMotherId;}

    /** returns pythia id of muon's grandmother (i.e. |id| < 6
	for quarks, id = 21 for gluon.
	Generated muon only.  */
    inline int genGrandmaId() const       {return theGrandmaId;}

    // some tests on the grandma id for ease of code reading
    inline bool grandmaIsqqbar() const    {return ((theGrandmaId > 0 && 
						    theGrandmaId <= 6) || 
						   (theGrandmaId < 0 && 
						    theGrandmaId >= -6));}
    inline bool grandmaIsgg() const       {return theGrandmaId == 21;}
    inline int grandmaType() const {
      if (grandmaIsqqbar()) return 0; 
      else if (grandmaIsgg()) return 1;
      else return -1; 
    }

    /** returns L1 track qulaity. L1 muons only. */
    inline int l1Quality() const          {return theQuality;}

    /** returns number of muon hits. */
    inline int nmuHits() const            {
      return theNrecHits - theNpixHits - theNsilHits;
    }

    /** returns number of pixel detector hits. */
    inline int npixHits() const           {return theNpixHits;}

    /** returns number of silicon tracker hits. */
    inline int nsilHits() const           {return theNsilHits;}

    /** returns total number of hits. */
    inline int nrecHits() const           {return theNrecHits;}

    /** returns index of the seed. */
    inline int seedIndex() const          {return theSeedIndex;}

    /** returns id of the closest muon at rec. level = level. */
    int closestId(const int level) const;

    /** returns id of the same seed muon at rec. level = level. */
    int sameSeedId(const int level) const;

    /** returns id of the matched muon at rec. level = level. */
    int matchedId(const int level) const;

    /** returns number of dof's. */
    inline int dof() const                {return theDof;}

    /** returns chi-squared. */
    inline double chi2() const            {return theChi2;}

    /** returns chi-squared of the backward fit. */
    inline double backChi2() const        {return theBackChi2;}

    /** returns chi-squared of the tracker only fit. */
    inline double trackerChi2() const     {return theTrackerChi2;}

    /** returns chi-squared of the muon only fit. */
    inline double muonFitChi2() const     {return theMuonFitChi2;}

    /** returns difference of chi-squares from the forward and the backward
	fits in the tracker. */
    inline double trackerChi2Diff() const {return theTrackerChi2Diff;}

    /** returns difference of chi-squares from the forward and the backward
	fits in the muon chambers. */
    inline double muonFitChi2Diff() const {return theMuonFitChi2Diff;}

    /** returns Pt of the forward fit. */
    inline double forwardPt() const       {return theForwardPt;}

    /** returns Pt of the fit in the tracker only. */
    inline double trackerPt() const       {return theTrackerPt;}

    /** returns Pt of the fit in the muon chambers only. */
    inline double muonFitPt() const       {return theMuonFitPt;}

    /** returns error on 1/Pt. */
    inline double errInvPt() const        {return theEInvPt;}

    /** returns error on 1/Pt from the forward fit. */
    inline double errForwardInvPt() const {return theEForwardInvPt;}

    /** returns error on 1/Pt from the fit in the tracker only. */
    inline double errTrackerInvPt() const {return theETrackerInvPt;}

    /** returns error on 1/Pt from the fit in the muon chambers only. */
    inline double errMuonFitInvPt() const {return theEMuonFitInvPt;}

    /** returns error on 1/P. */
    inline double errInvP() const         {return theEInvP;}

    /** returns X, Y or Z position of the vertex. */
    double vertexXYZ(const int index) const;

    math::XYZPoint vertex() const;

    /** returns X, Y or Z position of the track at the inner surface of
	the tracker. */
    double trackerXYZ(const int index) const;

    /** returns 4-momentum of the closest photon candidate. */
    inline TLorentzVector closestPhoton() const {return thePhoton;}

    // Derived quantities
    /** returns error on Pt. */
    double errPt() const;

    /** returns error on Pt from the forward fit. */
    double errForwardPt() const;

    /** returns error on Pt from the tracker only fit. */
    double errTrackerPt() const;

    /** returns error on Pt from the fit in the muon chambers only. */
    double errMuonFitPt() const;

    /** returns error on P. */
    double errP() const;

    // Set functions:
    /** sets everything to zero */
    void clear();

    /** sets basic class members from separate variables. */
    void fill(const bool valid, const int id,     const int level,
	      const int charge, const double phi, const double eta,
	      const double pt,  const double p);

    /** sets mother line number.  Generated muon only. */
    inline void setMotherLine(const int theValue)   {theMotherLine = theValue;}

    /** sets mother id. Generated muon only. */
    inline void setMotherId(const int theValue)     {theMotherId = theValue;}

    /** sets grandma id. Generated muon only. */
    inline void setGrandmaId(const int theValue)    {theGrandmaId = theValue;}

    /** sets L1 track quality.  L1 muons only. */
    inline void setQuality(const int theValue)      {theQuality = theValue;}

    /** sets number of pixel, silicon and total hits in the track. */
    void setHits(const int pixel, const int silicon, const int total);

    /** sets seed index of the track. */
    inline void setSeed(const int theValue)         {theSeedIndex = theValue;}

    /** sets id of the closest muon at rec. level = level. */
    void setClosestId(const int level, const int closest_id);

    /** sets id of the same seed muon at rec. level = level. */
    void setSameSeedId(const int level, const int sameseed_id);

    /** sets number of dof's. */
    inline void setDof(const int theValue)          {theDof = theValue;}

    /** sets chi-squared. */
    inline void setChi2(const double theValue)      {theChi2 = theValue;}

    /** sets chi-squared of the backward fit. */
    inline void setBackChi2(const double theValue)  {theBackChi2 = theValue;}

    /** sets chi-squared of the fit in the tracker only, and the difference
	between the forward and the backward fits. */
    void setTrackerChi2(const double trackerchi2, const double trackerchi2diff);

    /** sets chi-squared of the fit in the muon chambers only, and the difference
	between the forward and the backward fits. */
    void setMuonFitChi2(const double muonfitchi2, const double muonfitchi2diff);

    /** sets Pt from the forward fit. */
    inline void setForwardPt(const double theValue) {theForwardPt = theValue;}

    /** sets Pt from the fit in the tracker only. */
    inline void setTrackerPt(const double theValue) {theTrackerPt = theValue;}

    /** sets Pt from the fit in the muon  chambers only. */
    inline void setMuonFitPt(const double theValue) {theMuonFitPt = theValue;}

    /** sets error on 1/Pt. */
    inline void setEInvPt(const double theValue)    {theEInvPt = theValue;}

    /** sets error on 1/Pt from the forward fit. */
    inline void setEForwardInvPt(const double theValue)
    {theEForwardInvPt = theValue;}

    /** sets error on 1/Pt from the fit in the tracker only. */
    inline void setETrackerInvPt(const double theValue)
    {theETrackerInvPt = theValue;}

    /** sets error on 1/Pt from the fit in the muon chambers only. */
    inline void setEMuonFitInvPt(const double theValue)
    {theEMuonFitInvPt = theValue;}

    /** sets error on 1/P. */
    inline void setEInvP(const double theValue)  {theEInvP = theValue;}

    /** sets array containing X, Y and Z position of the vertex. */
    void setVertexXYZ(const double xv, const double yv, const double zv);

    /** sets array containing X, Y and Z position of the track at the inner
	surface of the tracker. */
    void setTrackerXYZ(const double xt, const double yt, const double zt);

    /** sets 4-momentum of the closest photon candidate. */
    inline void setClosestPhoton(const TLorentzVector & ph) {thePhoton = ph;}

    // Overloaded operators:
    Muon& operator = (const Muon&);
    /* 
       bool operator >  (const Muon&) const;
       bool operator <  (const Muon&) const;
       bool operator >= (const Muon&) const;
       bool operator <= (const Muon&) const;
       bool operator == (const Muon&) const;
       bool operator != (const Muon&) const;
    */

    friend std::ostream& operator << (std::ostream&, const Muon&);

  protected:
    int    theMotherLine;  // generated muons only
    int    theMotherId;
    int    theGrandmaId;

    int    theQuality;     // L1 muons only
 
    int    theNpixHits;    // number of pixel detector hits in the track
    int    theNsilHits;    // number of silicon tracker hits in the track
    int    theNrecHits;    // total number of hits in the track
    //  (pixel, silicon and muon)
    int    theSeedIndex;
    // id's of closest and same seed muons
    int    theClosestId[REC_LEVELS];  // at all rec. levels
    int    theSameSeedId[REC_LEVELS];

    int    theDof;
    double theChi2;
    double theBackChi2;    // chi-squared from the backward fit
    double theTrackerChi2;
    double theMuonFitChi2;
    double theTrackerChi2Diff;
    double theMuonFitChi2Diff;

    double theForwardPt;   // three different Pt's
    double theTrackerPt;
    double theMuonFitPt;
    double theEInvPt;      // errors on 1/Pt
    double theEForwardInvPt;
    double theETrackerInvPt;
    double theEMuonFitInvPt;
    double theEInvP;       // error on 1/P

    double theVertexXYZ[3];  // vertex coordinate (cm)
    double theTrackerXYZ[3]; // coordinate at the inner surface of the
    // tracker (cm)

    TLorentzVector thePhoton; // 4-momentum of the closest photon candidate
  };
}

#endif // ZP2MUMUON_H
