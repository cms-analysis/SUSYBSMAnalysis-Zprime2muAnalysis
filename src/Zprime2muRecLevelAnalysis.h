#ifndef ZP2MURLANALYSIS_H
#define ZP2MURLANALYSIS_H

#include <iosfwd>
#include <vector>
#include <string>

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

class Zprime2muRecLevelAnalysis : public Zprime2muAnalysis {
 public:
  explicit Zprime2muRecLevelAnalysis(const edm::ParameterSet& config);
  virtual ~Zprime2muRecLevelAnalysis() {}

  virtual void beginJob(const edm::EventSetup&) {}
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}

  // Dump the event, printing out the specified information at each
  // level of lepton reconstruction.
  virtual void dumpEvent(const bool trigOnly=false) const;

 protected:
  // Object to help with per-rec-level information.
  RecLevelHelper recLevelHelper;

  // The lepton collections, stored as CandidateBaseRefs.
  reco::CandidateBaseRefVector allLeptons[MAX_LEVELS];

  // The dilepton collections, stored as CompositeCandidates.
  reco::CompositeCandidateCollection allDileptons[MAX_LEVELS];

  // The collections of dileptons + brem photons.
  reco::CompositeCandidateCollection allResonances[MAX_LEVELS];
  
  ////////////////////////////////////////////////////////////////////
  // General utility
  ////////////////////////////////////////////////////////////////////
  
  // Return whether the rec level might be empty, based on the
  // useGen/Sim/Reco/etc. flags. (Useful especially in
  // Zprime2muResolution, to suppress the creation of some pages in
  // the huge postscript file.)
  bool skipRecLevel(const int level) const;

  ////////////////////////////////////////////////////////////////////
  // Print-outs
  ////////////////////////////////////////////////////////////////////

  // Print out all the relevant information about the lepton; but this
  // method is just as useful as documentation on how to access this
  // information.
  void dumpLepton(std::ostream& output, reco::CandidateBaseRef cand) const;

  // Print out all the relevant information about the dilepton, and
  // call dumpLepton on each of its daughters if dumpLeptons is true.
  void dumpDilepton(std::ostream& output,
		    const reco::CompositeCandidate& cand,
		    bool dumpLeptons=false) const;

};

#endif