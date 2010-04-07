#ifndef ZP2MUANALYSIS_H
#define ZP2MUANALYSIS_H

#include <iosfwd>
#include <vector>
#include <string>

#include "TLorentzVector.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/CutHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneratorUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerDecision.h"

// Verbosity levels.
enum Verbosity { VERBOSITY_NONE, VERBOSITY_SIMPLE,
		 VERBOSITY_LOTS, VERBOSITY_TOOMUCH };

class Zprime2muAnalysis : public edm::EDAnalyzer {
 public:
  explicit Zprime2muAnalysis(const edm::ParameterSet&);
  virtual ~Zprime2muAnalysis() {}

  virtual void beginJob() {}
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}

  // Dump the event, printing out the specified information at each
  // level of lepton reconstruction.
  virtual void dumpEvent(const bool trigOnly=false) const;

 protected:
  ////////////////////////////////////////////////////////////////////
  // Parameters read or determined from the config file:
  ////////////////////////////////////////////////////////////////////

  // verbosity controls the amount of debugging information printed;
  // levels are defined using the VERBOSITY_* codes above
  Verbosity verbosity;

  // Keep this many highest-invariant-mass dileptons.
  unsigned maxDileptons;

  // whether we are looking at electrons instead of muons;
  bool doingElectrons;

  // whether to look at generator-level information;
  bool useGen;

  // whether to look at GEANT tracks;
  bool useSim;

  // whether to look at trigger quantities;
  bool useTrigger;

  // whether to look at reconstructed quantities;
  bool useReco;

  // whether to include the extra muon reconstructors (FMS, PMR, etc);
  bool useOtherMuonRecos;

  // if the input file is only AOD, then don't use EDProducts that are
  // included only in the full RECO tier;
  bool usingAODOnly;

  // which rec level is to be taken as the "best" one (e.g. used by
  // analyses);
  int lBest;

  // basic quantities for the chosen lepton (muon or electron);
  unsigned int leptonFlavor; // PDG ID
  double       leptonMass;   // in GeV/c^2

  ////////////////////////////////////////////////////////////////////
  // Event data
  ////////////////////////////////////////////////////////////////////
  
  // Keep track of the event number for printing out.
  int eventNum;

  // Keep track of how many events total have been processed.
  int eventsDone;

  // A handle to the TFileService for convenience.
  edm::Service<TFileService> fs;

  // Helper object to extract the important generator-level objects:
  // the resonance and its decay products (including brem. photons).
  HardInteraction hardInteraction;

  // Helper object for extracting the trigger decision from the paths
  // we care about.
  TriggerDecision trigDecision;

  // Object to help with per-rec-level information.
  RecLevelHelper recLevelHelper;

  // Helper object for lepton/dilepton cuts.
  CutHelper cutHelper;

  // The lepton collections.
  reco::CandidateBaseRefVector allLeptons[MAX_LEVELS];

  // The dilepton collections.
  reco::CompositeCandidateCollection allDileptons[MAX_LEVELS];

  ////////////////////////////////////////////////////////////////////
  // General utility
  ////////////////////////////////////////////////////////////////////
  
  // Return whether the rec level might be empty, based on the
  // useGen/Sim/Reco/etc. flags.
  bool skipRecLevel(const int level) const;

  // Get the invariant mass of the two leptons plus their closest
  // photons (what we've in the past called the "resonance").
  double resonanceMass(const reco::CompositeCandidate& dil) const;

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

#endif // ZP2MUANALYSIS_H
