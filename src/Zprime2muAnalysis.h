#ifndef ZP2MUANALYSIS_H
#define ZP2MUANALYSIS_H

#include <iosfwd>
#include <vector>
#include <string>

#include "TLorentzVector.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneratorUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TeVMuHelper.h"
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

  virtual void beginJob(const edm::EventSetup&) {}
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}

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

  // Helper object for extracting the trigger decision from the paths
  // we care about.
  TriggerDecision trigDecision;
  
  // Helper object for TeV dimuon analysis selection.
  TeVMuHelper tevMuHelper;

  // The main dilepton collections: gen, hlt, default offline, and
  // "best" offline.
  edm::InputTag genDils, hltDils, recDils, bestDils;
  edm::Handle<reco::CompositeCandidateCollection>
    genDileptons, hltDileptons, recDileptons, bestDileptons;
};

#endif // ZP2MUANALYSIS_H
