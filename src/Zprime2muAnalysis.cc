//
// Authors: Jason Mumford, Jordan Tucker, Slava Valuev, UCLA
//

#include <algorithm>

#include "TH1.h"
#include "TStyle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
//#include "PhysicsTools/UtilAlgos/interface/AnySelector.h"
//#include "PhysicsTools/UtilAlgos/interface/AnyPairSelector.h"
//#include "PhysicsTools/CandUtils/interface/CandCombiner.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"
//#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/tdrstyle.h"

using namespace std;

//-----------------------------------------------------------------------------
//                         Hard-coded parameters
//-----------------------------------------------------------------------------

// minimum pT value for all muons
const double Zprime2muAnalysis::PTMIN = 1e-3;
// maximum number of dileptons to ever keep
const unsigned int Zprime2muAnalysis::MAX_DILEPTONS = 2;
// max value of eta to consider muons
const double Zprime2muAnalysis::ETA_CUT = 2.4;
const double Zprime2muAnalysis::ENDCAP_BARREL_CUT = 1.04;
// max value of eta to trigger at the lower rec levels
const double Zprime2muAnalysis::TRIGGER_ETA_CUT[NUM_REC_LEVELS] = {
  999., 2.1, 2.1, 2.1};
// Activate quality cuts on tracks
const bool Zprime2muAnalysis::DO_QCUTS  = false;
//const int    Zprime2muAnalysis::QSEL = 0; // Bob's cuts
const int Zprime2muAnalysis::QSEL = 3; // Optimized cuts

// Quality cuts on L2 quantities.
const double Zprime2muAnalysis::L2QCUT[NUM_Q_SETS][NUM_L2_CUTS] = 
// Nrhit |FNrm   |L3 L2 Sigma Diff |FMS NrmChi2|
 {{-999., -999., -999., -999.},
//{-999., -999., -999., -999.},
//{   4.,  100., -999., -999.},  // Bob's cuts
//{   1., -999., -999., -999.},  // Jason's cuts 1
//{   4., -999.,    3., -999.},  // Jason's cuts 2
//{   2., -999.,   5.6, -999.},  // Jason's cuts 3
  {-999., -999., -999.,  10.},
  {-999., -999., -999.,  20.},
  {-999., -999., -999.,  30.},
  {-999., -999., -999.,  40.},
  {-999., -999., -999.,  50.},
  {-999., -999., -999.,  60.},
  {-999., -999., -999.,  70.}};

// Quality cuts on L3 quantities.
const double Zprime2muAnalysis::L3QCUT[NUM_Q_SETS][NUM_L3_CUTS] =
// Npix |Nsil  |FNrm  |BNrm  |F-BNrm|F-M-T-TDff |B-M-T-TDff |p*Ep  |EpT*EpT
 {{-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.},
//{-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.},
//{-999.,    0., -999., -999., -999.,    -999.,      -999.,      .4,  -999.},
//{-999.,    5.,    8., -999.,   30.,    -999.,      -999.,   -999.,  -999.},
//{   0.,    3.,   14.,   13.,   14.,      15.,        15.,     14.,    .35},
//{-999., -999.,   29., -999.,    2.,    -999.,      -999.,   -999.,  -999.},
  {-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.},
  {-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.},
  {-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.},
  {-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.},
  {-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.},
  {-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.},
  {-999., -999., -999., -999., -999.,    -999.,      -999.,   -999.,  -999.}};

// Quality cuts on Tracker only fit quantities.
const double Zprime2muAnalysis::TRACKERQCUT[NUM_Q_SETS][NUM_TRACKER_CUTS] =
 // FNrm | F-B  |p*Ep  |pT*EpT
 {{-999., -999., -999., -999.},
//{-999., -999., -999., -999.},
//{-999., -999., -999., -999.},
//{-999., -999., -999., -999.},
//{   6.,   30.,    .3,   .3 },
//{  32.,   55., -999.,   .25},
  {-999., -999., -999., -999.},
  {-999., -999., -999., -999.},
  {-999., -999., -999., -999.},
  {-999., -999., -999., -999.},
  {-999., -999., -999., -999.},
  {-999., -999., -999., -999.},
  {-999., -999., -999., -999.}};

// Activate Trigger cuts at levels:                gen    L1     L2     L3
bool Zprime2muAnalysis::cutTrig[NUM_REC_LEVELS] = {true,  true,  true,  true};

////////////////////////////////////////////////////////////////////
// Configuration and driver methods
////////////////////////////////////////////////////////////////////

Zprime2muAnalysis::Zprime2muAnalysis(const edm::ParameterSet& config) 
  : eventNum(-1) {
  recLevelHelper.init(config, true);

  verbosity = VERBOSITY(config.getUntrackedParameter<int>("verbosity", 0));
  dateHistograms = config.getUntrackedParameter<bool>("dateHistograms");
  doingElectrons = config.getParameter<bool>("doingElectrons");
  doingHiggs = config.getParameter<bool>("doingHiggs");
  constructGenDil = config.getParameter<bool>("constructGenDil");
  generatedOnly = config.getParameter<bool>("generatedOnly");
  doingGeant4 = config.getParameter<bool>("doingGeant4");
  reconstructedOnly = config.getParameter<bool>("reconstructedOnly");
  useOtherMuonRecos = config.getParameter<bool>("useOtherMuonRecos");
  usingAODOnly = config.getParameter<bool>("usingAODOnly");
  useTriggerInfo = config.getParameter<bool>("useTriggerInfo");

  // no generator-level information or other TeV muon reconstructors
  // available in the AOD (yet)
  if (usingAODOnly) {
    reconstructedOnly = true;
    useOtherMuonRecos = false;
    generatedOnly = false;
  }

  if (doingElectrons)
    useOtherMuonRecos = false;

  if (generatedOnly || reconstructedOnly)
    doingGeant4 = false;

  // input tags for the reco collections we need
  l1ParticleMap   = config.getParameter<edm::InputTag>("l1ParticleMap");
  hltResults      = config.getParameter<edm::InputTag>("hltResults");

  if (doingElectrons) {
    leptonFlavor = 11;
    leptonMass = 0.000511;
  }
  else {
    leptonFlavor = 13;
    leptonMass = 0.10566;

    // Level-1 paths we want to use for the trigger decision.
    l1paths.push_back(l1extra::L1ParticleMap::kSingleMu7);
    l1paths.push_back(l1extra::L1ParticleMap::kDoubleMu3);

    // Level-2 paths (actually, the names of the modules ran)
    hltModules[0].push_back("SingleMuNoIsoL2PreFiltered");
    hltModules[0].push_back("DiMuonNoIsoL2PreFiltered");
    // Level-3 paths (module names)
    hltModules[1].push_back("SingleMuNoIsoL3PreFiltered");
    hltModules[1].push_back("DiMuonNoIsoL3PreFiltered");
    
    // HLT paths (the logical ANDs of L2 and L3 single/dimuon paths
    // above)
    hltPaths.push_back("HLT1MuonNonIso");
    hltPaths.push_back("HLT2MuonNonIso");
  }

  // Our preferred style.
  InitROOT();
  // Physics TDR style.
  // TH1::AddDirectory(false);
  // setTDRStyle();
}

void Zprime2muAnalysis::analyze(const edm::Event& event,
				const edm::EventSetup& eSetup) {
  clearValues();

  // could store the whole event/eventsetup object if needed
  eventNum = event.id().event();
  
  if (!doingElectrons) {
    storeL1Decision(event);
    storeHLTDecision(event);
  }

  // Store a reference to the generator-level particle collection.
  edm::Handle<reco::CandidateCollection> genp;
  event.getByLabel("genParticleCandidates", genp);
  genParticles = &*genp;
  
  // Store the particles from the resonant interaction especially.
  if (!doingHiggs)
    storeInteractionParticles(*genParticles, eventNum, intParticles);

  // Store per-rec-level information from the event (includes getting
  // all the match maps from the event).
  recLevelHelper.initEvent(event);

  // At each rec level (including the "best" leptons), store the
  // lepton CandidateBaseRefs and make dileptons, adding in photons
  // to make the resonance four-vector.
  for (int irec = 0; irec <= MAX_LEVELS; irec++) {
    storeLeptons(event, irec);
    makeDileptons(irec);
    addBremCandidates(irec);
  }

  // Retrieve the vector of original rec levels for the best leptons.
  edm::Handle<vector<int> > bestRecs;
  event.getByLabel("bestMuons", bestRecs);
  bestRecLevels = *bestRecs;
  
  // Also store the true generator-level resonance.
  addTrueResonance(event);

  // Dump the event if appropriate.
  if (verbosity >= VERBOSITY_SIMPLE || eventIsInteresting()) {
    dumpEvent(true,true,true,true,true,true);
    dumpDileptonMasses();
  }

  // Compare official and homemade trigger decisions.
  if (!doingElectrons)
    compareTrigDecision(event);
}

////////////////////////////////////////////////////////////////////
// Initialization
////////////////////////////////////////////////////////////////////

void Zprime2muAnalysis::InitROOT() {
  gROOT->SetStyle("Plain");
  gStyle->SetFillColor(0);
  if (dateHistograms)
    gStyle->SetOptDate();
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetMarkerSize(.1);
  gStyle->SetMarkerStyle(8);
  gStyle->SetGridStyle(3);
  gStyle->SetPaperSize(TStyle::kA4);
  gStyle->SetStatW(0.25); // width of statistics box; default is 0.19
  gStyle->SetStatFormat("6.4g"); // leave default format for now
  gStyle->SetTitleFont(52,"XY"); // italic font for axis
  gStyle->SetLabelFont(52,"XY"); // italic font for axis labels
  gStyle->SetStatFont(52);       // italic font for stat. box
}

void Zprime2muAnalysis::clearValues() {
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    allLeptons[i_rec].clear();
    allDileptons[i_rec].clear();
    dileptonResonances[i_rec].clear();
  }

  bestLeptons.clear();
  bestDileptons.clear();
  dileptonResonances[MAX_LEVELS].clear();

  for (int i_rec = 0; i_rec < NUM_REC_LEVELS; i_rec++) {
    passTrig[i_rec] = true;
    trigWord[i_rec] = 0;
  }
}

////////////////////////////////////////////////////////////////////
// Storing leptons/dileptons
////////////////////////////////////////////////////////////////////

bool Zprime2muAnalysis::storeLeptons(const edm::Event& event,
				     const int rec) {
  edm::View<reco::Candidate> lepCandView;
  if (!recLevelHelper.getView(event, rec, lepCandView))
    return false;

  const unsigned nlep = lepCandView.size();
  unsigned ilep;

  // Sort the refs by descending momentum magnitude. (RefVector does
  // not have a sort() that would delegate function calls to its
  // product, nor does it have an operator[] that works as an lvalue,
  // so sort using an auxillary array of indices.)
  vector<unsigned> sortIndices;

  // Don't bother sorting generator-level leptons.
  if (rec != lgen) {
    // Initialize the current order of the indices.
    for (ilep = 0; ilep < nlep; ilep++)
      sortIndices.push_back(ilep);

    // Insertion-sort the indices according to the momentum magnitude of
    // the corresponding lepton.
    for (int i = 1; i < int(nlep); i++) {
      int ndx = i;
      double p = lepCandView[i].p();
      int j = i - 1;
      while (j >= 0 && p > lepCandView[sortIndices[j]].p()) {
	sortIndices[j+1] = sortIndices[j];
	j--;
      }
      sortIndices[j+1] = ndx;
    }
  }

  // Store CandidateBaseRefs in the order we just found.
  LeptonRefVector& leps = rec == MAX_LEVELS ? bestLeptons : allLeptons[rec];
  for (ilep = 0; ilep < nlep; ilep++) {
    unsigned jlep = rec == lgen ? ilep : sortIndices[ilep];
    // don't cut on pt at gen level
    if (lepCandView[jlep].pt() > PTMIN)
      leps.push_back(lepCandView.refAt(jlep));
  }

  return true;
}

reco::CompositeCandidate*
Zprime2muAnalysis::newDilepton(const reco::CandidateBaseRef& dau1,
			       const reco::CandidateBaseRef& dau2) {
  reco::CompositeCandidate* dil = new reco::CompositeCandidate;
  dil->addDaughter(reco::ShallowCloneCandidate(reco::CandidateBaseRef(dau1)));
  dil->addDaughter(reco::ShallowCloneCandidate(reco::CandidateBaseRef(dau2)));
  AddFourMomenta addP4;
  addP4.set(*dil);
  //dil->setP4(dau1->polarP4() + dau2->polarP4());
  //dil->setCharge(dau1->charge() + dau2->charge());
  return dil;
}

void
Zprime2muAnalysis::removeDileptonOverlap(reco::CandidateCollection& dileptons) {
  reco::CandidateCollection::iterator pdi, qdi;
  for (pdi = dileptons.begin(); pdi != dileptons.end() - 1;) {
    for (qdi = pdi+1; qdi != dileptons.end(); qdi++) {
      const reco::CandidateBaseRef pm = dileptonDaughterByCharge(*pdi, -1);
      const reco::CandidateBaseRef pp = dileptonDaughterByCharge(*pdi, +1);
      const reco::CandidateBaseRef qm = dileptonDaughterByCharge(*qdi, -1);
      const reco::CandidateBaseRef qp = dileptonDaughterByCharge(*qdi, +1);
      if (id(pm) == id(qm) || id(pp) == id(qp)) {
	// If either lepton is shared, remove the second dilepton
	// (i.e. the one with lower invariant mass since we have
	// sorted the vector already), reset pointers and restart.
	dileptons.erase(qdi);
	pdi = dileptons.begin();
	break;
      }
      else {
	pdi++;
      }
    }
  }
}

void Zprime2muAnalysis::makeDileptons(const int rec) {
  const static bool debug = verbosity >= VERBOSITY_SIMPLE;
  
  const LeptonRefVector& leptons = getLeptons(rec);

  if (leptons.size() < 2) {
    if (debug)
      LogTrace("makeDileptons")
	<< "Only " << leptons.size() << " at level " << rec
	<< ", refusing to try making dileptons!";
    return;
  }

  // Don't use getDileptons() since we want non-const access.
  reco::CandidateCollection& dileptons
    = rec == MAX_LEVELS ? bestDileptons : allDileptons[rec];

  // Dileptons search in reconstructed muons.  First consider all
  // possible combinations of mu+ and mu-.
  const int totalCharge = 0;
  LeptonRefVector::const_iterator plep, qlep;
  for (plep = leptons.begin(); plep != leptons.end() - 1; plep++)
    for (qlep = plep+1; qlep != leptons.end(); qlep++)
      if ((*qlep)->charge() + (*plep)->charge() == totalCharge) {
	if (rec == lgen) {
	  // For generator-level dileptons, unless we are told to do
	  // otherwise by the "constructGenDil" parameter, make
	  // sure we pick up the actual leptons from the resonance --
	  // leptons with the same mother that has a PDG id of one of
	  // the resonances we care about.
	  if (!constructGenDil) {
	    const reco::Candidate* mom = sameMother(*plep, *qlep);
	    if (!mom || !isResonance(mom->pdgId()))
	      continue;
	  }
	}
	// Apply cuts at levels of reconstruction above
	// generator-level.
	else if (leptonIsCut(*plep) || leptonIsCut(*qlep))
	  continue;

	dileptons.push_back(newDilepton(*plep, *qlep));
      }

  if (debug)
    LogTrace("makeDileptons") << "Reconstructed " << dileptons.size()
			      << " dileptons at rec level " << rec;

  // If H -> ZZ* -> 4mu, sort dileptons to look for two highest-mass ones.
  if (dileptons.size() > 1)
    dileptons.sort(reverse_mass_sort());

  if (dileptons.size() > 1) {
    removeDileptonOverlap(dileptons);

    // Store only 1 dilepton for Z' and G*, and two highest mass dileptons
    // for H -> ZZ* -> 4 mu.
    const unsigned maxKept = doingHiggs ? MAX_DILEPTONS : 1;
    if (dileptons.size() > maxKept)
      dileptons.erase(dileptons.begin() + maxKept, dileptons.end());
  }
}

void Zprime2muAnalysis::addBremCandidates(const int rec) {
  const static bool debug = verbosity >= VERBOSITY_LOTS;

  // Do nothing for generator-level dimuons; already done in
  // addTrueResonance().
  if (rec == lgen)
    return;

  const double dRmax = 0.1;

  ostringstream out;
  if (debug) out << "addBremCandidates, rec level " << rec << ":\n";

  const reco::CandidateCollection& dileptons = getDileptons(rec);

  for (unsigned idi = 0; idi != dileptons.size(); idi++) {
    const reco::Candidate& dil = dileptons[idi];
    LorentzVector p4 = dil.p4();
    if (debug)
      out << " Dilepton #" << idi << " p4 = " << p4 << endl;

    // No way to improve measurement at Level-1 and Level-2
    if (rec != l1 && rec != l2) {
      vector<LorentzVector> photonP4s;

      for (unsigned idau = 0; idau < dil.numberOfDaughters(); idau++) {
	const reco::CandidateBaseRef& dau = dileptonDaughter(dil, idau);
	LorentzVector photonP4 = closestPhoton(dau);
	if (debug) out << "  Considering photon " << photonP4
		       << " from daughter " << dau->p4();
	
	if (photonP4.energy() > 0) {
	  double dR = deltaR(dau->p4(), photonP4);
	  if (debug) out << "; its dR = " << dR;

	  if (dR < dRmax) {
	    bool addedAlready = false;
	    // Make sure we don't add the same photon twice.
	    for (unsigned iph = 0; iph < photonP4s.size(); iph++) {
	      const LorentzVector diff = photonP4s[iph] - photonP4;
	      if (diff.P() < 0.001) {
		addedAlready = true;
		break;
	      }
	    }
	    
	    if (!addedAlready) {
	      if (debug) out << "; adding it";
	      photonP4s.push_back(photonP4);
	      p4 += photonP4;
	    }
	    else {
	      if (debug) out << "; already added, skipping it";
	    }
	  }
	}
	else {
	  if (debug) out << "; no photon";
	}

	if (debug) out << endl;
      }	      
    }

    dileptonResonances[rec].push_back(p4);

    if (debug) 
      out << " Include brem candidate(s) at level " << rec
	  << ": inv. mass w/o photons = " << dil.mass()
	  << "; w/  = " << p4.mass();
  }

  if (debug) LogTrace("addBremCandidates") << out.str();
}

void Zprime2muAnalysis::addTrueResonance(const edm::Event& event) {
  const static bool debug = verbosity >= VERBOSITY_LOTS;
  
  unsigned ndil = allDileptons[lgen].size();
  if (ndil == 0) {
    edm::LogWarning("addTrueResonance")
      << "No generated resonance found in the event!";
    return;
  }

  ostringstream out;
  out << "addTrueResonance:\n";

  for (unsigned idi = 0; idi < ndil; idi++) {
    const reco::Candidate& dil = allDileptons[lgen][idi];
    const reco::Candidate* mom = sameMother(dileptonDaughter(dil, 0),
					    dileptonDaughter(dil, 1));

    LorentzVector p4;

    // Take true resonance from PYTHIA only if generated mu+ and mu-
    // originated from the same resonance, otherwise just use the
    // dilepton p4.
    // SV: same approach as in ORCA, but perhaps needs some thinking.
    if (mom && isResonance(mom->pdgId())) {
      if (debug) out << "  Found mother resonance; taking its four-vector...\n";
      p4 = mom->p4();
    }
    else if (constructGenDil) {
      p4 = dil.p4();
      // If we are allowed to construct the generated dilepton, go and
      // add in the bremmed final-state photons. We obviously cannot
      // use here the reconstructed photons as in addBremCandidates().
      // We could also just add the four-vectors of the documentation
      // leptons. Stick with the former for now.
      if (debug) out << "  Allowed to construct generator-level dilepton,"
		     << " finding final-state brem photons...\n";
      for (unsigned dau = 0; dau < dil.numberOfDaughters(); dau++) {
	 // Look at the documentation lines for the leptons, and find
	 // their daughter photons if any.
	const reco::Candidate* docLepton
	  = dileptonDaughter(dil, dau)->mother();
	out << "   Doc lepton: pdgId: " << setw(3) << docLepton->pdgId() 
	    << " status: " << docLepton->status()
	    << " ndau: " << docLepton->numberOfDaughters()
	    << " p4: " << docLepton->p4() << endl;
	// If docLepton isn't really a documentation line, or it only
	// has one daughter (the final-state lepton we just came
	// from), never mind.
	unsigned nsis = docLepton->numberOfDaughters();
	if (docLepton->status() == 3 && nsis > 1) {
	  // docMom's daughters are potentially our lepton's sister photons.
	  for (unsigned isis = 0; isis < nsis; isis++) {
	    const reco::Candidate* sis = docLepton->daughter(isis);
	    out << "   Sister: pdgId: " << sis->pdgId()
		<< " status: " << sis->status()
		<< " p4: " << sis->p4();
	    if (sis->status() == 1 && sis->pdgId() == 22) {
	      out << ". Found sister photon!";
	      p4 += sis->p4();
	    }
	    out << endl;
	  }
	}
      }
    }
    else
      p4 = dil.p4();
    
    dileptonResonances[lgen].push_back(p4);

    if (debug) {
      out << " MC resonance: dimuon inv. mass = " << dil.mass()
	  << "; true inv. mass = " << p4.mass();
      LogTrace("addTrueResonance") << out.str();
    }
  }
}

////////////////////////////////////////////////////////////////////
// Using trigger info
////////////////////////////////////////////////////////////////////

void Zprime2muAnalysis::storeL1Decision(const edm::Event& event) {
  const static bool debug = (verbosity >= VERBOSITY_SIMPLE);

  // Get Level-1 decisions for trigger paths we are interested in.
  edm::Handle<l1extra::L1ParticleMapCollection> l1MapColl;
  event.getByLabel(l1ParticleMap, l1MapColl);

  if (!l1MapColl.isValid()) {
    edm::LogWarning("storeL1Decision")
      << "L1ParticleMapCollection with label [" << l1ParticleMap.encode()
      << "] not found!" << endl;
    return;
  }

  if (debug) LogTrace("storeL1Decision") << "storeL1Decision:";

  // Loop over chosen paths, check trigger decisions, and save them into
  // "trigbits".
  unsigned int trigbits = 0;
  int nl1paths = l1paths.size();
  for (int ipath = 0; ipath < nl1paths; ipath++) {
    const l1extra::L1ParticleMap& thisMap = (*l1MapColl)[l1paths[ipath]];
    bool fired = thisMap.triggerDecision();
    if (fired) trigbits = trigbits | (1 << ipath);
    if (debug) LogTrace("storeL1Decision")
      << "  " << thisMap.triggerName() << " (index " << l1paths[ipath]
      << "): decision " << fired;
  }

  if (debug)
    LogTrace("storeL1Decision") << " L1 official trigbits: " << trigbits;

  trigWord[l1] = trigbits;
  passTrig[l1] = (trigbits != 0);
}

void Zprime2muAnalysis::storeHLTDecision(const edm::Event& event) {
  const static bool debug = (verbosity >= VERBOSITY_SIMPLE);
  ostringstream out;

  // Getting the result of the entire HLT path is done easily with the
  // TriggerResults object, however to get at the separate decisions
  // for L2 and L3 we have to do a little hacky magic.

  // Try to get the HLT TriggerResults object now, before
  // trying to get the HLTFilterObjectWithRefs below
  // so that if there is no HLT information in the file, 
  // getByLabel will go ahead and throw an exception.
  edm::Handle<edm::TriggerResults> hltRes;
  event.getByLabel(hltResults, hltRes);
  edm::TriggerNames hltTrigNames;
  hltTrigNames.init(*hltRes);

  if (debug) out << "storeHLTDecision:\n";

  // Extract L2 and L3 decisions by seeing if the corresponding
  // HLTFilterObjectWithRefs exists and seeing how many muons it holds.
  for (unsigned int lvl = 0; lvl < 2; lvl++) {
    unsigned int l = l2 + lvl;
    unsigned int trigbits = 0;
    for (unsigned int ipath = 0; ipath < hltModules[lvl].size(); ipath++) {
      const string& trigName = hltModules[lvl][ipath];
      edm::Handle<reco::HLTFilterObjectWithRefs> hltFilterObjs;

      bool fired = true;
      try {
	event.getByLabel(trigName, hltFilterObjs);
      }
      catch (const cms::Exception& e) {
	fired = false;
      }

      // There may be an HLTFilterObject in the event even if the
      // trigger did not accept; the real check is to make sure that
      // the minimum number of muons for the trigger was met.
      const unsigned int minNMuons = ipath + 1;
      if (fired && hltFilterObjs->size() < minNMuons)
	fired = false;	  
	
      if (debug)
	out << "  " << trigName
	    << ": decision = " << fired << endl;

      if (fired) {
	trigbits = trigbits | (1 << ipath);

	if (debug) {
	  out << "  " << trigName << " filter result muons:\n";
	  reco::HLTFilterObjectWithRefs::const_iterator muItr;
	  int imu = 0;
	  for (muItr = hltFilterObjs->begin(); muItr != hltFilterObjs->end();
	       muItr++) {
	    out << "    #" << imu++ << " q = " << muItr->charge()
		<< " p = (" << muItr->px() << ", " << muItr->py()
		<< ", " << muItr->pz() << ", " << muItr->energy() << ")\n"
		<< "     pt = " << muItr->pt() << " eta = " << muItr->eta()
		<< " phi = " << muItr->phi() << endl;
	  }
	}
      }
    }

    trigWord[l] = trigbits;
    passTrig[l] = (trigbits != 0);
    out << "  trigWord[l" << l << "]: " << trigWord[l] << endl;
  }
    
  // Check that the official full HLT path decisions agree with
  // what we extracted for the official L3 decision.
  unsigned int hlt_trigbits = 0;
  for (unsigned int i = 0; i < hltPaths.size(); i++) {
    int ndx = hltTrigNames.triggerIndex(hltPaths[i]);
    bool fired = hltRes->accept(ndx);
    if (debug)
      out << " HLT path #" << ndx << ": " << hltPaths[i]
	  << " decision = " << fired << endl;
    if (fired) hlt_trigbits |= 1 << i;
  }
  if (hlt_trigbits != trigWord[l3]) {
    edm::LogWarning("storeHLTDecision")
      << "+++ Warning: official HLT"
      << " decision disagrees with extracted L3 decision:"
      << " official HLT: " << hlt_trigbits
      << " extracted L3: " << trigWord[l3] << " +++";
      
    if (verbosity == VERBOSITY_NONE)
      dumpEvent(false, true, true, true, false);
  }

  if (debug) LogTrace("storeHLTDecision") << out.str();
}

bool Zprime2muAnalysis::TriggerTranslator(const string& algo,
					  const unsigned int lvl,
					  const unsigned int nmu) const {
  // See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonHLT
  // https://twiki.cern.ch/twiki/bin/view/CMS/L1TriggerTableHLTExercise
  // L1 Quality codes are still the same, but we are suggested to use the 
  // methods useInSingleMuonTrigger() etc. instead of using the quality
  // code directly:
  // https://twiki.cern.ch/twiki/bin/view/CMS/GMTEmulator

  // Arrays are first-indexed by trigger level then by the number of
  // muons required.
  const unsigned int l = lvl - 1;
  const unsigned int n = nmu - 1;

  // these are for L = 10^{32}
  const double ptMin[3][2] = {
    { 7, 3},
    {16, 3},
    {16, 3}
  };
  const double etaMax[3] = { 2.5, 2.5, 2.5 };
  // 28/10/2007: in CMSSW_1_6_0, requirement is just that 
  // the vertex be < 200 microns from the beam axis, i.e. 
  // sqrt(vx**2+vy**2) < 0.02
  const double dxy2Max[3] = { 9999, 9999, 0.02*0.02 };
  const double dzMax[3] = { 9999, 9999, 9999 };
  const double nSigmaPt[3] = { 0, 3.9, 2.2 };
  // 12/05/2004: a muon must have more than 5 pixel plus silicon hits in
  // the tracker (DAQ TDR, p. 306).
  //const int minHits[3] = { 0, 0, 6 };
  // 28/10/2007: for CMSSW_1_6_0, min hits is 0
  // JMTBAD no hit counts in the AOD anyway...
  const int minHits[3] = { 0, 0, 0 };

  // Check if the algorithm string passed in is one we recognized.
  // String value is not currently used, but will be used in the
  // future to differentiate between, e.g., Iso and NonIso triggers.
  if (algo != "L1_SingleMu7" &&
      algo != "L1_DoubleMu3" &&
      algo != "SingleMuNoIsoL2PreFiltered" && 
      algo != "SingleMuNoIsoL3PreFiltered" &&
      algo != "DiMuonNoIsoL2PreFiltered" && 
      algo != "DiMuonNoIsoL3PreFiltered") {
    edm::LogWarning("TriggerTranslator")
      << "+++ unrecognized algorithm " << algo << "! +++";
    return false;
  }

  const static bool debug = verbosity >= VERBOSITY_SIMPLE;
  ostringstream out;

  if (debug) out << "TriggerTranslator, " << algo << ":\n";

  unsigned int muonsPass = 0;
  //vector<zp2mu::Muon>::const_iterator pmu; //, pmu_prev;

  //vector<double> zvtx;
  //vector<double>::const_iterator pvtx;
  for (unsigned imu = 0; imu < allLeptons[lvl].size(); imu++) {
    const reco::CandidateBaseRef& mu = allLeptons[lvl][imu];
    reco::TrackRef tk;
    if (lvl < l3)
      tk = mu->get<reco::TrackRef>();
    else
      tk = mu->get<reco::TrackRef, reco::CombinedMuonTag>();

    // JMTBAD HLTMuonPrefilter cuts not on pt but on what they call ptLx
    // nSigmaPt[l1] = 0, so ptLx(l1) = pt(l1), but safeguard against
    // accidentally setting nSigmaPt[l1] != 0:
    double ptLx;
    if (lvl == l1)
      ptLx = mu->pt();
    else
      ptLx = (1 + nSigmaPt[l]*invPError(tk)*mu->p())*mu->pt();
    if (debug)
      out << "  mu #" << id(mu)
	  << " pt: " << mu->pt() << " ptLx: " << ptLx << " (cut: " << ptMin[l][n] 
	  << ") eta: " << mu->eta() << " (cut: " << etaMax[l] << ")\n";
    bool pass = ptLx >= ptMin[l][n] && fabs(mu->eta()) <= etaMax[l];
    // don't bother evaluating the other constraints if they aren't set
    if (dxy2Max[l] != 9999) {
      double vx = mu->vx();
      double vy = mu->vy();
      pass = pass && vx*vx + vy*vy < dxy2Max[l];
    }
    if (dzMax[l] != 9999) {
      pass = pass && fabs(mu->vz()) < dzMax[l];
    }
    if (minHits[l] != 0) {
      // JMTBAD switch to reco::Track::numberOfValidHits()
      pass = pass && nHits(mu, HITS_TRK) > minHits[l];
    }
    if (pass) {
      // here check extra stuff depending on the trigger algorithm
      bool passExtra = true;

      if (lvl == l1) {
	int quality
	  = toConcrete<l1extra::L1MuonParticle>(mu).gmtMuonCand().quality();
	// In single muon trigger, GMT uses only muons with quality > 3.
	// In dimuon trigger, GMT uses muons with quality = 3 and 5-7.
	if ((nmu == 1 && quality < 4) ||
	    (nmu == 2 && (quality < 3 || quality == 4)))
	  passExtra = false;
      }
      /*
      // JMTBAD disabled for now, not currently in HLT
      else if (lvl == l3 && nmu == 2) {
        for (unsigned jmu = 0; jmu < allLeptons[lvl].size(); jmu++) {
	  const reco::CandidateBaseRef& mu_prev = allLeptons[lvl][jmu];
          // Skip ghost tracks (see p. 308 of DAQ TDR)
	  if (fabs(mu->eta() - mu_prev->eta()) < 0.01 &&
	      fabs(mu->phi() - mu_prev->phi()) < 0.05 &&
	      fabs(mu->pt()  - mu_prev->pt())  < 0.1) passExtra = false;
	}
      }
      */

      if (passExtra)
	muonsPass++;
    }

    if (muonsPass == nmu)
      break;
  }

  bool result;
  if (debug) out << "  TriggerTranslator result for " << algo << ": ";
  if (muonsPass >= nmu) {
    out << "pass!";
    result = true;
  }
  else {
    out << "fail!";
    result = false;
  }
  if (debug) LogTrace("TriggerTranslator") << out.str();

  return result;
}

void Zprime2muAnalysis::compareTrigDecision(const edm::Event& event,
					    bool old) const {
  const static bool debug = (verbosity >= VERBOSITY_SIMPLE);
  ostringstream out;

  if (debug) out << "compareTrigDecision:\n";

  edm::Handle<l1extra::L1ParticleMapCollection> l1MapColl;
  event.getByLabel(l1ParticleMap, l1MapColl);

  for (unsigned int lvl = l1; lvl <= l3; lvl++) {
    unsigned homemade_trigbits = 0;
    unsigned npaths = lvl == l1 ? l1paths.size() : hltModules[lvl-l2].size();
    for (unsigned ipath = 0; ipath < npaths; ipath++) {
      const string& trigName = lvl > l1 ? hltModules[lvl-l2][ipath]
	: (*l1MapColl)[l1paths[ipath]].triggerName();

      // Try to emulate HLT algorithms.
      bool fired = TriggerTranslator(trigName, lvl, ipath+1);
      
      // If the event passes, set the corresponding bit in trigbits.
      if (fired) homemade_trigbits = homemade_trigbits | (1 << ipath);
      if (debug)
	out << "  " << trigName << " (homemade): decision = " << fired << endl;
    }

    // "Official" muon HLTs are calculated only when corresponding
    // previous levels gave OK, while we calculate the decision for a
    // given level regardless of previous levels' decisions.
    if (lvl >= l2)
      homemade_trigbits &= trigWord[l1];
    if (lvl >= l3)
      homemade_trigbits &= trigWord[l2];
    // Compare official and homemade decisions.
    if (homemade_trigbits != trigWord[lvl]) {
      edm::LogWarning("compareTrigDecision")
	<< "+++ Warning: official L" << lvl
	<< " decision disagrees with homemade decision:"
	<< " official: " << trigWord[lvl]
	<< " homemade: " << homemade_trigbits << " +++";
      
      if (verbosity == VERBOSITY_NONE)
	dumpEvent(false, true, true, true, false);
    }
  }

  if (debug) LogTrace("compareTrigDecision") << out.str();
}

////////////////////////////////////////////////////////////////////
// Lepton/dilepton rec level utility functions
////////////////////////////////////////////////////////////////////

int Zprime2muAnalysis::id(const reco::CandidateRef& cand) const {
  if (cand.isNull())
    return -999;
  return cand.index();
}

int Zprime2muAnalysis::id(const reco::CandidateBaseRef& cand) const {
  if (cand.isNull())
    return -999;
  return cand.key();
}

int Zprime2muAnalysis::recLevel(const reco::CandidateBaseRef& cand) const {
  int level = recLevelHelper.recLevel(cand);
  if (level == lbest) {
    int ndx = id(cand);
    if (ndx >= 0 && ndx < int(bestRecLevels.size()))
      level = bestRecLevels[ndx];
  }
  return level;
}

int Zprime2muAnalysis::recLevel(const reco::Candidate& dil) const {
  int rec = -1;
  for (unsigned ilep = 0; ilep < dil.numberOfDaughters(); ilep++) {
    int r = recLevel(dileptonDaughter(dil, ilep));
    if (rec >= 0) {
      if (r != rec)
	return lbest;
    }
    else
      rec = r;
  }
  return rec;
}  

const Zprime2muAnalysis::LeptonRefVector&
Zprime2muAnalysis::getLeptons(const int rec) const {
  recLevelHelper.checkRecLevel(rec, "getLeptons");
  return rec < MAX_LEVELS ? allLeptons[rec] : bestLeptons;
}

const reco::CandidateCollection&
Zprime2muAnalysis::getDileptons(const int rec) const {
  recLevelHelper.checkRecLevel(rec, "getDileptons");
  return rec < MAX_LEVELS ? allDileptons[rec] : bestDileptons;
}

////////////////////////////////////////////////////////////////////
// Lepton/dilepton matching
////////////////////////////////////////////////////////////////////

bool Zprime2muAnalysis::matchDilepton(const reco::Candidate& dil,
				      const int level,
				      const reco::Candidate* newdil) const {
  const reco::CandidateBaseRef& lepp
    = matchedLepton(dileptonDaughterByCharge(dil, +1), level);
  const reco::CandidateBaseRef& lepm
    = matchedLepton(dileptonDaughterByCharge(dil, -1), level);
  int idp = id(lepp), idm = id(lepm);
  
  // look for a dilepton at that level which has the same two leptons
  const reco::CandidateCollection& dileptons = getDileptons(level);
  for (unsigned idi = 0; idi < dileptons.size(); idi++)
    if (id(dileptonDaughterByCharge(dileptons[idi], +1)) == idp &&
	id(dileptonDaughterByCharge(dileptons[idi], -1)) == idm) {
      newdil = &dileptons[idi];
      return true;
    }
  return false;
}

////////////////////////////////////////////////////////////////////
// Generator-level utility functions
////////////////////////////////////////////////////////////////////

void Zprime2muAnalysis::SetP4M(reco::Particle::LorentzVector& v,
	    double pt, double phi, double p, double theta, double m) const {
  v.SetCoordinates(pt*cos(phi), pt*sin(phi), p*cos(theta), sqrt(p*p + m*m));
}

const reco::Candidate*
Zprime2muAnalysis::mother(const reco::CandidateBaseRef& cand) const {
  int pId = cand->pdgId();
  const reco::Candidate* mom = cand->mother();
  while (mom != 0 && mom->pdgId() == pId)
    mom = mom->mother();
  return mom;
}
  
int Zprime2muAnalysis::motherId(const reco::CandidateBaseRef& cand) const {
  const reco::Candidate* mom = mother(cand);
  if (mom != 0) return mom->pdgId();
  else return 0;
  /*    throw cms::Exception("motherId")
        << "+++ cannot find mother for particle! pdgId: " << pId
        << " status: " << cand->status() << " +++\n";*/
}

int Zprime2muAnalysis::grandmotherId(const reco::CandidateBaseRef& cand) const {
  const reco::Candidate *mom = mother(cand);
  if (mom == 0) return 0;
  const reco::Candidate *gm = mom->mother();
  if (gm != 0) return gm->pdgId();
  else return 0;
  /* throw cms::Exception("grandmotherId")
     << "+++ cannot find grandmother for particle! pdgId: " << p->pdgId()
     << " status: " << cand->status() 
     << " mother pdgId: " << mId << " status: " << mom->status() << " +++\n";*/
}
 
const reco::Candidate*
Zprime2muAnalysis::sameMother(const reco::CandidateBaseRef& cand1,
			      const reco::CandidateBaseRef& cand2) const {
  const reco::Candidate* m1 = mother(cand1);
  const reco::Candidate* m2 = mother(cand2);

  bool sameMother = m1 != 0 && m2 != 0 && m1 == m2;
  //m1->pdgId() == m2->pdgId() &&
  //m1->status() == m2->status() && m1->p4() == m2->p4();

  if (sameMother) return m1;
  else return 0;
}

bool Zprime2muAnalysis::isResonance(int pdgId) const {
  // Z', Z0/gamma*, G, or G*
  return pdgId == 32 || pdgId == 23 || pdgId == 39 || pdgId == 5000039;
}

// Store the generator-level momenta for the particles in the resonant
// interaction from the MC record (assuming the event was of the form
// q qbar -> resonance -> l+l-); return success as a bool. The returned
// momenta include the true resonance, the final-state (after-brem)
// l-l+, the muons before bremsstrahlung, and the quark that entered
// the hard interaction.
bool Zprime2muAnalysis::storeInteractionParticles(
	                        const reco::CandidateCollection& genParticles, 
				int eventNum,
				InteractionParticles& ip) const{
  static const bool debug = verbosity >= VERBOSITY_TOOMUCH;
  ostringstream out;

  // Zero out the pointers to trap errors later.
  ip.genQuark = 0; 
  ip.genResonance = 0;
  ip.genLepPlus = 0;
  ip.genLepMinus = 0;
  ip.genLepPlusNoIB = 0;
  ip.genLepMinusNoIB = 0;

  // Count the total number of quarks entering the hard interaction.
  int iq = 0;

  // Loop over all particles stored in the current event.
  if (debug) out << "Looking for interaction quark:\n";
  reco::CandidateCollection::const_iterator p = genParticles.begin();
  for ( ; p != genParticles.end(); p++) {
    if (debug)
      out << "id: " << p->pdgId() << " status: " << p->status() << endl;

    if (p->status() == 3) {
      // Documentation lines.  We look for partons that enter the hard
      // interaction.  Their "mothers" are the partons that initiate parton
      // showers; the "mothers" of the "mothers" are primary protons, always
      // on lines 1 and 2.
      const reco::Candidate* mother = p->mother();
      if (mother != 0) {
	const reco::Candidate* grandmother = mother->mother();
	if (grandmother != 0 && grandmother->pdgId() == 2212) {
	  int id = p->pdgId();
	  if (id >= 1 && id <= 6) {
	    if (debug) out << "  quark in hard interaction found!\n";
	    iq++;
	    ip.genQuark = &*p;
	  } 
	}
	else {
	  if (debug) out << "  grandmother invalid!\n";
        }
      }
    }
  } 

  if (debug) LogDebug("storeInteractionParticles") << out.str();
  out.str("");

  if (iq != 1) { 
    edm::LogWarning("storeInteractionParticles")
      << "+++ Found " << iq << " quarks in the hard interaction;"
      << " skip the event # " << eventNum << " +++\n";
    return false;
  }

  // Find the resonance and its decay products.
  p = genParticles.begin();
  for ( ; p != genParticles.end(); p++) {
    // Look for resonance in documentation lines since the one in the
    // main body of event listing has no end vertex and no
    // descendants.
    int pid = p->pdgId();
    if (p->status() == 3 && isResonance(pid)) {
      if (debug)
	out << "Resonance found, id: " << p->pdgId() << endl;

      if (ip.genResonance == 0)
	ip.genResonance = &*p;
      else
	throw cms::Exception("storeInteractionParticles")
	  << "+++ Second generated resonance found in event # " << eventNum
	  << "! +++\n";

      if (debug) out << ip.genResonance ;

      // Resonance children are leptons before bremsstrahlung and the
      // resonance itself (??).  Use them later to get leptons before
      // brem.
      const reco::GenParticleCandidate* ppp =
	dynamic_cast<const reco::GenParticleCandidate*>(&*p);
      size_t ResChildren_size = ppp->numberOfDaughters();

      if (debug) {
	out << "Children:\n";
	for (unsigned int ic = 0; ic < ResChildren_size; ic++)
	  out << " " <<  ppp->daughter(ic) << endl;
	out << "Leptons before brem:\n";
      }

      // Look for the pre-brem leptons, (hepmc says these have status
      // 3).
      for (unsigned int ic = 0; ic < ResChildren_size; ic++) {
	if (ppp->daughter(ic)->status() == 3) {
	  if (ppp->daughter(ic)->pdgId() == int(leptonFlavor)) {
	    if (ip.genLepMinusNoIB == 0)
	      ip.genLepMinusNoIB = ppp->daughter(ic);
	    else
	      throw cms::Exception("storeInteractionParticles")
		<< "+++ Second before-brem l- found in resonance decay "
		<< "in event # " << eventNum << "! +++\n";

	    if (debug) out << ip.genLepMinusNoIB << endl;
	  }
	  else if (ppp->daughter(ic)->pdgId() == -int(leptonFlavor)) {
	    if (ip.genLepPlusNoIB == 0)
	      ip.genLepPlusNoIB = ppp->daughter(ic);
	    else
	      throw cms::Exception("storeInteractionParticles")
		<< "+++ Second before-brem l+ found in resonance decay "
		<< "in event # " << eventNum << "! +++\n";

	    if (debug) out << ip.genLepPlusNoIB << endl;
	  }
	}
      }

      // Z' descendants also include final-state muons and photons produced
      // by initial-state muons.
      // CL: add up to granddaughter, enough?
      typedef vector<const reco::Candidate*> CandidateVector;
      CandidateVector ResDescendants;
      CandidateVector ResGrandChildren;
      for (unsigned int ida = 0; ida != p->numberOfDaughters(); ida++)
	ResDescendants.push_back(p->daughter(ida));

      CandidateVector::const_iterator id = ResDescendants.begin();
      for ( ; id != ResDescendants.end(); id++) {
        if ((*id)->numberOfDaughters() == 0)
	  continue;
        for (unsigned igdd = 0; igdd < (*id)->numberOfDaughters(); igdd++)
	  ResGrandChildren.push_back((*id)->daughter(igdd));
      }
      ResDescendants.insert(ResDescendants.end(),
			    ResGrandChildren.begin(), ResGrandChildren.end());

      if (debug) {
	out << "Descendants:\n";
	for (id = ResDescendants.begin(); id != ResDescendants.end(); id++)
	  out << *id << endl;

	out << "Stable leptons:\n";
      }

      for (id = ResDescendants.begin(); id != ResDescendants.end(); id++) {
	if ((*id)->status() == 1) {
	  if ((*id)->pdgId() == int(leptonFlavor)) {
	    if (ip.genLepMinus == 0)
	      ip.genLepMinus = *id;
	    else
	      throw cms::Exception("storeInteractionParticles")
		<< "+++ Second mu- found in resonance decay in event # "
		<< eventNum << "! +++\n";

	    if (debug) out << ip.genLepMinus << endl;
	  }
	  else if ((*id)->pdgId() == -int(leptonFlavor)) {
	    if (ip.genLepPlus == 0)
	      ip.genLepPlus = *id;
	    else
	      throw cms::Exception("storeInteractionParticles")
		<< "+++ Second mu+ found in resonance decay in event # "
		<< eventNum << "! +++\n";

	    if (debug) out << ip.genLepPlus;
	  }
	}
      }
      // Do not break for now...
    }
  }

  if (debug) LogDebug("storeInteractionParticles") << out.str();

  if (ip.genResonance == 0) {
    edm::LogWarning("storeInteractionParticles")
      << " +++ There is no resonance generated in event # "
      << eventNum << "! +++\n";
    return false;
  }
  if (ip.genLepPlus == 0 || ip.genLepMinus == 0 ||
      ip.genLepMinusNoIB == 0 || ip.genLepPlusNoIB == 0) {
    edm::LogWarning("storeInteractionParticles")
      << " +++ At least one lepton from resonance decay is not found in event"
      << " #" << eventNum << "! +++\n";
    return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////////
// Track utility functions
////////////////////////////////////////////////////////////////////

reco::TrackBaseRef
Zprime2muAnalysis::getMainTrack(const reco::CandidateBaseRef& cand) const {
  const int rec = recLevel(cand);
  if (doingElectrons)
    return reco::TrackBaseRef(cand->get<reco::GsfTrackRef>());
  else if (rec == l2)
    return reco::TrackBaseRef(cand->get<reco::TrackRef>());
  else if (rec >= l3)
    return reco::TrackBaseRef(cand->get<reco::TrackRef,
			      reco::CombinedMuonTag>());
  else
    return reco::TrackBaseRef(); // an invalid ref
}

template <typename TrackType>
double Zprime2muAnalysis::pError(const TrackType& track) const {
  // dumb identity:
  // p = q / qoverp
  // dp_dqoverp = - q/qoverp^2

  const double qoverp = track->qoverp();
  const double q = track->charge();
  
  double sig_p = track->qoverpError() * q/qoverp/qoverp;
  if (sig_p < 0) sig_p = -sig_p;
  return sig_p;
}

template <typename TrackType>
double Zprime2muAnalysis::ptError(const TrackType& track) const {
  return track->ptError();
}

template <typename TrackType>
double Zprime2muAnalysis::invPtError(const TrackType& tk) const {
  return invError(tk->pt(), tk->ptError());
}
  
template <typename TrackType>
double Zprime2muAnalysis::invPError(const TrackType& tk) const {
  return invError(tk->p(), pError(tk));
}

double Zprime2muAnalysis::ptError(const reco::CandidateBaseRef& cand) const {
  return getMainTrack(cand)->ptError();
}

double Zprime2muAnalysis::pError(const reco::CandidateBaseRef& cand) const {
  return pError(getMainTrack(cand));
}

double Zprime2muAnalysis::invPtError(const reco::CandidateBaseRef& cand) const {
  return invPtError(getMainTrack(cand));
}

double Zprime2muAnalysis::invPError(const reco::CandidateBaseRef& cand) const {
  return invPError(getMainTrack(cand));
}

template <typename TrackType>
void Zprime2muAnalysis::numSubDetHits(const TrackType& theTrack,
				      int& nPixHits, int& nSilHits,
				      int& nOtherHits,
				      DetId::Detector otherDet) const {
  nPixHits = nSilHits = nOtherHits = 0;

  if (theTrack.isNull()) {
    edm::LogError("numSubDetHits")
      << "theTrack is a null reference!";
    return;
  }
    
  int nRecHits = theTrack->recHitsSize();
  // loop over all the rechits and count how many are in the pixels and
  // how many are in the strips
  for (trackingRecHit_iterator trh = theTrack->recHitsBegin();
       trh != theTrack->recHitsEnd(); trh++) {
    DetId did = (*trh)->geographicalId();
    DetId::Detector det = did.det();
    if (det == DetId::Tracker) {
      int subdetId = did.subdetId();
      // trying not to use magic numbers -- but pixel subdets are 1 and 2,
      // strip subdets are 3-6
      switch (subdetId) {
      case StripSubdetector::TIB:
      case StripSubdetector::TOB:
      case StripSubdetector::TID:
      case StripSubdetector::TEC:
	nSilHits++;
	break;
      case PixelSubdetector::PixelBarrel:
      case PixelSubdetector::PixelEndcap:
	nPixHits++;
	break;
      default:
	edm::LogWarning("numSubDetHits")
	  << "unknown subid for TrackingRecHit";
      }
    }
    else if (det == otherDet) // DetId::Muon, DetId::Ecal
      nOtherHits++;
  }
  // redundant check, hopefully
  if (nPixHits + nSilHits + nOtherHits != nRecHits)
    edm::LogWarning("numSubDetHits")
      << "TrackingRecHits from outside the tracker or 'other' systems!";
}

int Zprime2muAnalysis::nHits(const reco::CandidateBaseRef& cand, 
			     const int type) const {
  const int rec = recLevel(cand);
  const reco::TrackBaseRef tk = getMainTrack(cand);
  if (type == HITS_ALL || (rec == l2 && type == HITS_OTH))
    return tk->recHitsSize();
  else if (rec >= l3) {
    int pix, sil, oth;
    if (doingElectrons)
      numSubDetHits(tk, pix, sil, oth, DetId::Ecal);
    else
      numSubDetHits(tk, pix, sil, oth, DetId::Muon);
    if      (type == HITS_OTH) return oth;
    else if (type == HITS_TRK) return pix + sil;
    else if (type == HITS_PIX) return pix;
    else if (type == HITS_SIL) return sil;
  }
  // return nonsensical number of hits if the above logic failed
  return -1;
}

bool Zprime2muAnalysis::matchTracks(const reco::CandidateBaseRef& cand1,
				    const reco::CandidateBaseRef& cand2) const {
  return fabs(deltaPhi(cand1->phi(), cand2->phi())) < .5 &&
    fabs(cand1->eta() - cand2->eta()) < .5;
}

double
Zprime2muAnalysis::ptLx(const reco::TrackRef& theTrack, const int rec) const {
  static const double nSigma[3] = { 0, 3.9, 2.2 };
  if (rec < l1 || rec > l3) {
    edm::LogError("ptLx") << "cannot compute ptLx for generated/offline muons!";
    return 0;
  }
  double ptLx = theTrack->pt();
  double err0 = theTrack->error(0);
  double abspar0 = fabs(theTrack->parameter(0));
  if (abspar0>0) ptLx += nSigma[rec-1]*err0/abspar0*ptLx;
  return ptLx;
}

////////////////////////////////////////////////////////////////////
// Dilepton utility functions
////////////////////////////////////////////////////////////////////

const Zprime2muAnalysis::LorentzVector&
Zprime2muAnalysis::resV(const int rec, const int idil) const {
  recLevelHelper.checkRecLevel(rec, "resV");
  return dileptonResonances[rec][idil];
}

int Zprime2muAnalysis::numDaughtersInAcc(const reco::Candidate& dil,
					 const double etaCut) const {
  int count = 0;
  for (unsigned ilep = 0; ilep < dil.numberOfDaughters(); ilep++)
    if (fabs(dil.daughter(ilep)->eta()) < etaCut)
      count++;
  return count;
}

const reco::CandidateBaseRef
Zprime2muAnalysis::dileptonDaughter(const reco::Candidate& dil,
				    const unsigned i) const {
  if (i < 0 || i >= dil.numberOfDaughters())
    return reco::CandidateBaseRef(); // an invalid reference
  return dil.daughter(i)->masterClone();
}

const reco::CandidateBaseRef
Zprime2muAnalysis::dileptonDaughterByCharge(const reco::Candidate& dil,
					    const int charge) const {
  for (unsigned ilep = 0; ilep < dil.numberOfDaughters(); ilep++) {
    if (dil.daughter(ilep)->charge() == charge)
      return dileptonDaughter(dil, ilep);
  }
  return reco::CandidateBaseRef(); // an invalid reference
}

////////////////////////////////////////////////////////////////////
// Print-outs
////////////////////////////////////////////////////////////////////

void Zprime2muAnalysis::dumpLepton(ostream& output,
				   reco::CandidateBaseRef cand) const {
  // Make sure we're looking at the master ref (to handle shallow
  // clones which are daughters of dileptons).
  if (cand->hasMasterClone())
    cand = cand->masterClone().castTo<reco::CandidateBaseRef>();

  const int level = recLevel(cand);

  output << recLevelHelper.levelName(level, true)
	 << " #"          << setw(1) << id(cand)
	 << "  Q: "       << setw(2) << cand->charge();

  if (level == lgen) {
    //const reco::GenParticleCandidate& genmu
    //  = toConcrete<reco::GenParticleCandidate>(cand);
    // JMTBAD fake motherline for tkdiff
    output << " Origin: " << setw(4) << motherId(cand)
	   << "/"         << setw(4) << 0 //(candEx.id() < 2 ? 6 : 0) // rhs.genMotherLine()
	   << " Phi: "    << setw(7) << setprecision(4) << cand->phi()
	   << " Eta: "    << setw(7) << setprecision(4) << cand->eta()
	   << " Pt: "     << setw(7) << setprecision(5) << cand->pt()
	   << " P: "      << setw(7) << setprecision(5) << cand->p()
	   << endl;
  }
  else if (level == l1) {
    const l1extra::L1MuonParticle& l1mu
      = toConcrete<l1extra::L1MuonParticle>(cand);
    output << " Quality: "<< setw(2) << l1mu.gmtMuonCand().quality()
	   << " Phi: "    << setw(7) << setprecision(4) << cand->phi()// + 0.0218
	   << " Eta: "    << setw(7) << setprecision(4) << cand->eta()
	   << " Pt: "     << setw(7) << setprecision(5) << cand->pt()
	   << " P: "      << setw(7) << setprecision(5) << cand->p()
	   << endl;
  }
  else if (level == l2) {
    const reco::TrackBaseRef& tk = getMainTrack(cand);
    output << " Phi: "      << setw(8) << setprecision(4) << cand->phi()
	   << "      Eta: " << setw(8) << setprecision(4) << cand->eta()
	   << endl;
    output << "   Muhits: "    << setw(3) << nHits(cand, HITS_OTH)
	   << "   Chi2/Ndof: " << setw(8) << setprecision(4) << tk->chi2()
	   << "/"              << setw(2) << tk->ndof() << endl;
    output << "   P:     " << setw(7) << setprecision(5) << cand->p()
	   << " +/- " << pError(tk)
	   << "   "
	   << "   Pt:    " << setw(7) << setprecision(5) << cand->pt()
	   << " +/- " << ptError(tk)
	   << endl;
  }
  else {
    const reco::RecoCandidate& lep = toConcrete<reco::RecoCandidate>(cand);
    const reco::TrackBaseRef glbtrk = getMainTrack(cand);
    reco::TrackRef tktrk, statrk;
    double sumptr03 = 0;
    if (!doingElectrons) {
      tktrk = lep.track();
      statrk = lep.standAloneMuon();
      const reco::Muon& mu = toConcrete<reco::Muon>(cand);
      if (mu.isIsolationValid()) sumptr03 = mu.getIsolationR03().sumPt;
    }
    
    output << " Phi: "      << setw(8) << setprecision(4) << cand->phi()
	   << "      Eta: " << setw(8) << setprecision(4) << cand->eta()
	   << endl;
    output << "   Pixhits: " << setw(3) << nHits(cand, HITS_PIX)
	   << " Silhits: "   << setw(3) << nHits(cand, HITS_SIL)
	   << " Rechits: "   << setw(3) << nHits(cand, HITS_ALL)
	   << " Seed: "      << setw(3) << seedIndex(cand)
	   << endl;
    output << "   Closest  :";
    for (int i = 0; i < MAX_LEVELS; i++) {
      output << setw(5) << id(recLevelHelper.closestLepton(cand, i))
	     << "(" << recLevelHelper.levelName(i, true) << ")";
    }
    output << endl;
    output << "   Same-seed:";
    for (int i = 0; i < MAX_LEVELS; i++) {
      output << setw(5) << id(recLevelHelper.sameSeedLepton(cand, i))
	     << "(" << recLevelHelper.levelName(i, true) << ")";
    }
    output << endl;
    output << "   P:          " << setw(7) << setprecision(5) << cand->p();
    if (!glbtrk.isNull()) 
      output << " +/- "    << pError(glbtrk);
    output << endl;
    output << "   Pt:         " << setw(7) << cand->pt();
    if (!glbtrk.isNull())
      output << " +/- "    << ptError(glbtrk)
	     << "   Chi2/Ndof: "  << setw(11) << glbtrk->chi2()
	     << "/"        << setw(2) << glbtrk->ndof();
    output << endl;
    output << "   Forward Pt: " << setw(7) << 0 //rhs.forwardPt()
	   << " +/- "    << 0 //rhs.errForwardPt()
	   << "   Back Chi2: "  << setw(11) << 0 //rhs.backChi2()
	   << endl;
    double tkpt, tkpterr, tkchi2;
    if (tktrk.isNull() || doingElectrons)
      tkpt = tkpterr = tkchi2 = 0;
    else {
      tkpt = tktrk->pt();
      tkpterr = ptError(tktrk);
      tkchi2 = tktrk->chi2();
    }
    output << "   Tracker Pt: " << setw(7) << tkpt << " +/- " << tkpterr
	   << "   Tracker Chi2: " << setw(8) << tkchi2 << endl;
    double stapt, stapterr, stachi2;
    if (statrk.isNull() || doingElectrons)
      stapt = stapterr = stachi2 = 0;
    else {
      stapt = statrk->pt();
      stapterr = ptError(statrk);
      stachi2 = statrk->chi2();
    }
    output << "   MuonFit Pt: " << setw(7) << stapt
	   << " +/- " << stapterr
	   << "   MuonFit Chi2: " << setw(8) << stachi2 << endl;
    output << "   Tracker Chi2 diff.: " << setw(7) << 0 // rhs.trackerChi2Diff()
	   << "   MuonFit Chi2 diff.: " << setw(7) << 0 // rhs.muonFitChi2Diff()
	   << endl;
    if (!glbtrk.isNull()) {
      output << "   Vertex position: " << setw(11) << cand->vx() 
	     << " "                    << setw(11) << cand->vy() 
	     << " "                    << setw(11) << cand->vz() << endl;
      output << "   Track. position: " << setw(11) << glbtrk->innerPosition().X()
	     << " "                    << setw(11) << glbtrk->innerPosition().Y()
	     << " "                    << setw(11) << glbtrk->innerPosition().Z()
	     << endl;
    }

    output << "   Closest photon: " << closestPhoton(cand)
	   << "   Sum pT (dR<0.3): " << sumptr03 << endl;
    if (!doingElectrons)
      output << "   p4 (p, E):                    " << cand->p4() << endl
	     << "   Combined track: charge: " << setw(2) << glbtrk->charge()
	     << " p: " << glbtrk->momentum() << endl
	     << "   Standalone mu : charge: " << setw(2) << statrk->charge()
	     << " p: " << statrk->momentum() << endl
	     << "   Tracker track : charge: " << setw(2) << tktrk->charge()
	     << " p: " << tktrk->momentum() << endl;
    output << "   Cut code: " << leptonIsCut(cand) << endl;
  }
}

void
Zprime2muAnalysis::dumpDilepton(ostream& output,
				const reco::CompositeCandidate& cand,
				bool dumpLeptons) const {
  output << "Dilepton: charge: " << cand.charge()
	 << " pt: " << cand.pt() << " eta: " << cand.eta()
	 << " phi: " << cand.phi() << " mass: " << cand.mass() << endl;

  int larger = dileptonDaughter(cand, 0)->p() > dileptonDaughter(cand, 1)->p() ? 0 : 1;
  int smaller = larger == 0 ? 1 : 0;
  const reco::CandidateBaseRef& cand1 = dileptonDaughter(cand, larger);
  const reco::CandidateBaseRef& cand2 = dileptonDaughter(cand, smaller);

  if (dumpLeptons) {
    output << "Higher momentum daughter:\n";
    dumpLepton(output, cand1);
    output << "Lower momentum daughter:\n";
    dumpLepton(output, cand2);
  }
  else
    output << "  Higher-p daughter: " << id(cand1)
	   << " lower-p daughter: " << id(cand2) << endl;
}

void Zprime2muAnalysis::dumpDileptonMasses() const {
  ostringstream out;

  out << "Dilepton masses at levels  ";
  for (int i_rec = 0; i_rec <= MAX_LEVELS; i_rec++) {
    out << " " << i_rec;
    if (i_rec < MAX_LEVELS) out << "      ";
  }

  out << setw(6) << setprecision(5);
  for (unsigned int i_dil = 0; i_dil < MAX_DILEPTONS; i_dil++) {
    out << "\n Dilepton # " << i_dil << "; dimu mass: ";
    for (int i_rec = 0; i_rec <= MAX_LEVELS; i_rec++) {
      if (getDileptons(i_rec).size() > i_dil)
	out << getDileptons(i_rec)[i_dil].mass();
      else
	out << " ---  ";
      if (i_rec < MAX_LEVELS) out << ", ";
    }
    out << "\n";

    out << "                res mass: ";
    for (int i_rec = 0; i_rec <= MAX_LEVELS; i_rec++) {
      if (dileptonResonances[i_rec].size() > i_dil)
	out << dileptonResonances[i_rec][i_dil].mass();
      else
	out << " ---  ";
      if (i_rec < MAX_LEVELS) out << ", ";
    }
  }

  /*
    // JMTBAD before 17X LorentzVectors are stored in XYZT instead of
    // PtEtaPhiM rep; this leads to slightly different numbers in the
    // above (and below, where the mass should be exactly that of the
    // muon mass)

  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++)
    if (allDileptons[i_rec].size() >= 1)
      out << "\n 4vecs: " << allDileptons[i_rec][0].p4() << allDileptons[i_rec][0].p4().mag()
	  << " " << dileptonDaughterByCharge(allDileptons[i_rec][0], -1)->p4()	  << " " << dileptonDaughterByCharge(allDileptons[i_rec][0], -1)->p4().mag()
	  << " " << dileptonDaughterByCharge(allDileptons[i_rec][0], +1)->p4()	  << " " << dileptonDaughterByCharge(allDileptons[i_rec][0], +1)->p4().mag() << endl;
  */

  LogTrace("dumpDileptonMasses") << out.str();
}

void Zprime2muAnalysis::dumpDilQuality() const {
  // JMTBAD not implemented
}

void Zprime2muAnalysis::dumpEvent(const bool printGen, const bool printL1,
				  const bool printL2, const bool printL3,
				  const bool printBest,
				  const bool printDileptons) const {
  unsigned imu, idil;
  int irec;
  ostringstream out;

  out << "\n******************************** Event " << eventNum << endl;

  for (irec = 0; irec < MAX_LEVELS; irec++) {
    if (irec == 0 && !printGen ||
	irec == 1 && !printL1 ||
	irec == 2 && !printL2 ||
	irec >= 3 && !printL3) continue;
    if (irec >= l3)
      out << endl;
    for (imu = 0; imu < allLeptons[irec].size(); imu++)
      dumpLepton(out, allLeptons[irec][imu]);
    if (printDileptons)
      for (idil = 0; idil < allDileptons[irec].size(); idil++)
	dumpDilepton(out, allDileptons[irec][idil]);

  }
  if (printBest) {
    out << "\nBest off-line muons: \n";
    for (imu = 0; imu < bestLeptons.size(); imu++)
      dumpLepton(out, bestLeptons[imu]);
    if (printDileptons)
      for (idil = 0; idil < bestDileptons.size(); idil++)
	dumpDilepton(out, bestDileptons[idil]);
  }
  LogTrace("") << out.str();
}

////////////////////////////////////////////////////////////////////
// Quality cuts
////////////////////////////////////////////////////////////////////

bool Zprime2muAnalysis::TrackQCheck(const reco::CandidateBaseRef& lepton,
				    const int qsel, int& ncut) const {
  // JMTBAD not implemented
  return true;
}

bool Zprime2muAnalysis::dilQCheck(const reco::Candidate& dilepton,
				  const int qsel, int& ncut_mum,
				  int& ncut_mup) const {
  // JMTBAD not implemented
  return true;
}

////////////////////////////////////////////////////////////////////
// Trigger decision queries
////////////////////////////////////////////////////////////////////

bool Zprime2muAnalysis::passTrigger(const int irec) const {
  if (irec < 0 || irec > l3) {
    throw cms::Exception("Zprime2muAnalysis")
      << "+++ passTrigger error: L" << irec << " trigger is unknown +++\n";
  }

  return passTrig[irec];
}

bool Zprime2muAnalysis::passTrigger() const {
  unsigned int decision = 1;

  for (int itrig = l1; itrig <= l3; itrig++) {
    decision *= trigWord[itrig];
  }
  return (decision != 0);
}

////////////////////////////////////////////////////////////////////
// Lepton/dilepton acceptance
////////////////////////////////////////////////////////////////////

Zprime2muAnalysis::WhereLepton
Zprime2muAnalysis::whereIsLepton(const reco::CandidateBaseRef& lepton) {
  double abseta = fabs(lepton->eta());
  if (doingElectrons) {
    // JMTBAD are these numbers taken from ECAL TDR applicable/useful?
    if      (abseta < 1.48)  return W_BARREL;
    else if (abseta < 3.0)   return W_ENDCAP;
    else                     return W_OUTSIDE;
  }
  else {
    if      (abseta < 0.9) return W_BARREL;
    else if (abseta < 1.2) return W_OVERLAP;
    else if (abseta < 2.4) return W_ENDCAP;
    else                   return W_OUTSIDE;
  }
}

Zprime2muAnalysis::WhereDilepton
Zprime2muAnalysis::whereIsDilepton(const reco::Candidate& dil) {
  WhereLepton wlep1 = whereIsLepton(dileptonDaughter(dil, 0));
  WhereLepton wlep2 = whereIsLepton(dileptonDaughter(dil, 1));
  // get a unique ordering so testing below is easier
  if (int(wlep1) > int(wlep2)) {
    WhereLepton tmp = wlep1;
    wlep1 = wlep2;
    wlep2 = tmp;
  }

  if      (wlep1 == W_BARREL  && wlep2 == W_BARREL)   return W_BARRELBARREL;
  else if (wlep1 == W_BARREL  && wlep2 == W_OVERLAP)  return W_BARRELOVERLAP;
  else if (wlep1 == W_BARREL  && wlep2 == W_ENDCAP)   return W_BARRELENDCAP;
  else if (wlep1 == W_BARREL  && wlep2 == W_OUTSIDE)  return W_BARRELOUTSIDE;
  else if (wlep1 == W_OVERLAP && wlep2 == W_OVERLAP)  return W_OVERLAPOVERLAP;
  else if (wlep1 == W_OVERLAP && wlep2 == W_ENDCAP)   return W_OVERLAPENDCAP;
  else if (wlep1 == W_OVERLAP && wlep2 == W_OUTSIDE)  return W_OVERLAPOUTSIDE;
  else if (wlep1 == W_ENDCAP  && wlep2 == W_ENDCAP)   return W_ENDCAPENDCAP;
  else if (wlep1 == W_ENDCAP  && wlep2 == W_OUTSIDE)  return W_ENDCAPOUTSIDE;
  else //if (wlep1 == W_OUTSIDE && wlep2 == W_OUTSIDE)
    return W_OUTSIDEOUTSIDE;
}

////////////////////////////////////////////////////////////////////
// Analysis level function (cuts, etc.)
////////////////////////////////////////////////////////////////////

bool Zprime2muAnalysis::eventIsInteresting() {
  return false;
  bool result = false;
  if (bestDileptons.size() > 0) {
    const reco::Candidate& dil = bestDileptons[0];
    if (dil.mass() > 800)
      result = true;
    const reco::CandidateBaseRef mum = dileptonDaughterByCharge(dil, -1);
    const reco::CandidateBaseRef mup = dileptonDaughterByCharge(dil, +1);
    if (fabs(mum->phi() - mup->phi()) < .4)
      result = true;
  }

  if (result)
    edm::LogInfo("Zprime2muAnalysis") << "'Interesting' event found!";

  return result;
}

unsigned Zprime2muAnalysis::leptonIsCut(const reco::CandidateBaseRef& lepton) const {
  unsigned retval = 0;
  
  if (doingElectrons) {
    // JMTBAD no cuts yet implemented for electrons
  }
  else {
    // pT cut of 20 GeV
    if (lepton->pt() <= 20.)
      retval |= 0x01;

    // sum pT in cone of dR=0.3 cut of 10 GeV
    if (recLevel(lepton) >= lgmr) {
      const reco::Muon& muon = toConcrete<reco::Muon>(lepton);
      if (muon.isIsolationValid() && muon.getIsolationR03().sumPt > 10)
	retval |= 0x02;
    }
  }

  return retval;
}

/*
ostream& operator<<(ostream& out, const TLorentzVector& vect) {
  out << setw(7) << vect.Eta() << " | " << setw(7) << vect.Phi() << " | "
      << setw(7) << vect.P()   << " | " << setw(7) << vect.Pt()  << " | " 
      << setw(7) << vect.Pz()  << " | " << setw(7) << vect.Rapidity() << " | "
      << setw(6) << vect.M();
  return out;
}
*/

ostream& operator<<(ostream& out, const TLorentzVector& vect) {
  out << "(" << vect.X() << "," << vect.Y() << "," << vect.Z() << "," << vect.E() << ")";
  return out;
}

ostream& operator<<(ostream& out, const reco::Candidate* par) {
  out << "GenParticle: " << " id:" << par->pdgId()
      << " p4: " << par->p4() << " status:" << par->status()
      << " vertex:" << par->vertex();
  return out;
}

#if 0

// JMTBAD methods still to port?

void Zprime2muAnalysis::matchStudy(const zp2mu::Muon& muon) {
  // Test routine, which compares match by seed with match in phi and eta.

  if (muon.isValid() == false) {
    edm::LogError("matchStudy")
      << "+++ muon does not exist! +++\n";
    return;
  }

  int muonInd = muon.id();
  int seedInd = muon.sameSeedId(ltk);
  int closInd = muon.closestId(ltk);

  if (seedInd != closInd) {
    if (seedInd != -999 && closInd == -999) {
      LogTrace("matchStudy") << "Only match by seed is found";
      LogTrace("matchStudy")
	<< "  index = " << muonInd << " seedInd = " << seedInd;
    }
    else if (seedInd == -999 && closInd != -999) {
      LogTrace("matchStudy") << "Only match in eta and phi is found";
      LogTrace("matchStudy")
	<< "  index = " << muonInd << " closInd = " << closInd;
    }
    else {
      LogTrace("matchStudy")
	<< "Difference in proximity and seed matches";
      LogTrace("matchStudy")
	<< "  index = " << muonInd << " seedInd = " << seedInd
	<< " closInd = " << closInd;
    }

    dumpEvent(true, false, false, true, false);
    
  }
}

int Zprime2muAnalysis::findMatchedDiMuonId(const int rec, 
					   const zp2mu::DiMuon& dimuon,
					   const bool debug) {
  int matched_id = -999;
  
  if (rec < 0 || rec >= MAX_LEVELS) {
    throw cms::Exception("findMatchedDiMuonId")
      << "+++ Unknown rec. level " << rec << " +++\n";
  }

  if (dimuon.recLevel() == rec)
    return dimuon.id();

  zp2mu::Muon mum = dimuon.muMinus();
  zp2mu::Muon mup = dimuon.muPlus();
  int mum_matched = mum.matchedId(rec);
  int mup_matched = mup.matchedId(rec);
  if (mum_matched != -999 && mum_matched != -999) {
    vector<zp2mu::DiMuon>::const_iterator pdi;
    for (pdi = allDiMuons[rec].begin(); pdi != allDiMuons[rec].end(); pdi++) {
      if (pdi->muMinus().id() == mum_matched &&
	  pdi->muPlus().id()  == mup_matched) {
	matched_id = pdi->id();
	break;
      }
    }
  }
  return matched_id;
}

//-----------------------------------------------------------------------------
//                             Track quality
//-----------------------------------------------------------------------------
bool Zprime2muAnalysis::TrackQCheck(const zp2mu::Muon& muon, const int qsel,
				    int& ncut, const bool debug) {
  //returns true if track ok quality
  //input is track index usuable in ntuple
  bool returnValue = false;

  int l2_id, fms_id;
  zp2mu::Muon l2mu, fmsmu;
  int l2_cut_index = 0, l3_cut_index = 0, tracker_cut_index = 0;

  if (muon.isValid() == false) {
    throw cms::Exception("TrackQCheck")
      << "+++ muon does not exist! +++\n";
  }

  // Keep track of cut number in case of failure
  ncut = 0;
  if (L3QCUT[qsel][l3_cut_index] > -998.) {
    if (muon.npixHits() <= L3QCUT[qsel][l3_cut_index]) {
      if (debug)
	LogTrace("TrackQCheck")
	  << "TrackQCheck 1 false, Npix = " << muon.npixHits();
      return false;
    }
  }
  ncut++; l3_cut_index++;

  if (L3QCUT[qsel][l3_cut_index] > -998.) {
    if (muon.nsilHits() <= L3QCUT[qsel][l3_cut_index]) {
      if (debug)
	LogTrace("TrackQCheck")
	  << "TrackQCheck 2 false, Nsil = " << muon.nsilHits();
      return false;
    }
  }
  ncut++; l3_cut_index++;

  // L2 cuts and initial cuts on L3 or ABCM L3 cuts.
  if (muon.recLevel() == l3 || muon.recLevel() == lgmr ||
      muon.recLevel() == lfms) { // TO BE CHECKED!
    l2_id = muon.matchedId(l2); // index of corresponding L2 track
    if (l2_id != -999) l2mu = muonRef(l2, l2_id);
    if (L2QCUT[qsel][l2_cut_index] > -998.) {
      if (l2_id != -999 && l2mu.nmuHits() <= L2QCUT[qsel][l2_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 3 false, L2 Nrhit = " << l2mu.nmuHits();
	return false;
      }
    }
    ncut++; l2_cut_index++;

    if (L2QCUT[qsel][l2_cut_index] > -998.) {
      if (l2_id != -999 && (l2mu.chi2()-l2mu.dof())/sqrt(2.*l2mu.dof()) >
	  L2QCUT[qsel][l2_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 4 false, L2 (Chi2-Dof)/sqrt(2.*Dof) = " 
	    << (l2mu.chi2()-l2mu.dof())/sqrt(2.*l2mu.dof());
	return false;
      }
    }
    ncut++; l2_cut_index++;

    if (L2QCUT[qsel][l2_cut_index] > -998.) {
      double sigma_tot = sqrt(muon.errTrackerPt()*muon.errTrackerPt() +
			      muon.errMuonFitPt()*muon.errMuonFitPt());
      if (abs(muon.trackerPt() - muon.muonFitPt()) > 
	  (sigma_tot*L2QCUT[qsel][l2_cut_index])) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 5 false, |trackerPt - muonPt| = |" 
	    << muon.trackerPt() << " - " << muon.muonFitPt() << "| = " 
	    << abs(muon.trackerPt()-muon.muonFitPt())
	    << ", sigma = " << sigma_tot;
	return false;
      }
    }
    ncut++; l2_cut_index++;

    if (L2QCUT[qsel][l2_cut_index] > -998.) {
      fms_id = muon.matchedId(lfms); // index of corresponding fms track
      if (fms_id != -999) {
	int    trackerDof;
	double trackerNChi2, fmsNChi2, deltaNChi2;
	double fmsProb, trackerProb;
	fmsmu = muonRef(lfms, fms_id);
	trackerDof   = 2*(fmsmu.npixHits() + fmsmu.nsilHits()) - 5;
	trackerNChi2 = (fmsmu.trackerChi2() - trackerDof)/sqrt(2.*trackerDof);
	fmsNChi2     = (fmsmu.chi2() - fmsmu.dof())/sqrt(2.*fmsmu.dof());
	deltaNChi2   = abs(fmsNChi2) - abs(trackerNChi2);
	fmsProb     = -log(TMath::Prob(fmsmu.chi2(), fmsmu.dof()));
	trackerProb = -log(TMath::Prob(fmsmu.trackerChi2(), trackerDof));

	/* double l3Prob[2], trackerProb[2], abcmProb[2];
	for (int i_muon = 0; i_muon<numMuons[3]; i_muon++) {
	trackerId = OldClosestBySeed(3, 4, i_muon);
	abcmId    = OldClosestBySeed(3, 5, i_muon);
	l3Prob[0] = TMath::Prob(Chi2[1][i_muon],     Dof[1][i_muon]);
	l3Prob[1] = TMath::Prob(backChi2[0][i_muon], Dof[1][i_muon]);
	LogTrace("TrackQCheck")
          << "@@@@@@@@@ OptimizeL3Index Debug @@@@@@@@@@@";
	LogTrace("TrackQCheck")
          << "For L3 Prob         = "   << Chi2[1][i_muon] << "/" 
	  << Dof[1][i_muon] << " = "    << l3Prob[0]       << ", Back = " 
	  << backChi2[0][i_muon] << "/" << Dof[1][i_muon]  << " = "
          << l3Prob[1];
	if (trackerId != -999) {
	  trackerProb[0] = TMath::Prob(Chi2[2][trackerId],
				       Dof[2][trackerId]);
	  trackerProb[1] = TMath::Prob(backChi2[1][trackerId],
				       Dof[2][trackerId]);
	  LogTrace("TrackQCheck")
            << "For Tracker Prob    = " << Chi2[2][trackerId] << "/" 
	    << Dof[2][trackerId]        << " = " << trackerProb[0] 
	    <<", Back = " << backChi2[1][trackerId]   << "/" 
	    << Dof[2][trackerId] << " = " << trackerProb[1];
	}
	if (abcmId != -999) {
	  abcmProb[0] = TMath::Prob(Chi2[3][abcmId],     Dof[3][abcmId]);
	  abcmProb[1] = TMath::Prob(backChi2[2][abcmId], Dof[3][abcmId]);
	  LogTrace("TrackQCheck")
            << "For ABCM L3 Prob    = " << Chi2[3][abcmId] << "/" 
	    << Dof[3][abcmId] << " = " << abcmProb[0] <<", Back = " 
	    << backChi2[2][abcmId] << "/" << Dof[3][abcmId] << " = "
            << abcmProb[1];
	}
	LogTrace("TrackQCheck")
          << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
	} */

	if (abs(fmsProb-trackerProb) > L2QCUT[qsel][l2_cut_index]) {
	//if (abs(deltaNChi2) > L2QCUT[qsel][l2_cut_index]) {
	  if (debug) {
	    //if (qsel==QSEL) {
	    ostringstream out;
	    out << "TrackQCheck 6 false: " << endl
		<< "   trackerDof = " << trackerDof
		<< ", trackerChi2 = " << fmsmu.trackerChi2() << endl
		<< "   fmsDof = "     << fmsmu.dof()
		<< ", fmsChi2 = "     << fmsmu.chi2() << endl
		<< "   trackerNChi2 = abs(" << fmsmu.trackerChi2() << "-"
		<< trackerDof  << ")/sqrt(2*" << trackerDof << ") = " 
		<< trackerNChi2 << endl
		<< "   fmsNChi2     = abs(" << fmsmu.chi2() << "-"
		<< fmsmu.dof() << ")/sqrt(2*" << fmsmu.dof() << ") = "
		<< fmsNChi2 << endl
		<< "   fmsNChi2 - trackerNChi2 = " << abs(deltaNChi2) << endl
		<< "   tracker Prob = "
		<< TMath::Prob(fmsmu.trackerChi2(),trackerDof)
		<< ", dof Prob = " << TMath::Prob(fmsmu.chi2(), fmsmu.dof())
		<< endl << "   -ln(fms Prob) = "    << fmsProb
		<< ", -ln(tracker Prob) = " << trackerProb
		<< ", diff = " << abs(fmsProb-trackerProb);
	    LogTrace("TrackQCheck") << out.str();
	  }
	  return false;
	}
      }
    }
    ncut++; l2_cut_index++;

    int dof = muon.dof();
    if (L3QCUT[qsel][l3_cut_index] > -998.) {
      double normChi2 = (muon.chi2() - dof)/sqrt(2.*dof);
      if (normChi2 > L3QCUT[qsel][l3_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 7 false, (Chi2-dof)/sqrt(2.*dof) = " << normChi2;
	return false;
      }
    }
    ncut++; l3_cut_index++;

    if (L3QCUT[qsel][l3_cut_index] > -998.) {
      double normBackChi2 = (muon.backChi2() - dof)/sqrt(2.*dof);
      if (normBackChi2 > L3QCUT[qsel][l3_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 8 false, Back Chi2 Norm = " << normBackChi2;
	return false;
      }
      ncut++; l3_cut_index++;

      if (L3QCUT[qsel][l3_cut_index] > -998.) {
	double ndiffChi2 = abs(muon.chi2() - muon.backChi2())/dof;
	if (ndiffChi2 > L3QCUT[qsel][l3_cut_index]) {
	  if (debug)
	    LogTrace("TrackQCheck")
	      << "TrackQCheck 9 false, |Chi2 - Back Chi2|/dof = " << ndiffChi2;
	  return false;
	}
      }
      ncut++; l3_cut_index++;

      if (L3QCUT[qsel][l3_cut_index] > -998.) {
	double deltaChi2 = muon.chi2() - muon.muonFitChi2() -
                             muon.trackerChi2() - muon.trackerChi2Diff();
	if (deltaChi2 > L3QCUT[qsel][l3_cut_index]) {
	  if (debug)
	    LogTrace("TrackQCheck")
	      <<"TrackQCheck 10 false, For-Muon-Tracker-TrackerChi2Diff = "
	      << deltaChi2;
	  return false;
	}
      }
      ncut++; l3_cut_index++;

      if (L3QCUT[qsel][l3_cut_index] > -998.) {
	double deltaChi2 = muon.backChi2() - muon.muonFitChi2() -
                             muon.trackerChi2() - muon.trackerChi2Diff();
	if (deltaChi2 > L3QCUT[qsel][l3_cut_index]) {
	  if (debug)
	    LogTrace("TrackQCheck")
	      <<"TrackQCheck 11 false, Back-Muon-Tracker-TrackerChi2Diff= " 
	      << deltaChi2;
	  return false;
	}
      }
      ncut++; l3_cut_index++;
    }
    else {
      ncut += 4; l3_cut_index +=4; // skip 4 cuts if no extra L3 info
    }

    if (L3QCUT[qsel][l3_cut_index] > -998.) {
      double dPoP = muon.errP()/muon.p()*(2500./(muon.p()+1000.));
      if (dPoP > L3QCUT[qsel][l3_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 12 false, (Ep/p)*(2500/(p+1000)) = " << dPoP;
	return false;
      }
    }
    ncut++; l3_cut_index++;

    if (L3QCUT[qsel][l3_cut_index] > -998.) {
      double dPtoPt =	muon.errPt()/muon.pt()*2500./(muon.pt()+1000.);
      if (dPtoPt > L3QCUT[qsel][l3_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 13 false, (EpT/pT)*(2500/(pT+1000)) = " << dPtoPt;
	return false;
      }
    }
    ncut++;
  }
  // Cuts on tracks that failed the first cut, and are replaced by the tracker
  else if (muon.recLevel() == ltk) { // TO BE CHECKED!
    ncut = NUM_L2_CUTS + NUM_L3_CUTS;
    if (TRACKERQCUT[qsel][tracker_cut_index] > -998.) {
      double normChi2 = (muon.chi2() - muon.dof())/sqrt(2.*muon.dof());
      if (normChi2 > TRACKERQCUT[qsel][tracker_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 14 false, Tracker (Chi2-Dof)/sqrt(2.*Dof) = "
	    << normChi2;
	return false;
      }
    }
    ncut++; tracker_cut_index++;

    if (TRACKERQCUT[qsel][tracker_cut_index] > -998.) {
      double diffChi2 = muon.chi2() - muon.backChi2();
      if (abs(diffChi2) > TRACKERQCUT[qsel][tracker_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 15 false, Tracker Chi2 - Back Chi2 = " << diffChi2;
	return false;
      }
    }
    ncut++; tracker_cut_index++;

    if (TRACKERQCUT[qsel][tracker_cut_index] > -998.) {
      if (muon.errP()/muon.p() > TRACKERQCUT[qsel][tracker_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 16 false, Tracker Ep/p = " << muon.errP()/muon.p();
	return false;
      }
    }
    ncut++; tracker_cut_index++;
    
    if (TRACKERQCUT[qsel][tracker_cut_index] > -998.) {
      if (muon.errPt()/muon.pt() > TRACKERQCUT[qsel][tracker_cut_index]) {
	if (debug)
	  LogTrace("TrackQCheck")
	    << "TrackQCheck 17 false, Tracker pT*EpT = "
	    << muon.errPt()/muon.pt();
	return false;
      }
    }
  }

  ncut = -999;
  returnValue = true;
  return returnValue;
}

bool Zprime2muAnalysis::dilQCheck(const zp2mu::DiMuon& dimuon, const int qsel,
				  int& ncut_mum, int& ncut_mup,
				  const bool debug) {
  //returns true if ok quality

  if (dimuon.isValid() == false) {
    throw cms::Exception("dilQCheck")
      << "+++ dimuon does not exist! +++\n";
  }
  // Only does anything for rec >= 3.
  // JMTBAD this skips lpmr... do we want to do that?
  if (dimuon.recLevel() < l3 || dimuon.recLevel() > lfms) {
    edm::LogWarning("dilQCheck")
      << "+++ dimuon recLevel = " << dimuon.recLevel() << " +++\n";
    return false;
  }

  return (TrackQCheck(dimuon.muMinus(), qsel, ncut_mum, debug) &&
	  TrackQCheck(dimuon.muPlus(),  qsel, ncut_mup, debug));
}

void Zprime2muAnalysis::dumpDilQuality() {
  // Dumps to help analyze discrepancies between generated and reconstructed
  // dileptons.  Dump shows differences between various fits.
  int    gendi_id, gen_id, l2_id;
  double residual, normChi2, scale, scaledError;
  string mus[2] = {"mu-", "mu+"};
  vector<zp2mu::DiMuon>::const_iterator pdi;
  zp2mu::Muon     muon, genmu, l2mu;

  ostringstream strstrm;

  for (pdi = allDiMuons[0].begin(); pdi != allDiMuons[0].end(); pdi++)
    strstrm << "Gen dimuon # "  << pdi->id()
	    << ", mass = "      << pdi->dimuV().M()
	    << ", Ptgen mu- = " << pdi->muMinus().pt()
	    << ", mu+ = "       << pdi->muPlus().pt() << endl
	    << "                                eta mu- = "
	    << pdi->muMinus().eta()
	    << ", mu+ = "       << pdi->muPlus().eta() << endl;

  for (pdi = allDiMuons[l3].begin(); pdi != allDiMuons[l3].end(); pdi++) {
    strstrm << "L3  dimuon # "  << pdi->id()
	    << ", mass = "      << pdi->dimuV().M()
	    << ", Ptrec mu- = " << pdi->muMinus().pt()
	    << ", mu+ = "       << pdi->muPlus().pt() << endl;
    gendi_id = findMatchedDiMuonId(0, *pdi);
    if (gendi_id != -999) {
      zp2mu::DiMuon gendi = dimuonRef(0, gendi_id);
      strstrm << "                mass res = "
	      << (pdi->dimuV().M()-gendi.dimuV().M())/gendi.dimuV().M()
	      << endl;
    }
    else
      strstrm << " matched generated dimuon is not found " << endl;
  }

  for (pdi = allDiMuons[ltk].begin(); pdi != allDiMuons[ltk].end(); pdi++) {
    strstrm << "TK  dimuon # "  << pdi->id()
	    << ", mass = "      << pdi->dimuV().M()
	    << ", Ptrec mu- = " << pdi->muMinus().pt()
	    << ", mu+ = "       << pdi->muPlus().pt() << endl;
    gendi_id = findMatchedDiMuonId(0, *pdi);
    if (gendi_id != -999) {
      zp2mu::DiMuon gendi = dimuonRef(0, gendi_id);
      strstrm << "                mass res = "
	      << (pdi->dimuV().M()-gendi.dimuV().M())/gendi.dimuV().M()
	      << endl;
    }
    else
      strstrm << ", matched generated dimuon is not found " << endl;
  }

  strstrm << endl << "                   pT           "
	  << "   Chi2/Dof    normChi2    1/pT Res   scaledError" << endl;
  for (int i_part = 0; i_part < 2; i_part++) {
    strstrm << mus[i_part] << endl;

    for (pdi = allDiMuons[l3].begin(); pdi != allDiMuons[l3].end(); pdi++) {
      if (i_part == 0) muon = pdi->muMinus();
      else             muon = pdi->muPlus();

      // Level 2
      l2_id = muon.matchedId(l2);
      if (l2_id != -999) {
	l2mu = muonRef(l2, l2_id);
	normChi2 = (l2mu.chi2()-l2mu.dof())/sqrt(2.*l2mu.dof());
	strstrm << "  L2      " << l2mu.pt() << " +/- "
		<< l2mu.errPt() << "    " << l2mu.chi2() << "/"
		<< l2mu.dof()   << "    " << normChi2 << endl;
      }

      // Level 3
      normChi2 = (muon.chi2()-muon.dof())/sqrt(2.*muon.dof());
      strstrm << "  GMR     " << muon.pt() <<" +/- "
	      << muon.errPt() << "    " << muon.chi2() << "/" 
	      << muon.dof()   << "    " << normChi2 << "    ";
      gen_id = muon.matchedId(0);
      if (gen_id != -999) {
	genmu = muonRef(0, gen_id);
	residual = (1./muon.pt() - 1./genmu.pt())/(1./genmu.pt());
	strstrm << residual;
      }
      scale = muon.chi2()/muon.dof();
      if (scale < 1.) scale = 1.;
      scaledError = (muon.errPt()/muon.pt())*scale;
      strstrm << "     " << scaledError << endl;
    }

    for (pdi = allDiMuons[ltk].begin(); pdi != allDiMuons[ltk].end(); pdi++) {
      if (i_part == 0) muon = pdi->muMinus();
      else             muon = pdi->muPlus();
      normChi2 = (muon.chi2()-muon.dof())/sqrt(2.*muon.dof());
      strstrm << "  Tracker " << muon.pt() <<" +/- "
	      << muon.errPt() << "    " << muon.chi2() << "/"
	      << muon.dof()   << "    " << normChi2 << "    ";
      gen_id = muon.matchedId(0);
      if (gen_id != -999) {
	genmu = muonRef(0, gen_id);
	residual = (1./muon.pt() - 1./genmu.pt())/(1./genmu.pt());
	strstrm << residual;
      }
      scale = muon.chi2()/muon.dof();
      if (scale < 1.) scale = 1.;
      scaledError = (muon.errPt()/muon.pt())*scale;
      strstrm << "     " << scaledError << endl;
    }

    for (pdi = allDiMuons[lfms].begin(); pdi != allDiMuons[lfms].end(); pdi++) {
      if (i_part == 0) muon = pdi->muMinus();
      else             muon = pdi->muPlus();
      normChi2 = (muon.chi2()-muon.dof())/sqrt(2.*muon.dof());
      strstrm << "  FMS     " << muon.pt() <<" +/- "
	      << muon.errPt() << "    " << muon.chi2() << "/"
	      << muon.dof()   << "    " << normChi2 << "    ";
      gen_id = muon.matchedId(0);
      if (gen_id != -999) {
	genmu = muonRef(0, gen_id);
	residual = (1./muon.pt() - 1./genmu.pt())/(1./genmu.pt());
	strstrm << residual;
      }
      scale = muon.chi2()/muon.dof();
      if (scale < 1.) scale = 1.;
      scaledError = (muon.errPt()/muon.pt())*scale;
      strstrm << "     " << scaledError << endl;
    }
  }
  strstrm << endl << "-------------------------";

  LogTrace("dumpDilQuality") << strstrm.str();
}

#endif
