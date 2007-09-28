//
// Authors: Jason Mumford, Jordan Tucker, Slava Valuev, UCLA
//

#include <algorithm>

#include "TH1.h"
#include "TStyle.h"

#include "CLHEP/GenericFunctions/CumulativeChiSquare.hh"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

using namespace std;

// muon mass in GeV/c^2
const double Zprime2muAnalysis::MUMASS = 0.10566;
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

Zprime2muAnalysis::Zprime2muAnalysis(const edm::ParameterSet& config) 
  : eventNum(-1) {
  // verbosity controls the amount of debugging information printed;
  // levels are defined in an enum
  verbosity = VERBOSITY(config.getUntrackedParameter<int>("verbosity", 0));
  // JMTBAD probably a better way to check this
  doingHiggs = config.getParameter<bool>("doingHiggs");
  // whether only to look at generator-level muons (i.e. don't bother
  // trying to store globalMuons, standAloneMuons, etc)
  generatedOnly = config.getParameter<bool>("generatedOnly");

  // input tags for the reco collections we need
  standAloneMuons = config.getParameter<edm::InputTag>("standAloneMuons");
  genMuons = config.getParameter<edm::InputTag>("genMuons");
  globalMuons = config.getParameter<edm::InputTag>("globalMuons");
  globalMuonsFMS = config.getParameter<edm::InputTag>("globalMuonsFMS");
  globalMuonsPMR = config.getParameter<edm::InputTag>("globalMuonsPMR");
  photons = config.getParameter<edm::InputTag>("photons");

  InitROOT();
}

Zprime2muAnalysis::~Zprime2muAnalysis() {
}

void Zprime2muAnalysis::InitROOT() {
  TH1::AddDirectory(false);
  gROOT->SetStyle("Plain");
  gStyle->SetFillColor(0);
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
  for (int i_rec = 0; i_rec < NUM_REC_LEVELS; i_rec++) {
    passTrig[i_rec] = true;
    calcTrig[i_rec] = false;
    trigWord[i_rec] = 0;
  }
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    searchedDileptons[i_rec] = false;
    allMuons[i_rec].clear();
    allDiMuons[i_rec].clear();
  }
  bestMuons.clear();
  bestDiMuons.clear();
}

// Store generated muons.
// JMTBAD we store info in our own Muon
// classes, since for now, we want to store things not available in
// some of the reco:: classes (e.g. an "id" number of our own, etc.) 
// also, because of the deriving hierarchy, it is hard to lump
// gracefully reco::GenParticleCandidates, reco::Muons, l1Muons, etc
// so for now, hack out and reuse the Muon, DiMuon classes this has
// the side-effect of making the dilepton reconstruction under our
// control... but: There are existing modules that can do all sorts of
// selection on the HepMCProduct -- picking only muons, filtering by
// minimum pT, etc. -- so that we don't have to iterate over all
// candidates ourselves.  The facility also exists to make candidate
// dimuons... useful later?  see
// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCandidateModules
// and
// https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookGenParticleCandidate
void Zprime2muAnalysis::storeGeneratedMuons(const edm::Event& event) {
  static bool debug = verbosity >= VERBOSITY_TOOMUCH;
  int imu = 0;

  // JMTBAD as far as I can see, GenParticleCandidateProducer doesn't
  // preserve anything like a PYTHIA line index for easy
  // distinguishing of, say in the case of H->ZZ->4mu, which Z the
  // muon came from (is there an == operator?). use "raw" HepMCProduct
  // for now, which has this in the form of a decay vertex id
  // basic info on how to use at 
  // https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookGeneration and
  // http://lcgapp.cern.ch/project/simu/HepMC/refman/HepMC.2.00.02/annotated.html

  // Generated particles (i.e. PYTHIA).
  edm::Handle<reco::CandidateCollection> genParticles;
  event.getByLabel("genParticleCandidates", genParticles);

  // Simulated tracks (i.e. GEANT particles).
  edm::Handle<edm::SimTrackContainer> simtracks;
  if (!generatedOnly)
    event.getByLabel("g4SimHits", simtracks);

  if (debug)
    LogTrace("storeGeneratedMuons") << "\n";  

  // First, we loop over generated (PYTHIA) particles and take all muons.
  // Next, we loop over GEANT particles and add extra muons produced in
  // secondary interactions and decays (not in PYTHIA).  The reason for
  // doing it this way instead of getting *all* muons from GEANT list
  // is that apparently only particles within a certain eta range
  // (|eta| < 5 or so) are saved in GEANT container (why?), and a
  // fraction of muons from high-mass resonance decays would not be found.
  int idx = -1;
  for (reco::CandidateCollection::const_iterator ipl = genParticles->begin();
       ipl != genParticles->end(); ipl++) {
    idx++;
    if ((ipl)->status() == 1) { // skip documentation lines
      int id = (ipl)->pdgId();
      if (abs(id) != 13) continue;
      const reco::Particle::Vector& pmu = (ipl)->momentum();
      int charge = -id/abs(id);
      if (debug)
	LogTrace("storeGeneratedMuons")
	  << "Generated muon: pt = " << pmu.Rho() << " eta = " << pmu.Eta()
	  << " q = " << charge << " ind = " << idx;

      // Find its mother and grandmother.
      const reco::Candidate* genmuon = (&*ipl);
      int motherInd = 0;
      int motherId  = 0, grandmaId = 0;
      while (motherInd == 0) {
	const reco::Candidate* mother = genmuon->mother();
	if (mother == 0) {
	  edm::LogWarning("storeGeneratedMuons")
	    << "+++ empty mother for particle #" << idx << " +++\n";
	  break;
	}
	else {
	  int mId = mother->pdgId();

	  // This loop recovers the functionality of barcode() method.
	  int mother_idx = -1;
	  for (unsigned int j = 0; j < genParticles->size(); ++j) {
	    const reco::Candidate* ref = &((*genParticles)[j]);
	    if (mother->px() == ref->px() && mother->py() == ref->py() &&
		mother->pz() == ref->pz() &&
		mother->status() == ref->status() && mId == ref->pdgId()) {
	      mother_idx = j;
	      break;
	    }
	  }

	  if (debug) 
	    LogTrace("storeGeneratedMuons")
	      << "Mother: id = " << mId << " pt = "
	      << mother->momentum().Rho() << " M = " << mother->mass()
	      << " line = " << mother_idx;
	  // If muon's mother is a muon, go up the decay chain.  This is
	  // mainly to avoid stopping at muons in documentation lines,
	  // which are declared to be ancestors of muons produced in
	  // hard interaction (i.e., ancestors of themselves).
	  if (abs(mId) == 13) {
	    genmuon = mother;
	  }
	  else {
	    motherInd = mother_idx;
	    motherId  = mId;

	    // Also look for the muon's grandmother so we can
	    // sort by type (gg/qqbar) for graviton stuff.
	    std::vector<const reco::Candidate*> grandmas;
            for (unsigned int igm = 0; igm != mother->numberOfMothers(); igm++) 
	      grandmas.push_back(mother->mother(igm)); 

	    // Could loop over list of parents, but assume
	    // that using the first grandma (q/qbar or g) works
	    // just as well as the other.
	    const reco::Candidate* grandma = grandmas[0];
	    if (grandma != 0) {
	      grandmaId = grandma->pdgId();
	    }
	    break;
	  }
	}
      }

      if (motherInd == 0) {
	edm::LogWarning("storeGeneratedMuons")
	  << "+++ cannot find motherId for particle #" << idx << " +++\n";

	for (reco::CandidateCollection::const_iterator ipl = genParticles->begin();
	     ipl != genParticles->end(); ipl++) {
	  edm::LogWarning("storeGeneratedMuons")
	    << "Status = "  << (ipl)->status()
	    << " id = "     << (ipl)->pdgId()
	    << " line = "   << idx
	    << " pT = "     << (ipl)->momentum().Rho()
	    << " eta = "    << (ipl)->momentum().Eta()
	    << " M = "      << (ipl)->mass();
	}
      }

      zp2mu::Muon thisMu;
      thisMu.fill(1, imu, lgen, charge,
		  pmu.phi(), pmu.eta(), pmu.rho(),  pmu.R());
      thisMu.setMotherId(motherId);
      thisMu.setMotherLine(motherInd);
      thisMu.setGrandmaId(grandmaId);
      thisMu.setVertexXYZ((*ipl).vx(), (*ipl).vy(), (*ipl).vz()); 
      allMuons[lgen].push_back(thisMu);
      imu++;
    }
  }

  if (!generatedOnly) {
    // Add GEANT muons which were not in PYTHIA list.
    for (edm::SimTrackContainer::const_iterator isimtrk = simtracks->begin();
	 isimtrk != simtracks->end(); ++isimtrk) {
      // Store only muons above a certain pT threshold.
      if (abs(isimtrk->type()) == 13) {
	//HepMC::FourVector pmu = (*cand)->momentum();
	const CLHEP::HepLorentzVector& pmu = isimtrk->momentum();
	float pt  = pmu.perp();
	if (pt > PTMIN) {
	  int pythiaInd = isimtrk->genpartIndex();

	  if (pythiaInd != -1) continue; // Skip PYTHIA muons

	  float eta = pmu.pseudoRapidity();
	  int charge = static_cast<int>(isimtrk->charge());
	  if (debug)
	    LogTrace("storeGeneratedMuons")
	      << "Simulated muon: pt = " << pt << " eta = " << eta
	      << " q = " << charge << " ind = " << pythiaInd;

	  zp2mu::Muon thisMu;
	  thisMu.fill(1, imu, lgen, charge,
		      pmu.phi(), pmu.eta(), pmu.perp(), pmu.rho());
	  allMuons[lgen].push_back(thisMu);
	  imu++;
	}
	else {
	  edm::LogWarning("storeGeneratedMuons")
	    << "skipping Gen muon with pT = " << pt << " p = " << pmu.rho()
	    << " eta " << pmu.eta();
	}
      }
    }
  }
}

// JMTBAD ptError has disappeared? 
// http://cms.cern.ch/iCMS/jsp/openfile.jsp?type=NOTE&year=2006&files=NOTE2006_001.pdf
// but read:
// https://hypernews.cern.ch/HyperNews/CMS/get/recoTracking/150.html?inline=-1
// (incorrect formula that doesn't include covariance that I found in the lxr, too)

// so: utility function to calculate error on 1/pt from the covariance matrix
// actually stored in the track info
double Zprime2muAnalysis::invPtError(const reco::TrackRef& track) {
  // pt = sqrt(px^2 + py^2)
  //    = cos(lambda) * q / qoverp
  // with qoverp = q/|p|

  double qoverp = track->qoverp();
  double lambda = track->lambda();
  double q = track->charge();

  // need the derivatives to propagate error
  double dpt_dlambda = -sin(lambda)*q/qoverp;
  double dpt_dqoverp = -cos(lambda)*q/qoverp/qoverp;

  const int i_lambda = reco::TrackBase::i_lambda;
  const int i_qoverp = reco::TrackBase::i_qoverp;

  // use standard error propagation formula
  double var_pt =
    track->covariance(i_lambda, i_lambda)   * dpt_dlambda * dpt_dlambda +
    track->covariance(i_qoverp, i_qoverp)   * dpt_dqoverp * dpt_dqoverp +
    2*track->covariance(i_lambda, i_qoverp) * dpt_dlambda * dpt_dqoverp;

  double sig_pt = sqrt(var_pt);
  double pt = track->pt();

  // further propagate the error to 1/pt
  double sig_invpt = 1/pt/pt*sig_pt;
  return sig_invpt;
}

// JMTBAD calculate error on p ourselves, too
double Zprime2muAnalysis::invPError(const reco::TrackRef& track) {
  // dumb identity:
  // p = q / qoverp
  // dp_dqoverp = - q/qoverp^2

  double qoverp = track->qoverp();
  double q = track->charge();
  
  double sig_p = track->qoverpError() * q/qoverp/qoverp;
  if (sig_p < 0) sig_p = -sig_p;

  // further propagate the error to 1/p
  double p = track->p();
  double sig_invp = 1/p/p*sig_p;
  return sig_invp;
}

// Return the index into the standAloneMuons collection staTracks that
// this Muon's standAloneMuon (passed in as track) most closely matches 
// -- use to mock up seedIndex
// JMTBAD there's probably a better way to do this, since the collections for
// GMR, TK, FMS, PMR all seem sorted in the same way
int Zprime2muAnalysis::matchStandAloneMuon(
                         const edm::Handle<reco::TrackCollection> staTracks,
                         const reco::TrackRef& track) {
  int closest_ndx = -1;
  double min_mag2 = 1e99;
  const double MAX_MAG2 = 1; // more than enough for a comparison between doubles
  int ndx = 0;
  for (reco::TrackCollection::const_iterator stamu = staTracks->begin();
       stamu != staTracks->end(); stamu++) {
    if (stamu->charge() == track->charge()) {
      // calling mag2 results in complaints from the linker, even though I'm linking
      // against ROOT... grr
      double mag2 = (stamu->momentum() - track->momentum()).mag2();
      if (mag2 < min_mag2) {
	closest_ndx = ndx;
	min_mag2 = mag2;
      }
    }
    ndx++;
  }
  if (min_mag2 > MAX_MAG2) {
    edm::LogWarning("matchStandAloneMuon")
      << "could not match stand-alone muon!\n"
      << "muon 3-momentum p = " << track->momentum();
    return -999;
  }
  else
    return closest_ndx;
}

// utility function to match muons from same vertex
//
bool Zprime2muAnalysis::haveSameVertex(const zp2mu::Muon& lhs, const zp2mu::Muon& rhs) const {
  return ( (lhs.vertex() - rhs.vertex()).R() < 1e-6 );

}

bool Zprime2muAnalysis::haveSameVertex(const zp2mu::Muon& lhs, const math::XYZPoint& rhs) const {
  return ( (lhs.vertex() - rhs).R() < 1e-6 );

}

// check if a GenParticleCandidate is a mother of 2 zp2mu::Muons 
// temporary fix to work as checking barcode

bool Zprime2muAnalysis::isMotherOf(const reco::GenParticleCandidate& particle, const zp2mu::Muon& mup, const zp2mu::Muon& mum) const {
  if ( particle.numberOfDaughters() < 2 ) return false;
  if ( particle.charge() != 0 ) return false; 

  bool mupOK = false;
  bool mumOK = false;

  for (unsigned int igd = 0; igd != particle.numberOfDaughters(); igd++) {
      if ( abs(particle.daughter(igd)->pdgId()) != 13 ) continue;
      if ( particle.daughter(igd)->pdgId() == 13 ) {
         math::XYZVector mom = particle.daughter(igd)->momentum();
         if ( (fabs(mom.Eta() - mum.eta()) < 1e-5)
            &&(fabs(mom.Phi() - mum.phi()) < 1e-5)
            &&(fabs(mom.Rho() - mum.pt()) < 1e-5) ) mumOK = true;
      }
      if ( particle.daughter(igd)->pdgId() != -13 ) {
         math::XYZVector mom = particle.daughter(igd)->momentum();
         if ( (fabs(mom.Eta() - mup.eta()) < 1e-5)
            &&(fabs(mom.Phi() - mup.phi()) < 1e-5)
            &&(fabs(mom.Rho() - mup.pt()) < 1e-5) ) mumOK = true;
      }
  }
  return ( mupOK && mumOK );
}


// Fill a vector of Level-2 muons.  Use offline standalone muons for now
// instead of true Level-2 muons (seeded by Level-1 tracks).
void Zprime2muAnalysis::storeL2Muons(const edm::Event& event) {

  // Null reference to tracker-only tracks.
  reco::TrackRef tkTrack;

  // Standalone muons.
  edm::Handle<reco::TrackCollection> staTracks;
  event.getByLabel(standAloneMuons, staTracks);

  // Empty collection of photons since they are not applicable.
  reco::PhotonCollection photonCollection;

  int imu = 0;
  // Loop over tracks.
  for (unsigned int ista = 0; ista < staTracks->size(); ++ista) {
    // Build the reference to the standalone track.
    reco::TrackRef muTrack = reco::TrackRef(staTracks, ista);

    // Do the rest only if muTrack is not empty.
    if (muTrack.isNonnull()) {
      bool isStored = storeOfflineMuon(imu, l2, muTrack, tkTrack, muTrack,
				       ista, photonCollection);
      if (isStored) imu++;
    }
  }
}

// Fill a vector of off-line muons.
//  whichMuons is one of the InputTags we stored earlier, e.g. globalMuonsFMS
//    for the muons fit using the tracker+first muon station
//  irec is our rec level (e.g. lfms)
//  trackerOnly dictates whether to use the track reconstructed with tracker
//    information only (default false)
void Zprime2muAnalysis::storeOfflineMuons(const edm::Event& event, 
					  const edm::InputTag& whichMuons,
					  RECLEVEL irec, bool trackerOnly) {
  edm::Handle<reco::MuonTrackLinksCollection> muons;
  event.getByLabel(whichMuons, muons);

  // we want to look at the standalone tracks as seeds for the different
  // reconstructed muons -- can use the index into the collection as a 
  // seedIndex?
  edm::Handle<reco::TrackCollection> staTracks;
  event.getByLabel(standAloneMuons, staTracks);

  // Get collection of "corrected" (calibrated?) photons.
  // Should we move this block one level up, to avoid calling if for every
  // muon collection?
  edm::Handle<reco::PhotonCollection> thePhotons;
  event.getByLabel(photons, thePhotons);
  const reco::PhotonCollection photonCollection = *(thePhotons.product());
  //LogTrace("storeOfflineMuons")
  //  << "Found " << photonCollection.size() << " photons";

  int imu = 0;
  for (reco::MuonTrackLinksCollection::const_iterator muon = muons->begin();
       muon != muons->end(); muon++) {
    // Null references to combined, tracker-only, and standalone muon tracks.
    reco::TrackRef theTrack;
    reco::TrackRef tkTrack;
    reco::TrackRef muTrack;

    // Fill (some of) them depending on the flag.
    if (trackerOnly) {
      theTrack = muon->trackerTrack();
      tkTrack  = theTrack;
    }
    else {
      theTrack = muon->globalTrack();
      tkTrack  = muon->trackerTrack();
      muTrack  = muon->standAloneTrack();
    }

    // Do the rest only if theTrack is not empty.
    if (theTrack.isNonnull()) {
      // seedIndex is now the index into the standAloneMuon track collection
      int seedIndex = matchStandAloneMuon(staTracks, muon->standAloneTrack());

      bool isStored = storeOfflineMuon(imu, irec, theTrack, tkTrack, muTrack,
				       seedIndex, photonCollection);
      if (isStored) imu++;
    }
  }
}

bool Zprime2muAnalysis::storeOfflineMuon(const int imu, const RECLEVEL irec,
					 const reco::TrackRef& theTrack,
					 const reco::TrackRef& tkTrack,
					 const reco::TrackRef& muTrack,
					 const int seedIndex,
			 const reco::PhotonCollection& photonCollection) {
  bool isStored = false;

  // Skip muons with non-positive number of dof.
  // JMTBAD is this still relevant (does GMR keep muons with dof < 0?)
  if (theTrack->ndof() <= 0) {
    edm::LogWarning("storeOfflineMuon")
      << "muon's dof is " << theTrack->ndof() << "; skipping it";
    return isStored;
  }

  // Use momentum measured at vertex (taken from the interaction point).
  // JMTBAD assume this is what the reco people did in storing the momentum
  // in the muon class
  double pt = theTrack->pt();
  if (pt > PTMIN) {
    zp2mu::Muon thisMu;
    thisMu.fill(1, imu, irec, theTrack->charge(), theTrack->phi(),
		theTrack->eta(), pt, theTrack->p());

    int nRecHits = theTrack->recHitsSize();
    int nPixHits = 0, nSilHits = 0, nMuonHits = 0;
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
	  edm::LogWarning("storeOfflineMuon")
	    << "unknown subid for TrackingRecHit";
	}
      }
      else if (det == DetId::Muon)
	nMuonHits++;
    }
    // redundant check, hopefully
    if (nPixHits + nSilHits + nMuonHits != nRecHits)
      edm::LogWarning("storeOfflineMuon")
	<< "TrackingRecHits from outside the tracker or muon systems!";
    // last arg is sum of pixel, sil. and muon hits
    thisMu.setHits(nPixHits, nSilHits, nRecHits);

    // JMTBAD ndof for TrackBase is a double...
    thisMu.setDof(int(theTrack->ndof()));
    thisMu.setChi2(theTrack->chi2());

    // seedIndex is now the index into the standAloneMuon track collection
    thisMu.setSeed(seedIndex);

    // JMTBAD backward fit?
    //thisMu.setBackChi2((*il3).BackChi2());

    if (tkTrack.isNonnull()) {
      // JMTBAD difference between forward and backward?
      thisMu.setTrackerChi2(tkTrack->chi2(), 0);
      thisMu.setTrackerPt(tkTrack->pt());
      // store eInvPt also for the muon reconstructed in the tracker alone
      thisMu.setETrackerInvPt(invPtError(tkTrack));
    }

    if (muTrack.isNonnull()) {
      thisMu.setMuonFitChi2(muTrack->chi2(), 0);
      thisMu.setMuonFitPt(muTrack->pt());
      // store eInvPt also for the standalone muon (just muon systems)
      thisMu.setEMuonFitInvPt(invPtError(muTrack));
    }

    // JMTBAD forward fit
    //thisMu.setForwardPt((*il3).ForwardPt());
    //thisMu.setEForwardInvPt((*il3).ForwardEInvPt());

    thisMu.setEInvPt(invPtError(theTrack));
    thisMu.setEInvP(invPError(theTrack));

    // set the vertex found by the track fit
    thisMu.setVertexXYZ(theTrack->vx(), theTrack->vy(), theTrack->vz());
    thisMu.setTrackerXYZ(theTrack->innerPosition().X(),
			 theTrack->innerPosition().Y(),
			 theTrack->innerPosition().Z());

    if (verbosity >= VERBOSITY_TOOMUCH)
      LogTrace(" ") << "irec = " << irec;
    TLorentzVector photon = findClosestPhoton(theTrack, photonCollection);
    thisMu.setClosestPhoton(photon);

    allMuons[irec].push_back(thisMu);
    isStored = true;
  }
  else {
    edm::LogWarning("storeOfflineMuon")
      << "skipping muon at rec level " << irec << " with pT = " << pt
      << " p = " << theTrack->p() << " eta " << theTrack->eta();
  }

  return isStored;
}

TLorentzVector Zprime2muAnalysis::findClosestPhoton(
                              const reco::TrackRef& theTrack,
			      const reco::PhotonCollection& photonCollection) {
  static bool debug = verbosity >= VERBOSITY_TOOMUCH;
  double phdist = 999.0;
  reco::Particle::LorentzVector phclos; // 4-momentum of the closest PhotonCandidate

  double muon_eta = theTrack->eta();
  double muon_phi = theTrack->phi();

  for (reco::PhotonCollection::const_iterator iph = photonCollection.begin();
       iph != photonCollection.end(); iph++) {

    reco::Photon thePhoton(*iph);

    // Photon direction can be improved by taking into account precise
    // position of the primary vertex - leave for later.
    // thePhoton.setVertex(vtx);

    double phot_eta = thePhoton.eta();
    double phot_phi = thePhoton.phi();
    double dist = deltaR(muon_eta, muon_phi, phot_eta, phot_phi);

    if (debug) LogTrace("findClosestPhoton")
      << "Photon: eta = " << phot_eta << " phi = " << phot_phi
      << " energy = " << thePhoton.energy() << " dR = " << dist;

    if (dist < phdist) {
      phclos = thePhoton.p4();
      phdist = dist;
    }
  }

  if (debug) LogTrace("findClosestPhoton")
    << "Closest photon: eta = " << phclos.eta() << " phi = " << phclos.phi()
    << " px = " << phclos.px() << " py = " << phclos.py()
    << " pz = " << phclos.pz() << " energy = " << phclos.e()
    << " dR  = " << phdist;

  TLorentzVector photon(phclos.px(), phclos.py(), phclos.pz(), phclos.e());
  return photon;
}

double Zprime2muAnalysis::deltaR(const double eta1, const double phi1,
				 const double eta2, const double phi2) const {
  double dEta = eta1 - eta2;
  double dPhi = fabs(phi1 - phi2);
  while (dPhi > M_PI) dPhi -= 2*M_PI;
  return sqrt(dEta*dEta + dPhi*dPhi);
}

// Take the info from the event and store the generated muons and
// muons at various levels of reconstruction in our own vectors for
// sorting and producing dileptons
// does the duty of what was done in MuExAnalysis::initialize()
void Zprime2muAnalysis::storeMuons(const edm::Event& event) {
  clearValues();
  // Generated particles: GEANT and PYTHIA info
  storeGeneratedMuons(event);
  if (!generatedOnly) {
    // Standalone muons; should be replaced by Level-2 muons
    storeL2Muons(event);
    // off-line muons reconstructed with GlobalMuonProducer
    storeOfflineMuons(event, globalMuons, lgmr);
    // store the muons using the tracker-only tracks
    storeOfflineMuons(event, globalMuons, ltk, true);
    // muons using tracker+first muon station
    storeOfflineMuons(event, globalMuonsFMS, lfms);
    // muons using the PickyMuonReconstructor
    storeOfflineMuons(event, globalMuonsPMR, lpmr);
  }

  // Sort the list of found muons using overloaded criteria.
  // (by momentum magnitude)
  for (int irec = 1; irec < MAX_LEVELS; irec++) {
    if (allMuons[irec].size() > 1) {
      sort(allMuons[irec].begin(), allMuons[irec].end(),
	   greater<zp2mu::BaseMuon>());
    }
  }

  // Find closest (and same seed, if available) muons among the generated
  // muons and at other rec. levels.
  matchAllMuons();

  // Select the "best"-reconstructed muons from the four available options
  // and sort them.
  bestMuons = findBestMuons();
  if (bestMuons.size() > 1) {
    sort(bestMuons.begin(), bestMuons.end(), greater<zp2mu::BaseMuon>());
  }

  // Dump all muon vectors
  bool debug = verbosity >= VERBOSITY_TOOMUCH;
  if (debug) dumpEvent(true, true, true, true, true);

  // Construct opposite-sign dileptons out of lists of muons.
  // Pass ref. to Event only to save true resonance; need to think of a better
  // solution?
  makeAllDileptons(event, debug);
}

// Cocktail of tracks from various reconstructors.  Two possibilities:
// Piotr's tune (default) or Norbert's tune.  Both are available in CMSSW
// starting with 1_4_0_pre2, so we could switch to official versions once
// we start using 1_4_0 (after making sure that they agree with our version).
vector<zp2mu::Muon> Zprime2muAnalysis::findBestMuons() {
  bool debug = verbosity >= VERBOSITY_TOOMUCH;
  vector<zp2mu::Muon> tmpV;
  vector<zp2mu::Muon>::const_iterator pmu;
  static int ntrk = 0, ngmr = 0, nfms = 0, npmr = 0, ngpr = 0, ntot = 0;

  // Start with tracker-only tracks, as Norbert does.
  for (pmu = allMuons[ltk].begin(); pmu != allMuons[ltk].end(); pmu++) {
    zp2mu::Muon trkmu, gmrmu, fmsmu, pmrmu;
    int gmr_id = pmu->sameSeedId(lgmr);
    int fms_id = pmu->sameSeedId(lfms);
    int pmr_id = pmu->sameSeedId(lpmr);
    trkmu = *pmu;
    if (gmr_id != -999) gmrmu = muonRef(lgmr, gmr_id);
    if (fms_id != -999) fmsmu = muonRef(lfms, fms_id);
    if (pmr_id != -999) pmrmu = muonRef(lpmr, pmr_id);
    //zp2mu::Muon bestMu = NorbertsCocktail(trkmu, gmrmu, fmsmu, pmrmu, debug);
    zp2mu::Muon bestMu = PiotrsCocktail(trkmu, fmsmu, pmrmu, debug);
    if (debug) {
      if (bestMu == trkmu) {
	ntrk++; LogTrace("findBestMuons") << "  --> select Tracker only";
      }
      else if (bestMu == gmrmu) {
	ngmr++; LogTrace("findBestMuons") << "  --> select GMR";
	if (gmrmu.isValid() && pmrmu.isValid() &&
	    fabs(gmrmu.phi() - pmrmu.phi()) < 0.001 &&
	    fabs(gmrmu.eta() - pmrmu.eta()) < 0.001 &&
	    fabs(gmrmu.pt()  - pmrmu.pt()) < 0.001) ngpr++;
      }
      else if (bestMu == fmsmu) {
	nfms++; LogTrace("findBestMuons") << "  --> select TPFMS";
      }
      else if (bestMu == pmrmu) {
	npmr++; LogTrace("findBestMuons") << "  --> select PMR";
      }
      else if (!bestMu.isValid()) {
	ntot--; edm::LogWarning("findBestMuons") << "  --> reject muon\n";
      }
      else {
	throw cms::Exception("findBestMuons") << "+++ Unknown outcome! +++\n";
      }
      ntot++;
    }
    if (bestMu.isValid()) {
      tmpV.push_back(bestMu);
    }
  }

  // Fractions of events from different reconstructors.
  if (debug && (eventNum%1000 == 0 && eventNum > 0)) {
    LogTrace("findBestMuons")
      << "\nfindBestMuons summary: "
      << "TO: "    << ntrk/double(ntot) << "; GMR: " << ngmr/double(ntot)
      << "; FMS: " << nfms/double(ntot) << "; PMR: " << npmr/double(ntot);
    LogTrace("findBestMuons")
      << "GMR=PMR: " << ((ngmr > 0) ? ngpr/double(ngmr) : 0);
  }

  return tmpV;
}

zp2mu::Muon Zprime2muAnalysis::NorbertsCocktail(const zp2mu::Muon& trk,
						const zp2mu::Muon& gmr,
						const zp2mu::Muon& fms,
						const zp2mu::Muon& pmr,
						const bool debug) const {
  zp2mu::Muon result;

  double prob0 = 0.0, prob1 = 0.0, prob2 = 0.0, prob3 = 0.0;
  if (trk.isValid()) {
    Genfun::CumulativeChiSquare myCumulativeChiSquare(trk.dof());
    prob0 = -log(1.-myCumulativeChiSquare(trk.chi2()));
  }
  if (gmr.isValid()) {
    Genfun::CumulativeChiSquare myCumulativeChiSquare(gmr.dof());
    prob1 = -log(1.-myCumulativeChiSquare(gmr.chi2()));
  }
  if (fms.isValid()) {
    Genfun::CumulativeChiSquare myCumulativeChiSquare(fms.dof());
    prob2 = -log(1.-myCumulativeChiSquare(fms.chi2()));
  }
  if (pmr.isValid()) {
    Genfun::CumulativeChiSquare myCumulativeChiSquare(pmr.dof());
    prob3 = -log(1.-myCumulativeChiSquare(pmr.chi2()));
  }

  if (debug) {
    LogTrace("NorbertsCocktail")
      << " Event " << eventNum << " Probabilities: trk = " << prob0
      << "; gmr = " << prob1 << "; fms = " << prob2 << "; pmr = " << prob3;
  }

  if (gmr.isValid()) result = gmr;
  if (!gmr.isValid() && pmr.isValid()) result = pmr;
  
  if (gmr.isValid() && pmr.isValid() && ((prob1 - prob3) > 0.05)) result = pmr;

  if (trk.isValid() && !gmr.isValid() && !fms.isValid() && !pmr.isValid()) {
    result = trk;
    return result;
  }
  if (trk.isValid() && fms.isValid() && fabs(prob2 - prob0) > 30.) {
    result = trk;
    return result;
  }

  if (!gmr.isValid() && !pmr.isValid() && fms.isValid()) result = fms;

  zp2mu::Muon tmin;
  double probmin = 0.0;
  if (gmr.isValid() && pmr.isValid()) {
    probmin = prob3; tmin = pmr;
    if (prob1 < prob3) { probmin = prob1; tmin = gmr; }
  }
  else if (!pmr.isValid() && gmr.isValid()) { 
    probmin = prob1; tmin = gmr; 
  }
  else if (!gmr.isValid() && pmr.isValid()) {
    probmin = prob3; tmin = pmr; 
  }
  if (tmin.isValid() && fms.isValid() && ((probmin - prob2) > 3.5)) {
    result = fms;
  }

  return result;
}

// Piotr's version.
zp2mu::Muon Zprime2muAnalysis::PiotrsCocktail(const zp2mu::Muon& trk,
					      const zp2mu::Muon& fms,
					      const zp2mu::Muon& pmr,
					      const bool debug) const {
  zp2mu::Muon result;
  double prob_trk = 1000.0, prob_fms = 1000.0, prob_pmr = 1000.0;

  if (trk.isValid()) {
    Genfun::CumulativeChiSquare myCumulativeChiSquare(trk.dof());
    prob_trk = -log(1.-myCumulativeChiSquare(trk.chi2()));
  }
  if (fms.isValid()) {
    Genfun::CumulativeChiSquare myCumulativeChiSquare(fms.dof());
    prob_fms = -log(1.-myCumulativeChiSquare(fms.chi2()));
  }
  if (pmr.isValid()) {
    Genfun::CumulativeChiSquare myCumulativeChiSquare(pmr.dof());
    prob_pmr = -log(1.-myCumulativeChiSquare(pmr.chi2()));
  }
  if (debug) {
    LogTrace("PiotrsCocktail")
      << " Event " << eventNum << " Probabilities: trk = " << prob_trk
      << "; fms = " << prob_fms << "; pmr = " << prob_pmr;
  }

  if (fms.isValid() && !pmr.isValid()) {
    //if (trk.isValid() && prob_fms - prob_trk > 30.) result = trk;
    //else                                            result = fms;
    result = fms;
  }
  else if (!fms.isValid() && pmr.isValid()) {
    result = pmr;
    //if (trk.isValid() && prob_pmr - prob_trk > 30.) result = trk;
  }
  else if (fms.isValid() && pmr.isValid()) {
    if (prob_fms - prob_pmr > 0.9) result = pmr;
    //else {
    //  if (trk.isValid() && prob_fms - prob_trk > 30.) result = trk;
    //  else                                            result = fms;
    //}
    else result = fms;
  }

  return result;
}


//-----------------------------------------------------------------------------
//          Proximity match for generated and reconstructed tracks
//-----------------------------------------------------------------------------
bool Zprime2muAnalysis::matchEta(const double eta1, const double eta2, 
				 double *eta_diff) {

  bool debug = verbosity >= VERBOSITY_TOOMUCH;

  *eta_diff = abs(eta1-eta2);
  // if (debug)
  // LogTrace("Zprime2muAnalysis") << "Eta1-Eta2: " << *eta_diff;

  // See if the difference in eta is less than or equal to .5
  if (*eta_diff <= .5) {
    if (debug)
      LogTrace("matchEta")
	<< "Eta Match!  Difference is: " << *eta_diff;
    return true;
  }
  else
    return false;
}

bool Zprime2muAnalysis::matchPhi(const double phi1, const double phi2, 
				 double *phi_diff) {
  // See if difference in phi is less than or equal to .5 rad.

  bool debug = verbosity >= VERBOSITY_TOOMUCH;

  *phi_diff = abs(phi1-phi2);
  if (debug)
    LogTrace("matchPhi") << "phi1 : " << phi1 << " phi2 : " << phi2;

  // Compensate for two measurements near 0;
  if (*phi_diff > TMath::Pi())
    *phi_diff = abs((2.*TMath::Pi())-*phi_diff);

  if (*phi_diff <= .5) {
    if (debug)
      LogTrace("matchPhi") << "Phi Match!  Diff: " << *phi_diff;
    return true;
  }
  else
    return false;
}

bool Zprime2muAnalysis::matchEtaAndPhi(const double eta1, const double eta2,
				       const double phi1, const double phi2,
				       double *match_diff) {
  // See if eta and phi match for the 2 tracks.
  bool debug = false;
  double eta_diff = 999., phi_diff = 999.;

  bool isMatch = matchEta(eta1, eta2, &eta_diff) &&
                   matchPhi(phi1, phi2, &phi_diff);

  // Compute this value to give an estimate of how close 2 tracks are
  *match_diff = sqrt((phi_diff*phi_diff)+(eta_diff*eta_diff));

  if (isMatch) {
    if (debug)
      LogTrace("matchEtaAndPhi") << "Match!  Diff: " << *match_diff;
    return true;
  }
  else {
    if (debug)
      LogTrace("matchEtaAndPhi")
	<< "No Match! ------   Diff: " << *match_diff;
    return false;
  }
}

bool Zprime2muAnalysis::matchTracks(const zp2mu::Muon& muon1, 
				    const zp2mu::Muon& muon2,
				    const bool debug) {
  // See if tracks corresponding to mu+ and mu- at different reconstruction
  // recs correspond to each other by making an eta and phi match.
  bool   returnValue = false;
  double match_diff = 999.;

  if (muon1.isValid() == false) {
    edm::LogError("matchTracks") 
      << "+++ muon1 does not exist! +++\n";
    return returnValue;
  }
  if (muon2.isValid() == false) {
    edm::LogError("matchTracks")
      << "+++ muon2 does not exist! +++\n";
    return returnValue;
  }

  if (debug) {
    ostringstream out;
    out << "For muons " << muon1.id() << " (rec = " << muon1.recLevel()
	<< ") and "     << muon2.id() << " (rec = " << muon2.recLevel()
	<< ")" << endl;
    out << "muon1 eta: " << setw(7) << muon1.eta()
	<< "      phi: " << setw(7) << muon1.phi() << endl;
    out << "muon2 eta: " << setw(7) << muon2.eta()
	<< "      phi: " << setw(7) << muon2.phi();
    LogTrace("matchTracks") << out.str();
  }

  if (matchEtaAndPhi(muon1.eta(), muon2.eta(),
		     muon1.phi(), muon2.phi(), &match_diff)) {
    if (debug)
      LogTrace("matchTracks") << "Match is a success!";
    returnValue = true;
  }
  else {
    if (debug)
      LogTrace("matchTracks") << "Match is NOT a success!";
  }
  return returnValue;
}

bool Zprime2muAnalysis::findClosestId(const int rec, const zp2mu::Muon& muon,
				      int *closest_id, const bool debug) {
  // rec:   reconstruction level of muon to be matched to a given muon.
  // muon:  reference to a given muon.
  // returns: true/false and the index of the rec muon which is the closest
  // in eta and phi to a given muon.

  bool returnValue = false;
  double  match_diff = 999., match_mindiff = 999.;
  *closest_id = -999;

  if (rec < lgen || rec >= MAX_LEVELS) {
    edm::LogError("findClosestId")
      << "+++ Unknown rec. level = " << rec << " +++\n";
    return returnValue;
  }
  if (muon.isValid() == false) {
    edm::LogError("findClosestId") 
      << "+++ muon does not exist! +++\n";
    return returnValue;
  }
  if (muon.closestId(rec) >= 0) {
    edm::LogError("findClosestId")
      << "+++ this muon already has a matched muon #"
      << muon.closestId(rec) << " at level " << rec << "! +++\n";
    return returnValue;
  } 
  if (muon.recLevel() == lgen && muon.pt() < 1.) {
    // pT of a generated muon is too low to try to match
    return returnValue;
  } 

  vector<zp2mu::Muon>::const_iterator pmu;
  for (pmu = allMuons[rec].begin(); pmu != allMuons[rec].end(); pmu++) {
    if (rec == lgen && pmu->pt() < 1.) continue;
    if (debug) {
      ostringstream out;
      out << "For muon " << muon.id() << "----------------------------\n";
      out << "rec_a eta: " << muon.eta() << "   phi: " << muon.phi() << endl;
      out << "rec_b eta: " << pmu->eta() << "   phi: " << pmu->phi();
      LogTrace("findClosestId") << out.str();
    }

    // If the tracks lie within a certain region, store the value
    // of the match difference (eta_diff * phi_diff).
    if (matchEtaAndPhi(muon.eta(), pmu->eta(),
		       muon.phi(), pmu->phi(), &match_diff)) {
      if (debug) LogTrace("findClosestId") << "Match is a success!";
      returnValue = true;
      // If the match is closer than before, store the id of a muon.
      if (match_diff < match_mindiff) {
	*closest_id = pmu->id();
	match_mindiff = match_diff;
	if (debug) {
	  LogTrace("findClosestId")
	    << "match_diff: " << match_diff << " for muon #" << *closest_id;
	}
      }
    }
    else {
      if (debug)
	LogTrace("findClosestId") << "Match is NOT a success!";
    }
  }
  return returnValue;
}

// Find the seed of the alternate fit which matches the seed of the
// principle fit.  This function can only compare different high level
// fits (not Gen, L1, or L2).
bool Zprime2muAnalysis::findSameSeedId(const int rec, const zp2mu::Muon& muon, 
				       int *sameseed_id, const bool debug) {
  // rec:   reconstruction level of muon to be matched to a given muon.
  // muon:  reference to a given muon.
  // returns: true/false and the index of the rec muon which has the same
  // seed as a given muon.

  bool returnValue = false;
  *sameseed_id = -999;

  if (rec < l3 || rec >= MAX_LEVELS) {
    edm::LogError("findSameSeedId")
      << "+++ Invalid rec. level = " << rec << " +++\n";
    return returnValue;
  }
  if (muon.isValid() == false) {
    edm::LogError("findSameSeedId") << "+++ muon does not exist! +++\n";
    return returnValue;
  }
  if (muon.recLevel() < l3) {
    edm::LogError("findSameSeedId")
      << "+++ Invalid rec. level = " << muon.recLevel() << " of muon +++\n";
    return returnValue;
  }
  if (muon.sameSeedId(rec) >= 0) {
    edm::LogError("findSameSeedId")
      << "+++ this muon already has a matched muon #"
      << muon.sameSeedId(rec) << " at level " << rec << "! +++\n";
    return returnValue;
  } 

  // Get seed of principle muon.
  int princ_seed = muon.seedIndex();

  // Loop over all seeds of alternative fit and look for matching seed.
  vector<zp2mu::Muon>::const_iterator pmu;
  for (pmu = allMuons[rec].begin(); pmu != allMuons[rec].end(); pmu++) {
    if (pmu->seedIndex() == princ_seed) {
      returnValue = true;
      *sameseed_id = pmu->id();
      if (debug) {
	LogTrace("findSameSeedId")
	  << "Principle muon index = " << muon.id()
	  << ", same seed muon index = " << *sameseed_id
	  << ", their seeds = " << princ_seed;
      }
      return returnValue;
    }
  }
  if (debug)
    LogTrace("findSameSeedId") << "No seed match!";
  return returnValue;
}

void Zprime2muAnalysis::matchAllMuons(const bool debug) {
  // For every muon at every rec. level, look for the same seed and
  // closest muons at other rec. levels, and assign the id's of found
  // muons to the appropriate members of the Muon class.
  // The seed info is not available for generated, L1 and L2 muons, so only
  // look for a closest muon.  For muons at other levels, look for both the
  // same seed and the closest muon.
  int closest_id, sameseed_id;
  vector<zp2mu::Muon>::iterator pmu;

  for (int irec = MAX_LEVELS-1; irec >= lgen; irec--) {
    for (pmu = allMuons[irec].begin(); pmu != allMuons[irec].end(); pmu++) {
      for (int jrec = lgen; jrec < MAX_LEVELS; jrec++) {
	if (jrec != irec) {
	  if (findClosestId(jrec, *pmu, &closest_id, debug)) {
	    pmu->setClosestId(jrec, closest_id);
	  }
	  if (irec > l2 && jrec > l2 &&
	      findSameSeedId(jrec, *pmu, &sameseed_id, debug)) {
	    pmu->setSameSeedId(jrec, sameseed_id);
	  }
	}
	else {
	  pmu->setClosestId(jrec,  pmu->id()); // its own id
	  if (irec > l2 && jrec > l2) pmu->setSameSeedId(jrec, pmu->id());
	}
      }
    }
  }
}

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
//        Construct opposite-sign dileptons out of a list of muons
//-----------------------------------------------------------------------------
void Zprime2muAnalysis::makeAllDileptons(const edm::Event& event,
					 const bool debug) {
  // Try to make dileptons at all levels of reconstruction.
  for (int rec = 0; rec < MAX_LEVELS; rec++) {
    makeDileptons(rec, debug);
  }

  // Try to construct "best" dileptons from the best-reconstructed muons.
  if (!bestMuons.empty()) bestDiMuons = makeDileptons(lbest, bestMuons, debug);

  // Save true MC resonance for comparisons.
  addTrueResonance(event, allDiMuons[lgen]);

  // Include 4-momentum of close-by photon candidates.
  for (int rec = l1; rec <= MAX_LEVELS; rec++) {
    if (rec < MAX_LEVELS)       addBremCandidates(allDiMuons[rec], debug);
    else if (rec == MAX_LEVELS) addBremCandidates(bestDiMuons, debug);
  }

  if (debug) dumpDiMuonMasses();
}

void Zprime2muAnalysis::makeDileptons(const int rec, const bool debug) {
  if (rec < 0 || rec >= MAX_LEVELS) {
    throw cms::Exception("makeDileptons")
      << "+++ Unknown rec. level = " << rec << " +++\n";
  }

  if (searchedDileptons[rec]) {
    edm::LogWarning("makeDileptons")
      << "+++ dileptons for rec. level " << rec
      << " were already looked for +++\n";
    return;
  }

  vector<zp2mu::Muon>& muons = allMuons[rec];
  allDiMuons[rec] = makeDileptons(rec, muons, debug);

  // Reset dilepton id's.
  for (unsigned int i_dil = 0; i_dil < allDiMuons[rec].size(); i_dil++) {
    allDiMuons[rec][i_dil].setId(i_dil);
  }

  // Keep track of the fact that we have done the search already.
  searchedDileptons[rec] = true;
}

vector<zp2mu::DiMuon> 
Zprime2muAnalysis::makeDileptons(const int rec,
				 const vector<zp2mu::Muon>& muons,
				 const bool debug) {
  unsigned int i_dil, n_dil = 0;
  TLorentzVector vmu1, vmu2, vdil;
  vector<zp2mu::DiMuon> diMuons;
  vector<zp2mu::Muon>::const_iterator pmu, qmu;
  vector<zp2mu::DiMuon>::iterator pdi;

  if (debug) {
    LogTrace("makeDileptons")
      << "=======================================================";
    LogTrace("makeDileptons") << "Dilepton search for rec: " << rec;
  }

  // In the special case of the generated dilepton, we will use the actual
  // muons that were produced by the resonance.  In order to do this, we make
  // sure that the muons have a mother ID equal to the ID of the resonance.
  if (rec == 0 && muons.size() > 1) {
   // int idx[MAX_DILEPTONS]; //CL: determine by vertex pos instead of barcode

    math::XYZPoint vertexs[MAX_DILEPTONS];
    for (pmu = muons.begin(); pmu != muons.end(); pmu++) {
      if (pmu->genMotherId() == 32) {
	bool isnewdil = true;
	int  idxdil = -999;
	for (i_dil = 0; i_dil < n_dil; i_dil++) {
//	  if (pmu->genMotherLine() == idx[i_dil]) { //FIXME
	  if ( haveSameVertex((*pmu), vertexs[i_dil])) { 
	    idxdil = i_dil;
	    isnewdil = false;
	    break;
	  }
	}
	// Encounter this mother Id for the first time
	if (isnewdil && n_dil < MAX_DILEPTONS) {
	  zp2mu::DiMuon thisDiMu;
	  thisDiMu.fill(true, n_dil, 0);
	  diMuons.push_back(thisDiMu);
	  idxdil = n_dil;
//	  idx[idxdil] = pmu->genMotherLine();
	  vertexs[idxdil] = pmu->vertex();
	  n_dil++;
	}

	if (idxdil > -999) {
	  if (pmu->charge() == 1 &&
	      diMuons[idxdil].muPlus().isValid() == false) {
	    diMuons[idxdil].setMuPlus(*pmu);
	  }
	  else if (pmu->charge() == -1 &&
		   diMuons[idxdil].muMinus().isValid() == false) {
	    diMuons[idxdil].setMuMinus(*pmu);
	  }
	  else {
	    throw cms::Exception("makeDileptons")
	      << "+++" << "index_pos[0][" << idxdil << "] = "
	      << diMuons[idxdil].muPlus().isValid()
	      << "index_neg[0][" << idxdil << "] = "
	      << diMuons[idxdil].muMinus().isValid()
	      << "  charge = " << pmu->charge() << " +++\n";
	  }
	}
      }
    }

    for (pdi = diMuons.begin(); pdi != diMuons.end();) {
      if (pdi->muPlus().isValid()  == false ||
	  pdi->muMinus().isValid() == false) {
	edm::LogWarning("makeDileptons")
	  << "+++ dilepton # " << pdi->id()
	  << " at level " << rec << " is made of "
	  << (pdi->muPlus().isValid()  ? "valid" : "invalid") << " mu+ and " 
	  << (pdi->muMinus().isValid() ? "valid" : "invalid") << " mu-;"
	  << " removing it +++\n";
	if (debug) {
	  dumpEvent(true, false, false, false, false);
	}
	pdi = diMuons.erase(pdi); // points to the next element
      }
      else {
	zp2mu::Muon mup = pdi->muPlus();
	zp2mu::Muon mum = pdi->muMinus();
	vmu1.SetPtEtaPhiM(mup.pt(), mup.eta(), mup.phi(), MUMASS);
	vmu2.SetPtEtaPhiM(mum.pt(), mum.eta(), mum.phi(), MUMASS);
	vdil = vmu1 + vmu2;
	pdi->setDimuV(vdil);
	pdi++;
      }
    }
  }

  // Dileptons search in reconstructed muons.  First consider all possible
  // combinations of mu+ and mu-.
  if (rec > 0 && muons.size() > 1) {
    for (pmu = muons.begin(); pmu != muons.end()-1; pmu++) {

      vmu1.SetPtEtaPhiM(pmu->pt(), pmu->eta(), pmu->phi(), MUMASS);

      for (qmu = pmu+1; qmu != muons.end(); qmu++) {
	if (qmu->charge() != pmu->charge()) {
	  vmu2.SetPtEtaPhiM(qmu->pt(), qmu->eta(), qmu->phi(), MUMASS);
	  vdil = vmu1 + vmu2;

	  zp2mu::DiMuon thisDiMu;
	  thisDiMu.fill(true, n_dil, rec);
	  if (pmu->charge() == 1) {
	    thisDiMu.setMuPlus(*pmu);
	    thisDiMu.setMuMinus(*qmu);
	  }
	  else {
	    thisDiMu.setMuMinus(*pmu);
	    thisDiMu.setMuPlus(*qmu);
	  }
	  thisDiMu.setDimuV(vdil);
	  diMuons.push_back(thisDiMu);
	  n_dil++;
	}
      }
    }
  }

  // If H -> ZZ* -> 4mu, sort dileptons to look for two highest-mass ones.
  // If Z' or G*, do nothing because all muons were sorted from the highest
  // P to the lowest, and therefore the very first dimuon is already made
  // of two highest P muons.
  if (doingHiggs) {
    // Sort the list of found dileptons using overloaded criteria.
    if (diMuons.size() > 1) {
      sort(diMuons.begin(), diMuons.end(), greater<zp2mu::DiMuon>());
    }
  }

  if (rec > 0 && diMuons.size() > 1) {
    // If more than one dilepton was formed, we accept only those containing
    // distinct muons.  Again, the preference is given to higher mass
    // dilepton over lower mass dilepton.
    for (pdi = diMuons.begin(); pdi != diMuons.end()-1;) {
      vector<zp2mu::DiMuon>::iterator qdi;
      for (qdi = pdi+1; qdi != diMuons.end(); qdi++) {
	if (qdi->muMinus().id() == pdi->muMinus().id() ||
	    qdi->muPlus().id()  == pdi->muPlus().id()) {
	  //if ((*qdi) == (*pdi)) { // use overloaded ==
	  diMuons.erase(qdi);
	  pdi = diMuons.begin(); // reset pointers and restart
	  break;
	}
	else {
	  pdi++;
	}
      }
    }

    // Store only 1 dilepton for Z' and G*, and two highest mass dileptons
    // for H -> ZZ* -> 4 mu.
    if (doingHiggs) {
      if (diMuons.size() > MAX_DILEPTONS) {
	diMuons.erase(diMuons.begin() + MAX_DILEPTONS, diMuons.end());
      }
    }
    else {
      if (diMuons.size() > 1) {
	diMuons.erase(diMuons.begin() + 1, diMuons.end());
      }
    }
  }

  if (debug) {
    for (pdi = diMuons.begin(); pdi != diMuons.end(); pdi++) {
      LogTrace("makeDileptons")
	<< "Dilepton #" << pdi->id()
	<< " is found for muon+ #" << pdi->muPlus().id()
	<< " and muon- #"          << pdi->muMinus().id()
	<< "; its mass: "          << setprecision(5) << pdi->dimuV().M();
    }
  }

  return diMuons;
}

void Zprime2muAnalysis::addTrueResonance(const edm::Event& event,
					 vector<zp2mu::DiMuon>& diMuons) {
  bool debug = verbosity >= VERBOSITY_TOOMUCH;

  // Generated particles (i.e. PYTHIA).
  edm::Handle<reco::CandidateCollection> genParticles;
  event.getByLabel("genParticleCandidates", genParticles);

  for (vector<zp2mu::DiMuon>::iterator pdi = diMuons.begin();
       pdi != diMuons.end(); pdi++) {

    // Check that the dimuon is indeed the generated one.
    int rec = pdi->recLevel();
    if (rec != lgen)
      throw cms::Exception("addTrueResonance")
	<< "+++ Invalid rec. level = " << rec << " +++\n";

    // Take true resonance from PYTHIA only if generated mu+ and mu-
    // originated from the same resonance.
    // SV: same approach as in ORCA, but perhaps needs some thinking.
    if (pdi->muPlus().genMotherLine() == pdi->muMinus().genMotherLine()) {
      //if ( haveSameVertex(pdi->muPlus(), pdi->muMinus()) ) {

      int idx = -1;
      for (reco::CandidateCollection::const_iterator ipl = genParticles->begin();
	   ipl != genParticles->end(); ipl++) {
	idx++;

	if (idx == pdi->muPlus().genMotherLine()) {
	//if ( isMotherOf((*(dynamic_cast<const reco::GenParticleCandidate*>(&*ipl))), pdi->muPlus(), pdi->muMinus()) ) { //FIXME
	  pdi->setResV(TLorentzVector((ipl)->p4().x(),
				      (ipl)->p4().y(),
				      (ipl)->p4().z(),
				      (ipl)->p4().t()));
	  if (debug) LogTrace("addTrueResonance")
	    << "Resonance: id = " << (ipl)->pdgId() << " line = " << idx
	    << " M = " << (ipl)->mass();
	  break;
	}
      }
    }
    else {
      // Just use dimuon mass.
      pdi->setResV(pdi->dimuV());
    }

    if (debug) {
      LogTrace("addTrueResonance")
	<< "MC resonance: dimuon inv. mass = " << pdi->dimuV().M()
	<< "; true inv. mass = " << pdi->resV().M();
    }
  }
}

// If there is a photon candidate at the phi-eta distance dRmax from 
// either of the muons, combine its 4-momentum with that of the dimuon
// and put the result into dimuon's TLorentzVector theResV.
void Zprime2muAnalysis::addBremCandidates(vector<zp2mu::DiMuon>& diMuons,
					  const bool debug) {
  const double dRmax = 0.1;
  TLorentzVector vmu1, vmu2;

  vector<zp2mu::DiMuon>::iterator pdi;
  for (pdi = diMuons.begin(); pdi != diMuons.end(); pdi++) {
    int rec = pdi->recLevel();

    if (rec == lgen) {
      // Do nothing; work done in addTrueResonance()
      return;
    }
    else if (rec == l1 || rec == l2) {
      // No way to improve measurement at Level-1 and Level-2
      pdi->setResV(pdi->dimuV());
    }
    else {
      zp2mu::Muon mup = pdi->muPlus();
      zp2mu::Muon mum = pdi->muMinus();
      vmu1.SetPtEtaPhiM(mup.pt(), mup.eta(), mup.phi(), MUMASS);
      vmu2.SetPtEtaPhiM(mum.pt(), mum.eta(), mum.phi(), MUMASS);

      TLorentzVector vph1 = mup.closestPhoton();
      double dr1 = 999.;
      if (vph1.E() > 0.) dr1 = vmu1.DeltaR(vph1);

      TLorentzVector vph2 = mum.closestPhoton();
      double dr2 = 999.;
      if (vph2.E() > 0.) dr2 = vmu2.DeltaR(vph2);

      bool add_ph1 = false, add_ph2 = false;
      if (dr1 < dRmax) {add_ph1 = true;}

      if (dr2 < dRmax &&
	  (add_ph1 == false ||
	   (fabs(vph2.Px()-vph1.Px()) > 0.001 &&
	    fabs(vph2.Py()-vph1.Py()) > 0.001 &&
	    fabs(vph2.Pz()-vph1.Pz()) > 0.001))) {add_ph2 = true;}

      TLorentzVector vres = pdi->dimuV();
      if (add_ph1) {vres += vph1;}
      if (add_ph2) {vres += vph2;}
      pdi->setResV(vres);
    }

    if (debug) {
      LogTrace("addBremCandidates")
	<< "Include brem candidate(s) at level " << rec
	<< ": inv. mass w/o photons = " << pdi->dimuV().M()
	<< "; w/  = " << pdi->resV().M();
    }
  }
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

//-----------------------------------------------------------------------------
//                           Print-out methods
//-----------------------------------------------------------------------------
void Zprime2muAnalysis::dumpEvent(const bool printGen, const bool printL1,
				  const bool printL2,  const bool printL3,
				  const bool printBest) const {
  vector<zp2mu::Muon>::const_iterator pmu;

  ostringstream out;
  out << "\n******************************** Event " << eventNum << "\n";
  if (printGen)
    for (pmu = allMuons[0].begin(); pmu != allMuons[0].end(); pmu++) {
      out << *pmu;
    }
  if (printL1)
    for (pmu = allMuons[l1].begin(); pmu != allMuons[1].end(); pmu++) {
      out << *pmu;
    }
  if (printL2)
    for (pmu = allMuons[l2].begin(); pmu != allMuons[2].end(); pmu++) {
      out << *pmu;
    }
  if (printL3) {
    for (int rec = l3; rec < MAX_LEVELS; rec++)
      for (pmu = allMuons[rec].begin(); pmu != allMuons[rec].end(); pmu++) {
	out << *pmu;
      }
  }
  if (printBest) {
    out << "\nBest off-line muons: \n";
    for (pmu = bestMuons.begin(); pmu != bestMuons.end(); pmu++) {
      out << *pmu;
    }
  }
  LogTrace("dumpEvent") << out.str();
}

void Zprime2muAnalysis::dumpDiMuonMasses() const {
  ostringstream out;
  out << "Dilepton masses at levels  ";
  for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
    out << " " << i_rec;
    if (i_rec < MAX_LEVELS-1) out << "      ";
  }

  out << setw(6) << setprecision(5);
  for (unsigned int i_dil = 0; i_dil < MAX_DILEPTONS; i_dil++) {
    out << "\n Dilepton # " << i_dil << "; dimu mass: ";
    for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
      if (allDiMuons[i_rec].size() > i_dil)
	out << allDiMuons[i_rec][i_dil].dimuV().M();
      else
	out << " ---  ";
      if (i_rec < MAX_LEVELS-1) out << ", ";
    }
    out << "\n";

    out << "                res mass: ";
    for (int i_rec = 0; i_rec < MAX_LEVELS; i_rec++) {
      if (allDiMuons[i_rec].size() > i_dil)
	out << allDiMuons[i_rec][i_dil].resV().M();
      else
	out << " ---  ";
      if (i_rec < MAX_LEVELS-1) out << ", ";
    }
  }
  LogTrace("dumpDiMuonMasses") << out.str();
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

zp2mu::Muon& Zprime2muAnalysis::muonRef(const int rec, const int id) {
  // Access method which returns reference to the muon with muonId=id.

  if (rec < 0 || rec >= MAX_LEVELS) {
    throw cms::Exception("muonRef")
      << "+++ Unknown rec. level " << rec << " +++\n";
  }

  vector<zp2mu::Muon>::iterator pmu;
  for (pmu = allMuons[rec].begin(); pmu != allMuons[rec].end(); pmu++)
    if (pmu->id() == id)
      return (*pmu);

  // if execution gets to here, we haven't found the muon; throw
  throw cms::Exception("muonRef") << "+++ muonRef not found! +++\n";
}

zp2mu::DiMuon& Zprime2muAnalysis::dimuonRef(const int rec, const int id) {
  // Access method which returns reference to the dimuon with dimuonId=id.

  if (rec < 0 || rec >= MAX_LEVELS) {
    throw cms::Exception("dimuonRef")
      << "+++ Unknown rec. level " << rec << " +++\n";
  }

  vector<zp2mu::DiMuon>::iterator pdi;
  for (pdi = allDiMuons[rec].begin(); pdi != allDiMuons[rec].end(); pdi++)
    if (pdi->id() == id)
      return (*pdi);

  // if execution gets to here, we haven't found the dimuon; throw
  throw cms::Exception("dimuonRef") << "+++ dimuonRef not found! +++\n";
}

//-----------------------------------------------------------------------------
//                                 Trigger
//-----------------------------------------------------------------------------
bool Zprime2muAnalysis::passTrigger(const int irec) {

  if (calcTrig[irec])
    return passTrig[irec];

  // Keep track of the fact that you have checked the trigger already
  calcTrig[irec] = true;

  // JMTBAD all events pass trigger for now...
  trigWord[irec] = 1; 
  passTrig[irec] = true;
  return true;

  /*
  unsigned int trigbits = 0;
  switch(irec) {
  case 0:
    // Generated events always pass.
    trigbits = 1;
    break;
  case 1:
    trigbits = PassL1Trigger();
    break;
  case 2:
    trigbits = PassL2Trigger();
    break;
  case 3:
    trigbits = PassL3Trigger();
    break;
  default:
    edm::LogError("Zprime2muAnalysis")
      << "+++ passTrigger error: L" << irec << " trigger is unknown +++\n";
    break;
  }
  trigWord[irec] = trigbits;
  passTrig[irec] = (trigbits != 0);
  return passTrig[irec];
  */
}

bool Zprime2muAnalysis::passTrigger() {
  unsigned int decision = 1;

  for (int itrig = l1; itrig <= l3; itrig++) {
    if (!calcTrig[itrig]) {
      decision *= passTrigger(itrig);
    }
    else {
      decision *= trigWord[itrig];
    }
  }
  return (decision != 0);
}

void Zprime2muAnalysis::analyze(const edm::Event& event,
				const edm::EventSetup& eSetup) {
  // could store the whole event/eventsetup object if needed
  eventNum = event.id().event();
  storeMuons(event);
}

void Zprime2muAnalysis::beginJob(const edm::EventSetup& eSetup) {
}

void Zprime2muAnalysis::endJob() {
}

ostream& operator<<(ostream& out, const TLorentzVector& vect) {
  out << setw(7) << vect.Eta() << " | " << setw(7) << vect.Phi() << " | "
      << setw(7) << vect.P()   << " | " << setw(7) << vect.Pt()  << " | " 
      << setw(7) << vect.Pz()  << " | " << setw(7) << vect.Rapidity() << " | "
      << setw(6) << vect.M();
  return out;
}
