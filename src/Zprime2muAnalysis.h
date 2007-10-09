#ifndef ZP2MUANALYSIS_H
#define ZP2MUANALYSIS_H

#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Muon.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DiMuon.h"

// details about the number of reconstruction levels stored
const int NUM_REC_LEVELS = 4;
const int MAX_LEVELS = zp2mu::REC_LEVELS;
enum RECLEVEL { lgen, l1, l2, l3, lgmr, ltk, lfms, lpmr, lbest = 99 };
const std::string str_level[] = {
  "Gen", " L1", " L2", " L3", "GMR", "Tracker-only", "TPFMS", "PMR", "TMR"
};

// details about the number of quality cuts available
const int NUM_Q_SETS       = 8;
const int NUM_L2_CUTS      = 4;
const int NUM_L3_CUTS      = 9;
const int NUM_TRACKER_CUTS = 4;

namespace reco {class GenParticleCandidate;}

enum VERBOSITY { VERBOSITY_NONE, VERBOSITY_SIMPLE,
		 VERBOSITY_LOTS, VERBOSITY_TOOMUCH };

class Zprime2muAnalysis : public edm::EDAnalyzer {
 public:
  explicit Zprime2muAnalysis(const edm::ParameterSet&);
  virtual ~Zprime2muAnalysis();

  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  // keep track of the event number for printing out
  int eventNum;

  void InitROOT();

  double MUMASS;
  static const double PTMIN;
  static const unsigned int MAX_DILEPTONS;
  static const double ETA_CUT;
  static const double ENDCAP_BARREL_CUT;
  static const double TRIGGER_ETA_CUT[NUM_REC_LEVELS];
  static const bool DO_QCUTS;
  static const int QSEL;
  static const double L2QCUT[NUM_Q_SETS][NUM_L2_CUTS];
  static const double L3QCUT[NUM_Q_SETS][NUM_L3_CUTS];
  static const double TRACKERQCUT[NUM_Q_SETS][NUM_TRACKER_CUTS];
  static bool cutTrig[NUM_REC_LEVELS];

  // config file parameters
  bool doingHiggs; // determines whether one or two dimuons are kept
  bool generatedOnly; // whether only to look at generated muons
  bool doingElectrons; //determines whether to run on muons or electrons

  // Track quality studies and optimization.
  bool TrackQCheck(const zp2mu::Muon& muon, const int qsel,
		   int& ncut, const bool debug = false);
  bool dilQCheck(const zp2mu::DiMuon& dimuon, const int qsel,
		 int& ncut_mum, int& ncut_mup,
		 const bool debug = false);

  // General
  zp2mu::Muon& muonRef(const int rec, const int id);
  zp2mu::DiMuon& dimuonRef(const int rec, const int id);

  // Print-outs
  void dumpEvent(const bool printGen, const bool printL1 = false, 
		 const bool printL2 = false, const bool printL3 = false,
		 const bool printBest = false) const;
  void dumpDiMuonMasses() const;
  void dumpDilQuality();

  // Methods for matching muon tracks at different rec levels
  bool matchEta(const double eta1, const double eta2,
		double *eta_diff);
  bool matchPhi(const double phi1, const double phi2,
		double *phi_diff);
  bool matchEtaAndPhi(const double eta1, const double eta2,
		      const double phi1, const double phi2,
		      double *match_diff);
  bool matchTracks(const zp2mu::Muon& muon1, const zp2mu::Muon& muon2,
		   const bool debug = false);
  bool findClosestId(const int rec,  const zp2mu::Muon& muon, int *closest_id,
		     const bool debug = false);
  bool findSameSeedId(const int rec, const zp2mu::Muon& muon, int *sameseed_id,
		      const bool debug = false);
  void matchAllMuons(const bool debug = false);
  void matchStudy(const zp2mu::Muon& muon);
  int findMatchedDiMuonId(const int rec, const zp2mu::DiMuon& dimuon,
			  const bool debug = false);

  bool passTrigger(const int irec);
  bool passTrigger();

 private:
  VERBOSITY verbosity;

  void clearValues();
  void storeGeneratedMuons(const edm::Event&);
  void storeL2Muons(const edm::Event& event);
  void storeOfflineMuons(const edm::Event&, const edm::InputTag& whichMuons,
			 RECLEVEL irec, bool trackerOnly=false);
  bool storeOfflineMuon(const int imu, const RECLEVEL irec,
			const reco::TrackRef& theTrack,
			const reco::TrackRef& tkTrack,
			const reco::TrackRef& muTrack,
			const int seedIndex,
			const reco::PhotonCollection& photonCollection);
  bool storePixelMatchGsfElectron(const int imu, const RECLEVEL irec,
				  const reco::PixelMatchGsfElectron& theElectron,
				  const reco::PhotonCollection& photonCollection);
  void storePixelMatchGsfElectrons(const edm::Event&, const edm::InputTag& whichMuons,
				   RECLEVEL irec, bool trackerOnly=false);
  void storeMuons(const edm::Event&);

  std::vector<zp2mu::Muon> findBestMuons();
  zp2mu::Muon NorbertsCocktail(const zp2mu::Muon& trk, const zp2mu::Muon& gmr,
			       const zp2mu::Muon& fms, const zp2mu::Muon& pmr,
			       const bool debug) const;
  zp2mu::Muon PiotrsCocktail(const zp2mu::Muon& trk, const zp2mu::Muon& fms,
			     const zp2mu::Muon& pmr, const bool debug) const;

  TLorentzVector findClosestPhoton(const reco::TrackRef& muonTrack,
			       const reco::PhotonCollection& photonCollection);
  TLorentzVector findClosestPhoton(const reco::GsfTrackRef& muonTrack,
				   const reco::PhotonCollection& photonCollection);
  double deltaR(const double eta1, const double phi1,
		const double eta2, const double phi2) const;

  // utility functions needed to calculate error on 1/Pt, 1/P since such
  // methods do not exist in reco::Track as of now
  double invPtError(const reco::TrackRef& track);
  double invPtError(const reco::GsfTrackRef& track);
  double invPError(const reco::TrackRef& track);
  double invPError(const reco::GsfTrackRef& track);

  // utility function to match standalone muons (to match seeds)
  int matchStandAloneMuon(const edm::Handle<reco::TrackCollection> staTracks,
			  const reco::TrackRef& track);

  // utility function to match muons from same vertex
  bool haveSameVertex(const zp2mu::Muon& , const zp2mu::Muon&) const;

  // utility function to match muons from same vertex
  bool haveSameVertex(const zp2mu::Muon& , const math::XYZPoint&) const;

  bool isMotherOf(const reco::GenParticleCandidate&, const zp2mu::Muon&, const zp2mu::Muon&) const;

  edm::InputTag standAloneMuons;
  edm::InputTag genMuons;
  edm::InputTag globalMuons;
  edm::InputTag globalMuonsFMS;
  edm::InputTag globalMuonsPMR;
  edm::InputTag photons;
  edm::InputTag pixelMatchGsfElectrons;

  // Methods to build up dileptons from stored muons
  void makeAllDileptons(const edm::Event& event, const bool debug = false);
  void makeDileptons(const int rec, const bool debug = false);
  std::vector<zp2mu::DiMuon> makeDileptons(const int rec, 
				     const std::vector<zp2mu::Muon>& muons,
				     const bool debug = false);
  void addTrueResonance(const edm::Event& event, 
			std::vector<zp2mu::DiMuon>& diMuons);
  void addBremCandidates(std::vector<zp2mu::DiMuon>& diMuons, 
			 const bool debug = false);

 protected:  // the muons need to be accessible from our derived classes
  std::vector<zp2mu::Muon> allMuons[MAX_LEVELS];
  std::vector<zp2mu::Muon> bestMuons;

  bool searchedDileptons[MAX_LEVELS];
  std::vector<zp2mu::DiMuon> allDiMuons[MAX_LEVELS];
  std::vector<zp2mu::DiMuon> bestDiMuons;

  bool passTrig[NUM_REC_LEVELS];
  bool calcTrig[NUM_REC_LEVELS];
  unsigned int trigWord[NUM_REC_LEVELS];

  unsigned int leptonFlavor;

};

ostream& operator<<(ostream& out, const TLorentzVector& vect);

#endif // ZP2MUANALYSIS_H
