#include <cstdarg>
#include <ostream>
#include <boost/foreach.hpp>

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Dumpers.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"

int mlprintf(const char* category, const char* fmt, ...) {
  static const size_t bufsize = 10240; // big enough?
  static char buf[bufsize];
  va_list args;
  va_start(args, fmt);
  int ret = vsnprintf(buf, bufsize, fmt, args);
  va_end(args);
  edm::LogInfo(category) << buf;
  return ret;
}

int osprintf(std::ostream& out, const char* fmt, ...) {
  static const size_t bufsize = 10240; // big enough?
  static char buf[bufsize];
  va_list args;
  va_start(args, fmt);
  int ret = vsnprintf(buf, bufsize, fmt, args);
  va_end(args);
  out << buf;
  return ret;
}

std::ostream& operator<<(std::ostream& out, const reco::GenParticle& gen) {
  out << "pdgId: " << gen.pdgId() << " status: " << gen.status() << " q: " << gen.charge() << " pt: " << gen.pt()
      << " eta: " << gen.eta() << " phi: " << gen.phi() << " mass: " << gen.mass() << " vertex: " << gen.vertex();
  return out;
}

std::ostream& operator<<(std::ostream& out, const reco::HitPattern& hp) {
  out << "# hits: " << hp.numberOfHits()
      << "\n  # valid: tot: " << hp.numberOfValidHits() << " tk: " << hp.numberOfValidTrackerHits() << " pxb: " << hp.numberOfValidPixelBarrelHits() << " pxe: " << hp.numberOfValidPixelEndcapHits() << " tib: " << hp.numberOfValidStripTIBHits() << " tob: " << hp.numberOfValidStripTOBHits() << " tid: " << hp.numberOfValidStripTIDHits() << " tec: " << hp.numberOfValidStripTECHits() << " mu: " << hp.numberOfValidMuonHits() << " csc: " << hp.numberOfValidMuonCSCHits() << " dt: " << hp.numberOfValidMuonDTHits() << " rpc: " << hp.numberOfValidMuonRPCHits()
      << "\n  # lost: tot: " << hp.numberOfLostHits() << " tk: " << hp.numberOfLostTrackerHits() << " pxb: " << hp.numberOfLostPixelBarrelHits() << " pxe: " << hp.numberOfLostPixelEndcapHits() << " tib: " << hp.numberOfLostStripTIBHits() << " tob: " << hp.numberOfLostStripTOBHits() << " tid: " << hp.numberOfLostStripTIDHits() << " tec: " << hp.numberOfLostStripTECHits() << " mu: " << hp.numberOfLostMuonHits() << " csc: " << hp.numberOfLostMuonCSCHits() << " dt: " << hp.numberOfLostMuonDTHits() << " rpc: " << hp.numberOfLostMuonRPCHits()
      << "\n  # bad: tot: " << hp.numberOfBadHits() << " mu: " << hp.numberOfBadMuonHits()  << " csc: " << hp.numberOfBadMuonCSCHits()  << " dt: " << hp.numberOfBadMuonDTHits()  << " rpc: " << hp.numberOfBadMuonRPCHits()
      << "\n  # tk layers: with meas: " << hp.trackerLayersWithMeasurement() << " without: " << hp.trackerLayersWithoutMeasurement() << " totallyofforbad: " << hp.trackerLayersTotallyOffOrBad() << " null: " << hp.trackerLayersNull()
      << "\n  # px layers: with meas: " << hp.pixelLayersWithMeasurement() << " without: " << hp.pixelLayersWithoutMeasurement() << " totallyofforbad: " << hp.pixelLayersTotallyOffOrBad() << " null: " << hp.pixelLayersNull()
      << "\n  # si layers: with meas: " << hp.stripLayersWithMeasurement() << " without: " << hp.stripLayersWithoutMeasurement() << " totallyofforbad: " << hp.stripLayersTotallyOffOrBad() << " null: " << hp.stripLayersNull()
      << "\n  # dt with both views: " << hp.numberOfDTStationsWithBothViews() << " with rphi: " << hp.numberOfDTStationsWithRPhiView() << " with rz: " << hp.numberOfDTStationsWithRZView();
  return out;
}

std::ostream& operator<<(std::ostream& out, const reco::Track& tk) {
  out << "algo: " << tk.algoName() << " qualityMask: " << tk.qualityMask()
      << " q: " << tk.charge() << " p: " << tk.p()
      << " q/p error: " << tk.qoverpError() << " theta: " << tk.theta() << " theta error: " << tk.thetaError() << " phi error: " << tk.phiError()
      << " pt: " << tk.pt() << " pt error: " << tk.ptError() << " eta: " << tk.eta()
      << " phi: " << tk.phi() << " chi2: " << tk.chi2() << " dof: " << tk.ndof()
      << "\n  d0: " << tk.d0() << " d0 error: " << tk.d0Error()
      << " reference point: " << tk.referencePoint()
//    << "\n  innerPosition: " << tk.innerPosition() << " outerPosition: " << tk.outerPosition()
//    << "\n  innerMomentum: " << tk.innerMomentum() << " outerMomentum: " << tk.outerMomentum()
      << "\nhitpattern: " << tk.hitPattern();
  return out;
}

std::ostream& operator<<(std::ostream& out, const pat::Muon& mu) {
  out << "pt: " << mu.pt() << " eta: " << mu.eta() << " phi: " << mu.phi() << " p: " << mu.p() << "\nisGlobal: " << mu.isGlobalMuon() << " isTracker: " << mu.isTrackerMuon() << " isStandAlone: " << mu.isStandAloneMuon();
  
  if (mu.genParticle())
    out << "\nMC match: " << *mu.genParticle();

  if (mu.hasUserInt("trackUsedForMomentum")) 
    out << "\nTrack used for momentum: " << mu.userInt("trackUsedForMomentum");
  
  out << "\nTeV refit values:";
  osprintf(out, "\n%20s%20s%20s%20s%20s", "refit", "pt", "eta", "phi", "chi2/dof");
  for (size_t i = 0; i < patmuon::nTrackTypes; ++i) {
    reco::TrackRef tk = patmuon::trackByType(mu, patmuon::TrackType(i));
    osprintf(out, "\n%20s%20.1f%20.1f%20.1f%20.1f", patmuon::track_names[i].c_str(), tk->pt(), tk->eta(), tk->phi(), tk->normalizedChi2());
  }
    
  if (mu.isTrackerMuon())
    out << "\nTM number of matches (type=SegmentAndTrackArbitration): " << mu.numberOfMatches();

  if (mu.innerTrack().isNull())
    out << "\nTracker track ref is null!\n";
  else
    out << "\nTracker track:\n" << *mu.innerTrack();
  
  if (mu.outerTrack().isNull())
    out << "\nStand-alone track ref is null!\n";
  else
    out << "\nStand-alone track:\n" << *mu.outerTrack();
  
  if (mu.globalTrack().isNull())
    out << "\nGlobal track ref is null!\n";
  else
    out << "\nGlobal track:\n" << *mu.outerTrack();
  
  out << "\nTrigger match info: # trig obj matches: " << mu.triggerObjectMatches().size();
  int itosa = 0;
  BOOST_FOREACH(const pat::TriggerObjectStandAlone& tosa, mu.triggerObjectMatches()) {
    out << "\nTriggerObjectStandAlone #" << itosa++ << "\npaths:";
    BOOST_FOREACH(const std::string& s, tosa.pathNames())
      out << "\n" << s;
    out << "\nfilters:";
    BOOST_FOREACH(const std::string& s, tosa.filterLabels())
      out << "\n" << s;
  }

  static const char* id_algos[] = {
    "TrackerMuonArbitrated",
    "AllArbitrated",
    "GlobalMuonPromptTight",
    "TMLastStationLoose",
    "TMLastStationTight",
    "TM2DCompatibilityLoose",
    "TM2DCompatibilityTight",
    "TMOneStationLoose",
    "TMOneStationTight",
    "TMLastStationOptimizedLowPtLoose",
    "TMLastStationOptimizedLowPtTight",
    "GMTkChiCompatibility",
    "GMStaChiCompatibility",
    "GMTkKinkTight",
    "TMLastStationAngLoose",
    "TMLastStationAngTight",
    "TMOneStationAngLoose",
    "TMOneStationAngTight",
    "TMLastStationOptimizedBarrelLowPtLoose",
    "TMLastStationOptimizedBarrelLowPtTight",
    0
  };

  out << "\nMuon id algorithm results:";
  for (size_t i = 0; id_algos[i] != 0; ++i)
    out << "\n  " << std::setw(40) << id_algos[i] << ": " << muon::isGoodMuon(mu, muon::selectionTypeFromString(id_algos[i]));

  if (mu.isTimeValid()) {
    const reco::MuonTime& mt = mu.time();
    out << "\nTiming info: direction: " << mt.direction() << " nDof: " << mt.nDof << " timeAtIpInOut: " << mt.timeAtIpInOut << " +/- " << mt.timeAtIpInOutErr << " timeAtIpOutIn: " << mt.timeAtIpOutIn << " +/- " << mt.timeAtIpOutInErr;
  }
  else
    out << "\nMuonTime structure unavailable!";

  if (mu.isQualityValid()) {
    const reco::MuonQuality& mq = mu.combinedQuality();
    out << "\nQuality info: glbKink: " << mq.glbKink << " pos: " << mq.glbKink_position << " trkKink: " << mq.trkKink << " pos: " << mq.tkKink_position << " staRelChi2: " << mq.staRelChi2 << " trkRelChi2: " << mq.trkRelChi2;
  }
  else
    out << "\nMuonQuality structure unavailable!";

  return out;
}

std::ostream& operator<<(std::ostream& out, const pat::Electron& el) {
  out << "et: " << el.et() << " eta: " << el.eta() << " phi: " << el.phi() << " energy: " << el.energy();
  if (el.genParticle())
    out << "\nMC match: " << *el.genParticle();
  return out;
}
  
std::ostream& operator<<(std::ostream& out, const reco::CandidateBaseRef& cbr) {
  const pat::Muon* mu = dynamic_cast<const pat::Muon*>(&*cbr);
  if (mu)
    out << "pat::Muon: " << *mu;
  else {
    const pat::Electron* el = dynamic_cast<const pat::Electron*>(&*cbr);
    if (el)
      out << "pat::Electron: " << *el;
    else
      out << "dunno!";
  }

  return out;
}

std::ostream& operator<<(std::ostream& out, const pat::CompositeCandidate& dil) {
  out << "mass: " << dil.mass() << " pt: " << dil.pt() << " rapidity: " << dil.rapidity() << " number of daughters: " << dil.numberOfDaughters();
  for (size_t i = 0; i < dil.numberOfDaughters(); ++i)
    out << "\ndaughter " << i << ": pdgId: " << dil.daughter(i)->pdgId()
	<< "\n" << dil.daughter(i)->masterClone();
  
  return out;
}

std::ostream& operator<<(std::ostream& out, const edm::TriggerNames& t) {
  const size_t n = t.size();
  for (size_t i = 0; i < n; ++i) {
    out << i << "\t" << t.triggerName(i);
    if (i < n-1)
      out << "\n";
  }
  return out;
}
