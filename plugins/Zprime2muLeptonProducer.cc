#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"

class Zprime2muLeptonProducer : public edm::EDProducer {
public:
  explicit Zprime2muLeptonProducer(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  pat::Muon* cloneAndSwitchMuonTrack(const pat::Muon&) const;

  std::pair<pat::Electron*, int> doLepton(const pat::Electron&, const reco::CandidateBaseRef&) const;
  std::pair<pat::Muon*,     int> doLepton(const pat::Muon&,     const reco::CandidateBaseRef&) const;

  template <typename T> void doLeptons(edm::Event&, const edm::InputTag&, const std::string&) const;

  edm::InputTag muon_src;
  edm::InputTag electron_src;
  StringCutObjectSelector<pat::Muon> muon_selector;
  StringCutObjectSelector<pat::Electron> electron_selector;
  std::string muon_track_for_momentum;
  std::vector<std::string> muon_tracks_for_momentum;
  edm::InputTag muon_photon_match_src;
  edm::Handle<reco::CandViewMatchMap> muon_photon_match_map;
};

Zprime2muLeptonProducer::Zprime2muLeptonProducer(const edm::ParameterSet& cfg)
  : muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    electron_src(cfg.getParameter<edm::InputTag>("electron_src")),
    muon_selector(cfg.getParameter<std::string>("muon_cuts")),
    electron_selector(cfg.getParameter<std::string>("electron_cuts")),
    muon_track_for_momentum(cfg.getParameter<std::string>("muon_track_for_momentum")),
    muon_photon_match_src(cfg.getParameter<edm::InputTag>("muon_photon_match_src"))
{
  if (cfg.existsAs<std::vector<std::string> >("muon_tracks_for_momentum"))
    muon_tracks_for_momentum = cfg.getParameter<std::vector<std::string> >("muon_tracks_for_momentum");

  if (muon_tracks_for_momentum.size())
    for (size_t i = 0; i < muon_tracks_for_momentum.size(); ++i)
      produces<pat::MuonCollection>(muon_tracks_for_momentum[i]);
  else
    produces<pat::MuonCollection>("muons");

  produces<pat::ElectronCollection>("electrons");
}

pat::Muon* Zprime2muLeptonProducer::cloneAndSwitchMuonTrack(const pat::Muon& muon) const {
  // Muon mass to make a four-vector out of the new track.
  static const double mass = 0.10566;

  // Start with null track/invalid type before we find the right one.
  reco::TrackRef newTrack;
  patmuon::TrackType type = patmuon::nTrackTypes;

  // If the muon has the track embedded using the UserData mechanism,
  // take it from there first. Otherwise, try to get the track the
  // standard way.
  if (muon.hasUserData(muon_track_for_momentum))
    newTrack = patmuon::userDataTrack(muon, muon_track_for_momentum);
  else {
    type = patmuon::trackNameToType(muon_track_for_momentum);
    newTrack = patmuon::trackByType(muon, type);
  }
  
  // If we didn't find the appropriate track, indicate failure by a
  // null pointer.
  if (newTrack.isNull())
    return 0;

  // Make up a real Muon from the track so found.
  reco::Particle::Point vtx(newTrack->vx(), newTrack->vy(), newTrack->vz());
  reco::Particle::LorentzVector p4;
  const double p = newTrack->p();
  p4.SetXYZT(newTrack->px(), newTrack->py(), newTrack->pz(), sqrt(p*p + mass*mass));

  // The caller will own this pointer and is responsible for deleting
  // it.
  pat::Muon* mu = muon.clone();
  mu->setCharge(newTrack->charge());
  mu->setP4(p4);
  mu->setVertex(vtx);

  // Store the type code for the track used in the pat::Muon so it can
  // be easily recovered later.
  mu->addUserInt("trackUsedForMomentum", type);

  return mu;
}

std::pair<pat::Electron*,int> Zprime2muLeptonProducer::doLepton(const pat::Electron& el, const reco::CandidateBaseRef&) const {
  // Caller will own this pointer.
  pat::Electron* new_el = el.clone();
  int cutFor = new_el == 0 ? -1 : electron_selector(*new_el) ? 0 : 1; // JMTBAD todo heep flags also in this bitword
  return std::make_pair(new_el, cutFor);
}

std::pair<pat::Muon*,int> Zprime2muLeptonProducer::doLepton(const pat::Muon& mu, const reco::CandidateBaseRef& cand) const {
  // Failure is indicated by a null pointer as the first member of the
  // pair.

  // To use one of the refit tracks, we have to have a global muon.
  if (!mu.isGlobalMuon())
    return std::make_pair((pat::Muon*)(0), -1);

  // Copy the input muon, and switch its p4/charge/vtx out for that of
  // the selected refit track.
  pat::Muon* new_mu = cloneAndSwitchMuonTrack(mu);

  if (new_mu == 0)
    return std::make_pair(new_mu, -1);

  // Simply store the photon four-vector for now in the muon as a
  // userData.
  if (muon_photon_match_map.isValid()) {
    const reco::CandViewMatchMap& mm = *muon_photon_match_map;
    if (mm.find(cand) != mm.end()) {
      new_mu->addUserData<reco::Particle::LorentzVector>("photon_p4", mm[cand]->p4());
      new_mu->addUserInt("photon_index", mm[cand].key());
    }
  }

  // Evaluate cuts here with string object selector, and any code that
  // cannot be done in the string object selector (none so far).
  int cutFor = muon_selector(*new_mu) ? 0 : 1;

  return std::make_pair(new_mu, cutFor);
}

template <typename T>
void Zprime2muLeptonProducer::doLeptons(edm::Event& event, const edm::InputTag& src, const std::string& instance_label) const {
  typedef std::vector<T> TCollection;
  edm::Handle<TCollection> leptons; 
  event.getByLabel(src, leptons); 

  if (!leptons.isValid()) {
    edm::LogWarning("LeptonsNotFound") << src << " for " << instance_label << " not found, not producing anything.";
    return;
  }

  edm::Handle<reco::CandidateView> lepton_view;
  event.getByLabel(src, lepton_view);

  std::auto_ptr<TCollection> new_leptons(new TCollection);

  for (size_t i = 0; i < leptons->size(); ++i) {
    std::pair<T*,int> res = doLepton(leptons->at(i), lepton_view->refAt(i));
    if (res.first == 0)
      continue;
    res.first->addUserInt("cutFor", res.second);
    new_leptons->push_back(*res.first);
    delete res.first;
  }

  event.put(new_leptons, instance_label);
}

void Zprime2muLeptonProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  event.getByLabel(muon_photon_match_src, muon_photon_match_map);

  if (muon_tracks_for_momentum.size()) {
    for (size_t i = 0; i < muon_tracks_for_momentum.size(); ++i) {
      muon_track_for_momentum = muon_tracks_for_momentum[i];
      doLeptons<pat::Muon>(event, muon_src, muon_track_for_momentum);
    }
  }
  else
    doLeptons<pat::Muon>(event, muon_src, "muons");

  doLeptons<pat::Electron>(event, electron_src, "electrons");
}

DEFINE_FWK_MODULE(Zprime2muLeptonProducer);
