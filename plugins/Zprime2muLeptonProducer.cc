#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/ConsumesCollector.h>
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerUtilities.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "TRandom.h"

class Zprime2muLeptonProducer : public edm::EDProducer {
public:
  explicit Zprime2muLeptonProducer(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  pat::Electron* cloneAndSwitchElectronEnergy(const pat::Electron&) const;
  pat::Muon*     cloneAndSwitchMuonTrack     (const pat::Muon&)     const;
 

  void embedTriggerMatch(pat::Muon*, const std::string&, const trigger::TriggerObjectCollection&, std::vector<int>&);
  void embedTriggerMatch_or(pat::Muon*, const std::string&, const trigger::TriggerObjectCollection&, const trigger::TriggerObjectCollection&, std::vector<int>&, std::vector<int>&);

  std::pair<pat::Electron*, int> doLepton(const edm::Event&, const pat::Electron&, const reco::CandidateBaseRef&);
  std::pair<pat::Muon*,     int> doLepton(const edm::Event&, const pat::Muon&,     const reco::CandidateBaseRef&);

  template <typename T> edm::OrphanHandle<std::vector<T> > doLeptons(edm::Event&, const edm::InputTag&, const edm::InputTag&, const std::string&);

  edm::InputTag muon_src;
  edm::InputTag electron_src;
  //two different sources need to be defined in order to have the function doleptons working. Because of consumes, each input tag need to point a specific class.
  edm::InputTag muon_srcSecond;
  edm::InputTag electron_srcSecond;
  StringCutObjectSelector<pat::Muon> muon_selector;
  StringCutObjectSelector<pat::Electron> electron_selector;
  std::string muon_track_for_momentum;
  std::string muon_track_for_momentum_CSC;
  std::string muon_track_for_momentum_primary;
  std::vector<std::string> muon_tracks_for_momentum;
  std::vector<std::string> muon_tracks_for_momentum_CSC;
  edm::InputTag muon_photon_match_src;
  edm::Handle<reco::CandViewMatchMap> muon_photon_match_map;
  double electron_muon_veto_dR;
  std::vector<std::pair<float,float> > muon_eta_phis;
  edm::InputTag trigger_summary_src;
  double trigger_match_max_dR;
  trigger::TriggerObjectCollection L3_muons;
  trigger::TriggerObjectCollection L3_muons_2;
  std::vector<int> L3_muons_matched;
  std::vector<int> L3_muons_matched_2;
  trigger::TriggerObjectCollection prescaled_L3_muons;
  std::vector<int> prescaled_L3_muons_matched;

  
};

Zprime2muLeptonProducer::Zprime2muLeptonProducer(const edm::ParameterSet& cfg)
  : muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    electron_src(cfg.getParameter<edm::InputTag>("electron_src")),
    muon_srcSecond(cfg.getParameter<edm::InputTag>("muon_srcSecond")),
    electron_srcSecond(cfg.getParameter<edm::InputTag>("electron_srcSecond")),
    muon_selector(cfg.getParameter<std::string>("muon_cuts")),
    electron_selector(cfg.getParameter<std::string>("electron_cuts")),
    muon_track_for_momentum(cfg.getParameter<std::string>("muon_track_for_momentum")),
    muon_track_for_momentum_CSC(cfg.getParameter<std::string>("muon_track_for_momentum_CSC")),
    muon_track_for_momentum_primary(muon_track_for_momentum),
    muon_photon_match_src(cfg.getParameter<edm::InputTag>("muon_photon_match_src")),
    electron_muon_veto_dR(cfg.getParameter<double>("electron_muon_veto_dR")),
    trigger_summary_src(cfg.getParameter<edm::InputTag>("trigger_summary_src")),
    trigger_match_max_dR(cfg.getParameter<double>("trigger_match_max_dR"))
{

  // need to consumes all inputTag in 7_6 in order to get access to the information
  consumes<pat::MuonCollection>(muon_src);
  consumes<pat::ElectronCollection>(electron_src);
  consumes<reco::CandViewMatchMap >(muon_photon_match_src);
  consumes<trigger::TriggerEvent>(trigger_summary_src);
  consumes<reco::CandidateView>(muon_srcSecond);
  consumes<reco::CandidateView>(electron_srcSecond);
  

  if (cfg.existsAs<std::vector<std::string> >("muon_tracks_for_momentum"))
    muon_tracks_for_momentum = cfg.getParameter<std::vector<std::string> >("muon_tracks_for_momentum");

  if (cfg.existsAs<std::vector<std::string> >("muon_tracks_for_momentum_CSC"))
    muon_tracks_for_momentum_CSC = cfg.getParameter<std::vector<std::string> >("muon_tracks_for_momentum_CSC");

  if (muon_tracks_for_momentum.size())
    for (size_t i = 0; i < muon_tracks_for_momentum.size(); ++i)
      produces<pat::MuonCollection>(muon_tracks_for_momentum[i]);

  if (muon_tracks_for_momentum_CSC.size())
   for (size_t i = 0; i < muon_tracks_for_momentum_CSC.size(); ++i)
     produces<pat::MuonCollection>(muon_tracks_for_momentum_CSC[i]);

  produces<pat::MuonCollection>("muons");
  produces<pat::ElectronCollection>("electrons");


}

pat::Electron* Zprime2muLeptonProducer::cloneAndSwitchElectronEnergy(const pat::Electron& electron) const {
  // HEEP recommendation is to use always the calorimeter energy
  // instead of the GsfElectron/pat::Electron default, which uses a
  // weighted combination of the calorimeter and track-fit energy. See
  // the section "General Comments" at
  // https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronID.
  pat::Electron* el = electron.clone();
  el->setP4(electron.p4() * (electron.caloEnergy() / electron.energy()));
  return el;
}

pat::Muon* Zprime2muLeptonProducer::cloneAndSwitchMuonTrack(const pat::Muon& muon) const {
  // Muon mass to make a four-vector out of the new track.
  static const double mass = 0.10566;
  /*double a(0.);
  double b(0.);
  double c(0.);
  double sigma(0.);*/

  pat::Muon* mu = muon.clone();
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
  
  mu->setCharge(newTrack->charge());

  mu->setP4(p4);
  mu->setVertex(vtx);

  // Store the type code for the track used in the pat::Muon so it can
  // be easily recovered later.
  mu->addUserInt("trackUsedForMomentum", type);
  return mu;
}

void Zprime2muLeptonProducer::embedTriggerMatch(pat::Muon* new_mu, const std::string& ex, const trigger::TriggerObjectCollection& L3, std::vector<int>& L3_matched) {
  int best = -1;
  float defaultpTvalue = 20.;
  float best_dR = trigger_match_max_dR;
  //std::cout<<"embedded trigger function"<<std::endl;
  for (size_t i = 0; i < L3.size(); ++i) {
    // std::cout<<"embedded trigger function 1"<<std::endl;
    // Skip those already used.
    if (L3_matched[i])
      continue;
    // std::cout<<"embedded trigger function 2"<<std::endl;
    const float dR = reco::deltaR(L3[i], *new_mu);
    if (dR < best_dR) {
      best = int(i);
      best_dR = dR;
    }
  }

  // std::cout<<"best value ="<<best<<std::endl;
  //if (best < 0)
  //return;
   
     if (best >= 0) {
       // std::cout<<"embedded trigger function 3"<<std::endl;
  const trigger::TriggerObject& L3_mu = L3[best];
  L3_matched[best] = 1;
 
  int id = L3_mu.id();
  new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
  new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
  new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
  new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
  //std::cout<<"ex + trigger match = "<<ex<<"...TriggerMatchPt = "<<L3_mu.pt()<<std::endl;
  //std::cout<<"TriggerMatchPt muon producer = "<<new_mu->hasUserFloat(ex + "TriggerMatchPt")<<std::endl;
}
     else {
       // std::cout<<"embedded trigger function 4"<<std::endl;
       new_mu->addUserFloat(ex + "TriggerMatchPt",    defaultpTvalue);
       
       //std::cout<<"TriggerMatchPt muon producer own = "<<new_mu->hasUserFloat(ex + "TriggerMatchPt")<<std::endl;

    }
}

void Zprime2muLeptonProducer::embedTriggerMatch_or(pat::Muon* new_mu, const std::string& ex, const trigger::TriggerObjectCollection& L3, const trigger::TriggerObjectCollection& L3_or, std::vector<int>& L3_matched, std::vector<int>& L3_matched_2) {
    int best_1 = -1;
    int best_2 = -1;
    float defaultpTvalue = 20.;
    float best_dR_1 = trigger_match_max_dR;
    float best_dR_2 = trigger_match_max_dR;
//    std::cout<<"embedded trigger function"<<std::endl;
    for (size_t i = 0; i < L3.size(); ++i) {
        // Skip those already used.
        if (L3_matched[i])
        continue;
        
        const float dR = reco::deltaR(L3[i], *new_mu);
        if (dR < best_dR_1) {
            best_1 = int(i);
            best_dR_1 = dR;
        }
    }
    //std::cout<<"filtro2"<<std::endl;
    for (size_t i = 0; i < L3_or.size(); ++i) {
        // Skip those already used.
        if (L3_matched_2[i])
        continue;
        
        const float dR = reco::deltaR(L3_or[i], *new_mu);
        if (dR < best_dR_2) {
            best_2 = int(i);
            best_dR_2 = dR;
        }
    }
//    std::cout<<" best1 "<<best_1<<" best2 "<<best_2<<" best_dR_1 "<<best_dR_1<<" best_dR_2 "<<best_dR_2<<std::endl;
//    if (best_1 < 0 && best_2 < 0 )
//    return;
    
    if (best_2 <0 && best_1 >= 0){
        const trigger::TriggerObject& L3_mu = L3[best_1];
        L3_matched[best_1] = 1;
        
        int id = L3_mu.id();
        new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
        new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
        new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
        new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
        std::cout<<"ex + trigger match = "<<ex<<"...TriggerMatchPt = "<<L3_mu.pt()<<std::endl;
        std::cout<<"TriggerMatchPt muon producer = "<<new_mu->hasUserFloat(ex + "TriggerMatchPt")<<std::endl;
    }
    else if (best_1 <0 && best_2 >= 0){
        const trigger::TriggerObject& L3_mu = L3_or[best_2];
        L3_matched_2[best_2] = 1;
        
        int id = L3_mu.id();
        new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
        new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
        new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
        new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
    }
    else if (best_1 >=0 && best_2 >=0 && best_dR_1 <= best_dR_2){
        const trigger::TriggerObject& L3_mu = L3[best_1];
        L3_matched[best_1] = 1;
        
        int id = L3_mu.id();
        new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
        new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
        new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
        new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
    }
    else if (best_1 >=0 && best_2 >=0 && best_dR_1 > best_dR_2){
        const trigger::TriggerObject& L3_mu = L3_or[best_2];
        L3_matched_2[best_2] = 1;
        
        int id = L3_mu.id();
        new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
        new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
        new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
        new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
    }
    else {
        // std::cout<<"embedded trigger function 4"<<std::endl;
        new_mu->addUserFloat(ex + "TriggerMatchPt",    defaultpTvalue);
        //std::cout<<"TriggerMatchPt muon producer own = "<<new_mu->hasUserFloat(ex + "TriggerMatchPt")<<std::endl;
    }
}

std::pair<pat::Electron*,int> Zprime2muLeptonProducer::doLepton(const edm::Event& event, const pat::Electron& el, const reco::CandidateBaseRef&) {
  // Electrons can be faked by muons leaving energy in the ECAL. Don't
  // keep the electron if it's within the dR specified of any global
  // muon (a la HEEP). Also keep track of the minimum such dR for
  // embedding in the saved electron.
  double min_muon_dR = 1e99;
  for (std::vector<std::pair<float,float> >::const_iterator mu = muon_eta_phis.begin(), end = muon_eta_phis.end(); mu != end; ++mu) {
    double muon_dR = reco::deltaR(mu->first, mu->second, el.eta(), el.phi());
    if (muon_dR < min_muon_dR)
      min_muon_dR = muon_dR;

    if (muon_dR < electron_muon_veto_dR)
      return std::make_pair((pat::Electron*)(0), -1);
  }

  // Caller will own this pointer.
  pat::Electron* new_el = cloneAndSwitchElectronEnergy(el);

  if (new_el == 0)
    return std::make_pair(new_el, -1);

  // Store the minimum muon dR for later use.
  new_el->addUserFloat("min_muon_dR", float(min_muon_dR));

  // Evaluate cuts here with string object selector, and any code that
  // cannot be done in the string object selector (so far, just the
  // minimum muon dR above).
  int cutFor = electron_selector(*new_el) ? 0 : 1;

  return std::make_pair(new_el, cutFor);
}

std::pair<pat::Muon*,int> Zprime2muLeptonProducer::doLepton(const edm::Event& event, const pat::Muon& mu, const reco::CandidateBaseRef& cand) {
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

  // std::cout<<"doLepton"<<std::endl;

  // Do our own trigger matching and embed the results. After the next
  // pair of function calls, there will be new user floats:
  // {TriggerMatch, prescaledTriggerMatch} x {Pt, Eta, Phi,
  // Charge}. (Maybe embed whole candidates later.)
//  embedTriggerMatch(new_mu, "",          L3_muons,           L3_muons_matched);
  embedTriggerMatch_or(new_mu, "",         L3_muons, L3_muons_2,        L3_muons_matched, L3_muons_matched_2);
  embedTriggerMatch(new_mu, "prescaled", prescaled_L3_muons, prescaled_L3_muons_matched);

  // Evaluate cuts here with string object selector, and any code that
  // cannot be done in the string object selector (none so far).
  int cutFor = muon_selector(*new_mu) ? 0 : 1;
  return std::make_pair(new_mu, cutFor);
}

template <typename T>
edm::OrphanHandle<std::vector<T> > Zprime2muLeptonProducer::doLeptons(edm::Event& event, const edm::InputTag& src, const edm::InputTag& srcSecond,const std::string& instance_label) {
  typedef std::vector<T> TCollection;
  edm::Handle<TCollection> leptons; 
  event.getByLabel(src, leptons); 

  static std::map<std::string, bool> warned;
  if (!leptons.isValid()) {
    if (!warned[instance_label]) {
      edm::LogWarning("LeptonsNotFound") << src << " for " << instance_label << " not found, not producing anything -- not warning any more either.";
      warned[instance_label] = true;
    }
    return edm::OrphanHandle<std::vector<T> >();
  }

  edm::Handle<reco::CandidateView> lepton_view;
  event.getByLabel(srcSecond, lepton_view);

  std::auto_ptr<TCollection> new_leptons(new TCollection);

  for (size_t i = 0; i < leptons->size(); ++i) {
    std::pair<T*,int> res = doLepton(event, leptons->at(i), lepton_view->refAt(i));
    if (res.first == 0)
      continue;
    res.first->addUserInt("cutFor", res.second);
    new_leptons->push_back(*res.first);
    delete res.first;
  }

  return event.put(new_leptons, instance_label);
}

void Zprime2muLeptonProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  // Grab the match map between PAT photons and PAT muons so we can
  // embed the photon candidates later.
  event.getByLabel(muon_photon_match_src, muon_photon_match_map);
  static bool warned = false;
  if (!warned && !muon_photon_match_map.isValid()) {
    edm::LogWarning("PhotonsNotFound") << muon_photon_match_src << " for photons not found, not trying to embed their matches to muons -- not warning any more either.";
    warned = true;
  }

  // Store the L3 muons for trigger match embedding, and clear the
  // vector of matched flags that indicate whether the particular L3
  // muon has been used in a match already. This means matching
  // ambiguities are resolved by original sort order of the
  // candidates; no attempt is done to find a global best
  // matching. (This is how it was done in our configuration of the
  // PATTrigger matcher previously, so why not.) We do this for both
  // the main path and the prescaled path.

 
  Zprime2muTriggerPathsAndFilters pandf(event);
  if (!pandf.valid)
    throw cms::Exception("Zprime2muLeptonProducer") << "could not determine the HLT path and filter names for this event\n";
  L3_muons           = get_L3_muons(event, pandf.filter,           trigger_summary_src);
  L3_muons_2         = get_L3_muons(event, pandf.filter_2,         trigger_summary_src);
  prescaled_L3_muons = get_L3_muons(event, pandf.prescaled_filter, trigger_summary_src);
//    std::cout<<"filter "<<pandf.filter<<std::endl;
//    std::cout<<"L3_muons.size()"<<L3_muons.size()<<std::endl;
//    std::cout<<"filter_2 "<<pandf.filter_2<<std::endl;
//    std::cout<<"L3_muons_2.size()"<<L3_muons_2.size()<<std::endl;

//    std::cout<<"prescaled filter "<<pandf.prescaled_filter<<std::endl;
//    std::cout<<"prescaled_L3_muons.size()"<<prescaled_L3_muons.size()<<std::endl;
  L3_muons_matched.clear();
  L3_muons_matched.resize(L3_muons.size(), 0);
  L3_muons_matched_2.clear();
  L3_muons_matched_2.resize(L3_muons_2.size(), 0);
  prescaled_L3_muons_matched.clear();
  prescaled_L3_muons_matched.resize(prescaled_L3_muons.size(), 0);
  // std::cout<<"quel trigger path = "<<pandf.filter<<std::endl;
  //std::cout<<"nombre de muon that triggered HLT50 = "<<L3_muons.size()<<std::endl;
  // std::cout<<"nombre de muon that triggered prescale = "<<prescaled_L3_muons.size()<<std::endl;
 

  // Using the main choice for momentum assignment, make the primary
  // collection of muons, which will have branch name
  // e.g. leptons:muons.
  muon_track_for_momentum = muon_track_for_momentum_primary;
  edm::OrphanHandle<pat::MuonCollection> muons = doLeptons<pat::Muon>(event, muon_src, muon_srcSecond, "muons");

  // If requested, prepare list of eta,phi coordinates of muons for
  // vetoing electrons near them.
  muon_eta_phis.clear();
  if (electron_muon_veto_dR > 0) {
    if (!muons.isValid())
      throw cms::Exception("Zprime2muLeptonProducer") << "requested to veto electrons near muons, but main muon collection is invalid!\n";

    for (pat::MuonCollection::const_iterator mu = muons->begin(), end = muons->end(); mu != end; ++mu)
      muon_eta_phis.push_back(std::make_pair(mu->eta(), mu->phi()));
  }

  // Now make secondary collections of muons using the momentum
  // assignments specified. They will come out as e.g. leptons:tpfms,
  // leptons:picky, ...
 
 
  for (size_t i = 0; i < muon_tracks_for_momentum.size(); ++i) {
    // Reset the flags so the matching can be redone.
    L3_muons_matched.clear();
    L3_muons_matched.resize(L3_muons.size(), 0);
    L3_muons_matched_2.clear();
    L3_muons_matched_2.resize(L3_muons_2.size(), 0);
    prescaled_L3_muons_matched.clear();
    prescaled_L3_muons_matched.resize(prescaled_L3_muons.size(), 0);
  

 
    muon_track_for_momentum = muon_tracks_for_momentum[i];
    doLeptons<pat::Muon>(event, muon_src, muon_srcSecond, muon_track_for_momentum);
   
      
 }

  // And now make the HEEP electron collection, which will be
  // e.g. leptons:electrons.
  doLeptons<pat::Electron>(event, electron_src, electron_srcSecond, "electrons");
}

DEFINE_FWK_MODULE(Zprime2muLeptonProducer);
