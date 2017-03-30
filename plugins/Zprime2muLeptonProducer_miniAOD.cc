#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerUtilities.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralizedEndpoint.h"
#include "TLorentzVector.h"


class Zprime2muLeptonProducer_miniAOD : public edm::EDProducer {
public:
  explicit Zprime2muLeptonProducer_miniAOD(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  pat::Electron* cloneAndSwitchElectronEnergy(const pat::Electron&) const;
  pat::Muon*     cloneAndSwitchMuonTrack     (const pat::Muon&, const edm::Event& event)     const;

  void embedTriggerMatch(pat::Muon*, const std::string&, const pat::TriggerObjectStandAloneCollection&, std::vector<int>&);
  void embedTriggerMatch_or(pat::Muon*, const std::string&, const pat::TriggerObjectStandAloneCollection&, const pat::TriggerObjectStandAloneCollection&, std::vector<int>&, std::vector<int>&);
    
  // std::pair<pat::Electron*, int> doLepton(const edm::Event&, const pat::Electron&, const reco::CandidateBaseRef&);
  std::pair<pat::Electron*, int> doLepton(const edm::Event&, const pat::Electron&);
  std::pair<pat::Muon*,     int> doLepton(const edm::Event&, const pat::Muon&,     const reco::CandidateBaseRef&);

  template <typename T> edm::OrphanHandle<std::vector<T> > doLeptons(edm::Event&, const edm::InputTag&, const edm::InputTag&, const std::string&);

  template <typename T> edm::OrphanHandle<std::vector<T> > doLeptons(edm::Event&, edm::Handle<edm::ValueMap<bool> >&, edm::Handle<edm::View<pat::Electron> >&, const std::string&);

  edm::InputTag muon_src;
  edm::InputTag muon_view_src;
  StringCutObjectSelector<pat::Muon> muon_selector;
  std::string muon_track_for_momentum;
  std::string muon_track_for_momentum_primary;
  std::vector<std::string> muon_tracks_for_momentum;
  edm::InputTag muon_photon_match_src;
  edm::Handle<reco::CandViewMatchMap> muon_photon_match_map;
  double electron_muon_veto_dR;
  std::vector<std::pair<float,float> > muon_eta_phis;
  double trigger_match_max_dR;
  
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigger_summary_src_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  
  
  pat::TriggerObjectStandAloneCollection L3_muons;
  pat::TriggerObjectStandAloneCollection L3_muons_2;
  pat::TriggerObjectStandAloneCollection prescaled_L3_muons;
  std::vector<int> L3_muons_matched;
  std::vector<int> L3_muons_matched_2;
  std::vector<int> prescaled_L3_muons_matched;

  edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_;  
  edm::EDGetTokenT<edm::ValueMap<bool> > vidToken_;
  
  
};

Zprime2muLeptonProducer_miniAOD::Zprime2muLeptonProducer_miniAOD(const edm::ParameterSet& cfg)
  : muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    muon_view_src(cfg.getParameter<edm::InputTag>("muon_src")),
    muon_selector(cfg.getParameter<std::string>("muon_cuts")),
    muon_track_for_momentum(cfg.getParameter<std::string>("muon_track_for_momentum")),
    muon_track_for_momentum_primary(muon_track_for_momentum),
    muon_photon_match_src(cfg.getParameter<edm::InputTag>("muon_photon_match_src")),
    electron_muon_veto_dR(cfg.getParameter<double>("electron_muon_veto_dR")),
    //trigger_summary_src(cfg.getParameter<edm::InputTag>("trigger_summary_src")),
    trigger_match_max_dR(cfg.getParameter<double>("trigger_match_max_dR")),
    
    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),
    trigger_summary_src_(consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("trigger_summary"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("prescales"))),
    electronToken_(consumes<edm::View<pat::Electron> >(cfg.getParameter<edm::InputTag>("electron_src"))),
    vidToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("electron_id")))
{
  consumes<edm::View<reco::Candidate>>(muon_view_src);
  consumes<pat::MuonCollection>(muon_src);
  consumes<reco::CandViewMatchMap >(muon_photon_match_src);

  if (cfg.existsAs<std::vector<std::string> >("muon_tracks_for_momentum"))
    muon_tracks_for_momentum = cfg.getParameter<std::vector<std::string> >("muon_tracks_for_momentum");

  if (muon_tracks_for_momentum.size())
    for (size_t i = 0; i < muon_tracks_for_momentum.size(); ++i)
      produces<pat::MuonCollection>(muon_tracks_for_momentum[i]);

  produces<pat::MuonCollection>("muons");
  produces<pat::ElectronCollection>("electrons");
}


pat::Electron* Zprime2muLeptonProducer_miniAOD::cloneAndSwitchElectronEnergy(const pat::Electron& electron) const {
  // HEEP recommendation is to use always the calorimeter energy
  // instead of the GsfElectron/pat::Electron default, which uses a
  // weighted combination of the calorimeter and track-fit energy. See
  // the section "General Comments" at
  // https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronID.
  pat::Electron* el = electron.clone();
  el->setP4(electron.p4() * (electron.caloEnergy() / electron.energy()));
  return el;
}

pat::Muon* Zprime2muLeptonProducer_miniAOD::cloneAndSwitchMuonTrack(const pat::Muon& muon, const edm::Event& event) const {
  
  // Muon mass to make a four-vector out of the new track.
  
  
  
  pat::Muon* mu = muon.clone();
  
  
  
  // Start with null track/invalid type before we find the right one.
  reco::TrackRef newTrack;
  newTrack = muon.tunePMuonBestTrack();
  patmuon::TrackType type = patmuon::nTrackTypes;

  // If the muon has the track embedded using the UserData mechanism,
  // take it from there first. Otherwise, try to get the track the
  // standard way.


  if (muon.pt() > 100.) {
	mu->addUserInt("hasTeVMuons", 1);
  }
  else{

	 mu->addUserInt("hasTeVMuons", 0); 
  }
  if (muon.hasUserData(muon_track_for_momentum))
    newTrack = patmuon::userDataTrack(muon, muon_track_for_momentum);
  else {
 
	type = patmuon::trackNameToType(muon_track_for_momentum);
    	newTrack = patmuon::trackByType(*mu, type);
  }
  
  // If we didn't find the appropriate track, indicate failure by a
  // null pointer.
  
  if (newTrack.isNull()){
    std::cout << "No TuneP" << std::endl;
    return 0;
    
  }
  
  
  
  static const double mass = 0.10566;
  
  reco::Particle::Point vtx(newTrack->vx(), newTrack->vy(), newTrack->vz());
  reco::Particle::LorentzVector p4;

  //////////   Comment following lines to apply pt bias correction /////
   const double p = newTrack->p();  
   p4.SetXYZT(newTrack->px(), newTrack->py(), newTrack->pz(), sqrt(p*p + mass*mass));  
  //////////   Comment previous lines to apply pt bias correction ----->  Uncomment following lines /////


  
	///////// uncomment following lines to apply pt bias correction -----> comment previous lines /////////
//  float phi = newTrack->phi()*TMath::RadToDeg();

//  float mupt = GeneralizedEndpoint().GeneralizedEndpointPt(newTrack->pt(),newTrack->charge(),newTrack->eta(),phi,-1,1); //for DATA
//  float mupt = GeneralizedEndpoint().GeneralizedEndpointPt(newTrack->pt(),newTrack->charge(),newTrack->eta(),phi,0,1);  // for MC


//	float px = mupt*TMath::Cos(newTrack->phi());
//	float py = mupt*TMath::Sin(newTrack->phi());
//	float pz = mupt*TMath::SinH(newTrack->eta());
//	float p = mupt*TMath::CosH(newTrack->eta());
//	p4.SetXYZT(px, py, pz, sqrt(p*p + mass*mass));

// 	std::cout<<"my definition = "<<mupt<<std::endl;
	/////// uncomment previous lines to apply pt bias correction /////////


  mu->setP4(p4);  

  mu->setCharge(newTrack->charge());

  mu->setVertex(vtx);

  mu->addUserInt("trackUsedForMomentum", type);
  
  return mu;
}


void Zprime2muLeptonProducer_miniAOD::embedTriggerMatch(pat::Muon* new_mu, const std::string& ex, const pat::TriggerObjectStandAloneCollection& L3, std::vector<int>& L3_matched) {
  
  int best = -1;
  float defaultpTvalue = 20.;
  float best_dR = trigger_match_max_dR;
  //std::cout << "size of L3 collection: " << L3.size() << std::endl;
  for (size_t i = 0; i < L3.size(); ++i) {
    // Skip those already used.
    if (L3_matched[i])
      continue;

    const float dR = reco::deltaR(L3[i], *new_mu);
    if (dR < best_dR) {
      best = int(i);
      best_dR = dR;
    }
  }

//  if (best < 0)
//    return;
  if (best >= 0) {
  const pat::TriggerObjectStandAlone& L3_mu = L3[best];
  L3_matched[best] = 1;
  
  int id = L3_mu.pdgId();
  new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
  new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
  new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
  new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
  }
    else{
        new_mu->addUserFloat(ex + "TriggerMatchPt",    defaultpTvalue);
    }
  
  
}

void Zprime2muLeptonProducer_miniAOD::embedTriggerMatch_or(pat::Muon* new_mu, const std::string& ex, const pat::TriggerObjectStandAloneCollection& L3, const pat::TriggerObjectStandAloneCollection& L3_or, std::vector<int>& L3_matched, std::vector<int>& L3_matched_2) {
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
        const pat::TriggerObjectStandAlone& L3_mu = L3[best_1];
        L3_matched[best_1] = 1;
        
        int id = L3_mu.pdgId();
        new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
        new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
        new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
        new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
//        std::cout<<"ex + trigger match = "<<ex<<"...TriggerMatchPt = "<<L3_mu.pt()<<std::endl;
//        std::cout<<"TriggerMatchPt muon producer = "<<new_mu->hasUserFloat(ex + "TriggerMatchPt")<<std::endl;
    }
    else if (best_1 <0 && best_2 >= 0){
        const pat::TriggerObjectStandAlone& L3_mu = L3_or[best_2];
        L3_matched_2[best_2] = 1;
        
        int id = L3_mu.pdgId();
        new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
        new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
        new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
        new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
    }
    else if (best_1 >=0 && best_2 >=0 && best_dR_1 <= best_dR_2){
        const pat::TriggerObjectStandAlone& L3_mu = L3[best_1];
        L3_matched[best_1] = 1;
        
        int id = L3_mu.pdgId();
        new_mu->addUserFloat(ex + "TriggerMatchCharge", -id/abs(id));
        new_mu->addUserFloat(ex + "TriggerMatchPt",     L3_mu.pt());
        new_mu->addUserFloat(ex + "TriggerMatchEta",    L3_mu.eta());
        new_mu->addUserFloat(ex + "TriggerMatchPhi",    L3_mu.phi());
    }
    else if (best_1 >=0 && best_2 >=0 && best_dR_1 > best_dR_2){
        const pat::TriggerObjectStandAlone& L3_mu = L3_or[best_2];
        L3_matched_2[best_2] = 1;
        
        int id = L3_mu.pdgId();
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


std::pair<pat::Electron*,int> Zprime2muLeptonProducer_miniAOD::doLepton(const edm::Event& event, const pat::Electron& el) {

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
  // int cutFor = electron_selector(*new_el) ? 0 : 1;
  
  //legacy variable keept as 0, which means that the lepton passes the cuts
  int cutFor = 0; //selection moved to doLeptons where this function doLepton is called only for electrons passes EleID
  

  return std::make_pair(new_el, cutFor);
}

std::pair<pat::Muon*,int> Zprime2muLeptonProducer_miniAOD::doLepton(const edm::Event& event, const pat::Muon& mu, const reco::CandidateBaseRef& cand) {
  // Failure is indicated by a null pointer as the first member of the
  // pair.

  // To use one of the refit tracks, we have to have a global muon.
      
  if (!mu.isGlobalMuon())
    return std::make_pair((pat::Muon*)(0), -1);

  // Copy the input muon, and switch its p4/charge/vtx out for that of
  // the selected refit track.
  
  pat::Muon* new_mu = cloneAndSwitchMuonTrack(mu, event);

  if (new_mu == 0){
    return std::make_pair(new_mu, -1);
    //std::cout << "Warning" << std::endl;
  }  

  // Simply store the photon four-vector for now in the muon as a
  // userData.
  if (muon_photon_match_map.isValid()) {
    const reco::CandViewMatchMap& mm = *muon_photon_match_map;
    if (mm.find(cand) != mm.end()) {
      new_mu->addUserData<reco::Particle::LorentzVector>("photon_p4", mm[cand]->p4());
      new_mu->addUserInt("photon_index", mm[cand].key());
    }    
  }

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
edm::OrphanHandle<std::vector<T> > Zprime2muLeptonProducer_miniAOD::doLeptons(edm::Event& event, const edm::InputTag& src, const edm::InputTag& view_src, const std::string& instance_label) {
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
  event.getByLabel(view_src, lepton_view);

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


template <typename T>
edm::OrphanHandle<std::vector<T> > Zprime2muLeptonProducer_miniAOD::doLeptons(edm::Event& event, edm::Handle<edm::ValueMap<bool> > &vid, edm::Handle<edm::View<pat::Electron> > &patEles, const std::string& instance_label) {
 
  typedef std::vector<T> TCollection;
 
  static std::map<std::string, bool> warned;
  if(!patEles.isValid()){
    if (!warned[instance_label]) {
      edm::LogWarning("LeptonsNotFound") << patEles << " for " << instance_label << " not found, not producing anything -- not warning any more either.";
      warned[instance_label] = true;
    }
    return edm::OrphanHandle<std::vector<T> >();
  }
  
  
  std::auto_ptr<TCollection> new_leptons(new TCollection);
  
  if(patEles.isValid()){ 
    for(auto ele=patEles->begin();ele!=patEles->end();++ele){
      const edm::Ptr<pat::Electron> elePtr(patEles,ele-patEles->begin()); //value map is keyed of edm::Ptrs so we need to make one
     
      //Electron selection is done here using VID
      //A bool true is returned if electron passes ID
      //and only in this case a new lepton is created
      bool passID = (*vid)[elePtr];
      if(passID) {	
	const pat::Electron Electrons = *ele;
	std::pair<T*,int> res = doLepton(event, Electrons);
	if (res.first == 0)
	  continue;
	
	res.first->addUserInt("cutFor", res.second);
	new_leptons->push_back(*res.first);
	delete res.first;
      }
      
    }
  }

  return event.put(new_leptons, instance_label);
}



void Zprime2muLeptonProducer_miniAOD::produce(edm::Event& event, const edm::EventSetup& setup) {
  // Grab the match map between PAT photons and PAT muons so we can
  // embed the photon candidates later.
  //std::cout << event.id() << std::endl;
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
    throw cms::Exception("Zprime2muLeptonProducer_miniAOD") << "could not determine the HLT path and filter names for this event\n";
 

    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> trigger_summary_src;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

    event.getByToken(triggerBits_, triggerBits);
    event.getByToken(trigger_summary_src_, trigger_summary_src);
    event.getByToken(triggerPrescales_, triggerPrescales);
    
    const edm::TriggerNames &names = event.triggerNames(*triggerBits);
    

    int j = 0;
    
    L3_muons.clear();
    L3_muons_2.clear();
    prescaled_L3_muons.clear();
    for (pat::TriggerObjectStandAlone obj : *trigger_summary_src) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);
	
	//if (obj.collection() == "hltL3MuonCandidates::HLT"){
	for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {
	  
	  
	 //this should not be hard coded! 
	//std::cout << obj.filterLabels()[h] << std::endl;   
	if (obj.filterLabels()[h] == pandf.filter){ 
	    //FilterMatched[j] = 1;
	    L3_muons.push_back(obj);
    }
    if (obj.filterLabels()[h] == pandf.filter_2){
         //FilterMatched[j] = 1;
         L3_muons_2.push_back(obj);
    }
    
	  if (obj.filterLabels()[h] ==	pandf.prescaled_filter){
	    //FilterMatched[j] = 1;
	    prescaled_L3_muons.push_back(obj);
	  } 
	 }
       // }
	j++;
    }
//    std::cout<<"quel trigger path = "<<pandf.filter<<std::endl;
//    std::cout<<"nombre de muon that triggered HLT_50 = "<<L3_muons.size()<<std::endl;
//    std::cout<<"quel trigger path = "<<pandf.filter_2<<std::endl;
//    std::cout<<"nombre de muon that triggered HLT_Tk50 = "<<L3_muons_2.size()<<std::endl;
    
    L3_muons_matched.clear();
    L3_muons_matched.resize(L3_muons.size(), 0);
    L3_muons_matched_2.clear();
    L3_muons_matched_2.resize(L3_muons_2.size(), 0);
    prescaled_L3_muons_matched.clear();
    prescaled_L3_muons_matched.resize(prescaled_L3_muons.size(), 0);
//    std::cout<<"filter "<<pandf.filter<<std::endl;
//    std::cout<<"L3_muons.size()"<<L3_muons.size()<<std::endl;
//    std::cout<<"prescaled filter "<<pandf.prescaled_filter<<std::endl;
//    std::cout<<"prescaled_L3_muons.size()"<<prescaled_L3_muons.size()<<std::endl;
    
    // Using the main choice for momentum assignment, make the primary
    // collection of muons, which will have branch name
    // e.g. leptons:muons.
    muon_track_for_momentum = muon_track_for_momentum_primary;
    edm::OrphanHandle<pat::MuonCollection> muons = doLeptons<pat::Muon>(event, muon_src, muon_view_src, "muons");

    // If requested, prepare list of eta,phi coordinates of muons for
    // vetoing electrons near them.
    muon_eta_phis.clear();
    if (electron_muon_veto_dR > 0) {
      if (!muons.isValid())
	throw cms::Exception("Zprime2muLeptonProducer_miniAOD") << "requested to veto electrons near muons, but main muon collection is invalid!\n";
      
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
      doLeptons<pat::Muon>(event, muon_src, muon_view_src, muon_track_for_momentum);
    }
    
    // And now make the HEEP electron collection, which will be
    // e.g. leptons:electrons.
    // doLeptons<pat::Electron>(event, electron_src, electron_view_src, "electrons");
   
       
    edm::Handle<edm::ValueMap<bool> > vid;
    event.getByToken(vidToken_,vid);
        
    edm::Handle<edm::View<pat::Electron> > patEles;
    event.getByToken(electronToken_,patEles);

    doLeptons<pat::Electron>(event, vid, patEles,"electrons");
   
}

DEFINE_FWK_MODULE(Zprime2muLeptonProducer_miniAOD);
