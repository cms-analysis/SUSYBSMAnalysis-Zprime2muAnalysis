#include <memory>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GenEventTopology.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/EMuBackgroundsNtupleDumper.h"

template <typename CollType, typename LabelType>
bool getByLabel(const edm::Event& event, edm::Handle<CollType>& handle,
		LabelType label) {
  bool ok = true;
  try {
    event.getByLabel(label, handle);
  } catch (const cms::Exception& e) {
    edm::LogWarning("EMuBackgroundsNtupleDumper") << "getByLabel failed for " << label;
    ok = false;
  }

  return ok && handle.isValid() && !handle.failedToGet();
}

const int MAXLEPS = 50;

namespace jmt {
  enum {
    WEIGHTTYPE=1<<0,
    WEIGHTSELF=1<<1,
    MET=1<<2,
    TRIG=1<<3,
    JETS=1<<4
  };

  enum {
    LEPEXTRA=1<<0,
    LEPTRACK=1<<1
  };
}

struct jmt_event_t {
  unsigned char status;

  unsigned evt;
  unsigned run;

  float weight;
  unsigned short proc_id;

  unsigned char gen_status;
  short gen_pdg_id[2];
  float gen_pt[2];
  float gen_eta[2];
  float gen_phi[2];
  float gen_energy[2];

  float met_x;
  float met_y;
  unsigned short trig_bits;

  int nleptons; // has to be an int for ROOT to use it for variable size arrays

  unsigned short collection[MAXLEPS];
  unsigned char coll_status;
  unsigned short lep_status[MAXLEPS]; // can't be unsigned char for some reason, gets interpreted by the reader of the tree as a number

  short pdg_id[MAXLEPS];
  float pt[MAXLEPS];
  float eta[MAXLEPS];
  float phi[MAXLEPS];
  float energy[MAXLEPS];
    
  unsigned short iso_njets[MAXLEPS];
  unsigned short iso_ntracks[MAXLEPS];
  float iso_tk[MAXLEPS];
  float iso_ecal[MAXLEPS];
  float iso_hcal[MAXLEPS];

  float tk_pt[MAXLEPS];
  float tk_d0[MAXLEPS];
  float tk_dz[MAXLEPS];
  float tk_chi2dof[MAXLEPS];

  unsigned short tk_hits_px[MAXLEPS];
  unsigned short tk_hits_si[MAXLEPS];
  unsigned short tk_hits_mu[MAXLEPS];
  unsigned short tk_hits_px_lost[MAXLEPS];
  unsigned short tk_hits_si_lost[MAXLEPS];
  unsigned short tk_hits_mu_lost[MAXLEPS];

  short lepton_id[MAXLEPS]; // -1, 0, 1: 0 is unset

  float h_over_e[MAXLEPS];
  float delta_phi_in[MAXLEPS];
  float delta_eta_in[MAXLEPS];
  float sigma_eta_eta[MAXLEPS];
  float e_seed_over_p_in[MAXLEPS];

  int njets; // ditto above
  float jet_pt[MAXLEPS];
  float jet_eta[MAXLEPS];
  float jet_phi[MAXLEPS];
  float jet_energy[MAXLEPS];
};

void jmt_event_add_etc(jmt_event_t& ev, const edm::Event& event) {
  /*
  // Try to get the CSA07 info (won't be there in DY samples).
  edm::Handle<double> weightHandle;
  const double intLumi = 100;
  if (getByLabel(event, weightHandle, edm::InputTag("csa07EventWeightProducer:weight"))) {
    ev.weight = float(*weightHandle);
    ev.proc_id = (unsigned short)(csa07::csa07ProcessId(event, intLumi));
    ev.status |= jmt::WEIGHTTYPE;
  }
  */

  // MET
  edm::Handle<std::vector<pat::MET> > METs;
  if (getByLabel(event, METs, "selectedLayer1METs")) {
    const pat::MET& MET = METs->at(0);
    ev.met_x = MET.px();
    ev.met_y = MET.py();

    ev.status |= jmt::MET;
  }

  // trigger bits
  static const unsigned nPaths = 7;
  static const char* hltPaths[nPaths] = {
    "HLT1MuonNonIso",
    "HLT2MuonNonIso",
    "HLTXElectronMuonRelaxed",
    "HLT1ElectronRelaxed",
    "HLT2ElectronRelaxed",
    "HLT1EMHighEt",
    "HLT1EMVeryHighEt"
  };

  edm::Handle<edm::TriggerResults> hltRes;
  if (getByLabel(event, hltRes, edm::InputTag("TriggerResults::HLT"))) {
    const edm::TriggerNames& hltTrigNames = event.triggerNames(*hltRes);
    
    for (unsigned i = 0; i < nPaths; ++i) {
      int ndx = hltTrigNames.triggerIndex(hltPaths[i]);
      if (hltRes->accept(ndx)) 
	ev.trig_bits |= 1 << i;
    }

    ev.status |= jmt::TRIG;
  }

  edm::Handle<edm::View<reco::Candidate> > jets;
  if (getByLabel(event, jets, "selectedLayer1Jets")) {
    for (unsigned j = 0; j < jets->size(); j++) {
      const reco::Candidate& jet = jets->at(j);
      if (jet.pt() > 30 && fabs(jet.eta()) < 2.4) {
	unsigned nj = ev.njets++;

	ev.jet_pt[nj] 	  = jet.pt();
	ev.jet_eta[nj] 	  = jet.eta();
	ev.jet_phi[nj] 	  = jet.phi();
	ev.jet_energy[nj] = jet.energy();
      }
    }

    ev.status |= jmt::JETS;
  }
}

void jmt_event_add_lepton(jmt_event_t& ev,
			  const edm::Event& event,
			  const reco::Candidate* cand,
			  unsigned short coll) {
  unsigned n = ev.nleptons++;

  ev.collection[n] = coll;
  ev.pdg_id[n] 	  = cand->pdgId();
  ev.pt[n]     	  = cand->pt();
  ev.eta[n]    	  = cand->eta();
  ev.phi[n]    	  = cand->phi();
  ev.energy[n] 	  = cand->energy();

  const pat::Electron* patEl = dynamic_cast<const pat::Electron*>(cand);
  if (patEl != 0)
    ev.lepton_id[n] = patEl->electronID("tight") ? 1 : -1;

  const reco::Track* tk = 0;

  const reco::Muon* muon = dynamic_cast<const reco::Muon*>(cand);
  if (muon != 0) {
    const reco::MuonIsolation& iso = muon->isolationR03();
    ev.iso_njets[n]   = iso.nJets;
    ev.iso_ntracks[n] = iso.nTracks;
    ev.iso_tk[n]      = iso.sumPt;
    ev.iso_ecal[n]    = iso.emEt;
    ev.iso_hcal[n]    = iso.hadEt + iso.hoEt;

    ev.lep_status[n] |= jmt::LEPEXTRA;

    tk = muon->bestTrack();
  }

  if (muon == 0) {
    const reco::GsfElectron* electron = dynamic_cast<const reco::GsfElectron*>(cand);
    if (electron != 0) {
      // electron id variables
      ev.h_over_e[n]     = electron->hadronicOverEm();
      ev.delta_phi_in[n] = electron->deltaPhiSuperClusterTrackAtVtx();
      ev.delta_eta_in[n] = electron->deltaEtaSuperClusterTrackAtVtx();

      float eta = electron->eta();
      assert(0); // prevent trying to use this method until we figure
		 // out how to get covEtaEta in the new version
      //ev.sigma_eta_eta[n] = sqrt(getClusterShape(electron, event)->covEtaEta());
      if (eta >= 1.479) ev.sigma_eta_eta[n] = ev.sigma_eta_eta[n] - 0.02*(fabs(eta) - 2.3);

      ev.e_seed_over_p_in[n] = electron->superCluster()->seed()->energy()/electron->trackMomentumAtVtx().R();

      // electron isolation
      const reco::GsfTrackRef& elTk = electron->gsfTrack();
      edm::Handle<reco::TrackCollection> tracks;
      if (elTk.isNonnull() && getByLabel(event, tracks, "generalTracks")) { // JMTBAD need to not hardcode this and other labels...
	reco::TrackCollection::const_iterator tk = tracks->begin();
	for ( ; tk != tracks->end(); ++tk) {
	  double pt = tk->pt();
	  if (pt > 1 &&
	      fabs(tk->d0()) < 0.1 &&
	      tk->numberOfValidHits() >= 7 &&
	      fabs(elTk->vz() - tk->vz()) < 0.5) {
	  
	    double dR = reco::deltaR(*tk, *elTk);
	    if (dR > 0.01 && dR < 0.3) { // guess at a veto cone 
	      ev.iso_tk[n] += pt;
	      ev.iso_ntracks[n]++;
	    }
	  }
	}

	ev.lep_status[n] |= jmt::LEPEXTRA;
      }

      tk = electron->bestTrack();
    }
  }

  if (tk != 0) {
    ev.tk_pt[n]      = tk->pt();
    ev.tk_d0[n]      = tk->d0();
    ev.tk_dz[n]      = tk->dz();
    ev.tk_chi2dof[n] = tk->normalizedChi2();
      
    const reco::HitPattern& hp = tk->hitPattern();
    ev.tk_hits_px[n]      = hp.numberOfValidPixelHits();
    ev.tk_hits_si[n]      = hp.numberOfValidTrackerHits() - ev.tk_hits_px[n];
    ev.tk_hits_mu[n]      = hp.numberOfValidMuonHits();
    ev.tk_hits_px_lost[n] = hp.numberOfLostPixelHits();
    ev.tk_hits_si_lost[n] = hp.numberOfLostTrackerHits() - ev.tk_hits_px_lost[n];
    ev.tk_hits_mu_lost[n] = hp.numberOfLostMuonHits();

    ev.lep_status[n] |= jmt::LEPTRACK;
  }
}

jmt_event_t jmt_event;

EMuBackgroundsNtupleDumper::EMuBackgroundsNtupleDumper(const edm::ParameterSet& config) 
  : selfProcId(config.getUntrackedParameter<int>("selfProcId", 999))
{
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("events", "my event tree (g'bye edm)");
  
  tree->Branch("status", &jmt_event.status, "status/b");

  tree->Branch("evt", &jmt_event.evt, "evt/i");
  tree->Branch("run", &jmt_event.run, "run/i");

  tree->Branch("weight", &jmt_event.weight, "weight/F");
  tree->Branch("proc_id", &jmt_event.proc_id, "proc_id/s");

  tree->Branch("gen_status", &jmt_event.gen_status, "gen_status/b");
  tree->Branch("gen_pdg_id", jmt_event.gen_pdg_id, "gen_pdg_id[2]/S");
  tree->Branch("gen_pt", jmt_event.gen_pt, "gen_pt[2]/F");
  tree->Branch("gen_eta", jmt_event.gen_eta, "gen_eta[2]/F");
  tree->Branch("gen_phi", jmt_event.gen_phi, "gen_phi[2]/F");
  tree->Branch("gen_energy", jmt_event.gen_energy, "gen_energy[2]/F");

  tree->Branch("met_x", &jmt_event.met_x, "met_x/F");
  tree->Branch("met_y", &jmt_event.met_y, "met_y/F");
  tree->Branch("trig_bits", &jmt_event.trig_bits, "trig_bits/s");

  tree->Branch("nleptons", &jmt_event.nleptons, "nleptons/I");

  tree->Branch("coll_status", &jmt_event.coll_status, "coll_status/b");
  tree->Branch("collection", jmt_event.collection, "collection[nleptons]/s");
  tree->Branch("lep_status", jmt_event.lep_status, "lep_status[nleptons]/s");

  tree->Branch("pdg_id", jmt_event.pdg_id, "pdg_id[nleptons]/S");
  tree->Branch("pt", jmt_event.pt, "pt[nleptons]/F");
  tree->Branch("eta", jmt_event.eta, "eta[nleptons]/F");
  tree->Branch("phi", jmt_event.phi, "phi[nleptons]/F");
  tree->Branch("energy", jmt_event.energy, "energy[nleptons]/F");
    
  tree->Branch("iso_njets", jmt_event.iso_njets, "iso_njets[nleptons]/s");
  tree->Branch("iso_ntracks", jmt_event.iso_ntracks, "iso_ntracks[nleptons]/s");
  tree->Branch("iso_tk", jmt_event.iso_tk, "iso_tk[nleptons]/F");
  tree->Branch("iso_ecal", jmt_event.iso_ecal, "iso_ecal[nleptons]/F");
  tree->Branch("iso_hcal", jmt_event.iso_hcal, "iso_hcal[nleptons]/F");

  tree->Branch("tk_pt", jmt_event.tk_pt, "tk_pt[nleptons]/F");
  tree->Branch("tk_d0", jmt_event.tk_d0, "tk_d0[nleptons]/F");
  tree->Branch("tk_dz", jmt_event.tk_dz, "tk_dz[nleptons]/F");
  tree->Branch("tk_chi2dof", jmt_event.tk_chi2dof, "tk_chi2dof[nleptons]/F");

  tree->Branch("tk_hits_px", jmt_event.tk_hits_px, "tk_hits_px[nleptons]/s");
  tree->Branch("tk_hits_si", jmt_event.tk_hits_si, "tk_hits_si[nleptons]/s");
  tree->Branch("tk_hits_mu", jmt_event.tk_hits_mu, "tk_hits_mu[nleptons]/s");
  tree->Branch("tk_hits_px_lost", jmt_event.tk_hits_px_lost, "tk_hits_px_lost[nleptons]/s");
  tree->Branch("tk_hits_si_lost", jmt_event.tk_hits_si_lost, "tk_hits_si_lost[nleptons]/s");
  tree->Branch("tk_hits_mu_lost", jmt_event.tk_hits_mu_lost, "tk_hits_mu_lost[nleptons]/s");

  tree->Branch("lepton_id", jmt_event.lepton_id, "lepton_id[nleptons]/S");

  tree->Branch("h_over_e", jmt_event.h_over_e, "h_over_e[nleptons]/F");
  tree->Branch("delta_phi_in", jmt_event.delta_phi_in, "delta_phi_in[nleptons]/F");
  tree->Branch("delta_eta_in", jmt_event.delta_eta_in, "delta_eta_in[nleptons]/F");
  tree->Branch("sigma_eta_eta", jmt_event.sigma_eta_eta, "sigma_eta_eta[nleptons]/F");
  tree->Branch("e_seed_over_p_in", jmt_event.e_seed_over_p_in, "e_seed_over_p_in[nleptons]/F");

  tree->Branch("njets", &jmt_event.njets, "njets/I");
  tree->Branch("jet_pt", jmt_event.jet_pt, "jet_pt[njets]/F");
  tree->Branch("jet_eta", jmt_event.jet_eta, "jet_eta[njets]/F");
  tree->Branch("jet_phi", jmt_event.jet_phi, "jet_phi[njets]/F");
  tree->Branch("jet_energy", jmt_event.jet_energy, "jet_energy[njets]/F");
}

void EMuBackgroundsNtupleDumper::analyze(const edm::Event& event, const edm::EventSetup& eSetup) {
  memset(&jmt_event, 0, sizeof(jmt_event));

  jmt_event.evt = event.id().event();
  jmt_event.run = event.id().run();

  // Extract the generator-level information.
  int pdgId[2];
  reco::Particle::LorentzVector genP4[2];
  jmt_event.gen_status = GenEventTopology(event, pdgId, genP4) & 0xFF;
  for (unsigned i = 0; i < 2; i++) {
    jmt_event.gen_pdg_id[i] = pdgId[i];
    jmt_event.gen_pt[i]     = genP4[i].pt();
    jmt_event.gen_eta[i]    = genP4[i].eta();
    jmt_event.gen_phi[i]    = genP4[i].phi();
    jmt_event.gen_energy[i] = genP4[i].energy();
  }

  // Fill MET, jets, trigger bits.
  jmt_event_add_etc(jmt_event, event);

  if (!(jmt_event.status & jmt::WEIGHTTYPE)) {
    jmt_event.proc_id = (unsigned short)(selfProcId);
    jmt_event.status |= jmt::WEIGHTSELF;
  }

  // Flatten lepton collections.
  static const unsigned nColls = 5;
  static const char* collection[nColls] = {
     "muons",
     "bestMuons",
     "selectedLayer1Electrons",
     "heepSelectorVincent",
     "heepSelector"
  };

  for (unsigned icoll = 0; icoll < nColls; icoll++) {
    edm::Handle<edm::View<reco::Candidate> > leptons;
    if (getByLabel(event, leptons, collection[icoll])) {
      for (unsigned ilep = 0; ilep < leptons->size(); ilep++) {
	const reco::Candidate* lep = &leptons->at(ilep);
	jmt_event_add_lepton(jmt_event, event, lep, icoll);
      }

      jmt_event.coll_status |= 1 << icoll;
    }
  }

  tree->Fill();

  static unsigned eventsDone = 0;
  edm::LogInfo("EMuBackgroundsNtupleDumper")
    << "#done: " << ++eventsDone
    << " status: " << unsigned(jmt_event.status)
    << " weight: " << jmt_event.weight
    << " proc_id: " << jmt_event.proc_id << '\n'
    << " gen status: " << unsigned(jmt_event.gen_status)
    << " gen leptons pdg ids: " << jmt_event.gen_pdg_id[0] << " " << jmt_event.gen_pdg_id[1] << '\n'
    << " met x: " << jmt_event.met_x << " y: " << jmt_event.met_y << '\n'
    << " trigbits: " << jmt_event.trig_bits << '\n'
    << " nleptons: " << jmt_event.nleptons;
}

DEFINE_FWK_MODULE(EMuBackgroundsNtupleDumper);
