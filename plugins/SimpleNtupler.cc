#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

class SimpleNtupler : public edm::EDAnalyzer {
 public:
  explicit SimpleNtupler(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  struct tree_t {
    unsigned run;
    unsigned lumi;
    unsigned event;
    float dil_mass;
    float dil_pt;
    float dil_rap;
    float dil_eta;
    float dil_phi;
    float cos_angle;
    float vertex_chi2;
    float cos_cs;
    float phi_cs;
    int lep_id[2];
    float lep_pt[2];
    float lep_eta[2];
    float lep_phi[2];
    float lep_tk_pt[2];
    float lep_tk_eta[2];
    float lep_tk_phi[2];
    float lep_glb_pt[2];
    float lep_glb_eta[2];
    float lep_glb_phi[2];
    float lep_triggerMatchPt[2];
    float lep_chi2dof[2];
    float lep_dB[2];
    float lep_sumPt[2];
    float lep_emEt[2];
    float lep_hadEt[2];
    float lep_hoEt[2];
    short lep_tk_numberOfValidTrackerHits[2]; 
    short lep_tk_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidTrackerHits[2]; 
    short lep_glb_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidMuonHits[2];
    short lep_glb_muonStationsWithValidHits[2];
    short lep_numberOfMatches[2];
    bool lep_isGlobalMuon[2];
    bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool HLTPhysicsDeclared;
    bool GoodVtx;
    bool NoScraping;
    bool HLT_Single;
    bool HLT_Double;
  };

  tree_t t;
  TTree* tree;

  edm::InputTag hlt_src;
  edm::InputTag dimu_src;
};

SimpleNtupler::SimpleNtupler(const edm::ParameterSet& cfg)
  : hlt_src(cfg.getParameter<edm::InputTag>("hlt_src")),
    dimu_src(cfg.getParameter<edm::InputTag>("dimu_src"))
{
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("tt", &t, "run/i:lumi:event:dil_mass/F:dil_pt:dil_rap:dil_eta:dil_phi:cos_angle:vertex_chi2:cos_cs:phi_cs:lep_id[2]/I:lep_pt[2]/F:lep_eta[2]:lep_phi[2]:lep_tk_pt[2]:lep_tk_eta[2]:lep_tk_phi[2]:lep_glb_pt[2]:lep_glb_eta[2]:lep_glb_phi[2]:lep_triggerMatchPt[2]:lep_chi2dof[2]:lep_dB[2]:lep_sumPt[2]:lep_emEt[2]:lep_hadEt[2]:lep_hoEt[2]:lep_tk_numberOfValidTrackerHits[2]/S:lep_tk_numberOfValidPixelHits[2]:lep_glb_numberOfValidTrackerHits[2]:lep_glb_numberOfValidPixelHits[2]:lep_glb_numberOfValidMuonHits[2]:lep_glb_muonStationsWithValidHits[2]:lep_numberOfMatches[2]:lep_isGlobalMuon[2]/O:lep_isTrackerMuon[2]:GoodDataRan:HLTPhysicsDeclared:GoodVtx:NoScraping:HLT_Single:HLT_Double");

  tree->SetAlias("OurSel",
		 "("							\
		 "lep_isGlobalMuon[0] && "				\
		 "lep_pt[0] > 20 && "					\
		 "lep_tk_numberOfValidTrackerHits[0] >= 10 && "		\
		 "lep_sumPt[0] / lep_tk_pt[0] < 0.1"			\
		 ") && ("						\
		 "lep_isGlobalMuon[1] && "				\
		 "lep_pt[1] > 20 && "					\
		 "lep_tk_numberOfValidTrackerHits[1] >= 10 && "		\
		 "lep_sumPt[1] / lep_tk_pt[1] < 0.1"			\
		 ") && ( ("						\
		 "abs(lep_dB[0]) < 0.2 && "				\
		 "lep_chi2dof[0] < 10 && "				\
		 "lep_tk_numberOfValidPixelHits[0] >= 1 && "		\
		 "lep_glb_muonStationsWithValidHits[0] >= 2 && "	\
		 "lep_isTrackerMuon[0] && "				\
		 "lep_triggerMatchPt[0] >= 15"				\
		 ") || ("						\
		 "abs(lep_dB[1]) < 0.2 && "				\
		 "lep_chi2dof[1] < 10 && "				\
		 "lep_tk_numberOfValidPixelHits[1] >= 1 && "		\
		 "lep_glb_muonStationsWithValidHits[1] >= 2 && "	\
		 "lep_isTrackerMuon[1] && "				\
		 "lep_triggerMatchPt[1] >= 15"				\
		 ") ) && "						\
		 "lep_id[0] + lep_id[1] == 0 && "			\
		 "cos_angle > -0.9998 && "				\
		 "vertex_chi2 < 10 && "					\
		 "GoodDataRan && "					\
		 "HLTPhysicsDeclared && "				\
		 "NoScraping && "					\
		 "GoodVtx");

  tree->SetAlias("VBTFSel",
		 "lep_isGlobalMuon[0] && "				\
		 "lep_isTrackerMuon[0] && "				\
		 "lep_tk_pt[0] > 20 && "				\
		 "abs(lep_tk_eta[0]) < 2.1 && "				\
		 "abs(lep_dB[0]) < 0.2 && "				\
		 "(lep_sumPt[0] + lep_emEt[0] + lep_hadEt[0]) / lep_tk_pt[0] < 0.15 && " \
		 "lep_glb_numberOfValidTrackerHits[0] > 10 && "		\
		 "lep_glb_numberOfValidPixelHits[0] >= 1 && "		\
		 "lep_glb_numberOfValidMuonHits[0] > 0 && "		\
		 "lep_numberOfMatches[0] >= 2 && "			\
		 "lep_isGlobalMuon[1] && "				\
		 "lep_isTrackerMuon[1] && "				\
		 "lep_tk_pt[1] > 20 && "				\
		 "abs(lep_tk_eta[1]) < 2.1 && "				\
		 "abs(lep_dB[1]) < 0.2 && "				\
		 "(lep_sumPt[1] + lep_emEt[1] + lep_hadEt[1]) / lep_tk_pt[1] < 0.15 && " \
		 "lep_glb_numberOfValidTrackerHits[1] > 10 && "		\
		 "lep_glb_numberOfValidPixelHits[1] >= 1 && "		\
		 "lep_glb_numberOfValidMuonHits[1] > 0 && "		\
		 "lep_numberOfMatches[1] >= 2 && "			\
		 "(lep_triggerMatchPt[0] >= 15 || lep_triggerMatchPt[1] >= 15) && " \
		 "lep_id[0] + lep_id[1] == 0");

  tree->SetAlias("OurNewSel",
		 "lep_isGlobalMuon[0] && "				\
		 "lep_isTrackerMuon[0] && "				\
		 "lep_tk_pt[0] > 20 && "				\
		 "abs(lep_dB[0]) < 0.2 && "				\
		 "lep_sumPt[0] / lep_tk_pt[0] < 0.1 && "		\
		 "lep_glb_numberOfValidTrackerHits[0] > 10 && "		\
		 "lep_glb_numberOfValidPixelHits[0] >= 1 && "		\
		 "lep_glb_numberOfValidMuonHits[0] > 0 && "		\
		 "lep_glb_muonStationsWithValidHits[0] >= 2 && "	\
		 "lep_isGlobalMuon[1] && "				\
		 "lep_isTrackerMuon[1] && "				\
		 "lep_tk_pt[1] > 20 && "				\
		 "abs(lep_dB[1]) < 0.2 && "				\
		 "lep_sumPt[1] / lep_tk_pt[1] < 0.1 && "		\
		 "lep_glb_numberOfValidTrackerHits[1] > 10 && "		\
		 "lep_glb_numberOfValidPixelHits[1] >= 1 && "		\
		 "lep_glb_numberOfValidMuonHits[1] > 0 && "		\
		 "lep_glb_muonStationsWithValidHits[1] >= 2 && "	\
		 "(lep_triggerMatchPt[0] >= 15 || lep_triggerMatchPt[1] >= 15) && " \
		 "lep_id[0] + lep_id[1] == 0 && "			\
		 "cos_angle > -0.9998 && "				\
		 "vertex_chi2 < 10 && "					\
		 "GoodDataRan && "					\
		 "HLTPhysicsDeclared && "				\
		 "NoScraping && "					\
		 "GoodVtx");
}

void SimpleNtupler::analyze(const edm::Event& event, const edm::EventSetup&) {
  memset(&t, 0, sizeof(tree_t));

  t.run = event.id().run();
  t.lumi = event.luminosityBlock();
  t.event = event.id().event();

  edm::Handle<edm::TriggerResults> respat;
  event.getByLabel(edm::InputTag("TriggerResults", "", "PAT"), respat);
  const edm::TriggerNames& namespat = event.triggerNames(*respat);
  if (namespat.triggerIndex("goodDataHLTPhysicsDeclared") < respat->size()) {
    t.GoodDataRan = 1;
    t.HLTPhysicsDeclared = respat->accept(namespat.triggerIndex("goodDataHLTPhysicsDeclared"));
    t.GoodVtx = respat->accept(namespat.triggerIndex("goodDataPrimaryVertexFilter"));
    t.NoScraping = respat->accept(namespat.triggerIndex("goodDataNoScraping"));
  }

  edm::Handle<edm::TriggerResults> reshlt;
  event.getByLabel(hlt_src, reshlt);
  const edm::TriggerNames& nameshlt = event.triggerNames(*reshlt);

  const unsigned r = event.id().run();
  if (r <= 147119) {
    t.HLT_Single = reshlt->accept(nameshlt.triggerIndex("HLT_Mu9")); // changing this to pt > 15 is taken care of by the selection
    t.HLT_Double = reshlt->accept(nameshlt.triggerIndex("HLT_DoubleMu3"));
  }
  else {
    t.HLT_Single = reshlt->accept(nameshlt.triggerIndex("HLT_Mu15_v1"));
    t.HLT_Double = reshlt->accept(nameshlt.triggerIndex("HLT_DoubleMu3_v2"));
  }

  edm::Handle<pat::CompositeCandidateCollection> dils;
  event.getByLabel(dimu_src, dils);

  BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils) {
    t.dil_mass = dil.mass();
    t.dil_pt = dil.pt();
    t.dil_rap = dil.rapidity();
    t.dil_eta = dil.eta();
    t.dil_phi = dil.phi();

    // Only deal with dileptons composed of e,mu for now.
    assert(dil.numberOfDaughters() == 2);
    assert(abs(dil.daughter(0)->pdgId()) == 11 || abs(dil.daughter(0)->pdgId()) == 13);
    assert(abs(dil.daughter(1)->pdgId()) == 11 || abs(dil.daughter(1)->pdgId()) == 13);

    const bool opp_sign = dil.daughter(0)->charge() + dil.daughter(1)->charge() == 0;
    const bool diff_flavor = abs(dil.daughter(0)->pdgId()) != abs(dil.daughter(1)->pdgId());

    for (size_t i = 0; i < 2; ++i) {
      // For e-mu dileptons, put the muon first. For opposite-sign
      // dileptons, always put the negative lepton first. Otherwise
      // don't mess with the order.
      size_t w = i;
      if (diff_flavor)
	w = abs(dil.daughter(i)->pdgId()) == 13 ? 0 : 1;
      else if (opp_sign)
	w = dil.daughter(i)->charge() < 0 ? 0 : 1;

      t.lep_id[w] = dil.daughter(i)->pdgId();
      t.lep_pt[w] = dil.daughter(i)->pt();
      t.lep_eta[w] = dil.daughter(i)->eta();
      t.lep_phi[w] = dil.daughter(i)->phi();

      if (abs(t.lep_id[w]) != 13) {
	t.lep_tk_pt[w] = -999;
	t.lep_tk_eta[w] = -999;
	t.lep_tk_phi[w] = -999;
	t.lep_glb_pt[w] = -999;
	t.lep_glb_eta[w] = -999;
	t.lep_glb_phi[w] = -999;
	t.lep_triggerMatchPt[w] = -999;
	t.lep_chi2dof[w] = -999;
	t.lep_dB[w] = -999;
	t.lep_sumPt[w] = -999;
	t.lep_emEt[w] = -999;
	t.lep_hadEt[w] = -999;
	t.lep_hoEt[w] = -999;
	t.lep_tk_numberOfValidTrackerHits[w] = -999; 
	t.lep_tk_numberOfValidPixelHits[w] = -999;
	t.lep_glb_numberOfValidTrackerHits[w] = -999; 
	t.lep_glb_numberOfValidPixelHits[w] = -999;
	t.lep_glb_numberOfValidMuonHits[w] = -999;
	t.lep_glb_muonStationsWithValidHits[w] = -999;
	t.lep_numberOfMatches[w] = -999;
	t.lep_isGlobalMuon[w] = false;
	t.lep_isTrackerMuon[w] = false;
      }
      else {
	const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
	assert(mu);

	t.lep_tk_pt[w] = mu->innerTrack()->pt();
	t.lep_tk_eta[w] = mu->innerTrack()->eta();
	t.lep_tk_phi[w] = mu->innerTrack()->phi();
	t.lep_glb_pt[w] = mu->globalTrack()->pt();
	t.lep_glb_eta[w] = mu->globalTrack()->eta();
	t.lep_glb_phi[w] = mu->globalTrack()->phi();
	if (!mu->triggerObjectMatchesByPath("HLT_Mu15_v1").empty())
	  t.lep_triggerMatchPt[w] = mu->triggerObjectMatchesByPath("HLT_Mu15_v1").at(0).pt();
	else if (!mu->triggerObjectMatchesByPath("HLT_Mu9").empty())
	  t.lep_triggerMatchPt[w] = mu->triggerObjectMatchesByPath("HLT_Mu9").at(0).pt();
	else
	  t.lep_triggerMatchPt[w] = -999;
	t.lep_chi2dof[w] = mu->globalTrack()->normalizedChi2();
	t.lep_dB[w] = mu->dB();
	t.lep_sumPt[w] = mu->isolationR03().sumPt;
	t.lep_emEt[w] = mu->isolationR03().emEt;
	t.lep_hadEt[w] = mu->isolationR03().hadEt;
	t.lep_hoEt[w] = mu->isolationR03().hoEt;
	t.lep_tk_numberOfValidTrackerHits[w] = mu->innerTrack()->hitPattern().numberOfValidTrackerHits();
	t.lep_tk_numberOfValidPixelHits[w] = mu->innerTrack()->hitPattern().numberOfValidPixelHits();
	t.lep_glb_numberOfValidTrackerHits[w] = mu->globalTrack()->hitPattern().numberOfValidTrackerHits();
	t.lep_glb_numberOfValidPixelHits[w] = mu->globalTrack()->hitPattern().numberOfValidPixelHits();
	t.lep_glb_numberOfValidMuonHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
	t.lep_glb_muonStationsWithValidHits[w] = mu->globalTrack()->hitPattern().muonStationsWithValidHits();
	t.lep_numberOfMatches[w] = mu->numberOfMatches();
	t.lep_isGlobalMuon[w] = mu->isGlobalMuon();
	t.lep_isTrackerMuon[w] = mu->isTrackerMuon();
      }
    }

    t.cos_angle = dil.hasUserFloat("cos_angle") ? dil.userFloat("cos_angle") : 999;
    t.vertex_chi2 = dil.hasUserFloat("vertex_chi2") ? dil.userFloat("vertex_chi2") : -999;

    if (opp_sign) {
      const reco::CandidateBaseRef mum = dileptonDaughterByCharge(dil, -1);
      const reco::CandidateBaseRef mup = dileptonDaughterByCharge(dil, +1);
    
      t.cos_cs = calcCosThetaCSAnal(mum->pz(), mum->energy(), mup->pz(), mup->energy(), dil.pt(), dil.pz(), dil.mass());
      t.phi_cs = calcPhiCSAnal(mum->px(), mum->py(), mup->px(), mup->py(), dil.pt(), dil.eta(), dil.phi(), dil.mass(), true);
    }
    else {
      t.cos_cs = -999;
      t.phi_cs = -999;
    }

    tree->Fill();
  }
}


DEFINE_FWK_MODULE(SimpleNtupler);
