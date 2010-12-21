#include "TH1F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerDecision.h"

class EfficiencyFromMC : public edm::EDAnalyzer {
 public:
  explicit EfficiencyFromMC(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);
 private:
  const unsigned nbins;
  const double min_mass;
  const double max_mass;
  const bool use_resonance_mass;
  const bool use_resonance_mass_denom;
  const edm::InputTag dimuon_src;
  const edm::InputTag hlt_obj_src;
  const double hlt_single_min_pt;  
  const double acceptance_max_eta;
  const double acceptance_min_pt;
  TriggerDecision triggerDecision;
  HardInteraction hardInteraction;

  typedef std::pair<TH1F*, TH1F*> effhistos;

  effhistos accnopt;
  effhistos acceptance;
  effhistos recowrtacc;
  effhistos recowrtacctrig;
  effhistos totalreco;
  std::vector<effhistos> l1_path_effs;
  std::vector<effhistos> hlt_path_effs;
  effhistos l1_or_eff;
  effhistos hlt_or_eff;
  effhistos total_trig_eff;

  std::pair<TH1F*, TH1F*> make_eff_pair(TString name, TString title) {
    edm::Service<TFileService> fs;
    TH1F* a = fs->make<TH1F>("Num" + name, title, nbins, min_mass, max_mass);
    TH1F* b = fs->make<TH1F>("Den" + name, title, nbins, min_mass, max_mass);
    return std::make_pair(a, b);
  }
};

std::string join(const std::vector<std::string>& strings, const std::string& join) {
  std::string r;
  for (size_t i = 0; i < strings.size()-1; ++i) {
    r += strings[i].c_str();
    r += join;
  }
  r += strings.back().c_str();
  return r;
}

EfficiencyFromMC::EfficiencyFromMC(const edm::ParameterSet& cfg)
  : nbins(cfg.getParameter<unsigned>("nbins")),
    min_mass(cfg.getParameter<double>("min_mass")),
    max_mass(cfg.getParameter<double>("max_mass")),
    use_resonance_mass(cfg.getParameter<bool>("use_resonance_mass")),
    use_resonance_mass_denom(cfg.getParameter<bool>("use_resonance_mass_denom")),
    dimuon_src(cfg.getParameter<edm::InputTag>("dimuon_src")),
    hlt_obj_src(cfg.getParameter<edm::InputTag>("hlt_obj_src")),
    hlt_single_min_pt(cfg.getParameter<double>("hlt_single_min_pt")),
    acceptance_max_eta(cfg.getParameter<double>("acceptance_max_eta")),
    acceptance_min_pt(cfg.getParameter<double>("acceptance_min_pt")),
    hardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction"))    
{
  triggerDecision.init(cfg.getParameter<edm::ParameterSet>("triggerDecision"));

  edm::Service<TFileService> fs;

  accnopt = make_eff_pair("AccNoPt", "Acceptance (no p_{T} cut) vs. mass");
  acceptance = make_eff_pair("Acceptance", "Acceptance vs. mass");
  recowrtacc = make_eff_pair("RecoWrtAcc", "Offline dimuon efficiency for events in acceptance vs. mass");
  recowrtacctrig = make_eff_pair("RecoWrtAccTrig", "Offline dimuon efficiency for events in acceptance and firing the trigger vs. mass");
  totalreco = make_eff_pair("TotalReco", "Total dimuon efficiency vs. mass");

  for (size_t i = 0; i < triggerDecision.l1_paths().size(); ++i)
    l1_path_effs.push_back(make_eff_pair(TString::Format("L1Path_%i_%s", i, triggerDecision.l1_paths().at(i).c_str()), TString::Format("#varepsilon(%s) vs. mass", triggerDecision.l1_paths().at(i).c_str())));
  for (size_t i = 0; i < triggerDecision.hlt_paths().size(); ++i)
    hlt_path_effs.push_back(make_eff_pair(TString::Format("HLTPath_%i_%s", i, triggerDecision.hlt_paths().at(i).c_str()), TString::Format("#varepsilon(%s) vs. mass", triggerDecision.hlt_paths().at(i).c_str())));
  
  l1_or_eff   = make_eff_pair("L1OrEff",   TString::Format("#varepsilon(%s) vs. mass", join(triggerDecision.l1_paths(),  std::string(" || ")).c_str()));
  hlt_or_eff  = make_eff_pair("HLTOrEff",  TString::Format("#varepsilon(%s) vs. mass", join(triggerDecision.hlt_paths(), std::string(" || ")).c_str()));

  total_trig_eff = make_eff_pair("TotalTrigEff", TString::Format("#varepsilon((%s) && (%s)) vs. mass", join(triggerDecision.l1_paths(), std::string(" || ")).c_str(), join(triggerDecision.hlt_paths(), std::string(" || ")).c_str()));
}

void EfficiencyFromMC::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  hardInteraction.Fill(event);
  triggerDecision.initEvent(event);

  if (!hardInteraction.IsValid()) {
    edm::LogWarning("EfficiencyFromMC") << "!hardInteraction.isValid()";
    return;
  }

  const double m = use_resonance_mass ? hardInteraction.resonance->mass() : hardInteraction.dilepton().mass();
  const double m_denom = use_resonance_mass_denom ? hardInteraction.resonance->mass() : hardInteraction.dilepton().mass();
  accnopt.second->Fill(m_denom);
  acceptance.second->Fill(m_denom);
  totalreco.second->Fill(m_denom);

  if (fabs(hardInteraction.lepMinus->eta()) < acceptance_max_eta && fabs(hardInteraction.lepPlus->eta()) < acceptance_max_eta)
    accnopt.first->Fill(m);

  // both gen leptons in acceptance?
  if (fabs(hardInteraction.lepMinus->eta()) < acceptance_max_eta && fabs(hardInteraction.lepPlus->eta()) < acceptance_max_eta &&
      hardInteraction.lepMinus->pt() > acceptance_min_pt && hardInteraction.lepPlus->pt() > acceptance_min_pt)
    acceptance.first->Fill(m);
  else
    // Trigger efficiencies below are with respect to events where
    // both gen muons were in acceptance, so stop processing this
    // event now.
    return;

  recowrtacc.second->Fill(m);

  for (size_t i = 0; i < triggerDecision.l1_paths().size(); ++i)
    l1_path_effs[i].second->Fill(m);
  for (size_t i = 0; i < triggerDecision.hlt_paths().size(); ++i)
    hlt_path_effs[i].second->Fill(m);
  l1_or_eff.second->Fill(m);
  hlt_or_eff.second->Fill(m);
  total_trig_eff.second->Fill(m);
  
  bool l1_or = false, hlt_or = false;
  
  for (size_t i = 0; i < triggerDecision.l1_paths().size(); ++i) {
    if (triggerDecision.l1_path_pass(i)) {
      l1_path_effs[i].first->Fill(m);
      l1_or = true;
    }
  }

  bool hlt_pass_overridden_pt = true;
  if (hlt_single_min_pt > 0) {
    edm::Handle<std::vector<reco::RecoChargedCandidate> > hlt_objs;
    event.getByLabel(hlt_obj_src, hlt_objs);
    bool pass = false;
    for (std::vector<reco::RecoChargedCandidate>::const_iterator hlt_obj = hlt_objs->begin(), hoe = hlt_objs->end(); hlt_obj != hoe; ++hlt_obj) {
      if (hlt_obj->pt() > hlt_single_min_pt) {
	pass = true;
	break;
      }
    }
    if (!pass) hlt_pass_overridden_pt = false;
  }

  for (size_t i = 0; i < triggerDecision.hlt_paths().size(); ++i) {
    if (triggerDecision.hlt_path_pass(i)) {
      if (i == 0 && !hlt_pass_overridden_pt)
	continue;
      hlt_path_effs[i].first->Fill(m);
      hlt_or = true;
    }
  }
  
  if (l1_or) l1_or_eff.first->Fill(m);
  if (hlt_or) hlt_or_eff.first->Fill(m);
  if (l1_or && hlt_or) total_trig_eff.first->Fill(m);
  
  if (l1_or && hlt_or)
    recowrtacctrig.second->Fill(m);
      
  // Look for an offline reconstructed dimuon while we're here to
  // measure the total reco efficiency. Loose match in dR for the two
  // muons to the gen two muons before we count it.
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  event.getByLabel(dimuon_src, dimuons);
  static const double dRmax = 0.5;
  for (pat::CompositeCandidateCollection::const_iterator di = dimuons->begin(), die = dimuons->end(); di != die; ++di) {
    reco::CandidateBaseRef dau0 = dileptonDaughter(*di, 0);
    reco::CandidateBaseRef dau1 = dileptonDaughter(*di, 1);
    if ((reco::deltaR(*dau0, *hardInteraction.lepPlus)  < dRmax || reco::deltaR(*dau1, *hardInteraction.lepPlus)  < dRmax) &&
	(reco::deltaR(*dau0, *hardInteraction.lepMinus) < dRmax || reco::deltaR(*dau1, *hardInteraction.lepMinus) < dRmax)) {
      recowrtacc.first->Fill(m);
      if (l1_or && hlt_or) {
	recowrtacctrig.first->Fill(m);
	totalreco.first->Fill(m);
      }
      break;
    }
  }
}

DEFINE_FWK_MODULE(EfficiencyFromMC);
