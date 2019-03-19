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
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerUtilities.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


class EfficiencyFromMCMini : public edm::EDAnalyzer {
 public:
  explicit EfficiencyFromMCMini(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);
 private:
  const unsigned nbins;
  const double min_mass;
  const double max_mass;
  const bool use_resonance_mass;
  const bool use_resonance_mass_denom;
  const edm::InputTag dimuon_src;
  std::vector<std::string> trigger_filters;
  std::vector<std::string> trigger_path_names;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigger_summary_src_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::InputTag src;
  const double hlt_single_min_pt;
  const double hlt_single_max_eta;
  const double acceptance_max_eta_1;
  const double acceptance_max_eta_2;
  const double acceptance_min_pt;
  HardInteraction hardInteraction;

  //-- e.g. { {L3s passing Mu50}, {L3s passing OldMu100}, {L3s passing TkMu100} }
  // or 
  // { {L3s passing Mu27} }
  std::vector<pat::TriggerObjectStandAloneCollection> vec_L3_muons;

  TH1F *GenWeight;

  typedef std::pair<TH1F*, TH1F*> effhistos;

  effhistos accnopt;
  effhistos acceptance;
  effhistos recowrtacc;
  effhistos recowrtacctrig;
  effhistos totalreco;
  std::vector<effhistos> hlt_path_effs;
  effhistos hlt_or_eff;

  effhistos accnopt_bb;
  effhistos acceptance_bb;
  effhistos recowrtacc_bb;
  effhistos recowrtacctrig_bb;
  effhistos totalreco_bb;
  std::vector<effhistos> hlt_path_effs_bb;
  effhistos hlt_or_eff_bb;

  effhistos accnopt_e;
  effhistos acceptance_e;
  effhistos recowrtacc_e;
  effhistos recowrtacctrig_e;
  effhistos totalreco_e;
  std::vector<effhistos> hlt_path_effs_e;
  effhistos hlt_or_eff_e;

  std::pair<TH1F*, TH1F*> make_eff_pair(TString name, TString title) {
    edm::Service<TFileService> fs;
    TH1F* a = fs->make<TH1F>("Num" + name, title, nbins, min_mass, max_mass);  a->Sumw2();
    TH1F* b = fs->make<TH1F>("Den" + name, title, nbins, min_mass, max_mass);  b->Sumw2();
    return std::make_pair(a, b);
  }

  std::string join(const std::vector<std::string>& strings, const std::string& join) {
    std::string r;
    for (size_t i = 0; i < strings.size()-1; ++i) {
      r += strings[i].c_str();
      r += join;
    }
    r += strings.back().c_str();
    return r;
  }

  double _eventWeight;
  bool   _useMadgraphWeight;
  double _theWeight;

};


EfficiencyFromMCMini::EfficiencyFromMCMini(const edm::ParameterSet& cfg)
  : nbins(cfg.getParameter<unsigned>("nbins")),
    min_mass(cfg.getParameter<double>("min_mass")),
    max_mass(cfg.getParameter<double>("max_mass")),
    use_resonance_mass(cfg.getParameter<bool>("use_resonance_mass")),
    use_resonance_mass_denom(cfg.getParameter<bool>("use_resonance_mass_denom")),
    dimuon_src(cfg.getParameter<edm::InputTag>("dimuon_src")),
    trigger_filters(cfg.getParameter<std::vector<std::string>>("trigger_filters")),
    trigger_path_names(cfg.getParameter<std::vector<std::string>>("trigger_path_names")),
    trigger_summary_src_(consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("trigger_summary"))),
    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),
    hlt_single_min_pt(cfg.getParameter<double>("hlt_single_min_pt")),
    hlt_single_max_eta(cfg.getParameter<double>("hlt_single_max_eta")),
    acceptance_max_eta_1(cfg.getParameter<double>("acceptance_max_eta_1")),
    acceptance_max_eta_2(cfg.getParameter<double>("acceptance_max_eta_2")),
    acceptance_min_pt(cfg.getParameter<double>("acceptance_min_pt")),
    hardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")),
    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    _theWeight(1.)
{
  consumes<reco::GenParticleCollection>(hardInteraction.src);
  consumes<pat::CompositeCandidateCollection>(dimuon_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));

  if(trigger_filters.size()<1 || (trigger_filters.size() != trigger_path_names.size()) )
    edm::LogError("EfficiencyFromMC") << "check trigger_filters and trigger_path_names";


  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  GenWeight = fs->make<TH1F>("GenWeight", "weight per event", 4, -2,2);

  accnopt           = make_eff_pair("AccNoPt", "Acceptance (no p_{T} cut) vs. mass");
  acceptance        = make_eff_pair("Acceptance", "Acceptance vs. mass");
  recowrtacc        = make_eff_pair("RecoWrtAcc", "Offline dimuon efficiency for events in acceptance vs. mass");
  hlt_or_eff        = make_eff_pair("TrigWrtAcc",  TString::Format("#varepsilon(%s) vs. mass", join(trigger_path_names, std::string(" || ")).c_str()));
  recowrtacctrig    = make_eff_pair("RecoWrtAccTrig", "Offline dimuon efficiency for events in acceptance and firing the trigger vs. mass");
  totalreco         = make_eff_pair("TotalReco", "Total dimuon efficiency vs. mass");

  accnopt_bb        = make_eff_pair("AccNoPt_bb", "Acceptance (no p_{T} cut) vs. mass");
  acceptance_bb     = make_eff_pair("Acceptance_bb", "Acceptance vs. mass");
  recowrtacc_bb     = make_eff_pair("RecoWrtAcc_bb", "Offline dimuon efficiency for events in acceptance vs. mass");
  hlt_or_eff_bb     = make_eff_pair("TrigWrtAcc_bb",  TString::Format("#varepsilon(%s) vs. mass", join(trigger_path_names, std::string(" || ")).c_str()));
  recowrtacctrig_bb = make_eff_pair("RecoWrtAccTrig_bb", "Offline dimuon efficiency for events in acceptance and firing the trigger vs. mass");
  totalreco_bb      = make_eff_pair("TotalReco_bb", "Total dimuon efficiency vs. mass");

  accnopt_e         = make_eff_pair("AccNoPt_e", "Acceptance (no p_{T} cut) vs. mass");
  acceptance_e      = make_eff_pair("Acceptance_e", "Acceptance vs. mass");
  recowrtacc_e      = make_eff_pair("RecoWrtAcc_e", "Offline dimuon efficiency for events in acceptance vs. mass");
  hlt_or_eff_e      = make_eff_pair("TrigWrtAcc_e",  TString::Format("#varepsilon(%s) vs. mass", join(trigger_path_names, std::string(" || ")).c_str()));
  recowrtacctrig_e  = make_eff_pair("RecoWrtAccTrig_e", "Offline dimuon efficiency for events in acceptance and firing the trigger vs. mass");
  totalreco_e       = make_eff_pair("TotalReco_e", "Total dimuon efficiency vs. mass");

  //-- Efficiency for each HLT path
  for (size_t i = 0; i < trigger_path_names.size(); ++i) {
    hlt_path_effs.push_back(make_eff_pair(TString::Format("HLTPath_%s", trigger_path_names.at(i).c_str()), TString::Format("#varepsilon(%s) vs. mass", trigger_path_names.at(i).c_str())));
    hlt_path_effs_bb.push_back(make_eff_pair(TString::Format("HLTPath_%s_bb", trigger_path_names.at(i).c_str()), TString::Format("#varepsilon(%s) vs. mass", trigger_path_names.at(i).c_str())));
    hlt_path_effs_e.push_back(make_eff_pair(TString::Format("HLTPath_%s_e", trigger_path_names.at(i).c_str()), TString::Format("#varepsilon(%s) vs. mass", trigger_path_names.at(i).c_str())));
  }

}




void EfficiencyFromMCMini::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  hardInteraction.Fill(event);

  edm::Handle<pat::TriggerObjectStandAloneCollection> trigger_summary_src;
  edm::Handle<edm::TriggerResults> triggerBits;
  event.getByToken(trigger_summary_src_, trigger_summary_src);
  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);
  
  if (!hardInteraction.IsValid()) {
    edm::LogWarning("EfficiencyFromMC") << "!hardInteraction.isValid()";
    return;
  }

  //-- Generator weights
  if (_useMadgraphWeight) {
    _eventWeight = 1.;
    _theWeight = 1.;

    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid()){
      _eventWeight = gen_ev_info->weight();
      _theWeight = ( _eventWeight > 0 ) ? 1.0 : -1.0;
    }
    GenWeight->Fill( _theWeight );
  }

  const double m = use_resonance_mass ? hardInteraction.resonance->mass() : hardInteraction.dilepton().mass();
  const double m_denom = use_resonance_mass_denom ? hardInteraction.resonance->mass() : hardInteraction.dilepton().mass();
  
  accnopt.second->Fill( m_denom, _theWeight );
  accnopt_bb.second->Fill( m_denom, _theWeight );
  accnopt_e.second->Fill( m_denom, _theWeight );

  acceptance.second->Fill( m_denom, _theWeight );
  acceptance_bb.second->Fill( m_denom, _theWeight );
  acceptance_e.second->Fill( m_denom, _theWeight );

  totalreco.second->Fill( m_denom, _theWeight );
  totalreco_bb.second->Fill( m_denom, _theWeight );
  totalreco_e.second->Fill( m_denom, _theWeight );

  // both gen leptons in acceptance?

  const double aeta_minus = fabs(hardInteraction.lepMinus->eta());
  const double aeta_plus  = fabs(hardInteraction.lepPlus ->eta());
  const double smaller_aeta = aeta_minus < aeta_plus ? aeta_minus : aeta_plus;
  const double larger_aeta  = aeta_minus > aeta_plus ? aeta_minus : aeta_plus;

  const bool in_acceptance_no_pt = smaller_aeta < acceptance_max_eta_1 && larger_aeta < acceptance_max_eta_2;

  if (in_acceptance_no_pt)
  {
    accnopt.first->Fill( m, _theWeight );
    accnopt_bb.first->Fill( m, _theWeight );
    accnopt_e.first->Fill( m, _theWeight );
  } // No Eta categories for Acceptance

  if (in_acceptance_no_pt &&
      hardInteraction.lepMinus->pt() > acceptance_min_pt &&
      hardInteraction.lepPlus ->pt() > acceptance_min_pt)
  {
    acceptance.first->Fill( m, _theWeight );
    acceptance_bb.first->Fill( m, _theWeight );
    acceptance_e.first->Fill( m, _theWeight );
  }  // No Eta categories for Acceptance
  else
    // Trigger efficiencies below are with respect to events where
    // both gen muons were in acceptance, so stop processing this
    // event now.
    return;


  //--- Fill dimuon events within acceptance
  //--- Note that both numerator and denominator are within eta category
  recowrtacc.second->Fill( m, _theWeight );
  hlt_or_eff.second->Fill( m, _theWeight );
  for (size_t i = 0; i < trigger_path_names.size(); ++i)
    hlt_path_effs[i].second->Fill( m, _theWeight );

  if (aeta_minus <= 1.2 && aeta_plus <= 1.2) {
    recowrtacc_bb.second->Fill( m, _theWeight );
    hlt_or_eff_bb.second->Fill( m, _theWeight );
    for (size_t i = 0; i < trigger_path_names.size(); ++i)
      hlt_path_effs_bb[i].second->Fill( m, _theWeight );
  } // BB

  else {
    recowrtacc_e.second->Fill( m, _theWeight );
    hlt_or_eff_e.second->Fill( m, _theWeight );
    for (size_t i = 0; i < trigger_path_names.size(); ++i)
      hlt_path_effs_e[i].second->Fill( m, _theWeight );
  } // BE + EE


  //-- HLT part
  if (hlt_single_min_pt > 0 || hlt_single_max_eta > 0) {

    vec_L3_muons.clear();
    for(unsigned i_f=0; i_f<trigger_filters.size(); ++i_f) {
      vec_L3_muons.push_back({});
      vec_L3_muons.back().clear();
    }
    for(pat::TriggerObjectStandAlone obj : *trigger_summary_src) {
      obj.unpackPathNames(names);
      obj.unpackFilterLabels(event, *triggerBits); // for 2017~
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {
        for(unsigned i_f=0; i_f<trigger_filters.size(); ++i_f) {
          if (obj.filterLabels()[h] == trigger_filters[i_f]) {
            if(!( obj.pt() > hlt_single_min_pt && fabs(obj.eta()) < hlt_single_max_eta ))
              edm::LogError("EfficiencyFromMC") << "check hlt_single_min_pt and hlt_single_max_eta";
            vec_L3_muons[i_f].push_back(obj);
          }
        }
      }
    }

  }

  bool hlt_or = false;
  for (size_t i_f = 0; i_f < trigger_path_names.size(); ++i_f) {
    if (vec_L3_muons[i_f].size()>0) {
      hlt_path_effs[i_f].first->Fill( m, _theWeight );
      if (aeta_minus <= 1.2 && aeta_plus <= 1.2) {
        hlt_path_effs_bb[i_f].first->Fill( m, _theWeight );
      } // BB
      else {
        hlt_path_effs_e[i_f].first->Fill( m, _theWeight );
      } // BE + EE
      hlt_or = true;
    }
  }

  if (hlt_or) {
    hlt_or_eff.first->Fill( m, _theWeight );
    recowrtacctrig.second->Fill( m, _theWeight );

    if (aeta_minus <= 1.2 && aeta_plus <= 1.2) {
      hlt_or_eff_bb.first->Fill( m, _theWeight );
      recowrtacctrig_bb.second->Fill( m, _theWeight );
    } // BB

    else { 
      hlt_or_eff_e.first->Fill( m, _theWeight );
      recowrtacctrig_e.second->Fill( m, _theWeight );
    } // BE + EE

    // Look for an offline reconstructed dimuon while we're here to
    // measure the total reco efficiency. Loose match in dR for the two
    // muons to the gen two muons before we count it.
    edm::Handle<pat::CompositeCandidateCollection> dimuons;
    event.getByLabel(dimuon_src, dimuons);
    static const double dRmax = 0.25;
    for (pat::CompositeCandidateCollection::const_iterator di = dimuons->begin(), die = dimuons->end(); di != die; ++di) {
      reco::CandidateBaseRef dau0 = dileptonDaughter(*di, 0);
      reco::CandidateBaseRef dau1 = dileptonDaughter(*di, 1);
      if ((reco::deltaR(*dau0, *hardInteraction.lepPlus)  < dRmax || reco::deltaR(*dau1, *hardInteraction.lepPlus)  < dRmax) &&
          (reco::deltaR(*dau0, *hardInteraction.lepMinus) < dRmax || reco::deltaR(*dau1, *hardInteraction.lepMinus) < dRmax)) 
      {
        recowrtacc.first->Fill( m, _theWeight );
        recowrtacctrig.first->Fill( m, _theWeight );
        totalreco.first->Fill( m, _theWeight );

        if (aeta_minus <= 1.2 && aeta_plus <= 1.2) {
          recowrtacc_bb.first->Fill( m, _theWeight );
          recowrtacctrig_bb.first->Fill( m, _theWeight );
          totalreco_bb.first->Fill( m, _theWeight );
        } // BB

        else {
          recowrtacc_e.first->Fill( m, _theWeight );
          recowrtacctrig_e.first->Fill( m, _theWeight );
          totalreco_e.first->Fill( m, _theWeight );
        } // BE+EE

        break;
      }
    }


  }

}



DEFINE_FWK_MODULE(EfficiencyFromMCMini);
