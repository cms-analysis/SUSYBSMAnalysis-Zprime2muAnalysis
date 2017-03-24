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

class EfficiencyFromMCMini : public edm::EDAnalyzer {
 public:
  explicit EfficiencyFromMCMini(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);
 private:
  const unsigned nbins;
  const double min_mass;
  const double max_mass;
  const bool check_l1;
  const bool use_resonance_mass;
  const bool use_resonance_mass_denom;
  const edm::InputTag dimuon_src;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigger_summary_src_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::InputTag l1GtObjectMap;
  edm::InputTag hltResults;
  edm::InputTag src;
  const double hlt_single_min_pt;
  const double hlt_single_max_eta;
  const double checking_prescaled_path;
  const double acceptance_max_eta_1;
  const double acceptance_max_eta_2;
  const double acceptance_min_pt;
  TriggerDecision triggerDecision;
  HardInteraction hardInteraction;

  pat::TriggerObjectStandAloneCollection l3_mus;


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
  effhistos trig_match_eff;

  effhistos accnopt_bb;
  effhistos acceptance_bb;
  effhistos recowrtacc_bb;
  effhistos recowrtacctrig_bb;
  effhistos totalreco_bb;
//  std::vector<effhistos> l1_path_effs_bb;
//  std::vector<effhistos> hlt_path_effs_bb;
//  effhistos l1_or_eff_bb;
  effhistos hlt_or_eff_bb;
//  effhistos total_trig_eff_bb;
  effhistos trig_match_eff_bb;
    
  effhistos accnopt_e;
  effhistos acceptance_e;
  effhistos recowrtacc_e;
  effhistos recowrtacctrig_e;
  effhistos totalreco_e;
//  std::vector<effhistos> l1_path_effs_e;
//  std::vector<effhistos> hlt_path_effs_e;
//  effhistos l1_or_eff_e;
  effhistos hlt_or_eff_e;
//  effhistos total_trig_eff_e;
  effhistos trig_match_eff_e;

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

};





EfficiencyFromMCMini::EfficiencyFromMCMini(const edm::ParameterSet& cfg)
  : nbins(cfg.getParameter<unsigned>("nbins")),
    min_mass(cfg.getParameter<double>("min_mass")),
    max_mass(cfg.getParameter<double>("max_mass")),
    check_l1(cfg.getParameter<bool>("check_l1")),
    use_resonance_mass(cfg.getParameter<bool>("use_resonance_mass")),
    use_resonance_mass_denom(cfg.getParameter<bool>("use_resonance_mass_denom")),
    dimuon_src(cfg.getParameter<edm::InputTag>("dimuon_src")),
    trigger_summary_src_(consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("trigger_summary"))),
    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),
    hlt_single_min_pt(cfg.getParameter<double>("hlt_single_min_pt")),
    hlt_single_max_eta(cfg.getParameter<double>("hlt_single_max_eta")),
    checking_prescaled_path(cfg.getParameter<bool>("checking_prescaled_path")),
    acceptance_max_eta_1(cfg.getParameter<double>("acceptance_max_eta_1")),
    acceptance_max_eta_2(cfg.getParameter<double>("acceptance_max_eta_2")),
    acceptance_min_pt(cfg.getParameter<double>("acceptance_min_pt")),
    hardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction"))    
{
  triggerDecision.init(cfg.getParameter<edm::ParameterSet>("triggerDecision")); 
  consumes<edm::TriggerResults>(cfg.getParameter<edm::ParameterSet>("triggerDecision").getParameter<edm::InputTag>("hltResults"));
  consumes<L1GlobalTriggerObjectMapRecord>(cfg.getParameter<edm::ParameterSet>("triggerDecision").getParameter<edm::InputTag>("l1GtObjectMap"));
  consumes<reco::GenParticleCollection>(hardInteraction.src);
  consumes<pat::CompositeCandidateCollection>(dimuon_src);


  edm::Service<TFileService> fs;

  accnopt = make_eff_pair("AccNoPt", "Acceptance (no p_{T} cut) vs. mass");
  acceptance = make_eff_pair("Acceptance", "Acceptance vs. mass");
  recowrtacc = make_eff_pair("RecoWrtAcc", "Offline dimuon efficiency for events in acceptance vs. mass");
  recowrtacctrig = make_eff_pair("RecoWrtAccTrig", "Offline dimuon efficiency for events in acceptance and firing the trigger vs. mass");
  totalreco = make_eff_pair("TotalReco", "Total dimuon efficiency vs. mass");
  trig_match_eff = make_eff_pair("TrigMatch","N-1 Trigger Matching Efficiency vs. mass");

  accnopt_bb = make_eff_pair("AccNoPt_bb", "Acceptance (no p_{T} cut) vs. mass");
  acceptance_bb = make_eff_pair("Acceptance_bb", "Acceptance vs. mass");
  recowrtacc_bb = make_eff_pair("RecoWrtAcc_bb", "Offline dimuon efficiency for events in acceptance vs. mass");
  recowrtacctrig_bb = make_eff_pair("RecoWrtAccTrig_bb", "Offline dimuon efficiency for events in acceptance and firing the trigger vs. mass");
  totalreco_bb = make_eff_pair("TotalReco_bb", "Total dimuon efficiency vs. mass");
  trig_match_eff_bb = make_eff_pair("TrigMatch_bb","N-1 Trigger Matching Efficiency vs. mass");
  
  accnopt_e = make_eff_pair("AccNoPt_e", "Acceptance (no p_{T} cut) vs. mass");
  acceptance_e = make_eff_pair("Acceptance_e", "Acceptance vs. mass");
  recowrtacc_e = make_eff_pair("RecoWrtAcc_e", "Offline dimuon efficiency for events in acceptance vs. mass");
  recowrtacctrig_e = make_eff_pair("RecoWrtAccTrig_e", "Offline dimuon efficiency for events in acceptance and firing the trigger vs. mass");
  totalreco_e = make_eff_pair("TotalReco_e", "Total dimuon efficiency vs. mass");
  trig_match_eff_e = make_eff_pair("TrigMatch_e","N-1 Trigger Matching Efficiency vs. mass");

  for (size_t i = 0; i < triggerDecision.l1_paths().size(); ++i)
    l1_path_effs.push_back(make_eff_pair(TString::Format("L1Path_%u_%s", unsigned(i), triggerDecision.l1_paths().at(i).c_str()), TString::Format("#varepsilon(%s) vs. mass", triggerDecision.l1_paths().at(i).c_str())));
  for (size_t i = 0; i < triggerDecision.hlt_paths().size(); ++i)
    hlt_path_effs.push_back(make_eff_pair(TString::Format("HLTPath_%u_%s", unsigned(i), triggerDecision.hlt_paths().at(i).c_str()), TString::Format("#varepsilon(%s) vs. mass", triggerDecision.hlt_paths().at(i).c_str())));
  
  l1_or_eff   = make_eff_pair("L1OrEff",   TString::Format("#varepsilon(%s) vs. mass", join(triggerDecision.l1_paths(),  std::string(" || ")).c_str()));
  hlt_or_eff  = make_eff_pair("HLTOrEff",  TString::Format("#varepsilon(%s) vs. mass", join(triggerDecision.hlt_paths(), std::string(" || ")).c_str()));
  hlt_or_eff_bb  = make_eff_pair("HLTOrEff_bb",  TString::Format("#varepsilon(%s) vs. mass", join(triggerDecision.hlt_paths(), std::string(" || ")).c_str()));
  hlt_or_eff_e  = make_eff_pair("HLTOrEff_e",  TString::Format("#varepsilon(%s) vs. mass", join(triggerDecision.hlt_paths(), std::string(" || ")).c_str()));

  total_trig_eff = make_eff_pair("TotalTrigEff", TString::Format("#varepsilon((%s) && (%s)) vs. mass", join(triggerDecision.l1_paths(), std::string(" || ")).c_str(), join(triggerDecision.hlt_paths(), std::string(" || ")).c_str()));
}




void EfficiencyFromMCMini::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  hardInteraction.Fill(event);
  triggerDecision.initEvent(event);

  edm::Handle<pat::TriggerObjectStandAloneCollection> trigger_summary_src;
  edm::Handle<edm::TriggerResults> triggerBits;
  event.getByToken(trigger_summary_src_, trigger_summary_src);
  event.getByToken(triggerBits_, triggerBits);
  
  if (!hardInteraction.IsValid()) {
    edm::LogWarning("EfficiencyFromMC") << "!hardInteraction.isValid()";
    return;
  }

  const double m = use_resonance_mass ? hardInteraction.resonance->mass() : hardInteraction.dilepton().mass();
  const double m_denom = use_resonance_mass_denom ? hardInteraction.resonance->mass() : hardInteraction.dilepton().mass();
  
  accnopt.second->Fill(m_denom);
  accnopt_bb.second->Fill(m_denom);
  accnopt_e.second->Fill(m_denom);

  acceptance.second->Fill(m_denom);
  acceptance_bb.second->Fill(m_denom);
  acceptance_e.second->Fill(m_denom);

  totalreco.second->Fill(m_denom);
  totalreco_bb.second->Fill(m_denom);
  totalreco_e.second->Fill(m_denom);

  // both gen leptons in acceptance?

  const double aeta_minus = fabs(hardInteraction.lepMinus->eta());
  const double aeta_plus  = fabs(hardInteraction.lepPlus ->eta());
  const double smaller_aeta = aeta_minus < aeta_plus ? aeta_minus : aeta_plus;
  const double larger_aeta  = aeta_minus > aeta_plus ? aeta_minus : aeta_plus;

  const bool in_acceptance_no_pt = smaller_aeta < acceptance_max_eta_1 && larger_aeta < acceptance_max_eta_2;

  if (in_acceptance_no_pt)
  {
    accnopt.first->Fill(m);
    accnopt_bb.first->Fill(m);
    accnopt_e.first->Fill(m);
  } // No Eta categories for Acceptance

  if (in_acceptance_no_pt &&
      hardInteraction.lepMinus->pt() > acceptance_min_pt &&
      hardInteraction.lepPlus ->pt() > acceptance_min_pt)
  {
    acceptance.first->Fill(m);
    acceptance_bb.first->Fill(m);
    acceptance_e.first->Fill(m);
  }  // No Eta categories for Acceptance
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

  if (hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2) {
    recowrtacc_bb.second->Fill(m);
    hlt_or_eff_bb.second->Fill(m);
  } // BB

  if ( !(hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2) ) {
    recowrtacc_e.second->Fill(m);
    hlt_or_eff_e.second->Fill(m);
  } // BE + EE

  
  bool l1_or = false, hlt_or = false;

  if (check_l1) {
    for (size_t i = 0; i < triggerDecision.l1_paths().size(); ++i) {
      if (triggerDecision.l1_path_pass(i)) {
	    l1_path_effs[i].first->Fill(m);
	    l1_or = true;
      }
    }
  }
  else
    // If check_l1 is not set, rely on HLT decision without explicitly
    // checking L1 decision.
    l1_or = true;

  bool hlt_pass_overridden = false;
 
  if (hlt_single_min_pt > 0 || hlt_single_max_eta > 0) {
    Zprime2muTriggerPathsAndFilters pandf(event);
    if (!pandf.valid) throw cms::Exception("Zprime2muTriggerPathsAndFilters") << "could not determine the HLT path and filter names for this event\n";

    const edm::TriggerNames &names = event.triggerNames(*triggerBits);
    int j = 0;
    
    l3_mus.clear();

    for (pat::TriggerObjectStandAlone obj : *trigger_summary_src) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);

         for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {
           if (checking_prescaled_path && obj.filterLabels()[h] ==  pandf.prescaled_filter){
     	     //FilterMatched[j] = 1;
	         l3_mus.push_back(obj);
	       }
           if (!checking_prescaled_path && obj.filterLabels()[h] ==  pandf.filter){ 
	         //FilterMatched[j] = 1;
	         l3_mus.push_back(obj);
	       }
           //trigger::TriggerObjectCollection l3_mus = get_L3_muons(event, checking_prescaled_path ? pandf.prescaled_filter : pandf.filter, trigger_summary_src);
         }
      j++;
    
    }

 
    
    bool pass = false;
   
    for (pat::TriggerObjectStandAloneCollection::const_iterator l3_mu = l3_mus.begin(), hoe = l3_mus.end(); l3_mu != hoe; ++l3_mu) {
      if (l3_mu->pt() > hlt_single_min_pt && fabs(l3_mu->eta()) < hlt_single_max_eta) {
	    pass = true;
	    break;
      }
    }
    if (!pass) hlt_pass_overridden = true;
  }

  for (size_t i = 0; i < triggerDecision.hlt_paths().size(); ++i) {
    if (triggerDecision.hlt_path_pass(i)) {
      if (hlt_pass_overridden) edm::LogWarning("EfficiencyFromMC")
	    << "Trigger decision overruled by cuts on Level-3 muon!\n";
      if (i == 0 && hlt_pass_overridden)
	    continue;

      hlt_path_effs[i].first->Fill(m);
      hlt_or = true;
    }
  }

 
  if (l1_or) l1_or_eff.first->Fill(m);
  if (hlt_or) {
    hlt_or_eff.first->Fill(m);

    if (hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2)
    { 
      hlt_or_eff_bb.first->Fill(m);
    } // BB

    if ( !(hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2) )
    { 
      hlt_or_eff_e.first->Fill(m);
    } // BE + EE

  }

  if (l1_or && hlt_or) total_trig_eff.first->Fill(m);
 
  if (l1_or && hlt_or)
  {
    recowrtacctrig.second->Fill(m);

    if (hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2)
    {
      recowrtacctrig_bb.second->Fill(m);
    } // BB

    if ( !(hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2) )
    {
      recowrtacctrig_e.second->Fill(m);
    } // BE + EE 

  }   
  
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
      recowrtacc.first->Fill(m);
      if (l1_or && hlt_or) {
        recowrtacctrig.first->Fill(m);
        trig_match_eff.second->Fill(m);
      } //All
            
      if (hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2) {
        recowrtacc_bb.first->Fill(m);
        if (l1_or && hlt_or) {
          recowrtacctrig_bb.first->Fill(m);
          trig_match_eff_bb.second->Fill(m);
        }
      } //BB

      if ( !(hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2) ) {
        recowrtacc_e.first->Fill(m);
        if (l1_or && hlt_or) {
          recowrtacctrig_e.first->Fill(m);
          trig_match_eff_e.second->Fill(m);
        }
      } //BE+EE

      // Trigger matching
      std::vector< Double_t > l_dR;
      for (pat::TriggerObjectStandAloneCollection::const_iterator l3_Mu = l3_mus.begin(), hme = l3_mus.end(); l3_Mu <= hme; ++l3_Mu) {
        //std::cout << "dR1 : " << reco::deltaR(*l3_Mu,*dau0) << std::endl;
        //std::cout << "dR2 : " << reco::deltaR(*l3_Mu,*dau1) << std::endl;
        Double_t dR = ( reco::deltaR(*l3_Mu,*dau0) < reco::deltaR(*l3_Mu,*dau1) ) ? reco::deltaR(*l3_Mu,*dau0) : reco::deltaR(*l3_Mu,*dau1);
        //std::cout << "dR : " << dR << std::endl;
        l_dR.push_back(dR);
      }
      //std::cout << "Size of l_dR : " << l_dR.size() << std::endl;

      Double_t dR_min = *( std::min_element( l_dR.begin(), l_dR.end() ) );
      //std::cout << "dR_min : " << dR_min << std::endl;

      if ( l1_or && hlt_or && ( dR_min < 0.2 ) ) {
        trig_match_eff.first->Fill(m);
        totalreco.first->Fill(m);

        if (hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2) {
          trig_match_eff_bb.first->Fill(m);
          totalreco_bb.first->Fill(m);
        } // BB

        if ( !(hardInteraction.lepMinus->eta()<=1.2 && hardInteraction.lepPlus->eta()<=1.2 && hardInteraction.lepMinus->eta()>=-1.2 && hardInteraction.lepPlus->eta()>=-1.2) ) {
          trig_match_eff_e.first->Fill(m);
          totalreco_e.first->Fill(m);
        } // BE + EE
      }
      break;
    }
  }
}



DEFINE_FWK_MODULE(EfficiencyFromMCMini);
