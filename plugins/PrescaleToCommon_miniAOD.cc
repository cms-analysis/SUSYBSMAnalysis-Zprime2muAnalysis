#include "CLHEP/Random/RandFlat.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"


class PrescaleToCommon_miniAOD : public edm::EDFilter {
public:
  explicit PrescaleToCommon_miniAOD(const edm::ParameterSet&);

private:
  //virtual bool beginRun(edm::Run&, const edm::EventSetup&);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  const std::string hlt_process_name;
  edm::InputTag TriggerResults_src;
  edm::InputTag L1Prescale_max_src;
  edm::InputTag L1Prescale_min_src;
  edm::InputTag Prescale_src;
  const std::vector<std::string> trigger_paths;
  const int overall_prescale;
  const bool assume_simulation_has_prescale_1;
  const bool debugInfo;
  const bool no_common;
  HLTConfigProvider hlt_cfg;
  TH1F* randoms;
  bool _disable;
  HLTPrescaleProvider hltPrescaleProvider_;
  
};

PrescaleToCommon_miniAOD::PrescaleToCommon_miniAOD(const edm::ParameterSet& cfg)
  : hlt_process_name(cfg.getParameter<edm::InputTag>("hlt_src").process()),
    TriggerResults_src(cfg.getParameter<edm::InputTag>("TriggerResults_src")),
    L1Prescale_max_src(cfg.getParameter<edm::InputTag>("L1Prescale_max_src")),
    L1Prescale_min_src(cfg.getParameter<edm::InputTag>("L1Prescale_min_src")),
    Prescale_src(cfg.getParameter<edm::InputTag>("Prescale_src")),
    trigger_paths(cfg.getParameter<std::vector<std::string> >("trigger_paths")),
    overall_prescale(cfg.getParameter<int>("overall_prescale")),
    assume_simulation_has_prescale_1(cfg.getParameter<bool>("assume_simulation_has_prescale_1")),
    debugInfo(cfg.getParameter<bool>("debugInfo")),
    no_common(cfg.getParameter<bool>("no_common")),
    _disable(cfg.getUntrackedParameter<bool>("disable",false)),
    hltPrescaleProvider_(cfg, consumesCollector(), *this)
{
  edm::Service<TFileService> fs;
  randoms = fs->make<TH1F>("randoms", "", 100, 0, 1); 
  consumes<edm::TriggerResults>(TriggerResults_src);
  consumes<pat::PackedTriggerPrescales>(L1Prescale_max_src);
  consumes<pat::PackedTriggerPrescales>(L1Prescale_min_src);
  consumes<pat::PackedTriggerPrescales>(Prescale_src);
}

//bool PrescaleToCommon_miniAOD::beginRun(edm::Run& run, const edm::EventSetup& setup) {
//  bool changed = true;
//  if (!hlt_cfg.init(run, setup, hlt_process_name, changed))
//    throw cms::Exception("PrescaleToCommon_miniAOD") << "HLTConfigProvider::init failed with process name " << hlt_process_name << "\n";
//  return true;
//}

void PrescaleToCommon_miniAOD::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{
    bool changed = true;
    // if (!hlt_cfg.init(run, setup, hlt_process_name, changed))
        if (!hltPrescaleProvider_.init(run,setup,hlt_process_name,changed))
        throw cms::Exception("PrescaleToCommon_miniAOD") << "HLTConfigProvider::init failed with process name " << hlt_process_name << "\n";
}

bool PrescaleToCommon_miniAOD::filter(edm::Event& event, const edm::EventSetup& setup) {
  // JMTBAD move this into common code with CheckPrescale!
  if (_disable) return true;
 
  // The idea behind having trigger_paths be a vector is to handle
  // cases like HLT_Mu24_v1, v2, etc. In that case we expect exactly
  // one of the versions to be available, so check for this. This is
  // not meant to handle checking more-fundamentally different trigger
  // paths, e.g. HLT_Mu15_v2, HLT_Mu24_v1, but only different versions
  // of the "same" path.
  //
  // CJS - If different trigger types have the same exact prescales then this is OK
  // e.g. HLT_Mu27_v* and HLT_TkMu27_v*


   edm::Handle<pat::PackedTriggerPrescales> hltPrescales;
   edm::Handle<pat::PackedTriggerPrescales> L1Prescales_max;
   edm::Handle<pat::PackedTriggerPrescales> L1Prescales_min;

   event.getByLabel(Prescale_src, hltPrescales);
   event.getByLabel(L1Prescale_max_src, L1Prescales_max);
   event.getByLabel(L1Prescale_min_src, L1Prescales_min);




  HLTConfigProvider const& hlt_cfg = hltPrescaleProvider_.hltConfigProvider();
  //bool found = false;
  std::string trigger_path;
  //unsigned path_index = hlt_cfg.size();
  
 
  if (debugInfo) {
    const std::vector<std::string>& pathList = hlt_cfg.triggerNames();
    std::cout<<"path size "<<pathList.size()<<std::endl;
  }

  edm::Handle<edm::TriggerResults> hlt_results;
  event.getByLabel(TriggerResults_src, hlt_results);

  //edm::TriggerResultsByName trbn = event.triggerResultsByName(*hlt_results);

  //const edm::TriggerNames &names = event.triggerNames(*hlt_results);
  //hltPrescales->setTriggerNames(event.triggerNames(*hlt_results));
  //hltPrescales->setTriggerNames();

  std::vector<unsigned> found_index;
  std::vector<std::string> found_path;
  std::vector<bool> found_result;
  std::vector<double> found_l1_prescale;
  std::vector<double> found_hlt_prescale;
    
  for (std::vector<std::string>::const_iterator path = trigger_paths.begin(), end = trigger_paths.end(); path != end; ++path) {
      if (debugInfo) {
          std::cout<<" trigger_path "<<*path<<std::endl;
          std::cout<<" hlt_cfg.triggerIndex(*path) "<<hlt_cfg.triggerIndex(*path)<<" hlt_cfg.size() "<<hlt_cfg.size()<<std::endl;
      }

    unsigned ndx = hlt_cfg.triggerIndex(*path);
    if (ndx == hlt_cfg.size())
      continue;

    //if (found)
      //throw cms::Exception("PrescaleToCommon_miniAOD") << "a version of the trigger path " << *path << " was already found; probably you misconfigured.\n";
    
    //found = true;
    //trigger_path = *path;
    //path_index = ndx;
    found_path.push_back(*path);
    found_index.push_back(ndx);
    found_result.push_back(hlt_results->accept(ndx));
    //found_result.push_back(trbn.accept(ind));
    //found_result.push_back(hltPrescales->triggerResults().accept(*path));
    found_hlt_prescale.push_back(hltPrescales->getPrescaleForIndex(ndx));
    //found_hlt_prescale.push_back(hltPrescales->getPrescaleForName(*path));
    //L1Prescale_max = L1Prescales_max->getPrescaleForName(trigger_path);
    //found_l1_prescale.push_back(L1Prescales_min->getPrescaleForName(*path));
    found_l1_prescale.push_back(L1Prescales_min->getPrescaleForIndex(ndx));
  }

  /*
  if (!found)
    throw cms::Exception("PrescaleToCommon_miniAOD") << "none of the trigger paths specified were found!\n";
  */

  // If none of the triggers are in the event then throw an error
  if (found_path.size()==0 || /*found_index.size()==0 ||*/ found_result.size()==0)
    throw cms::Exception("PrescaleToCommon_miniAOD") << "none of the trigger paths specified were found!\n";
  // This should never happen
  //if (found_path.size()!=found_index.size()) {
  if (found_path.size()!=found_result.size()) {
    throw cms::Exception("PrescaleToCommon_miniAOD") << "path and index vectors not same size!?!\n";
  }

  if (debugInfo) {
      for (unsigned int i=0; i<found_path.size(); i++) {
          std::cout << "Path " << found_path[i] << /*" index " << found_index[i] << */" result " << found_result[i] 
              << " l1 prescale " << found_l1_prescale[i] << " hlt prescale " << found_hlt_prescale[i] << std::endl;
      }
  }

  // Find result of trigger paths
  bool result = false;
  if (std::find(found_result.begin(), found_result.end(), true) != found_result.end()){
      result = true;
  }

  // if we do not want a common prescale or event is MC or none of the triggers fired in the event, then return the trigger result
  if (no_common || !event.isRealData() || result==false) {
      return result;
  }

  // Everything here assumes that at least one of the trigger paths has fired
  
  // Check if prescales for triggers are idential
  if (std::adjacent_find(found_hlt_prescale.begin(), found_hlt_prescale.end(), std::not_equal_to<>()) != found_hlt_prescale.end() || 
          std::adjacent_find(found_l1_prescale.begin(), found_l1_prescale.end(), std::not_equal_to<>()) != found_l1_prescale.end()) {
    throw cms::Exception("PrescaleToCommon_miniAOD") << "prescales not identical for all paths\n";
  }

  double hltPrescale = found_hlt_prescale[0]; // OK to use first element since all are identical
  double L1Prescale_min = found_l1_prescale[0];
  //double L1Prescale_max = 1;

  /*
  // For MC samples, can assume the prescales are 1.
  if (event.isRealData() || !assume_simulation_has_prescale_1) {
    hltPrescale = hltPrescales->getPrescaleForIndex(path_index);
    //L1Prescale_max = L1Prescales_max->getPrescaleForName(trigger_path);
    L1Prescale_min = L1Prescales_min->getPrescaleForIndex(path_index);
  }
  else
    //prescales = std::make_pair(1,1);
    // Do not filter out MC events with prescales=1, apply the
    // appropriate weights later.
    return true;
  */

  //std::cout<<"------PRESCALES: "<<overall_prescale<<"\t"<<prescales.first<<"\t"<<prescales.second<<std::endl;


  const int total_prescale_already = hltPrescale*L1Prescale_min;
  if (debugInfo) {
      std::cout << "hltPrescale " << hltPrescale << " L1Prescale_min " << L1Prescale_min << " total prescale " << total_prescale_already << std::endl;
  }

  if (total_prescale_already > overall_prescale)
    throw cms::Exception("PrescaleToCommon_miniAOD") << "total_prescale_already = " << total_prescale_already << " but overall_prescale requested is " << overall_prescale << "!\n";

  // If we're already there, keep the event.
  if (total_prescale_already == overall_prescale)
    return true;
  
  // Given overall_prescale factor, chance to keep this event is
  // total_prescale_already/overall_prescale. So get uniform number r
  // in ]0,1[, and keep the event if r < chance.
  edm::Service<edm::RandomNumberGenerator> rng;
  if (!rng.isAvailable()) throw cms::Exception("PrescaleToCommon_miniAOD") << "RandomNumberGeneratorService not available!\n";
  CLHEP::RandFlat rand(rng->getEngine(event.streamID()));
  const double rnd = rand.fire();
  randoms->Fill(rnd);
  return rnd < double(total_prescale_already)/overall_prescale;
}

DEFINE_FWK_MODULE(PrescaleToCommon_miniAOD);
