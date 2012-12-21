#include "CLHEP/Random/RandFlat.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class PrescaleToCommon : public edm::EDFilter {
public:
  explicit PrescaleToCommon(const edm::ParameterSet&);

private:
  virtual bool beginRun(edm::Run&, const edm::EventSetup&);
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  const std::string hlt_process_name;
  const std::vector<std::string> trigger_paths;
  const int overall_prescale;
  const bool assume_simulation_has_prescale_1;
  HLTConfigProvider hlt_cfg;
  TH1F* randoms;
  bool _disable;
};

PrescaleToCommon::PrescaleToCommon(const edm::ParameterSet& cfg)
  : hlt_process_name(cfg.getParameter<edm::InputTag>("hlt_src").process()),
    trigger_paths(cfg.getParameter<std::vector<std::string> >("trigger_paths")),
    overall_prescale(cfg.getParameter<int>("overall_prescale")),
    assume_simulation_has_prescale_1(cfg.getParameter<bool>("assume_simulation_has_prescale_1")),
    _disable(cfg.getUntrackedParameter<bool>("disable",false))
{
  edm::Service<TFileService> fs;
  randoms = fs->make<TH1F>("randoms", "", 100, 0, 1);
}

bool PrescaleToCommon::beginRun(edm::Run& run, const edm::EventSetup& setup) {
  bool changed = true;
  if (!hlt_cfg.init(run, setup, hlt_process_name, changed))
    throw cms::Exception("PrescaleToCommon") << "HLTConfigProvider::init failed with process name " << hlt_process_name << "\n";
  
  return true;
}

bool PrescaleToCommon::filter(edm::Event& event, const edm::EventSetup& setup) {
  // JMTBAD move this into common code with CheckPrescale!
  if (_disable) return true;
 
  // The idea behind having trigger_paths be a vector is to handle
  // cases like HLT_Mu24_v1, v2, etc. In that case we expect exactly
  // one of the versions to be available, so check for this. This is
  // not meant to handle checking more-fundamentally different trigger
  // paths, e.g. HLT_Mu15_v2, HLT_Mu24_v1, but only different versions
  // of the "same" path.
  bool found = false;
  std::string trigger_path;
  unsigned path_index = hlt_cfg.size();
  
  for (std::vector<std::string>::const_iterator path = trigger_paths.begin(), end = trigger_paths.end(); path != end; ++path) {
    unsigned ndx = hlt_cfg.triggerIndex(*path);
    if (ndx == hlt_cfg.size())
      continue;

    if (found)
      throw cms::Exception("PrescaleToCommon") << "a version of the trigger path " << *path << " was already found; probably you misconfigured.\n";
    
    found = true;
    trigger_path = *path;
    path_index = ndx;
  }

  if (!found)
    throw cms::Exception("PrescaleToCommon") << "none of the trigger paths specified were found!\n";

  // If the trigger path didn't fire for whatever reason, then go
  // ahead and skip the event.
  edm::Handle<edm::TriggerResults> hlt_results;
  event.getByLabel(edm::InputTag("TriggerResults", "", hlt_process_name), hlt_results);
  if (!hlt_results->accept(path_index))
    return false;
  
  std::pair<int, int> prescales;
  // For MC samples, can assume the prescales are 1.
  if (event.isRealData() || !assume_simulation_has_prescale_1) 
    prescales = hlt_cfg.prescaleValues(event, setup, trigger_path);
  else
    //prescales = std::make_pair(1,1);
    // Do not filter out MC events with prescales=1, apply the
    // appropriate weights later.
    return true;

  //std::cout<<"------PRESCALES: "<<overall_prescale<<"\t"<<prescales.first<<"\t"<<prescales.second<<std::endl;

  const int total_prescale_already = prescales.second * prescales.first;

  if (total_prescale_already > overall_prescale)
    throw cms::Exception("PrescaleToCommon") << "total_prescale_already = " << total_prescale_already << " but overall_prescale requested is " << overall_prescale << "!\n";

  // If we're already there, keep the event.
  if (total_prescale_already == overall_prescale)
    return true;
  
  // Given overall_prescale factor, chance to keep this event is
  // total_prescale_already/overall_prescale. So get uniform number r
  // in ]0,1[, and keep the event if r < chance.
  edm::Service<edm::RandomNumberGenerator> rng;
  if (!rng.isAvailable()) throw cms::Exception("PrescaleToCommon") << "RandomNumberGeneratorService not available!\n";
  CLHEP::RandFlat rand(rng->getEngine());
  const double rnd = rand.fire();
  randoms->Fill(rnd);
  return rnd < double(total_prescale_already)/overall_prescale;
}

DEFINE_FWK_MODULE(PrescaleToCommon);
