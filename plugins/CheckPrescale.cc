#include <boost/foreach.hpp>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class CheckPrescale : public edm::EDAnalyzer {
public:
  explicit CheckPrescale(const edm::ParameterSet&);

private:
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  const std::string hlt_process_name;
  const std::vector<std::string> trigger_paths;
  HLTConfigProvider hlt_cfg;
};

CheckPrescale::CheckPrescale(const edm::ParameterSet& cfg)
  : hlt_process_name(cfg.getParameter<std::string>("hlt_process_name")),
    trigger_paths(cfg.getParameter<std::vector<std::string> >("trigger_paths"))
{
}

void CheckPrescale::beginRun(const edm::Run& run, const edm::EventSetup& setup) {
  bool changed = true;
  if (!hlt_cfg.init(run, setup, hlt_process_name, changed))
    throw cms::Exception("CheckPrescale") << "HLTConfigProvider::init failed with process name " << hlt_process_name << "\n";

  //if (changed) {
  //  printf("new hlt config in run %u!\ndatasets:", run.run());
  //  hlt_cfg.dump("Datasets");
  //}
}

void CheckPrescale::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  // The idea behind having trigger_paths be a vector is to handle
  // cases like HLT_Mu24_v1, v2, etc. In that case we expect exactly
  // one of the versions to be available, so check for this. This is
  // not meant to handle checking more-fundamentally different trigger
  // paths, e.g. HLT_Mu15_v2, HLT_Mu24_v1, but only different versions
  // of the "same" path.
  bool found = false;
  for (std::vector<std::string>::const_iterator trigger_path = trigger_paths.begin(), end = trigger_paths.end(); trigger_path != end; ++trigger_path) {
    if (hlt_cfg.triggerIndex(*trigger_path) != hlt_cfg.size()) {
      if (found)
	throw cms::Exception("CheckPrescale") << "a version of the trigger path " << *trigger_path << " was already found; probably you misconfigured.\n";
      found = true;
    }
    else
      continue;

    int prescale = hlt_cfg.prescaleValue(event, setup, *trigger_path);
    if (prescale != 1)
      throw cms::Exception("CheckPrescale") << "trigger path " << *trigger_path << " has prescale != 1! prescale = " << prescale << "\n";
  }
}

DEFINE_FWK_MODULE(CheckPrescale);
