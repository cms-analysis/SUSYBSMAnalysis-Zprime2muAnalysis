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
  const std::string trigger_path;
  HLTConfigProvider hlt_cfg;
};

CheckPrescale::CheckPrescale(const edm::ParameterSet& cfg)
  : hlt_process_name(cfg.getParameter<std::string>("hlt_process_name")),
    trigger_path(cfg.getParameter<std::string>("trigger_path"))
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
  int prescale = hlt_cfg.prescaleValue(event, setup, trigger_path);
  if (prescale != 1)
    throw cms::Exception("CheckPrescale") << "trigger path " << trigger_path << " has prescale != 1! prescale = " << prescale << "\n";
}

DEFINE_FWK_MODULE(CheckPrescale);
