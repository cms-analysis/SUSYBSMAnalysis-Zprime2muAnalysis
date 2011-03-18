#include <boost/foreach.hpp>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class VetoOtherDataset : public edm::EDFilter {
public:
  explicit VetoOtherDataset(const edm::ParameterSet&);

private:
  virtual bool beginRun(edm::Run&, const edm::EventSetup&);
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  const std::string hlt_process_name;
  const std::string dataset_to_veto;
  HLTConfigProvider hlt_cfg;
};

VetoOtherDataset::VetoOtherDataset(const edm::ParameterSet& cfg)
  : hlt_process_name(cfg.getParameter<std::string>("hlt_process_name")),
    dataset_to_veto(cfg.getParameter<std::string>("dataset_to_veto"))
{
}

bool VetoOtherDataset::beginRun(edm::Run& run, const edm::EventSetup& setup) {
  bool changed = true;
  if (!hlt_cfg.init(run, setup, hlt_process_name, changed))
    throw cms::Exception("VetoOtherDataset") << "HLTConfigProvider::init failed with process name " << hlt_process_name << "\n";

  //if (changed) {
  //  printf("new hlt config in run %u!\ndatasets:", run.run());
  //  hlt_cfg.dump("Datasets");
  //}

  return true;
}

bool VetoOtherDataset::filter(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<edm::TriggerResults> trigger_results;
  event.getByLabel(edm::InputTag("TriggerResults", "", hlt_process_name), trigger_results);
  const edm::TriggerNames& trigger_names = event.triggerNames(*trigger_results);
  const size_t num_paths = trigger_results->size();

  // If any one of the paths in the dataset-to-veto fired, filter out this event.
  BOOST_FOREACH(const std::string& path, hlt_cfg.datasetContent(dataset_to_veto)) {
    const size_t index = trigger_names.triggerIndex(path);
    if (index >= num_paths)
      throw cms::Exception("VetoOtherDataset") << "path to veto " << path << " is not found in the TriggerResults object!\n";
    if (trigger_results->accept(index))
      return false;
  }

  // None of the paths in the dataset-to-veto fired, so keep this event.
  return true;
}

DEFINE_FWK_MODULE(VetoOtherDataset);
