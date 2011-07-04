#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
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
  const bool dump_prescales;
  const bool throw_on_prescale;

  struct tree_t {
    unsigned run;
    unsigned lumi;
    int l1;
    int hlt;
  };

  tree_t last_t, t;
  TTree* tree;
};

CheckPrescale::CheckPrescale(const edm::ParameterSet& cfg)
  : hlt_process_name(cfg.getParameter<std::string>("hlt_process_name")),
    trigger_paths(cfg.getParameter<std::vector<std::string> >("trigger_paths")),
    dump_prescales(cfg.getUntrackedParameter<bool>("dump_prescales", false)),
    throw_on_prescale(cfg.getUntrackedParameter<bool>("throw_on_prescale", !dump_prescales))
{
  if (dump_prescales) {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("t", "");
    tree->Branch("tt", &t, "run/i:lumi:l1/I:hlt");
  }
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

    std::pair<int, int> prescales = hlt_cfg.prescaleValues(event, setup, *trigger_path);
    if (prescales.first != 1 || prescales.second != 1) {
      std::ostringstream out;
      out << "trigger path " << *trigger_path << " has prescale != 1! L1 prescale = " << prescales.first << "  HLT prescale: " << prescales.second << "\n";
      if (throw_on_prescale)
        throw cms::Exception("CheckPrescale") << out.str();
      else
        edm::LogWarning("CheckPrescale") << out.str();
    }

    if (dump_prescales) {
      memset(&t, 0, sizeof(tree_t));

      t.run  = event.id().run();
      t.lumi = event.luminosityBlock();
      t.l1   = prescales.first;
      t.hlt  = prescales.second;

      if (t.run == last_t.run && t.lumi == last_t.lumi)
        assert(t.l1 == last_t.l1 && t.hlt == last_t.hlt);
      else
        tree->Fill();
      
      last_t = t;
    }
  }

  if (!found)
    throw cms::Exception("CheckPrescale") << "none of the trigger paths specified were found!\n";
}

DEFINE_FWK_MODULE(CheckPrescale);
