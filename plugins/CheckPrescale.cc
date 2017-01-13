#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

class CheckPrescale : public edm::EDAnalyzer {
public:
  explicit CheckPrescale(const edm::ParameterSet&);

private:
  //virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  const std::string hlt_process_name;
  const std::vector<std::string> trigger_paths;
  HLTConfigProvider hlt_cfg;
  const bool dump_prescales;
  const bool throw_on_prescale;
  HLTPrescaleProvider hltPrescaleProvider_;

  struct tree_t {
    unsigned run;
    unsigned lumi;
    unsigned event;
    int l1;
    int hlt;
  };

  tree_t last_t, t;
  TTree* tree;
};

CheckPrescale::CheckPrescale(const edm::ParameterSet& cfg)
  : hlt_process_name(cfg.getParameter<edm::InputTag>("hlt_src").process()),
    trigger_paths(cfg.getParameter<std::vector<std::string> >("trigger_paths")),
    dump_prescales(cfg.getUntrackedParameter<bool>("dump_prescales", false)),
    throw_on_prescale(cfg.getUntrackedParameter<bool>("throw_on_prescale", !dump_prescales)),
    hltPrescaleProvider_(cfg, consumesCollector(), *this)
{
  if (dump_prescales) {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("t", "");
    tree->Branch("tt", &t, "run/i:lumi:event:l1/I:hlt");
  }
}

//void CheckPrescale::beginRun(const edm::Run& run, const edm::EventSetup& setup) {
void CheckPrescale::beginRun(edm::Run const& run, edm::EventSetup const& setup){

  bool changed = true;
  if (!hltPrescaleProvider_.init(run, setup, hlt_process_name, changed))
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
  HLTConfigProvider const& hlt_cfg = hltPrescaleProvider_.hltConfigProvider();
  bool found = false;
  for (std::vector<std::string>::const_iterator trigger_path = trigger_paths.begin(), end = trigger_paths.end(); trigger_path != end; ++trigger_path) {
    if (hlt_cfg.triggerIndex(*trigger_path) != hlt_cfg.size()) {
      if (found)
	throw cms::Exception("CheckPrescale") << "a version of the trigger path " << *trigger_path << " was already found; probably you misconfigured.\n";
      found = true;
    }
    else
      continue;

    //std::pair<int, int> prescales = hlt_cfg.prescaleValues(event, setup, *trigger_path);
      std::pair<int, int> prescales;
      prescales.first=1;
      std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail;
      prescalesInDetail = hltPrescaleProvider_.prescaleValuesInDetail(event, setup, *trigger_path);
    //std::cout << "trigger path " << *trigger_path << " has L1/HLT prescales: " << prescales.first << "/" << prescales.second << "   event r/l/e: " << event.id().run() << "/" << event.luminosityBlock() << "/" << event.id().event() << " orbitNum: " << event.orbitNumber() << " BX: " << event.bunchCrossing() << " seconds since the epoch: " << event.time().unixTime() << "\n";

     
      for (unsigned int i=0; i<prescalesInDetail.first.size(); ++i) {
          //std::cout<<" prescalesInDetail.first[i].first "<<prescalesInDetail.first[i].first<<" prescalesInDetail.first[i].second "<<prescalesInDetail.first[i].second<<std::endl;
          prescales.first *= prescalesInDetail.first[i].second;
          //std::cout<<" prescales.first "<<prescales.first<<std::endl;
      }
      
      //std::cout<<" prescalesInDetail.second "<<prescalesInDetail.second<<std::endl;
      prescales.second = prescalesInDetail.second;

    if (prescales.first != 1 || prescales.second != 1) {
      std::ostringstream out;
      out << "trigger path " << *trigger_path << " has prescale != 1! L1 prescale = " << prescales.first << "  HLT prescale: " << prescales.second;
      if (throw_on_prescale)
        throw cms::Exception("CheckPrescale") << out.str() << "\n";
      else
        edm::LogWarning("CheckPrescale") << out.str();
    }

    if (dump_prescales) {
      memset(&t, 0, sizeof(tree_t));

      t.run  = event.id().run();
      t.lumi = event.luminosityBlock();
      t.event = event.id().event();
      t.l1   = prescales.first;
      t.hlt  = prescales.second;

      //if (t.run == last_t.run && t.lumi == last_t.lumi) {
      //  if (t.l1 != last_t.l1 || t.hlt != last_t.hlt)
      //    throw cms::Exception("CheckPrescale") << "for trigger path " << *trigger_path << ", inside (run, lumi) = (" << t.run << ", " << t.lumi << ") going to event " << event.id().event() << " different L1/HLT prescales found: before " << last_t.l1 << "/" << last_t.hlt << "  now: " << t.l1 << "/" << t.hlt << "\n";
      //}
      //else
        tree->Fill();
      
      last_t = t;
    }
  }

  if (!found)
    throw cms::Exception("CheckPrescale") << "none of the trigger paths specified were found!\n";
}

DEFINE_FWK_MODULE(CheckPrescale);
