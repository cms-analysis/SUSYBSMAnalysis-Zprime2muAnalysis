#include <boost/foreach.hpp>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Dumpers.h"

class PrintEvent : public edm::EDAnalyzer {
 public:
  explicit PrintEvent(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  bool dump_trigger_names;
  edm::InputTag trigger_results_src;
  bool dump_dileptons;
  edm::InputTag dilepton_src;
};

PrintEvent::PrintEvent(const edm::ParameterSet& cfg) {
  std::ostringstream out;
  out << "configuration:\n";

  dump_trigger_names = cfg.existsAs<edm::InputTag>("trigger_results_src");
  out << "dump_trigger_names: " << dump_trigger_names << " ";
  if (dump_trigger_names) {
    trigger_results_src = cfg.getParameter<edm::InputTag>("trigger_results_src");
    out << trigger_results_src;
  }
  out << "\n";

  dump_dileptons = cfg.existsAs<edm::InputTag>("dilepton_src");
  out << "dump_dileptons: " << dump_dileptons << " ";
  if (dump_dileptons) {
    dilepton_src = cfg.getParameter<edm::InputTag>("dilepton_src");
    out << dilepton_src;
  }
  out << "\n";

  edm::LogInfo("PrintEvent") << out.str();
}

void PrintEvent::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  if (dump_trigger_names) {
    std::ostringstream out;
    edm::Handle<edm::TriggerResults> res;
    event.getByLabel(trigger_results_src, res);
    out << "HLT paths in event:\n" << event.triggerNames(*res) << "\n";
    edm::LogInfo("PrintEvent") << out.str();
  }

  if (dump_dileptons) {
    std::ostringstream out;

    edm::Handle<pat::CompositeCandidateCollection> dils;
    event.getByLabel(dilepton_src, dils);

    if (!dils.isValid())
      edm::LogInfo("PrintEvent") << "WARNING! tried to get dileptons using " << dilepton_src << " and failed!";
    else if (dils->size()) {
      edm::Handle<edm::TriggerResults> res;
      event.getByLabel(edm::InputTag("TriggerResults", "", "PAT"), res);
      const edm::TriggerNames& names = event.triggerNames(*res);
      for (size_t i = 0; i < res->size(); ++i)
	out << "TriggerResults::PAT path #" << i << " name " << names.triggerName(i) << " fired? " << res->accept(i) << "\n";
    
      edm::Handle<std::vector<reco::Vertex> > vtxs;
      event.getByLabel("offlinePrimaryVertices", vtxs);
      BOOST_FOREACH(const reco::Vertex& vtx, *vtxs)
	out << "vtx with z " << vtx.z() << "\n";
    
      out << "dileptons:\n";
      BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils)
	out << dil << "\n";
    
      edm::LogInfo("PrintEvent") << out.str();
    }
  }
}

DEFINE_FWK_MODULE(PrintEvent);
