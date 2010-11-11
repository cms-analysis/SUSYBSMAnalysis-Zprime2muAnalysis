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
  edm::InputTag dimuon_src;
  edm::InputTag hlt_src;
};

PrintEvent::PrintEvent(const edm::ParameterSet& cfg)
  : dimuon_src(cfg.getParameter<edm::InputTag>("dimuon_src")),
    hlt_src(cfg.getParameter<edm::InputTag>("hlt_src"))
{
}

void PrintEvent::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  std::ostringstream out;

  edm::Handle<edm::TriggerResults> res;
  event.getByLabel(hlt_src, res);
  out << "HLT paths in event:\n" << event.triggerNames(*res) << "\n";

  edm::LogInfo("PrintEvent") << out.str();
  out.str("");

  edm::Handle<pat::CompositeCandidateCollection> dils;
  event.getByLabel(dimuon_src, dils);

  if (!dils.isValid())
    edm::LogInfo("PrintEvent") << "WARNING! tried to get dils and failed!";
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

DEFINE_FWK_MODULE(PrintEvent);
