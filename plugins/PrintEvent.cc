#include <boost/foreach.hpp>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Dumpers.h"

class PrintEvent : public edm::EDAnalyzer {
 public:
  explicit PrintEvent(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  bool dump_trigger_results;
  edm::InputTag trigger_results_src;
  bool dump_muons;
  edm::InputTag muon_src;
  bool dump_dileptons;
  edm::InputTag dilepton_src;
};

PrintEvent::PrintEvent(const edm::ParameterSet& cfg) {
  std::ostringstream out;
  out << "configuration:\n";

  dump_trigger_results = cfg.existsAs<edm::InputTag>("trigger_results_src");
  out << "dump_trigger_results: " << dump_trigger_results << " ";
  if (dump_trigger_results) {
    trigger_results_src = cfg.getParameter<edm::InputTag>("trigger_results_src");
    out << trigger_results_src;
  }
  out << "\n";

  dump_muons = cfg.existsAs<edm::InputTag>("muon_src");
  out << "dump_muons: " << dump_muons << " ";
  if (dump_muons) {
    muon_src = cfg.getParameter<edm::InputTag>("muon_src");
    out << muon_src;
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
  if (dump_trigger_results) {
    std::ostringstream out;
    edm::Handle<edm::TriggerResults> res;
    event.getByLabel(trigger_results_src, res);
    const edm::TriggerNames& names = event.triggerNames(*res);
    out << "TriggerResults: paths and results in event for " << trigger_results_src << ":\n";
    for (size_t i = 0; i < res->size(); ++i)
      out << "path #" << i << " name " << names.triggerName(i) << " fired? " << res->accept(i) << "\n";
    out << "\n";
    edm::LogInfo("PrintEvent") << out.str();
  }

  if (dump_muons) {
    std::ostringstream out;
    edm::Handle<reco::BeamSpot> beamSpot;
    event.getByLabel("offlineBeamSpot", beamSpot);
    reco::TrackBase::Point bs(beamSpot->x0(), beamSpot->y0(), beamSpot->z0());
    osprintf(out, "beamspot: %f %f %f\n", bs.x(), bs.y(), bs.z());
    
    edm::Handle<pat::MuonCollection> muons;
    event.getByLabel(muon_src, muons);
    out << "muons (size: " << muons->size() << "):\n";
    BOOST_FOREACH(const pat::Muon& mu, *muons) {
      out << mu << "\n";
    }
    out << "\n";
    edm::LogInfo("PrintEvent") << out.str();
  }

  if (dump_dileptons) {
    std::ostringstream out;

    edm::Handle<pat::CompositeCandidateCollection> dils;
    event.getByLabel(dilepton_src, dils);

    if (!dils.isValid())
      edm::LogInfo("PrintEvent") << "WARNING! tried to get dileptons using " << dilepton_src << " and failed!";
    else if (dils->size()) {
      edm::Handle<std::vector<reco::Vertex> > vtxs;
      event.getByLabel("offlinePrimaryVertices", vtxs);
      BOOST_FOREACH(const reco::Vertex& vtx, *vtxs)
	out << "vtx with z " << vtx.z() << "\n";
    
      out << "dileptons:\n";
      BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils)
	out << dil << "\n";

      out << "\n";
    
      edm::LogInfo("PrintEvent") << out.str();
    }
  }
}

DEFINE_FWK_MODULE(PrintEvent);
