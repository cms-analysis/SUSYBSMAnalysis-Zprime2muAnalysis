#include <boost/foreach.hpp>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
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
  bool dump_met;
  std::vector<edm::InputTag> met_srcs;
  bool dump_jets;
  edm::InputTag jet_src;
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

  dump_met = cfg.existsAs<std::vector<edm::InputTag> >("met_srcs");
  out << "dump_met: " << dump_met << " ";
  if (dump_met) {
    met_srcs = cfg.getParameter<std::vector<edm::InputTag> >("met_srcs");
    BOOST_FOREACH(const edm::InputTag& met_src, met_srcs)
      out << met_src.encode() << " ";
  }
  out << "\n";

  dump_jets = cfg.existsAs<edm::InputTag>("jet_src");
  out << "dump_jets: " << dump_jets << " ";
  if (dump_jets) {
    jet_src = cfg.getParameter<edm::InputTag>("jet_src");
    out << jet_src;
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

  if (dump_met) {
    std::ostringstream out;
    out << "mets:\n";
    BOOST_FOREACH(const edm::InputTag& met_src, met_srcs) {
      edm::Handle<pat::METCollection> met;
      event.getByLabel(met_src, met);
      out << met_src.encode() << ": et: " << met->front().et() << " phi: " << met->front().phi() << "\n";
    }
    edm::LogInfo("PrintEvent") << out.str();
  }

  if (dump_jets) {
    std::ostringstream out;
    edm::Handle<pat::JetCollection> jets;
    event.getByLabel(jet_src, jets);
    out << "jets (size: " << jets->size() << "):\n";
    BOOST_FOREACH(const pat::Jet& jet, *jets)
      out << "p4: " << jet.p4() << "\n";
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
    BOOST_FOREACH(const pat::Muon& mu, *muons)
      ::operator<<(out, mu) << "\n"; // *vomits uncontrollably*
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
      out << "Offline primary vertices:\n";
      BOOST_FOREACH(const reco::Vertex& vtx, *vtxs)
	out << "vtx with x, y, z: " << vtx.x() << " " << vtx.y() << " " << vtx.z() << "\n";
    
      out << "\n Dileptons:\n";
      BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils) {
	out << dil << "\n";
	if (dil.hasUserFloat("vertexX"))
	  out << "Common dimuon vertex: (x, y, z): "
	      << dil.userFloat("vertexX") << ", " << dil.userFloat("vertexY") << ", " << dil.userFloat("vertexZ")
	      << " chi2/dof: " << dil.userFloat("vertex_chi2") << "\n";
	if (dil.hasUserFloat("vertexM"))
	  out << "Mass computed with the common-vertex constraint: "
	      << dil.userFloat("vertexM") << " +/- " << dil.userFloat("vertexMError") << "\n";
	if (dil.hasUserFloat("cos_angle"))
	  out << "Cos of angle between muons: " << dil.userFloat("cos_angle")
	      << " largest dpT/pT: " << dil.userFloat("dpt_over_pt") << "\n";
      }

      out << "\n";
    
      edm::LogInfo("PrintEvent") << out.str();
    }
  }
}

DEFINE_FWK_MODULE(PrintEvent);
