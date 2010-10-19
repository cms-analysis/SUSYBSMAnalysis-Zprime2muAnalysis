#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

class MuonsFromDimuons : public edm::EDProducer {
public:
  explicit MuonsFromDimuons(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag dimuon_src;
};

MuonsFromDimuons::MuonsFromDimuons(const edm::ParameterSet& cfg)
  : dimuon_src(cfg.getParameter<edm::InputTag>("dimuon_src"))
{
  produces<pat::MuonCollection>();
}

void MuonsFromDimuons::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  event.getByLabel(dimuon_src, dimuons);

  std::auto_ptr<pat::MuonCollection> muons(new pat::MuonCollection);

  for (pat::CompositeCandidateCollection::const_iterator di = dimuons->begin(), die = dimuons->end(); di != die; ++di) {
    muons->push_back(toConcrete<pat::Muon>(dileptonDaughter(*di, 0)));
    muons->push_back(toConcrete<pat::Muon>(dileptonDaughter(*di, 1)));
  }

  event.put(muons);
}

DEFINE_FWK_MODULE(MuonsFromDimuons);
