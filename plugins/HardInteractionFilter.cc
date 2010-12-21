#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

class HardInteractionFilter : public edm::EDFilter {
public:
  explicit HardInteractionFilter(const edm::ParameterSet&);

private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  HardInteraction hardInteraction;
  const double min_mass;
  const double max_mass;
  const bool use_resonance_mass;
  const double max_muon_eta;
  const double min_muon_pt;
};

HardInteractionFilter::HardInteractionFilter(const edm::ParameterSet& cfg)
  : hardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")),
    min_mass(cfg.getParameter<double>("min_mass")),
    max_mass(cfg.getParameter<double>("max_mass")),
    use_resonance_mass(cfg.getParameter<bool>("use_resonance_mass")),
    max_muon_eta(cfg.getParameter<double>("max_muon_eta")),
    min_muon_pt(cfg.getParameter<double>("min_muon_pt"))
{
}

bool HardInteractionFilter::filter(edm::Event& event, const edm::EventSetup&) {
  hardInteraction.Fill(event);
  const double m = use_resonance_mass ? hardInteraction.resonance->mass() : hardInteraction.dilepton().mass();
  return
    m > min_mass &&
    m < max_mass &&
    fabs(hardInteraction.lepMinus->eta()) < max_muon_eta && 
    fabs(hardInteraction.lepPlus ->eta()) < max_muon_eta &&
    hardInteraction.lepMinus->pt() > min_muon_pt &&
    hardInteraction.lepPlus ->pt() > min_muon_pt;
}

DEFINE_FWK_MODULE(HardInteractionFilter);
