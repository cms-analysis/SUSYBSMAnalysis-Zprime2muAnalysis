#ifndef Zp2mu_TriggerUtilities_h
#define Zp2mu_TriggerUtilities_h

#include "DataFormats/HLTReco/interface/TriggerObject.h"

namespace edm { class Event; class InputTag; }

struct Zprime2muTriggerPathsAndFilters {
  std::string path;
  std::string filter;
  std::string path_2;
  std::string filter_2;
  std::string prescaled_path;
  std::string prescaled_filter;
  std::string prescaled_path_2;
  std::string prescaled_filter_2;
  bool valid;
  
  Zprime2muTriggerPathsAndFilters() : valid(false) {}
  Zprime2muTriggerPathsAndFilters(const edm::Event&);
};

trigger::TriggerObjectCollection get_L3_muons(const edm::Event& event,
					      const std::string& filter, 
					      const edm::InputTag& trigger_summary_src = edm::InputTag("hltTriggerSummaryAOD", "", "HLT"),
					      const std::string& collection="hltL3MuonCandidates");

#endif

