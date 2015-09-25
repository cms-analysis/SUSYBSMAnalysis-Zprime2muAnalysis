// JMTBAD make this fill lepton histos from all leptons + those that
// make it into dileptons always, not just depending on flag.

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <iostream>

class WeightAnalyzer : public edm::EDAnalyzer {
 public:
  explicit WeightAnalyzer(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void extractweight(const edm::Event&);

  private:
  //edm::InputTag GenEventInfo;
  double EventWeight;
  double mc_weight = 0;
  Int_t positive =0;
  Int_t negative =0;
  Int_t total =0;
  edm::InputTag genEventInfo_;

  };

WeightAnalyzer::WeightAnalyzer(const edm::ParameterSet& iConfig)
: genEventInfo_(iConfig.getParameter<edm::InputTag>("genEventInfo"))

{
   TH1::SetDefaultSumw2(true);
}

void WeightAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
 extractweight(event);
}
void WeightAnalyzer::extractweight(const edm::Event& event)
{
   
 edm::Handle<GenEventInfoProduct> gen_ev_info;
//   event.getByLabel(edm::InputTag("generator"), gen_ev_info);
 event.getByLabel(genEventInfo_, gen_ev_info);
 EventWeight = 1.;
 EventWeight = gen_ev_info->weight();
 mc_weight = ( EventWeight > 0 ) ? 1 : -1;
 std::cout<<"weight="<<mc_weight<<std::endl;
    if (mc_weight==1) positive++;
    else negative++;
    std::cout<<"negative weight="<<negative<<std::endl;
    std::cout<<"total weight = "<<positive+negative<<std::endl;
    std::cout<<"positive weight = "<<positive<<std::endl;
  
}


DEFINE_FWK_MODULE(WeightAnalyzer);
