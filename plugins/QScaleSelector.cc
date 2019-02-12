#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <iostream>
class QScaleSelector : public edm::EDFilter {
 public:
  explicit QScaleSelector(const edm::ParameterSet&);

 private:
 virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag src;
  const double min_mass;
  const double max_mass;
};


QScaleSelector::QScaleSelector(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),  
    min_mass(cfg.getParameter<double>("min_mass")),
    max_mass(cfg.getParameter<double>("max_mass"))
{
  consumes<GenEventInfoProduct>(src);
 

}

bool QScaleSelector::filter(edm::Event& event, const edm::EventSetup&) {


    edm::Handle<GenEventInfoProduct> genEvtInfo;
    event.getByLabel( src , genEvtInfo );
    double qScale = genEvtInfo->pdf()->scalePDF;  // in case of Pythia6, this will be pypars/pari(23)
   return qScale >= min_mass && qScale <=max_mass;
}
 

DEFINE_FWK_MODULE(QScaleSelector);
