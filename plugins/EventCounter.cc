// -*- C++ -*-
//
// Package:    Histograms
// Class:      EventCounter
// 
/**\class EventCounter EventCounter.cc brot/EventCounter/src/EventCounter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>

**/

// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <iostream> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


//ROOT
#include "TH1.h"
//
// class decleration
//

class EventCounter : public edm::EDAnalyzer {
public:
  explicit EventCounter(const edm::ParameterSet&);
  ~EventCounter();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::EDGetTokenT<GenEventInfoProduct> GenInfoTag_;

  //histos
  std::map<std::string, TH1F* > count_; 
  bool debug;
};

// constructors and destructor
EventCounter::EventCounter(const edm::ParameterSet& iConfig):
  GenInfoTag_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfoTag")))
{
  debug = false;
 edm::Service<TFileService> fs;
 std::string name = "Events";
 count_[name] = fs->make<TH1F>(name.c_str() , name.c_str() , 2 , 0.5 ,  2.5);
 count_[name]->GetXaxis()->SetBinLabel(1, "Events");
 name = "weights";
 count_[name] = fs->make<TH1F>(name.c_str() , name.c_str() , 2 , -1.5 ,  1.5);
 count_[name]->GetXaxis()->SetBinLabel(1, "weights");
 
}



EventCounter::~EventCounter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// member functions
// ------------ method called to for each event  ------------
void
EventCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // count every event
  count_["Events"]->Fill(1);	
  //get generator weigts
  edm::Handle<GenEventInfoProduct> genInfoProduct;
  iEvent.getByToken(GenInfoTag_, genInfoProduct);

  if (genInfoProduct.isValid()){
   	if ((*genInfoProduct).weight() < 0.0){
   	
   		count_["weights"]->Fill(-1);
   	}
   	else{
   		count_["weights"]->Fill(1);
   	}
  }
  else{
 
 	count_["weights"]->Fill(1); 
  
  }  



}

// ------------ method called once each job just before starting event loop  ------------
void 
EventCounter::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventCounter::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventCounter);
