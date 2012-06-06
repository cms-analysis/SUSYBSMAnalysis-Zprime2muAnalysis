// -*- C++ -*-
//
// Package:    GetPrescaleL1HLT
// Class:      GetPrescaleL1HLT
// 
/**\class GetPrescaleL1HLT GetPrescaleL1HLT.cc UserArea/GetPrescaleL1HLT/src/GetPrescaleL1HLT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Theodore Kypreos
//         Created:  Tue Jun  5 15:45:09 CDT 2012
// $Id$
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"


//
// class declaration
//

class GetPrescaleL1HLT : public edm::EDProducer {
   public:
      explicit GetPrescaleL1HLT(const edm::ParameterSet&);
      ~GetPrescaleL1HLT();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        int _l1prescale;
        int _hltprescale;
        int _prescaleWeight;
//        std::string hltpath;


   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      HLTConfigProvider hlt_cfg;
        edm::InputTag _tagHlt;
        std::string hlt_process_name;
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
GetPrescaleL1HLT::GetPrescaleL1HLT(const edm::ParameterSet& iConfig)
:_tagHlt(iConfig.getParameter<edm::InputTag>("hlt_src"))
{

    produces<int>("HLTPrescale");
    produces<int>("L1Prescale");
    produces<int>("TotalPrescale");
    hlt_process_name= _tagHlt.process();
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


GetPrescaleL1HLT::~GetPrescaleL1HLT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GetPrescaleL1HLT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{



  bool found = false;
  std::string trigger_path = "HLT_Mu24_eta2p1_v3";
  unsigned int path_index = hlt_cfg.size();
  
    unsigned int ndx = hlt_cfg.triggerIndex(trigger_path);
    if (ndx < hlt_cfg.size()) {  
        found = true;
        path_index = ndx;
//    trigger_path = *path;
    }

//    if (!found){
//        std::cout<<"you don't seem to have this trigger. I'm getting out of here while I can..."<<std::endl;
//        return;
//    }

    edm::Handle<edm::TriggerResults> hlt_results;
    iEvent.getByLabel(_tagHlt,hlt_results);
    std::pair<int, int> prescales(1,1);

    if (found) {
        if (hlt_results->accept(path_index)) {
//            std::cout<<"path is accepted"<<std::endl;
            if (iEvent.isRealData()) prescales = hlt_cfg.prescaleValues(iEvent, iSetup, trigger_path);
          }
    }

   std::auto_ptr<int> l1Prescale(new int(prescales.first));
   std::auto_ptr<int> hltPrescale(new int(prescales.second));
   std::auto_ptr<int> totalPrescale(new int(prescales.first*prescales.second));


  iEvent.put(hltPrescale, "HLTPrescale");
  iEvent.put(l1Prescale, "L1Prescale");
  iEvent.put(totalPrescale, "TotalPrescale");

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
GetPrescaleL1HLT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GetPrescaleL1HLT::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
GetPrescaleL1HLT::beginRun(edm::Run& run, edm::EventSetup const& setup)
{
  bool changed = true;
  if (!hlt_cfg.init(run, setup, hlt_process_name, changed))
    throw cms::Exception("PrescaleToCommon") << "HLTConfigProvider::init failed with process name " << hlt_process_name << "\n";
  
//  return true;
}

// ------------ method called when ending the processing of a run  ------------
void 
GetPrescaleL1HLT::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GetPrescaleL1HLT::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GetPrescaleL1HLT::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GetPrescaleL1HLT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GetPrescaleL1HLT);
