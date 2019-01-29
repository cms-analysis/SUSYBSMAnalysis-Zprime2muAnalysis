// -*- C++ -*-
//
// Package:    DielectronPreselector
// Class:      DielectronPreselector
// 
/**\class DielectronPreselector DielectronPreselector.cc SUSYBSMAnalysis/Zprime2muAnalysis/DielectronPreselector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/



// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//
// class declaration
//

class DielectronPreselector : public edm::EDFilter {
public:
  explicit DielectronPreselector(const edm::ParameterSet&);
  ~DielectronPreselector();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::EDGetTokenT< std::vector< pat::Electron > >  electronToken_;
  int ptCut_;
  int multiplicityCut_;


  bool debug;
};

// constructors and destructor
DielectronPreselector::DielectronPreselector(const edm::ParameterSet& iConfig):
 electronToken_ (consumes< std::vector< pat::Electron > > (iConfig.getParameter<edm::InputTag>("electrons")))
{
  ptCut_    = iConfig.getParameter < double > ("ptCut");
  multiplicityCut_    = iConfig.getParameter < double > ("nElectrons");

  debug = false;
}

DielectronPreselector::~DielectronPreselector(){}


// member functions
// ------------ method called on each new Event  ------------
bool
DielectronPreselector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector< pat::Electron > > electrons;  
  iEvent.getByToken(electronToken_, electrons);
  bool filter = false;

  int nElectrons = 0;

  for (std::vector<pat::Electron>::const_iterator ele = electrons->begin(); ele != electrons->end(); ele++){

	if ((*ele).pt() > ptCut_) nElectrons++;


  }

  if (nElectrons >= multiplicityCut_) filter = true;
  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DielectronPreselector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DielectronPreselector::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DielectronPreselector);
