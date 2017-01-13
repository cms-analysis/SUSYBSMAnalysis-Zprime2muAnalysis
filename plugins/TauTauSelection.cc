#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class TauTauSelection : public edm::EDFilter {
 public:
  explicit TauTauSelection(const edm::ParameterSet&);

 private:
 virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag src;
  
};


TauTauSelection::TauTauSelection(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src"))
  { 
  consumes<reco::GenParticleCollection>(src);
 

}

bool TauTauSelection::filter(edm::Event& event, const edm::EventSetup&) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel(src, genParticles);

  reco::GenParticleCollection::const_iterator genp = genParticles->begin();
 
  bool lepton1 = false;
  bool lepton2 = false;
  
   for (; genp != genParticles->end(); genp++) {
   

  if (genp->pdgId() == 15){
    const reco::Candidate* m = genp->mother();
    //std::cout<<"lep Id ="<<genp->pdgId()<<"muon status = "<<genp->status()<<"...muon mother ="<<m->pdgId()<<std::endl;

    if(m->pdgId()==23){
      //std::cout<<"lep 1 found"<<std::endl;
      lepton1=true;
  
     }
  }

   if (genp->pdgId()==-15){
     const reco::Candidate* m2 = genp->mother();
     //std::cout<<"lep Id ="<<genp->pdgId()<<"...muon 2 status = "<<genp->status()<<"..muon mother ="<<m2->pdgId()<<std::endl;

    if(m2->pdgId()==23){
      //  std::cout<<"lep 2 found"<<std::endl;
       lepton2=true;

  }
  }
  }

    

   return 
     lepton1==true && lepton2==true;
 }
 

DEFINE_FWK_MODULE(TauTauSelection);
