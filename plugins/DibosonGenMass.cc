#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class DibosonGenMass : public edm::EDFilter {
 public:
  explicit DibosonGenMass(const edm::ParameterSet&);

 private:
 virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag src;
  const double min_mass;
  const double max_mass;
};


DibosonGenMass::DibosonGenMass(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),  
    min_mass(cfg.getParameter<double>("min_mass")),
    max_mass(cfg.getParameter<double>("max_mass"))
{
  consumes<reco::GenParticleCollection>(src);
 

}

bool DibosonGenMass::filter(edm::Event& event, const edm::EventSetup&) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel(src, genParticles);

  reco::GenParticleCollection::const_iterator genp = genParticles->begin();
 
  float ener1 = 0;
  float ener2 = 0;
  float px1 =0;
  float px2 =0;
  float py1 =0;
  float py2 =0;
  float pz1 =0;
  float pz2 =0;
  float massWW = 1.;


  

   for (; genp != genParticles->end(); genp++) {
   

 

  if (genp->pdgId() == 13 || genp->pdgId() == 11  || genp->pdgId() == 15){
    const reco::Candidate* m = genp->mother();
    //std::cout<<"lep Id ="<<genp->pdgId()<<"muon status = "<<genp->status()<<"...muon mother ="<<m->pdgId()<<std::endl;

    if(m->pdgId()==-24){
    px1 = genp->px();
    py1 = genp->py();
    pz1 = genp->pz();
    ener1 = sqrt(genp->mass()*genp->mass()+genp->px()*genp->px()+genp->py()*genp->py()+genp->pz()*genp->pz());
     }
  }

   if (genp->pdgId() == -13 || genp->pdgId() ==-11  || genp->pdgId() == -15){
     const reco::Candidate* m2 = genp->mother();
     // std::cout<<"lep Id ="<<genp->pdgId()<<"...muon 2 status = "<<genp->status()<<"..muon mother ="<<m2->pdgId()<<std::endl;

    if(m2->pdgId()==24){
    px2 = genp->px();
    py2 = genp->py();
    pz2 = genp->pz();
    ener2 = sqrt(genp->mass()*genp->mass()+genp->px()*genp->px()+genp->py()*genp->py()+genp->pz()*genp->pz());

  }
  }
  }

    massWW = sqrt((ener1+ener2)*(ener1+ener2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));

    //std::cout<<"WW mass = "<<massWW<<std::endl;


   return 
     massWW >= min_mass && massWW<=max_mass;
 }
 

DEFINE_FWK_MODULE(DibosonGenMass);
