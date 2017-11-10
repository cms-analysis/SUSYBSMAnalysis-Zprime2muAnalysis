#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class TTbarSelection : public edm::EDFilter {
 public:
  explicit TTbarSelection(const edm::ParameterSet&);

 private:
 virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag src;
  const double min_mass;
  const double max_mass;
};


TTbarSelection::TTbarSelection(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),  
    min_mass(cfg.getParameter<double>("min_mass")),
    max_mass(cfg.getParameter<double>("max_mass"))
{
  consumes<reco::GenParticleCollection>(src);
 

}

bool TTbarSelection::filter(edm::Event& event, const edm::EventSetup&) {

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
  Double_t phi1 = -999;
  Double_t phi2 = -999;
  float massTT = 1.;
  int count = 0;


  

   for (; genp != genParticles->end(); genp++) {
   

 

  if (genp->pdgId() == 13 || genp->pdgId() == 11  || genp->pdgId() == 15){
    const reco::Candidate* m = genp->mother();

    if(m->pdgId()==-24){
    	const reco::Candidate* mm = m->mother();
    	if(mm->pdgId()==-6 || mm->pdgId()==-24){
	    	px1 = genp->px();
	    	py1 = genp->py();
		pz1 = genp->pz();
    		ener1 = sqrt(genp->mass()*genp->mass()+genp->px()*genp->px()+genp->py()*genp->py()+genp->pz()*genp->pz());
// 		    std::cout<<"lep Id ="<<genp->pdgId()<<"\tmuon status = "<<genp->status()<<"\tmuon mother ="<<m->pdgId()<<"\tmuon grandmother = "<<mm->pdgId()<<std::endl;
    	}
     }
  }

   if (genp->pdgId() == -13 || genp->pdgId() ==-11  || genp->pdgId() == -15){
     const reco::Candidate* m2 = genp->mother();
//      std::cout<<"lep Id ="<<genp->pdgId()<<"...muon 2 status = "<<genp->status()<<"..muon mother ="<<m2->pdgId()<<std::endl;

    if(m2->pdgId()==24){
    	const reco::Candidate* mm2 = m2->mother();
    	if(mm2->pdgId()==6 || mm2->pdgId()==24){
		    px2 = genp->px();
		    py2 = genp->py();
		    pz2 = genp->pz();
		    ener2 = sqrt(genp->mass()*genp->mass()+genp->px()*genp->px()+genp->py()*genp->py()+genp->pz()*genp->pz());
// 		    std::cout<<"lep Id ="<<genp->pdgId()<<"\tmuon status = "<<genp->status()<<"\tmuon mother ="<<m2->pdgId()<<"\tmuon grandmother = "<<mm2->pdgId()<<std::endl;
		}

  }
  }
  }

    massTT = sqrt((ener1+ener2)*(ener1+ener2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));

   return 
     massTT >= min_mass && massTT<=max_mass;
//      std::cout<<"Total event = "<<count<<std::endl;
 }
 

DEFINE_FWK_MODULE(TTbarSelection);
