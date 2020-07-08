#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TH1F.h"
#include "TFile.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


class GenMassHistos : public edm::EDAnalyzer {
 public:
  explicit GenMassHistos(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);


 private:

  edm::InputTag src;
  HardInteraction* hardInteraction;
  TH1F* GenMass;
};


GenMassHistos::GenMassHistos(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
    hardInteraction(new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")))
 
{
  consumes<reco::GenParticleCollection>(src);

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true); 

  GenMass            = fs->make<TH1F>("GenMass", "GenMass", 20000, 0, 20000);
}

void GenMassHistos::analyze(const edm::Event& event, const edm::EventSetup& setup) {

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
  //Double_t phi1 = -999;
  //Double_t phi2 = -999;
  float massTT = 1.;
  //int count = 0;


  

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

   hardInteraction->Fill(event);
   if(hardInteraction->IsValidForRes()) 
	{
		massTT = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).mass();
	}
   else massTT = sqrt((ener1+ener2)*(ener1+ener2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));
   GenMass->Fill(massTT);
 }
 

DEFINE_FWK_MODULE(GenMassHistos);
