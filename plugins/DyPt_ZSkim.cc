#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"



class DyPt_ZSkim : public edm::EDFilter {
 public:
  explicit DyPt_ZSkim(const edm::ParameterSet&);

 private:
 virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag src;
  const double min_mass;
  const double max_mass;
};


DyPt_ZSkim::DyPt_ZSkim(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),  
    min_mass(cfg.getParameter<double>("min_mass")),
    max_mass(cfg.getParameter<double>("max_mass"))
{
  consumes<reco::GenParticleCollection>(src);
 

}

bool DyPt_ZSkim::filter(edm::Event& event, const edm::EventSetup&) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel(src, genParticles);

  reco::GenParticleCollection::const_iterator genp = genParticles->begin();
 
  float pt1 = 0;
  float pt2 = 0;
  float eta1 =0;
  float eta2 =0;
  float phi1 =0;
  float phi2 =0;
  TLorentzVector mu1, mu2;
  TLorentzVector Z;

  float ptW = 0.;


  

   for (; genp != genParticles->end(); genp++) {
   

 

  if (genp->pdgId() == 13){
    const reco::Candidate* m = genp->mother();
    //std::cout<<"lep Id ="<<genp->pdgId()<<"muon status = "<<genp->status()<<"...muon mother ="<<m->pdgId()<<std::endl;

    if(m->pdgId()==23){
	    eta1 = genp->eta();
    	phi1 = genp->phi();
	    pt1 = genp->pt();

		mu1.SetPtEtaPhiM(pt1, eta1, phi1, 0.105);
     }
  }

   if (genp->pdgId() == -13){
     const reco::Candidate* m2 = genp->mother();
     // std::cout<<"lep Id ="<<genp->pdgId()<<"...muon 2 status = "<<genp->status()<<"..muon mother ="<<m2->pdgId()<<std::endl;
     
     	if(m2->pdgId()==23){
	    	eta2 = genp->eta();
	    	phi2 = genp->phi();
		    pt2 = genp->pt();

			mu2.SetPtEtaPhiM(pt2, eta2, phi2, 0.105);  
		}
	}
	}

	Z = mu1+ mu2;

	ptW = Z.Pt();
	
    //std::cout<<"WW mass = "<<massWW<<std::endl;


   return 
     ptW >= min_mass && ptW <= max_mass;
 }
 

DEFINE_FWK_MODULE(DyPt_ZSkim);
