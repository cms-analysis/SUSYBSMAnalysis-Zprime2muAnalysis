#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GenEventTopology.h"

unsigned GenEventTopology(const edm::Event& event, int pdgId[2], reco::Particle::LorentzVector p4[2]) {
  unsigned status = 0;

  pdgId[0] = pdgId[1] = 0;
  p4[0] = reco::Particle::LorentzVector();
  p4[1] = reco::Particle::LorentzVector();

  const reco::Candidate* ts[2] = {0, 0};
  const reco::Candidate* Ws[2] = {0, 0};
  const reco::Candidate* Zs[2] = {0, 0};

  // Try to find t and tbar, and the Ws.
  edm::Handle<reco::CandidateCollection> genParticles;
  event.getByLabel("genParticleCandidates", genParticles);

  for (reco::CandidateCollection::const_iterator genp = genParticles->begin(); genp != genParticles->end(); genp++) {
    if (genp->status() != 3) continue;

    int id = genp->pdgId();
    if      (id ==   6) ts[0] = &*genp;
    else if (id ==  -6) ts[1] = &*genp;	
    else if (id ==  24) Ws[0] = &*genp;
    else if (id == -24) Ws[1] = &*genp;	
  }

  if (ts[0] != 0 && ts[1] != 0) status |= topo_ttbar;
  if (Ws[0] != 0) status |= topo_Wp;
  if (Ws[1] != 0) status |= topo_Wm;

  // Assume the t and tbar went into W+ and W-, respectively.  Try
  // to find the W+ and W- from the t and the tbar.
  //if (ts_found && !Ws_found)
  //  for (unsigned i = 0; i < 2; i++)
  //    for (unsigned j = 0; j < 2; j++)
  //      if (abs(ts[i]->daughter(j)->pdgId()) == 24)
  //        Ws[i] = ts[i]->daughter(j);
    
  // Encode what the Ws went into.
  for (unsigned i = 0; i < 2; i++) {
    if (Ws[i] == 0) continue;
    
    int id  = abs(Ws[i]->daughter(0)->pdgId());
    
    if (id >= 11 && id <= 16) {
      // W -> l nu. Figure out which was the nu and which was not.
      unsigned not_nu = ((id - 11) % 2 == 0) ? 0 : 1;
      const reco::Candidate* lep = Ws[i]->daughter(not_nu);
      
      pdgId[i] = lep->pdgId();
      p4[i] = lep->p4();
    }
  }

  // Look for a DY decay too.
  unsigned Zcount = 0;
  for (reco::CandidateCollection::const_iterator genp = genParticles->begin(); genp != genParticles->end(); genp++) {
    if (genp->status() == 3 && abs(genp->pdgId()) == 23) {
      Zs[Zcount++] = &*genp;

      int id = abs(genp->daughter(0)->pdgId());

      if      (id == 11) status |= topo_Zee;
      else if (id == 13) status |= topo_Zmm;
      else if (id == 15) status |= topo_Ztt;
            
      if (Zcount == 2) break;
    }
  }

  if (Zcount > 0) status |= topo_Z;
  if (Zcount > 1) status |= topo_Z2;

  edm::LogInfo("GenEventTopology") << "pointers: t: " << ts[0] << " tbar: " << ts[1]
				   << " W+: " << Ws[0] << " W-: " << Ws[1]
				   << " Z1: " << Zs[0] << " Z2: " << Zs[1];
  
  return status;
}
