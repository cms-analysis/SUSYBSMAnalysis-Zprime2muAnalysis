#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

using namespace std;
using namespace edm;
using namespace reco;

class DimuonBremRecoverer : public EDProducer {
public:
  explicit DimuonBremRecoverer(const ParameterSet&);
  ~DimuonBremRecoverer() {}

private:
  virtual void produce(Event&, const EventSetup&);

  bool debug;

  // Which dimuons and photons to get from the event.
  InputTag dimuons;
  InputTag photonMatchMap;

  // Only add in photons that are within this dR of the daughter muons.
  double dRmax;

  typedef Particle::LorentzVector LorentzVector;
};

DimuonBremRecoverer::DimuonBremRecoverer(const ParameterSet& cfg)
  : debug(cfg.getUntrackedParameter<int>("verbosity", 0) > 0),
    dimuons(cfg.getParameter<InputTag>("dimuons")),
    photonMatchMap(cfg.getParameter<InputTag>("photonMatchMap")),
    dRmax(cfg.getParameter<double>("dRmax"))
{
  produces<CompositeCandidateCollection>();
}

void DimuonBremRecoverer::produce(Event& event,
				  const EventSetup& eSetup) {
  ostringstream out;

  // Get the existing dimuons.
  Handle<CompositeCandidateCollection> dimus;
  event.getByLabel(dimuons, dimus);

  // Get the photon match map for this event.
  Handle<CandViewMatchMap> photons;
  event.getByLabel(photonMatchMap, photons);

  // Set up the output collection.
  auto_ptr<CompositeCandidateCollection>
    resonances(new CompositeCandidateCollection);

  if (dimus.failedToGet()) {
    edm::LogWarning("DimuonBremRecoverer")
      << "Dimuons " << dimuons << " not found in event;"
      << " producing empty collection.";
  }
  else {
    // If the photon match map isn't found, just copy the dimuons.
    if (photons.failedToGet()) {
      edm::LogWarning("DimuonBremRecoverer")
	<< "Photon match map " << photonMatchMap << " not found in event;"
	<< " copying existing dimuons.";
      CompositeCandidateCollection::const_iterator dimu;
      for (dimu = dimus->begin(); dimu != dimus->end(); dimu++)
	resonances->push_back(*dimu);
    }
    else {
      CompositeCandidateCollection::const_iterator dimu;
      unsigned idi = 0;
      for (dimu = dimus->begin(); dimu != dimus->end(); dimu++, idi++) {
	if (debug)
	  out << " Dilepton #" << idi << " p4 = " << dimu->p4() << endl;

	// Keep track of the photons we added already.
	vector<LorentzVector> photonP4s;

	// Make a new dimuon.
	CompositeCandidate newDimu;

	// First, add in all the dimuon's daughters (of which hopefully
	// there are exactly two) so that the new dimuon's daughter muons
	// are in the same positions as for the old dimuon.
	for (unsigned idau = 0; idau < dimu->numberOfDaughters(); idau++) {
	  // Make sure we save the master clone of the muon into the new
	  // dimuon so we can still use it in all our match maps, etc.
	  const CandidateBaseRef& dau = dimu->daughter(idau)->masterClone();
	  newDimu.addDaughter(ShallowCloneCandidate(dau));
	}

	// Now loop over the daughter muons again, and add their closest
	// photons to the new dimuon.
	for (unsigned idau = 0; idau < dimu->numberOfDaughters(); idau++) {
	  const CandidateBaseRef& dau = dimu->daughter(idau)->masterClone();

	  // Skip non-muons.
	  if (abs(dau->pdgId()) != 13) continue;
    
	  if (photons->find(dau) != photons->end()) {
	    const CandidateBaseRef& photon = (*photons)[dau];
	    LorentzVector photonP4 = photon->p4();
	    if (debug) out << "  Considering photon " << photonP4
			   << " from daughter " << dau->p4();
	
	    if (photonP4.energy() > 0) {
	      double dR = deltaR(dau->p4(), photonP4);
	      if (debug) out << "; its dR = " << dR;

	      if (dR < dRmax) {
		bool addedAlready = false;
		// Try not to add the same photon twice.
		for (unsigned iph = 0; iph < photonP4s.size(); iph++) {
		  const LorentzVector diff = photonP4s[iph] - photonP4;
		  if (diff.P() < 0.001) {
		    addedAlready = true;
		    break;
		  }
		}
	    
		if (!addedAlready) {
		  if (debug) out << "; adding it";
		  photonP4s.push_back(photonP4);
		  newDimu.addDaughter(ShallowCloneCandidate(photon));
		}
		else {
		  if (debug) out << "; already added, skipping it";
		}
	      }
	    }
	    else {
	      if (debug) out << "; E=0, not considering\n";
	    }

	    if (debug) out << endl;
	  }	      
	}

	AddFourMomenta addP4;
	addP4.set(newDimu);

	resonances->push_back(newDimu);

	if (debug) 
	  out << " Include brem candidate(s): inv. mass w/o photons = "
	      << dimu->mass() << "; w/  = " << newDimu.mass();
	
	if (debug) LogDebug("DimuonBremRecoverer") << out.str();
      }
    }
  }
  
  event.put(resonances);
}

DEFINE_FWK_MODULE(DimuonBremRecoverer);
