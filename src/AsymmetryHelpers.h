#ifndef AsymmetryHelpers_h
#define AsymmetryHelpers_h

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFitData.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// This routine calculates the quantities to be used in the multiple
// dimension fit. Variables used in the fit are stored in the
// structure "data" returned by reference.  The flag "internalBremOn",
// is for switching on and off the effect of internal bremsstrahlung.
bool computeFitQuantities(const reco::GenParticleCollection& genParticles,
			  int leptonFlavor, bool internalBremOn,
			  AsymFitData& data);

#endif
