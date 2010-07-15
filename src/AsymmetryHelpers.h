#ifndef AsymmetryHelpers_h
#define AsymmetryHelpers_h

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFitData.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// This routine calculates the quantities to be used in the multiple
// dimension fit. Variables used in the fit are stored in the
// structure "data" returned by reference.  The flag "internalBremOn",
// is for switching on and off the effect of internal bremsstrahlung.
bool computeFitQuantities(const reco::GenParticleCollection& genParticles,
			  const bool doingElectrons, const bool internalBremOn,
			  AsymFitData& data);

void calcAsymmetry(double f, double b, double& A_FB, double& e_A_FB);
void calcAsymmetry(const TH1F* h_cos, double& A_FB, double& e_A_FB);
void calcAsymmetry(const TH1F* IdF, const TH1F* IdB, double& A_FB, double& e_A_FB);
void calcAsymmetry(const TH1F* IdF, const TH1F* IdB, TH1F* IdA, double& A_FB, double& e_A_FB);

void fitCosTheta(std::ostream& out, TH1F* h_cos, int num_params=3, bool fix_norm=false, bool draw=true);

#endif
