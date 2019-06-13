#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/LRWeightProducer.h"
#include <TLorentzVector.h>

LRWeightProducer::LRWeightProducer(const edm::ParameterSet cfg)
  : doingElectrons(cfg.getParameter<bool>("doingElectrons")),
    lambda(cfg.getParameter<int>("Lambda")),
    interference(cfg.getParameter<int>("interference")),
    doingLR(cfg.getParameter<bool>("doingLR")),
    calculate(cfg.getParameter<bool>("calculate")),
    leptonFlavor(doingElectrons ? 11 : 13),
    leptonMass(doingElectrons ? 0.000511 : 0.10566)
{
}


LRWeightProducer::~LRWeightProducer() {
}




double LRWeightProducer::calculateWeight(const edm::Event& event, HardInteraction*& hardInteraction, double alpEM) {
  if (!calculate) return 1.;
//  //Begin calculating fracLR and fracRL                                     
  std::complex<double> I(0.0, 1.0);
 // Complex amplitudes.                                                    
  std::complex<double> meLL(0., 0.);
  std::complex<double> meRR(0., 0.);
  std::complex<double> meLR(0., 0.);
  std::complex<double> meRL(0., 0.);
  std::complex<double> meLR_SM(0., 0.);
  std::complex<double> meRL_SM(0., 0.);
 if (!hardInteraction->IsValidForRes()) return 1.;
 int quarkId=hardInteraction->quark->pdgId();
 int idAbs=quarkId;

  int leptonMinusId=hardInteraction->lepMinusNoIB->pdgId();
  int idNew=leptonMinusId;

  double qCLambda2 = lambda*lambda;
  int qCetaLR = interference;

  double sH = (hardInteraction->quark->p4() + hardInteraction->antiquark->p4()).mag2();
  double tH = (hardInteraction->quark->p4() - hardInteraction->lepMinusNoIB->p4()).mag2();
  double uH = (hardInteraction->quark->p4() - hardInteraction->lepPlusNoIB->p4()).mag2();



  //double sH2 = sH*sH;
  double tH2 = tH*tH;
  double uH2 = uH*uH;


  double tmPgvf = 0.25 * af[idAbs] - 4. * 0.2223 * ef[idAbs];
  double tmPgaf = 0.25 * af[idAbs];

  double tmPgvl = 0.25 * af[idNew] - 4. * 0.2223 * ef[idNew];
  double tmPgal = 0.25 * af[idNew];
  double tmPgLf = tmPgvf + tmPgaf;
  double tmPgLl = tmPgvl + tmPgal;
  double tmPgRf = tmPgvf - tmPgaf;
  double tmPgRl = tmPgvl - tmPgal;
// Kinematics  
//double qCmNew  = fMasterGen.particleData.m0(idNew);
 //double qCmNew2 = qCmNew * qCmNew;  
  double qCmZ    = 91.1876;
  double qCmZ2   = qCmZ * qCmZ;
  double qCGZ    = 2.4952;
  double qCGZ2 = qCGZ * qCGZ;
  // Necessary variables to ampitude                                        	
 // First term 
 // double alpEM =1./137; 
  double tmPe2QfQl = 4. * M_PI * alpEM * ef[idAbs] * ef[idNew];
  double qCPropGm = 1./sH;
  //Second term.Model depended variables are defined using incoming quark and outgoing fermion information                                                 
     double tmPe2s2c2 = 4. * M_PI * alpEM;
     double denomPropZ = pow((sH - qCmZ2), 2) + qCmZ2 * qCGZ2;
     double qCrePropZ  = (sH - qCmZ2) / denomPropZ;
     double qCimPropZ = -qCmZ * qCGZ / denomPropZ;
 //Third term:4. * M_PI * qCetaLR / qCLambda2;    
 // Amplitudes, M = gamma + Z + CI.
      meLL = tmPe2QfQl * qCPropGm
	+ tmPe2s2c2 * tmPgLf * tmPgLl * (qCrePropZ + I * qCimPropZ);
      meRR = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgRf * tmPgRl * (qCrePropZ + I * qCimPropZ);
      meLR = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgLf * tmPgRl * (qCrePropZ + I * qCimPropZ)
        + 4. * M_PI * qCetaLR / qCLambda2;
      meRL = tmPe2QfQl * qCPropGm
	+ tmPe2s2c2 * tmPgRf * tmPgLl * (qCrePropZ + I * qCimPropZ)
+ 4. * M_PI * qCetaLR / qCLambda2;

      // According to Steve's idea, add SM matrix elements for sigmaNew.        
      // // Define standard model matrix elements of RL and LR model   
      meLR_SM = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgLf * tmPgRl * (qCrePropZ + I * qCimPropZ);

      meRL_SM = tmPe2QfQl * qCPropGm
+ tmPe2s2c2 * tmPgRf * tmPgLl * (qCrePropZ + I * qCimPropZ);
 // Calculate weighting facror
      double sigma0 = 1.0;
      double sigmaOld = sigma0 * uH2 * std::real(meLL*std::conj(meLL));
      sigmaOld += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
      sigmaOld += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
      sigmaOld += sigma0 * tH2 * std::real(meRL*std::conj(meRL));

      double sigmaNewLR = sigma0 * uH2 *std:: real(meLL*std::conj(meLL));
      sigmaNewLR += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
      sigmaNewLR += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
// sigma += sigma0 * tH2 * std::real(meRL*std::conj(meRL)); 
      sigmaNewLR += sigma0 * tH2 * std::real(meRL_SM *std::conj(meRL_SM));
      double fracLR = sigmaNewLR / sigmaOld;

      double sigmaNewRL = sigma0 * uH2 *std:: real(meLL*std::conj(meLL));
      sigmaNewRL += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
//sigmaNew += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
      sigmaNewRL += sigma0 * tH2 * std::real(meRL*std::conj(meRL));
      sigmaNewRL += sigma0 * tH2 * std::real(meLR_SM*std::conj(meLR_SM));
      double fracRL = sigmaNewRL / sigmaOld;

      if (doingLR) return fracLR;
      else return fracRL; 
}
