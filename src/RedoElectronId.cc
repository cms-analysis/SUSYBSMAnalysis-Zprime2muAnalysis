#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RedoElectronId.h"

// Stolen from EgammaAnalysis/ElectronIDAlgos to make it more convenient to run from code.

using namespace reco;
using namespace edm;
using namespace std;

// cut value arrays of form {hoe, sigmaEtaEta, dPhiIn, dEtaIn}

const double robustEleIDCutsBarrel[4] = {0.115, 0.0140, 0.090, 0.0090};
const double robustEleIDCutsEndcap[4] = {0.150, 0.0275, 0.092, 0.0105};

// cut value arrays of form {barrel cat 0, barrel cat 1, barrel cat 2, barrel cat high E/p,
//                           endcap cat 0, endcap cat 1, endcap cat 2, endcap cat high E/p}

const double looseEleIDCuts_hOverE[8]       = {0.115,   0.10,   0.055,  0.,  
					       0.145,   0.12,   0.150,  0.};
const double looseEleIDCuts_sigmaEtaEta[8]  = {0.0140,  0.0120, 0.0115, 0.,
					       0.0275,  0.0265, 0.0265, 0.};
const double looseEleIDCuts_deltaPhiIn[8]   = {0.05,    0.025,  0.053,  0.09, 
					       0.07,    0.03,   0.092,  0.092};
const double looseEleIDCuts_deltaEtaIn[8]   = {0.009,   0.0045, 0.0085, 0.,
					       0.0105,  0.0068, 0.010,  0.};
const double looseEleIDCuts_eSeedOverPin[8] = {0.11,    0.91,   0.11,   0.,   
					       0.,      0.85,   0.,     0.};

const double tightEleIDCuts_hOverE[8]          = {0.05,    0.042,  0.045,  0.,  
						  0.055,   0.037,  0.050,  0.};
const double tightEleIDCuts_sigmaEtaEta[8]     = {0.0125,  0.011,  0.01,   0.,
						  0.0265,  0.0252, 0.026,  0.};
const double tightEleIDCuts_deltaPhiIn[8]      = {0.032,   0.016,  0.0525, 0.09, 
						  0.025,   0.035,  0.065,  0.092};
const double tightEleIDCuts_deltaEtaIn[8]      = {0.0055,  0.0030, 0.0065, 0.,
						  0.0060,  0.0055, 0.0075, 0.};
const double tightEleIDCuts_eSeedOverPinMin[8] = {0.24,    0.94,   0.11,   0.,   
						  0.32,    0.83,   0.,     0.};
const double tightEleIDCuts_eSeedOverPinMax[8] = {99999.,  99999., 99999., 99999.,   
						  99999.,  99999., 99999., 99999.};

const ClusterShapeRef& getClusterShape(const GsfElectron* electron, const Event& e) {
  InputTag barrelClusterShapeAssocProducer_("hybridSuperClusters:hybridShapeAssoc");
  InputTag endcapClusterShapeAssocProducer_("islandBasicClusters:islandEndcapShapeAssoc");

  // Get association maps linking BasicClusters to ClusterShape.
  Handle<BasicClusterShapeAssociationCollection> clusterShapeHandleBarrel;
  Handle<BasicClusterShapeAssociationCollection> clusterShapeHandleEndcap;
  e.getByLabel(barrelClusterShapeAssocProducer_, clusterShapeHandleBarrel);
  e.getByLabel(endcapClusterShapeAssocProducer_, clusterShapeHandleEndcap);

  // Find entry in map corresponding to seed BasicCluster of SuperCluster
  BasicClusterShapeAssociationCollection::const_iterator seedShpItr;

  if (electron->classification()<100) {
    seedShpItr=clusterShapeHandleBarrel->find(electron->superCluster()->seed());
    if (electron->classification()==40 && seedShpItr == clusterShapeHandleBarrel->end()) 
        seedShpItr=clusterShapeHandleEndcap->find(electron->superCluster()->seed());
  } else {
    seedShpItr=clusterShapeHandleEndcap->find(electron->superCluster()->seed());
  }

  return seedShpItr->val;
}

int classify(const GsfElectron* electron) {
  
  double eta = electron->p4().Eta();
  double eOverP = electron->eSuperClusterOverP();
  double pin  = electron->trackMomentumAtVtx().R(); 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  
  int cat;
  
  if((fabs(eta)<1.479 && fBrem<0.06) || (fabs(eta)>1.479 && fBrem<0.1)) 
    cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) 
    cat=0;
  else 
    cat=2;
  
  return cat;
}

bool result(const GsfElectron* electron,
	    const Event& e,
	    const string& quality_) { 
  
  double eta = fabs(electron->p4().Eta());
  const ClusterShapeRef& shapeRef = getClusterShape(electron, e);
  
  double eOverP = electron->eSuperClusterOverP();
  double eSeed = electron->superCluster()->seed()->energy();
  double pin  = electron->trackMomentumAtVtx().R();   
  double eSeedOverPin = eSeed/pin; 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  
  double hOverE = electron->hadronicOverEm();
  double sigmaee = sqrt(shapeRef->covEtaEta());
  double deltaPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();
  double deltaEtaIn = electron->deltaEtaSuperClusterTrackAtVtx();
  
  int eb;
  if (eta < 1.479) 
    eb = 0;
  else {
    eb = 1; 
    sigmaee = sigmaee - 0.02*(fabs(eta) - 2.3);   //correct sigmaetaeta dependence on eta in endcap
  }

  //std::vector<double> cut;
  const double* cut;
    
  // ROBUST Selection
  if (quality_ == "robust") {

    // hoe, sigmaEtaEta, dPhiIn, dEtaIn
    if (eta < 1.479)
      cut = robustEleIDCutsBarrel; //cuts_.getParameter<std::vector<double> >("barrel");
    else
      cut = robustEleIDCutsEndcap; //cuts_.getParameter<std::vector<double> >("endcap");

    if (hOverE > cut[0]) 
      return false;    

    if (sigmaee > cut[1]) 
      return false;    

    if (fabs(deltaPhiIn) > cut[2]) 
      return false;    

    if (fabs(deltaEtaIn) > cut[3]) 
      return false;    
    
    return true;
  }
  
  int cat = classify(electron);

  // TIGHT Selection
  if (quality_ == "tight") {
    
    if ((eOverP < 0.8) && (fBrem < 0.2)) 
      return false;
    
    if (eOverP < 0.9*(1-fBrem))
      return false;
    
    cut = tightEleIDCuts_hOverE; //cuts_.getParameter<std::vector<double> >("hOverE");
    if (hOverE > cut[cat+4*eb]) 
      return false;    
    
    cut = tightEleIDCuts_sigmaEtaEta; //cuts_.getParameter<std::vector<double> >("sigmaEtaEta");
    if (sigmaee > cut[cat+4*eb]) 
      return false;    
    
    cut = tightEleIDCuts_deltaPhiIn; //cuts_.getParameter<std::vector<double> >("deltaPhiIn");
    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cut[cat+4*eb]) 
        return false;    
    } else {
      if (fabs(deltaPhiIn) > cut[3+4*eb])
        return false;
    }
    
    cut = tightEleIDCuts_deltaEtaIn; //cuts_.getParameter<std::vector<double> >("deltaEtaIn");
    if (fabs(deltaEtaIn) > cut[cat+4*eb]) 
      return false;    
    
    cut = tightEleIDCuts_eSeedOverPinMin; //cuts_.getParameter<std::vector<double> >("eSeedOverPinMin");
    if (eSeedOverPin < cut[cat+4*eb]) 
      return false;  

    cut = tightEleIDCuts_eSeedOverPinMax; //cuts_.getParameter<std::vector<double> >("eSeedOverPinMax");
    if (eSeedOverPin > cut[cat+4*eb]) 
      return false;  

    return true;
  }
  
    // LOOSE Selection
  if (quality_ == "loose") {
    
    if ((eOverP < 0.8) && (fBrem < 0.2)) 
      return false;
    
    cut = looseEleIDCuts_hOverE; //cuts_.getParameter<std::vector<double> >("hOverE");
    if (hOverE > cut[cat+4*eb]) 
      return false;    
    
    cut = looseEleIDCuts_sigmaEtaEta; //cuts_.getParameter<std::vector<double> >("sigmaEtaEta");
    if (sigmaee > cut[cat+4*eb]) 
      return false;    
    
    cut = looseEleIDCuts_deltaPhiIn; //cuts_.getParameter<std::vector<double> >("deltaPhiIn");
    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cut[cat+4*eb]) 
        return false;    
    } else {
      if (fabs(deltaPhiIn) > cut[3+4*eb])
        return false;
    }
    
    cut = looseEleIDCuts_deltaEtaIn; //cuts_.getParameter<std::vector<double> >("deltaEtaIn");
    if (fabs(deltaEtaIn) > cut[cat+4*eb]) 
      return false;    
    
    cut = looseEleIDCuts_eSeedOverPin; //cuts_.getParameter<std::vector<double> >("eSeedOverPin");
    if (eSeedOverPin < cut[cat+4*eb]) 
      return false;  
    
    return true; 
  }
  
  return false;
}
