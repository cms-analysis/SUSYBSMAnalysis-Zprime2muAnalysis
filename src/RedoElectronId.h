#ifndef RedoElectronId_h
#define RedoElectronId_h

#include <string>

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "FWCore/Framework/interface/Event.h"

const reco::ClusterShapeRef& getClusterShape(const reco::PixelMatchGsfElectron* electron, const edm::Event& e);
int classify(const reco::PixelMatchGsfElectron* electron);
bool result(const reco::PixelMatchGsfElectron* electron, const edm::Event& e, const std::string& quality_);

#endif // RedoElectronId_h
