#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muExample.h"

using namespace std;

Zprime2muExample::Zprime2muExample(const edm::ParameterSet& config) 
  : Zprime2muAnalysis(config)
{
  // Book histos, get parameters from config here.
  int    nBins     = config.getParameter<int>   ("nBins");
  double lowerMass = config.getParameter<double>("lowerMass");
  double upperMass = config.getParameter<double>("upperMass");

  // By using the TFileService (instantiated in Zprime2muAnalysis as
  // "fs"), saving to a ROOT file and deleting the TH1F object are
  // taken care of for you.
  hDilMass = fs->make<TH1F>("hDilMass",
			    "Dilepton (#mu^{+}#mu^{-}) mass spectrum",
			    nBins, lowerMass, upperMass);
}

void Zprime2muExample::analyze(const edm::Event& event, 
			       const edm::EventSetup& eSetup) {
  Zprime2muAnalysis::analyze(event, eSetup);

  // Do whatever you want with the dileptons here.
  for (unsigned idil = 0; idil < allDileptons[lOP].size(); idil++)
    hDilMass->Fill(allDileptons[lOP].at(idil).mass());
}

DEFINE_FWK_MODULE(Zprime2muExample);
