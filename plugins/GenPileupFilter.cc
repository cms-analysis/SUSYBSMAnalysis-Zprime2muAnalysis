#include "TH1F.h"
#include "TH2F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class GenPileupFilter : public edm::EDFilter {
public:
  explicit GenPileupFilter(const edm::ParameterSet&);

private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  int min_intime;
  int max_intime;
  int min_late;
  int max_late;

  TH1F* pudist[3];
  TH2F* intime_vs_late;
};

GenPileupFilter::GenPileupFilter(const edm::ParameterSet& cfg)
  : min_intime(cfg.getParameter<int>("min_intime")),
    max_intime(cfg.getParameter<int>("max_intime")),
    min_late(cfg.getParameter<int>("min_late")),
    max_late(cfg.getParameter<int>("max_late"))
{
  edm::Service<TFileService> fs;
  for (int i = 0; i < 3; ++i)
    pudist[i] = fs->make<TH1F>(TString::Format("pudist%i", i), "", 50, 0, 50);
  intime_vs_late = fs->make<TH2F>("intime_vs_late", "", 51, -1, 50, 51, -1, 50);
}

bool GenPileupFilter::filter(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<std::vector<PileupSummaryInfo> > pileup;
  event.getByLabel("addPileupInfo", pileup);

  int intime = -1;
  int late = -1;
  for (std::vector<PileupSummaryInfo>::const_iterator it = pileup->begin(), end = pileup->end(); it != end; ++it) {
    int n = it->getPU_NumInteractions();
    int bx = it->getBunchCrossing();
    assert(bx >= -1 && bx <= 1);

    if (bx == 0) {
      assert(intime == -1);
      intime = n;
    }
    else if (bx == 1) {
      assert(late == -1);
      late = n;
    }
      
    pudist[bx+1]->Fill(it->getPU_NumInteractions());
  }
  
  intime_vs_late->Fill(late, intime);

  return intime >= min_intime && intime <= max_intime && late >= min_late && late <= max_late;
}

DEFINE_FWK_MODULE(GenPileupFilter);
