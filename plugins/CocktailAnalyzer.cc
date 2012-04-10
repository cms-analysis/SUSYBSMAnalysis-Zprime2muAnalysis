#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"

struct tmr_datum {
  double tkonly_res2;
  double tpfms_res2;
  double tpfms_tkonly_probdiff;
  tmr_datum(double tkonly_res, double tpfms_res, double tkonly_prob, double tpfms_prob)
    : tkonly_res2(tkonly_res * tkonly_res),
      tpfms_res2(tpfms_res * tpfms_res),
      tpfms_tkonly_probdiff(tpfms_prob - tkonly_prob) {}
};

class CocktailAnalyzer : public edm::EDAnalyzer {
 public:
  explicit CocktailAnalyzer(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

 private:
  edm::InputTag muon_src;
  double eta_min, eta_max;

  enum { global, tkonly, tpfms, picky, num_ingredients };
  enum { global_tkonly, global_tpfms, global_picky, tkonly_tpfms, tkonly_picky, tpfms_picky, num_pairs };
  
  TH1F* choice_labels(TH1F* choice) {
    TAxis* ax = choice->GetXaxis();
    ax->SetBinLabel(1 + global, "Global");
    ax->SetBinLabel(1 + tkonly, "TkOnly");
    ax->SetBinLabel(1 + tpfms,  "TPFMS");
    ax->SetBinLabel(1 + picky,  "Picky");
    return choice;
  }

  void fill_choice(TH1F* choice, reco::TrackRef tracks[num_ingredients], reco::TrackRef chosen) {
    int i = 0;
    for ( ; i < num_ingredients; ++i)
      if (chosen == tracks[i])
	break;
    choice->Fill(i);
  }

  TH1F* TMRCocktailChoice;
  TH1F* TunePCocktailChoice;
  TH1F* SigmaSwitchCocktailChoice;
  TH2F* TrackLnChi2TailProb[num_pairs];
  TH1F* TMRSelectedTPFMSResolution;
  TH1F* TMRSelectedTkOnlyResolution;
  TH1F* TMRRejectedTPFMSResolution;
  TH1F* TMRRejectedTkOnlyResolution;

  std::vector<tmr_datum> tmr_data;
};

CocktailAnalyzer::CocktailAnalyzer(const edm::ParameterSet& cfg)
  : muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    eta_min(cfg.getParameter<double>("eta_min")),
    eta_max(cfg.getParameter<double>("eta_max"))
{
  edm::Service<TFileService> fs;

  int ni = int(num_ingredients);
  TMRCocktailChoice         = choice_labels(fs->make<TH1F>("TMRCocktailChoice",         "TMR cocktail choice",          ni, 0, ni));
  TunePCocktailChoice       = choice_labels(fs->make<TH1F>("TunePCocktailChoice",       "Tune P cocktail choice",       ni, 0, ni));
  SigmaSwitchCocktailChoice = choice_labels(fs->make<TH1F>("SigmaSwitchCocktailChoice", "Sigma-switch cocktail choice", ni, 0, ni));

  TrackLnChi2TailProb[global_tkonly] = fs->make<TH2F>("TrackLnChi2TailProbGlobalVsTkOnly", "-ln(P), tracker-only vs. global", 500, 0, 100, 500, 0, 100);
  TrackLnChi2TailProb[global_tpfms]  = fs->make<TH2F>("TrackLnChi2TailProbGlobalVsTPFMS",  "-ln(P), TPFMS vs. global",        500, 0, 100, 500, 0, 100);
  TrackLnChi2TailProb[global_picky]  = fs->make<TH2F>("TrackLnChi2TailProbGlobalVsPicky",  "-ln(P), picky vs. global",        500, 0, 100, 500, 0, 100);
  TrackLnChi2TailProb[tkonly_tpfms]  = fs->make<TH2F>("TrackLnChi2TailProbTkOnlyVsTPFMS",  "-ln(P), TPFMS vs. tracker-only",  500, 0, 100, 500, 0, 100);
  TrackLnChi2TailProb[tkonly_picky]  = fs->make<TH2F>("TrackLnChi2TailProbTkOnlyVsPicky",  "-ln(P), picky vs. tracker-only",  500, 0, 100, 500, 0, 100);
  TrackLnChi2TailProb[tpfms_picky]   = fs->make<TH2F>("TrackLnChi2TailProbTPFMSVsPicky",   "-ln(P), picky vs. TPFMS",         500, 0, 100, 500, 0, 100);

  TMRSelectedTPFMSResolution  = fs->make<TH1F>("TMRSelectedTPFMSResolution",  "TPFMS inv pT res when selected by TMR",  100, -0.3, 0.3);
  TMRSelectedTkOnlyResolution = fs->make<TH1F>("TMRSelectedTkOnlyResolution", "TkOnly inv pT res when selected by TMR", 100, -0.3, 0.3);
  TMRRejectedTPFMSResolution  = fs->make<TH1F>("TMRRejectedTPFMSResolution",  "TPFMS inv pT res when rejected by TMR",  100, -0.3, 0.3);
  TMRRejectedTkOnlyResolution = fs->make<TH1F>("TMRRejectedTkOnlyResolution", "TkOnly inv pT res when rejected by TMR", 100, -0.3, 0.3);
}

void CocktailAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&) {
  edm::Handle<pat::MuonCollection> muons;
  event.getByLabel(muon_src, muons);

  for (pat::MuonCollection::const_iterator mu = muons->begin(), mue = muons->end(); mu != mue; ++mu) {
    double aeta = fabs(mu->eta());
    if (aeta < eta_min || aeta > eta_max)
      continue;

    reco::TrackRef tks[num_ingredients] = {
      mu->globalTrack(),
      mu->innerTrack(),
      mu->tpfmsMuon(),
      mu->pickyMuon()
    };

    const reco::TrackRef& tmr = muon::TMR(mu->innerTrack(), mu->tpfmsMuon()).first;

    // Histogram the chosen track type for each cocktail (e.g. for
    // TMR, the number of times TPFMS or tracker-only are chosen).
    fill_choice(TMRCocktailChoice,         tks, tmr);
    fill_choice(TunePCocktailChoice,       tks, muon::tevOptimized(*mu).first);
    fill_choice(SigmaSwitchCocktailChoice, tks, muon::sigmaSwitch(*mu).first);

    // For each unique pair of cocktail ingredients, scatterplot the
    // -log(chi2 tail probability) of each (the figure of merit for
    // the TMR and Tune P cocktails).
    double probs[num_ingredients] = {0.};
    for (size_t i = 0; i < num_ingredients; ++i)
      probs[i] = muon::trackProbability(tks[i]);

    size_t k = 0;
    for (size_t i = 0; i < num_ingredients; ++i)
      for (size_t j = i+1; j < num_ingredients; ++j, ++k)
	TrackLnChi2TailProb[k]->Fill(probs[i], probs[j]);

    // For each TMR choice, if we have MC truth, plot separately the
    // q/pt resolutions for the selected tracks and the rejected
    // tracks.
    const reco::GenParticle* gen = mu->genParticle();
    if (gen != 0) {
      const double gen_qoverpt = gen->charge()/gen->pt();
      const double tkonly_qoverpt = tks[tkonly]->charge() / tks[tkonly]->pt();
      const double tpfms_qoverpt  = tks[tpfms]->charge() / tks[tpfms]->pt();
      const double tkonly_res = tkonly_qoverpt/gen_qoverpt - 1;
      const double tpfms_res = tpfms_qoverpt/gen_qoverpt - 1;

      if (tmr == tks[tpfms]) {
	TMRSelectedTPFMSResolution ->Fill(tpfms_res);
	TMRRejectedTkOnlyResolution->Fill(tkonly_res);
      }
      else if (tmr == tks[tkonly]) {
	TMRSelectedTkOnlyResolution->Fill(tkonly_res);
	TMRRejectedTPFMSResolution ->Fill(tpfms_res);
      }

      // Store the resolutions and the prob-diffs for TPFMS and TkOnly
      // tracks so we can study at the end of the job the effect of
      // moving the TMR cut.
      tmr_datum d(tkonly_res, tpfms_res, probs[tkonly], probs[tpfms]);
      tmr_data.push_back(d);
    }
  }
}

double tmr_res(const std::vector<tmr_datum>& data, const double cut) {
  double f = 0; 
  std::vector<tmr_datum>::const_iterator d = data.begin(), e = data.end(); 
  for ( ; d != e; ++d) 
    f += d->tpfms_tkonly_probdiff > cut ? d->tkonly_res2 : d->tpfms_res2; 
  return f; 
} 

void CocktailAnalyzer::endJob() {
  if (tmr_data.size() == 0)
    return;

  const int M = 1000; 
  const double cmin = -10, cmax = 10, dc = (cmax-cmin)/M;
  double xx[M] = {0.}, yy[M] = {0.};
  for (int i = 0; i < M; ++i) {
    xx[i] = cmin + dc*i;
    yy[i] = tmr_res(tmr_data, xx[i]);
  } 

  edm::Service<TFileService> fs;
  TGraph* graph = fs->make<TGraph>(M, xx, yy); 
  graph->GetXaxis()->SetRangeUser(cmin, cmax);
}

DEFINE_FWK_MODULE(CocktailAnalyzer);
