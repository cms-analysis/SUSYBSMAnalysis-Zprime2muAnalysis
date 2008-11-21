#include "TCanvas.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TText.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muMuonIdStudy.h"

using namespace edm;
using namespace reco;
using namespace std;

Zprime2muMuonIdStudy::Zprime2muMuonIdStudy(const ParameterSet& config)
  : Zprime2muRecLevelAnalysis(config) {
  bookHistos();
}

void Zprime2muMuonIdStudy::endJob() {
  Zprime2muRecLevelAnalysis::endJob();
  
  drawHistos();
}

void Zprime2muMuonIdStudy::analyze(const Event& event,
				   const EventSetup& eSetup) {
  Zprime2muRecLevelAnalysis::analyze(event, eSetup);

  const unsigned OK = 3;

  for (int i = 0; i < 2; i++) {
    const int rec = ltk + i;
    for (int ilep = 0; ilep < int(allLeptons[rec].size()); ilep++) {
      const CandidateBaseRef& cand = allLeptons[rec][ilep];
      const Muon* muon    = toConcretePtr<Muon>(cand);
      const Muon* gmrMuon = toConcretePtr<Muon>(recLevelHelper.sameSeedLepton(cand, lgmr));
      
      unsigned ptrOK = (muon != 0) | ((gmrMuon != 0) << 1);
      h_OK->Fill(ptrOK | (i << 2));
      if (ptrOK != OK)
	continue;

      unsigned eOK = muon->isEnergyValid() | (gmrMuon->isEnergyValid() << 1);
      h_eOK->Fill(eOK | (i << 2));
      if (eOK == OK) {
	const MuonEnergy& calE = muon->calEnergy();
	const MuonEnergy& gmrCalE = gmrMuon->calEnergy();
	
	h_calE_em   [i]->Fill(calE.em    - gmrCalE.em);
	h_calE_emS9 [i]->Fill(calE.emS9  - gmrCalE.emS9);
	h_calE_had  [i]->Fill(calE.had   - gmrCalE.had);
	h_calE_hadS9[i]->Fill(calE.hadS9 - gmrCalE.hadS9);
	h_calE_ho   [i]->Fill(calE.ho    - gmrCalE.ho);
	h_calE_hoS9 [i]->Fill(calE.hoS9  - gmrCalE.hoS9);
      }

      unsigned isoOK = muon->isIsolationValid() | (gmrMuon->isIsolationValid() << 1);
      h_isoOK->Fill(isoOK | (i << 2));
      if (isoOK == OK) {
	const MuonIsolation& iso03 = muon->isolationR03();
	const MuonIsolation& gmrIso03 = gmrMuon->isolationR03();

	h_iso03_emEt   [i]->Fill(iso03.emEt    - gmrIso03.emEt);
	h_iso03_hadEt  [i]->Fill(iso03.hadEt   - gmrIso03.hadEt);
	h_iso03_hoEt   [i]->Fill(iso03.hoEt    - gmrIso03.hoEt);
	h_iso03_nJets  [i]->Fill(iso03.nJets   - gmrIso03.nJets);
	h_iso03_nTracks[i]->Fill(iso03.nTracks - gmrIso03.nTracks);
	h_iso03_sumPt  [i]->Fill(iso03.sumPt   - gmrIso03.sumPt);

	const MuonIsolation& iso05 = muon->isolationR05();
	const MuonIsolation& gmrIso05 = gmrMuon->isolationR05();

	h_iso05_emEt   [i]->Fill(iso05.emEt    - gmrIso05.emEt);
	h_iso05_hadEt  [i]->Fill(iso05.hadEt   - gmrIso05.hadEt);
	h_iso05_hoEt   [i]->Fill(iso05.hoEt    - gmrIso05.hoEt);
	h_iso05_nJets  [i]->Fill(iso05.nJets   - gmrIso05.nJets);
	h_iso05_nTracks[i]->Fill(iso05.nTracks - gmrIso05.nTracks);
	h_iso05_sumPt  [i]->Fill(iso05.sumPt   - gmrIso05.sumPt);
      }

      unsigned timeOK = muon->isTimeValid() | (gmrMuon->isTimeValid() << 1);
      h_timeOK->Fill(timeOK | (i << 2));
      if (timeOK == OK) {
	const MuonTime& time = muon->time();
	const MuonTime& gmrTime = gmrMuon->time();

	h_time_freeInverseBeta   [i]->Fill(time.freeInverseBeta    - gmrTime.freeInverseBeta);
	h_time_freeInverseBetaErr[i]->Fill(time.freeInverseBetaErr - gmrTime.freeInverseBetaErr);
	h_time_inverseBeta       [i]->Fill(time.inverseBeta        - gmrTime.inverseBeta);
	h_time_inverseBetaErr    [i]->Fill(time.inverseBetaErr     - gmrTime.inverseBetaErr);
	h_time_nStations         [i]->Fill(time.nStations          - gmrTime.nStations);
	h_time_timeAtIpInOut     [i]->Fill(time.timeAtIpInOut      - gmrTime.timeAtIpInOut);
	h_time_timeAtIpInOutErr  [i]->Fill(time.timeAtIpInOutErr   - gmrTime.timeAtIpInOutErr);
	h_time_timeAtIpOutIn     [i]->Fill(time.timeAtIpOutIn      - gmrTime.timeAtIpOutIn);
	h_time_timeAtIpOutInErr  [i]->Fill(time.timeAtIpOutInErr   - gmrTime.timeAtIpOutInErr);
      }
       
      unsigned matchOK = muon->isMatchesValid() | (gmrMuon->isMatchesValid() << 1);
      h_matchOK->Fill(matchOK | (i << 2));
      if (matchOK == OK) {
	// check matches?
      }
    }      
  }
}

TH1F* Zprime2muMuonIdStudy::makeHist(const char* baseName, int i,
				     int nbins, double min, double max) const {
  string name = nameHist(baseName, i);
  TH1F* h = fs->make<TH1F>(name.c_str(), name.c_str(), nbins, min, max);
  return h;
}

void Zprime2muMuonIdStudy::bookHistos() {
  h_OK      = fs->make<TH1F>("h_OK",      "h_OK",      8, 0, 8);
  h_eOK     = fs->make<TH1F>("h_eOK",     "h_eOK",     8, 0, 8);
  h_isoOK   = fs->make<TH1F>("h_isoOK",   "h_isoOK",   8, 0, 8);
  h_timeOK  = fs->make<TH1F>("h_timeOK",  "h_timeOK",  8, 0, 8);
  h_matchOK = fs->make<TH1F>("h_matchOK", "h_matchOK", 8, 0, 8);

  for (int i = 0; i < 2; i++) {
    h_calE_em   [i] = makeHist("h_calE_em",    i, 50, -5, 5);
    h_calE_emS9 [i] = makeHist("h_calE_emS9",  i, 50, -5, 5);
    h_calE_had  [i] = makeHist("h_calE_had",   i, 50, -5, 5);
    h_calE_hadS9[i] = makeHist("h_calE_hadS9", i, 50, -5, 5);
    h_calE_ho   [i] = makeHist("h_calE_ho",    i, 50, -5, 5);
    h_calE_hoS9 [i] = makeHist("h_calE_hoS9",  i, 50, -5, 5);

    h_iso03_emEt   [i] = makeHist("h_iso03_emEt",    i, 50, -5, 5);
    h_iso03_hadEt  [i] = makeHist("h_iso03_hadEt",   i, 50, -5, 5);
    h_iso03_hoEt   [i] = makeHist("h_iso03_hoEt",    i, 50, -5, 5);
    h_iso03_nJets  [i] = makeHist("h_iso03_nJets",   i, 50, -5, 5);
    h_iso03_nTracks[i] = makeHist("h_iso03_nTracks", i, 50, -5, 5);
    h_iso03_sumPt  [i] = makeHist("h_iso03_sumPt",   i, 50, -5, 5);

    h_iso05_emEt   [i] = makeHist("h_iso05_emEt",    i, 50, -5, 5);
    h_iso05_hadEt  [i] = makeHist("h_iso05_hadEt",   i, 50, -5, 5);
    h_iso05_hoEt   [i] = makeHist("h_iso05_hoEt",    i, 50, -5, 5);
    h_iso05_nJets  [i] = makeHist("h_iso05_nJets",   i, 50, -5, 5);
    h_iso05_nTracks[i] = makeHist("h_iso05_nTracks", i, 50, -5, 5);
    h_iso05_sumPt  [i] = makeHist("h_iso05_sumPt",   i, 50, -5, 5);

    h_time_freeInverseBeta   [i] = makeHist("h_time_freeInverseBeta",    i, 50, -5, 5);
    h_time_freeInverseBetaErr[i] = makeHist("h_time_freeInverseBetaErr", i, 50, -5, 5);
    h_time_inverseBeta       [i] = makeHist("h_time_inverseBeta",        i, 50, -5, 5);
    h_time_inverseBetaErr    [i] = makeHist("h_time_inverseBetaErr",     i, 50, -5, 5);
    h_time_nStations         [i] = makeHist("h_time_nStations",          i, 50, -5, 5);
    h_time_timeAtIpInOut     [i] = makeHist("h_time_timeAtIpInOut",      i, 50, -5, 5);
    h_time_timeAtIpInOutErr  [i] = makeHist("h_time_timeAtIpInOutErr",   i, 50, -5, 5);
    h_time_timeAtIpOutIn     [i] = makeHist("h_time_timeAtIpOutIn",      i, 50, -5, 5);
    h_time_timeAtIpOutInErr  [i] = makeHist("h_time_timeAtIpOutInErr",   i, 50, -5, 5);
  }
}

void Zprime2muMuonIdStudy::drawHistos() const {
  TText t;
  t.SetTextFont(12);
  t.SetTextSize(.03);

  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 700);
  TPostScript *ps = new TPostScript("muonidstudy.ps", 111);

  const int NUM_PAGES = 90;
  TPad *pad[NUM_PAGES];
  for (int i_page = 0; i_page <= NUM_PAGES; i_page++)
    pad[i_page] = new TPad("","", .05, .05, .95, .93);

  int page = 0;
  ostringstream strpage;
  string tit;
  TPaveLabel *title = 0;
  
  gStyle->SetOptLogy(1);

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title = new TPaveLabel(0.1,0.94,0.9,0.98, "Info block statuses");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);  h_OK->Draw();
  pad[page]->cd(2);  h_eOK->Draw();
  pad[page]->cd(3);  h_isoOK->Draw();
  pad[page]->cd(4);  h_timeOK->Draw();
  pad[page]->cd(5);  h_matchOK->Draw();
  c1->Update();

  //TString names[2] = {"FMS", "PMR"};
  TString names[2] = {"TK", "FMS"};

  for (int i = 0; i < 2; i++) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    TString titl = "MuonEnergy (" + names[i] + ")";
    title = new TPaveLabel(0.1,0.94,0.9,0.98, titl);
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    pad[page]->cd(1);  h_calE_em   [i]->Draw();
    pad[page]->cd(2);  h_calE_emS9 [i]->Draw();
    pad[page]->cd(3);  h_calE_had  [i]->Draw();
    pad[page]->cd(4);  h_calE_hadS9[i]->Draw();
    pad[page]->cd(5);  h_calE_ho   [i]->Draw();
    pad[page]->cd(6);  h_calE_hoS9 [i]->Draw();
    c1->Update();
  }

  for (int i = 0; i < 2; i++) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    TString titl = "MuonIsolation (R03) (" + names[i] + ")";
    title = new TPaveLabel(0.1,0.94,0.9,0.98, titl);
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    pad[page]->cd(1);  h_iso03_emEt   [i]->Draw();
    pad[page]->cd(2);  h_iso03_hadEt  [i]->Draw();
    pad[page]->cd(3);  h_iso03_hoEt   [i]->Draw();
    pad[page]->cd(4);  h_iso03_nJets  [i]->Draw();
    pad[page]->cd(5);  h_iso03_nTracks[i]->Draw();
    pad[page]->cd(6);  h_iso03_sumPt  [i]->Draw();
    c1->Update();
  }

  for (int i = 0; i < 2; i++) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    TString titl = "MuonIsolation (R05) (" + names[i] + ")";
    title = new TPaveLabel(0.1,0.94,0.9,0.98, titl);
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,3);
    pad[page]->cd(1);  h_iso05_emEt   [i]->Draw();
    pad[page]->cd(2);  h_iso05_hadEt  [i]->Draw();
    pad[page]->cd(3);  h_iso05_hoEt   [i]->Draw();
    pad[page]->cd(4);  h_iso05_nJets  [i]->Draw();
    pad[page]->cd(5);  h_iso05_nTracks[i]->Draw();
    pad[page]->cd(6);  h_iso05_sumPt  [i]->Draw();
    c1->Update();
  }

  for (int i = 0; i < 2; i++) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    TString titl = "MuonTime (" + names[i] + ")";
    title = new TPaveLabel(0.1,0.94,0.9,0.98, titl);
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(3,3);
    pad[page]->cd(1);  h_time_freeInverseBeta   [i]->Draw();
    pad[page]->cd(2);  h_time_freeInverseBetaErr[i]->Draw();
    pad[page]->cd(3);  h_time_inverseBeta       [i]->Draw();
    pad[page]->cd(4);  h_time_inverseBetaErr    [i]->Draw();
    pad[page]->cd(5);  h_time_nStations         [i]->Draw();
    pad[page]->cd(6);  h_time_timeAtIpInOut     [i]->Draw();
    pad[page]->cd(7);  h_time_timeAtIpInOutErr  [i]->Draw();
    pad[page]->cd(8);  h_time_timeAtIpOutIn     [i]->Draw();
    pad[page]->cd(9);  h_time_timeAtIpOutInErr  [i]->Draw();
		   
    c1->Update();
  }

  ps->Close();

  delete title;
  delete ps;
  delete c1;
}

DEFINE_FWK_MODULE(Zprime2muMuonIdStudy);
