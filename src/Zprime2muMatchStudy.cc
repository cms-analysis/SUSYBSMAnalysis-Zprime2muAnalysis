#include "TCanvas.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TText.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muMatchStudy.h"

using namespace edm;
using namespace reco;
using namespace std;

Zprime2muMatchStudy::Zprime2muMatchStudy(const ParameterSet& config)
  : Zprime2muRecLevelAnalysis(config) {
  bookHistos();
}

void Zprime2muMatchStudy::endJob() {
  Zprime2muRecLevelAnalysis::endJob();
  
  drawHistos();
}

void Zprime2muMatchStudy::analyze(const Event& event,
				   const EventSetup& eSetup) {
  Zprime2muRecLevelAnalysis::analyze(event, eSetup);

  for (int irec = 0; irec < MAX_LEVELS; irec++) {
    for (int ilep = 0; ilep < int(allLeptons[irec].size()); ilep++) {
      const CandidateBaseRef& icand = allLeptons[irec][ilep];
      int iId = recLevelHelper.id(icand);

      for (int jrec = 0; jrec < MAX_LEVELS-1; jrec++) {
	if (irec == jrec) continue;

	const CandidateBaseRef jcand[2] = {
	  recLevelHelper.closestLepton(icand, jrec),
	  recLevelHelper.sameSeedLepton(icand, jrec)
	};

	const int closInd = recLevelHelper.id(jcand[0]);
	const int seedInd = recLevelHelper.id(jcand[1]);
	const int index[2] = {closInd, seedInd};

	for (int t = 0; t < 2; t++) {
	  if (index[t] != -999) {
	    hDeltaPt[t][irec][jrec]->Fill(jcand[t]->pt()     - icand->pt());
	    hDeltaQ [t][irec][jrec]->Fill(jcand[t]->charge() - icand->charge());
	    hDeltaR [t][irec][jrec]->Fill(deltaR(*jcand[t], *icand));
	  }
	}

	if (seedInd != closInd) {
	  if (seedInd != -999 && closInd == -999)
	    LogTrace("matchStudy") << "Only match by seed is found:\n"
				   << "  index = " << iId << " seedInd = " << seedInd;
	  else if (seedInd == -999 && closInd != -999)
	    LogTrace("matchStudy") << "Only match in eta and phi is found:\n"
				   << "  index = " << iId << " closInd = " << closInd;
	  else
	    LogTrace("matchStudy")
	      << "Difference in proximity and seed matches:\n"
	      << "  index = " << iId << " seedInd = " << seedInd
	      << " closInd = " << closInd;
	}
      }
    }
  }
}

TH1F* Zprime2muMatchStudy::makeHist(const char* baseName, int type, int irec, int jrec,
				    const char* title, int nbins, double min, double max) const {
  if (irec == jrec || (type == 1 && irec < l3 && jrec < l3)) return 0;
  string name = nameHist(baseName, type, irec, jrec);
  TH1F* h = fs->make<TH1F>(name.c_str(), title, nbins, min, max);
  return h;
}

void Zprime2muMatchStudy::bookHistos() {
  for (int type = 0; type < 2; type++) {
    for (int irec = 0; irec < MAX_LEVELS; irec++) {
      for (int jrec = 0; jrec < MAX_LEVELS-1; jrec++) {
	hDeltaR[type][irec][jrec]  = makeHist("hDeltaR",  type, irec, jrec, "DeltaR",   50,    0,   1);
	hDeltaQ[type][irec][jrec]  = makeHist("hDeltaQ",  type, irec, jrec, "DeltaQ",    5,   -2,   2);
	hDeltaPt[type][irec][jrec] = makeHist("hDeltaPt", type, irec, jrec, "DeltaPt", 100, -500, 500);
      }
    }
  }
}

void Zprime2muMatchStudy::drawHistos() const {
  TText t;
  t.SetTextFont(12);
  t.SetTextSize(.03);

  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 700);
  TPostScript *ps = new TPostScript("matchstudy.ps", 111);

  const int NUM_PAGES = 90;
  TPad *pad[NUM_PAGES];
  for (int i_page = 0; i_page <= NUM_PAGES; i_page++)
    pad[i_page] = new TPad("","", .05, .05, .95, .93);

  int page = 0;
  ostringstream strpage;
  string tit;
  TPaveLabel *title;

  
  for (int type = 0; type < 2; type++) {
    for (int irec = 0; irec < MAX_LEVELS; irec++) {
      for (int jrec = 0; jrec < MAX_LEVELS-1; jrec++) {
	if (irec == jrec || (type == 1 && (irec < l3 || jrec < l3))) continue;

	ps->NewPage();
	c1->Clear();
	c1->cd(0);
	TString titl = type == 0 ? "Closest match, " : "Same-seed match, ";
	titl += irec;
	titl += " -> ";
	titl += jrec;
	title = new TPaveLabel(0.1,0.94,0.9,0.98,titl);
	title->SetFillColor(10);
	title->Draw();
	strpage << "- " << (++page) << " -";
	t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
	gStyle->SetOptLogy(1);
	pad[page]->Draw();
	pad[page]->Divide(2,2);
	pad[page]->cd(1);  hDeltaR [type][irec][jrec]->Draw();
	pad[page]->cd(2);  hDeltaQ [type][irec][jrec]->Draw();
	pad[page]->cd(3);  hDeltaPt[type][irec][jrec]->Draw();
	gStyle->SetOptLogy(0);
	c1->Update();
      }
    }
  }

  ps->Close();

  delete title;
  delete ps;
  delete c1;
}

DEFINE_FWK_MODULE(Zprime2muMatchStudy);
