#ifndef FakeMuonAnalyzer_h
#define FakeMuonAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "NtupleBase.h"
#include "FakeMuonGenerator.h"
class FakeMuonAnalyzer : public NtupleBase {
	public :

		FakeMuonAnalyzer(TTree *tree=0):NtupleBase(tree){};
		virtual ~FakeMuonAnalyzer();
		//virtual void Loop();
		virtual void  Loop(TString s_output="", Long64_t entryToStart=0, Long64_t entriesToAnalyze = -1);

		void Initialize(TString s_input);
		void GetObservedDistributions(Long64_t jentry, double EventWeight);
		void GetPredictedDistributions(Long64_t jentry, double EventWeight);
		void FillPredictedHists(vector<FakeMuonGenerator::fakeMuon> blahs, double EventWeight);

		// - - hists related
		void BookingFakeHists();
		void BookingObsvHists();
		void WriteFakeHists();
		void WriteObsvHists();

	private :
		TFile * outputFile_ ;
		FakeMuonGenerator* myGenerator;

		//- - - - booking histograms
		TH1F *h_fakeDeltaR, *h_fakeDeltaEta, *h_fakeDeltaPhi, *h_fakeInvMass,
			 *h_fakeMuonP, *h_fakeMuonEta, *h_fakeMuonPhi, *h_fakeMuonPt;
		TH1F *h_obsvDeltaR, *h_obsvDeltaEta, *h_obsvDeltaPhi, *h_obsvInvMass,
			 *h_obsvMuonP, *h_obsvMuonEta, *h_obsvMuonPhi, *h_obsvMuonPt;
		//	TH1F *h_obsvInvMassPassLooseCut, *h_obsvInvMassPassLikelihoodRatioCut;
};

#endif

