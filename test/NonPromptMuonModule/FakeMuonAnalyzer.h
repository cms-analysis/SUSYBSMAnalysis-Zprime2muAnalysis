#ifndef FakeMuonAnalyzer_h
#define FakeMuonAnalyzer_h

#include <TRandom.h>
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
		
		// ------ determine if this event is W event
		bool KeepEventsWithTMassOnWtmPeak(Long64_t jentry);

		// ------ determine if this event is Z event
		bool KeepEventsWithMassOnZPeak(Long64_t jentry);

		// ------ determine if this muon passes a set of cuts
		bool AcceptedMuon(Long64_t jentry, int muonIndex);

		// ------ get index of the muon which is on W transverse mass peak with met
		int GetPromptMuonFromWtmPeak(Long64_t jentry);

		// ------ get indexes of two muons whose invariant mass is closest to Z mass
		vector<int> GetPromptMuonsFromZPeak(Long64_t jentry);

		// ------ get a muon from collection randomly as prompt muon
		int GetRandomPromptMuon(Long64_t jentry);

		// ------ get the muon with maximum pt as prompt muon
		int GetMaxPtPromptMuon(Long64_t jentry, int firstMuon = -1);

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

		TRandom _gRandom;

};

#endif

