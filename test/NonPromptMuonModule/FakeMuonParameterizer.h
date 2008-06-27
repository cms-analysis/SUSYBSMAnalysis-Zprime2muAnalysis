#ifndef FakeMuonParameterizer_h
#define FakeMuonParameterizer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <iostream>
#include <vector>
#include <TString.h>
#include "NtupleBase.h"
class FakeMuonParameterizer : public NtupleBase {
	public :

		FakeMuonParameterizer(TTree *tree=0) :NtupleBase(tree){};
		~FakeMuonParameterizer();
		void Initialize(TString s_input="");
		virtual void  Loop(TString s_output="", Long64_t entryToStart=0, Long64_t entriesToAnalyze = -1);
		// - - hists related
		void BookingHists();
		void WriteHists();
	private:
		TFile *outputFile_;
		TH2F *h_allJetPEta;//p:eta
		TH2F *h_matchedJetPEta;//p:eta
		TH2F *h_muonPVsJetP;
		TH2F *h_muonEtaVsJetEta;
		TH2F *h_muonPhiVsJetPhi;
		TH2F *h_muonPtVsJetPt;	
		TH1F *h_dRMuonJet;

};

#endif

