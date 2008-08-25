#ifndef TTbarAnalysis_h
#define TTbarAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "NtupleBase.h"
class TTbarAnalysis : public NtupleBase {
	public :

		TTbarAnalysis(TTree *tree=0):NtupleBase(tree){};
		virtual ~TTbarAnalysis();
		//virtual void Loop();
		//virtual void  Loop(TString s_output="", Long64_t entryToStart=0, Long64_t entriesToAnalyze = -1);
		void  Loop(TString s_output, Long64_t entryToStart, Long64_t entriesToAnalyze);

		void Initialize(TString s_input);
		int GetObservedDistributions(Long64_t jentry, int process, double EventWeight);
		
		void GetJetProbDistribution(Long64_t jentry, double EventWeight);
		void GetJetPerformance();

		int GetNumberOfBQuarkInAcceptance(Long64_t jentry, double EventWeight);

		void PrintEventsCountingInformation();
		void PlotDiMuonInvMassWith1Band2B(double A1, double A2, double n1, double n2);

		// - - hists related
		void BookingObsvHists();
		void WriteObsvHists();

	private :
		TFile * outputFile_ ;

		//- - - - booking histograms
		TH1F *h_obsvDeltaR, *h_obsvDeltaEta, *h_obsvDeltaPhi, *h_obsvInvMass,
			 *h_obsvMuonP, *h_obsvMuonEta, *h_obsvMuonPhi, *h_obsvMuonPt;
		//	TH1F *h_obsvInvMassPassLooseCut, *h_obsvInvMassPassLikelihoodRatioCut;
		TH1F *h_obsvInvMass1B, *h_obsvInvMass2B;

		
		TH1F *h_jetProb[8];// btag probablity for different flavoring jets,  0-for non-flavor jets, 1-d, 2-u, 3-s, 4-c, 5-b, 6-g, 7-for non-b flavor jets
		TH1F *h_bTagEffVsProbCut[8]; // for b-jet, it's efficiency vs. probablity cut, for non-b jet, it's faking rate vs. probability cut
				

		// --- --- test for acceptance of b parton from TTbar
		TH1F *h_bEta;
		TH1F *h_bPt;	


		// - - - - -Saving numbers here , you can add new ENTRY/bin here
		// - - - - -    to calc&save errors automatically
		TH1F *h_numbers;
		enum ENUMCountingNumbers { NTTbarEvents =1, 
			NTTbarEventsWith1BPartonAccepted,
			NTTbarEventsWith2BPartonAccepted, 
			NEventsWith1BTagged, 
			NEventsWith2BTagged, 
			NEventsWithMBTagged,
			NEventsWith1BTagged_ZMassVeto, 
			NEventsWith2BTagged_ZMassVeto, 
			NEventsWithMBTagged_ZMassVeto,
			NTTbarEventsWith1BTagged,
			NTTbarEventsWith2BTagged, 
			NTTbarEventsWithMBTagged,
			NTTbarEventsWith1BTagged_ZMassVeto,
			NTTbarEventsWith2BTagged_ZMassVeto, 
			NTTbarEventsWithMBTagged_ZMassVeto,
			};



		Float_t ttbarEventsTotal;
		// -- only count the b  from t
		Float_t ttbarEvetns1bAccepted;
		Float_t ttbarEvetns2bAccepted;

		// -- count b-tagged jets // all
		Float_t bTagging1JetEventsPassed;
		Float_t bTagging2JetEventsPassed;
		Float_t bTaggingMJetEventsPassed;
		Float_t bTagging1JetEventsPassed_ZMassVeto;
		Float_t bTagging2JetEventsPassed_ZMassVeto;
		Float_t bTaggingMJetEventsPassed_ZMassVeto;

		// --- count b-tagged jets for several processes
		enum ENUMProcess { NonSoupProcess=0, GumboDiJets=1, ChowderWJets, ChowderZJets, ChowderTTbar, GumboPhotonJets, GumboMinbias, GumboHiggs, GumboZprime, StewBBbarToJPsi, StewDiJets, StewBottomonium, StewCharmonium, StewBBE, StewPPElectronX, StewPPMuonX };
		TString S_PrimaryProcess[20] ; 
		Float_t bTagging1JetEventsPassed_Proc[20] ;
		Float_t bTagging2JetEventsPassed_Proc[20] ;
		Float_t bTaggingMJetEventsPassed_Proc[20] ;
		Float_t bTagging1JetEventsPassed_ZMassVeto_Proc[20] ;
		Float_t bTagging2JetEventsPassed_ZMassVeto_Proc[20] ;
		Float_t bTaggingMJetEventsPassed_ZMassVeto_Proc[20] ;
		TH1F *h_InvariantMassByPrimaryProcess[20];
		TH1F *h_InvMass1BByPrimaryProcess[20], *h_InvMass2BByPrimaryProcess[20];

};

#endif

