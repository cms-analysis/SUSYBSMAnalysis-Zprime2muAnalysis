#define TTbarAnalysis_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH3.h>
#include "TH1F.h"
#include "TSystem.h"
#include <stdio.h>
#include <TLorentzVector.h>
#include "TSystem.h"
#include <TF1.h>

#include "TTbarAnalysis.h"
#include "NtupleUtilities.C"
#include "Configure.C"

#include "PlotUtils/Utils.C"

const double MUMASS = 0.10566;
void TTbarAnalysis::Initialize(TString s_input)
{
	if(s_input == "") s_input = _InputFileToDoAnalyzing;
	TTree *tree=0;
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(s_input);
		if (!f) {
			f = new TFile(s_input);
		}
		tree = (TTree*)gDirectory->Get("myTree");

	}
	Init(tree);

	 ttbarEventsTotal = 0;
	// -- only count the b  from t
	 ttbarEvetns1bAccepted = 0;
	 ttbarEvetns2bAccepted = 0;

	// -- count b-tagged jets // all
	 bTagging1JetEventsPassed = 0;
	 bTagging2JetEventsPassed = 0;
	 bTaggingMJetEventsPassed = 0;
	 bTagging1JetEventsPassed_ZMassVeto = 0;
	 bTagging2JetEventsPassed_ZMassVeto = 0;
	 bTaggingMJetEventsPassed_ZMassVeto = 0;
	 for(int i=0; i<20; i++ )
	 {
		 bTagging1JetEventsPassed_Proc[i] = 0 ;
		 bTagging2JetEventsPassed_Proc[i] = 0 ;
		 bTaggingMJetEventsPassed_Proc[i] = 0 ;
		 bTagging1JetEventsPassed_ZMassVeto_Proc[i] = 0 ;
		 bTagging2JetEventsPassed_ZMassVeto_Proc[i] = 0 ;
		 bTaggingMJetEventsPassed_ZMassVeto_Proc[i] = 0 ;
	 }

	 S_PrimaryProcess[NonSoupProcess] = "NonSoupProcess";
	 S_PrimaryProcess[GumboDiJets] = "GumboDiJets";
	 S_PrimaryProcess[GumboHiggs]="GumboHiggs";
	 S_PrimaryProcess[GumboZprime]="GumboZprime";
	 S_PrimaryProcess[GumboPhotonJets]="GumboPhotonJets";
	 S_PrimaryProcess[GumboMinbias]="GumboMinbias";
	 S_PrimaryProcess[ChowderTTbar]="ChowderTTbar";
	 S_PrimaryProcess[ChowderZJets]="ChowderZJets";
	 S_PrimaryProcess[ChowderWJets]="ChowderWJets";
	 S_PrimaryProcess[StewPPMuonX]="StewPPMuonX";
	 S_PrimaryProcess[StewPPElectronX]="StewPPElectronX";
	 S_PrimaryProcess[StewBBE]="StewBBE";
	 S_PrimaryProcess[StewBBbarToJPsi]="StewBBbarToJPsi";
	 S_PrimaryProcess[StewDiJets]="StewDiJets";
	 S_PrimaryProcess[StewBottomonium]="StewBottomonium";
	 S_PrimaryProcess[StewCharmonium]="StewCharmonium";


}
//void TTbarAnalysis::Loop()
void TTbarAnalysis::Loop(TString s_output, Long64_t entryToStart, Long64_t entriesToAnalyze)
	//void TTbarAnalysis::Loop(TString s_output, Long64_t entryToStart , Long64_t entriesToAnalyze)
{
	if (fChain == 0) return;

	//-- Booking hists and creating file
	if(s_output=="") s_output = _OutputFile;
	outputFile_  = new TFile( s_output, "RECREATE" ) ;
	outputFile_->cd();
	BookingObsvHists();

	cout<<"Intialize done..."<<endl;

	Long64_t nentries = fChain->GetEntriesFast();
	cout<<"All Entries_"<<nentries<<endl;
	if(entriesToAnalyze < 0 ) entriesToAnalyze = nentries;
	cout<<"Starting from "<<entryToStart<<"th entry"<<endl;
	cout<<"Entries to be analyzed ="<<entriesToAnalyze<<endl;
	Long64_t nbytes = 0, nb = 0;
	Long64_t jentry= entryToStart;
	vector<Double_t> vweight ;
	entriesToAnalyze += entryToStart;
	if(entriesToAnalyze > nentries ) entriesToAnalyze = nentries;
	for (; jentry<entriesToAnalyze;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry%10000==0) cout<<"Entry: "<<jentry<<endl;

		int process = 99 ; // signal sample, non-soup sample
		if(_ifRunOnSoup) process = GetProcessId(genEventProcID, genEventScale, genEventFilterEff, evtWeight, 100.);

		//- - - -  count the different weights 
		bool weightAdded = false;
		for(int i = 0; i<(int)vweight.size(); i++)
		{
			if(evtWeight == vweight[i])
			{
				weightAdded = true;
				break;
			}
		}
		if(weightAdded==false) vweight.push_back(evtWeight);

		// - - - do we want to take into accound evt weight ?
		if( _ifTakeIntoAccountWeight_Analyzing == 0 )evtWeight=1;
		Double_t EventWeight = evtWeight;
		if(!_ifRunOnSoup) EventWeight =1;


		if(process < 27 && process > 21 )  //ChowderTTbar
		{
			ttbarEventsTotal += 1*EventWeight; 
			h_numbers->Fill(NTTbarEvents, EventWeight);
			int numb = GetNumberOfBQuarkInAcceptance(jentry, _BQuarkAcceptanceEta);
			if (numb == 1 ) 
			{
				ttbarEvetns1bAccepted += 1*EventWeight;
				h_numbers->Fill(NTTbarEventsWith1BPartonAccepted, EventWeight);
			}
			if (numb == 2 ) 
			{
				ttbarEvetns2bAccepted += 1*EventWeight;
				h_numbers->Fill(NTTbarEventsWith2BPartonAccepted, EventWeight);
			}
		}

		// count the jets above 20 GeV/c
		Int_t numOfJetsWithPtLargerThan20GeV =0;
		for(int i =0; i<jetSize; i++){
			double jetPt_ = (*jetPt)[i];
			if (jetPt_<_MinJetPt) continue;
			if (fabs((*jetEta)[i]) > _MaxJetEta ) continue;
			numOfJetsWithPtLargerThan20GeV++;
		}
		if (numOfJetsWithPtLargerThan20GeV< _MinNumberOfJetsAboveThreshold) continue;

		GetObservedDistributions(jentry, process ,EventWeight);
		GetJetProbDistribution(jentry, EventWeight);

	}//end entries loop



	GetJetPerformance();

	cout<<"Print weights:"<<endl;
	for(int i=0; i<(int)vweight.size();i++)
	{
		cout<<i <<" th weight= "<<vweight[i]<<endl;
	}

	PrintEventsCountingInformation();



}//end Loop() 
void TTbarAnalysis::PrintEventsCountingInformation()
{
	//--------------------------------------------------
	//--- printout numbers .............................
	cout<<"|          *MC Counts*             |||"<<endl;
	Float_t A1 = ttbarEvetns1bAccepted/(float)ttbarEventsTotal;
	Float_t A2 = ttbarEvetns2bAccepted/(float)ttbarEventsTotal;
	printf("|Total TTbar Events =                               | %.2e ||\n", ttbarEventsTotal);
	printf("|TTbar with one b in acceptance ( fabs(eta)< %3.1f )  | %.2e |  A1 = %5.3f |\n", _BQuarkAcceptanceEta, ttbarEvetns1bAccepted , A1);
	printf("|TTbar with two b in acceptance ( fabs(eta)< %3.1f )  | %.2e |  A2 = %5.3f |\n", _BQuarkAcceptanceEta, ttbarEvetns2bAccepted , A2); 

	printf("\n\n|        *Selection Criteria*          ||| \n");
	printf("|Muon Pt > %3.1f (GeV), fabs(eta)< %3.1f ||| \n", _MinMuonPt, _MaxMuonEta);
	if(_ifApplyIsolationCut ) printf("|Iso03sumPt < %3.1f (GeV), ", _MaxIso03sumPt);
	else printf("|No Isolation Cut ,");
	if(_ifApplyImpactParCut ) printf("fabs(dxy)< %3.2f (cm) |||\n", _MaxDxy);
	else printf("No Impact Parameter Cut |||\n");
	printf("|Di-Muon Mass > %4.1f , DeltaR > %4.1f |||\n", _DiMuonMassCut, _DeltaRCutForTwoMuons );	
	printf("|Jet Pt > %4.1f , Jet fabs(eta) < %3.1f |||\n", _MinJetPt, _MaxJetEta );
	printf("|BTagging Discriminator = %3.1f|||\n", _BTagDiscriminator);

	printf("\n\n");
	printf("|Events with one b-tagged, n1 = %.2e |||\n", bTagging1JetEventsPassed);
	printf("|Events with two b-tagged, n2 = %.2e |||\n", bTagging2JetEventsPassed);
	printf("|Events with >2  b-tagged, nM = %.2e |||\n", bTaggingMJetEventsPassed);

	float b_eff = 0; 
	if(A2!=0 && bTagging2JetEventsPassed!=0)
		b_eff = (A1/A2+2)/(bTagging1JetEventsPassed/bTagging2JetEventsPassed+2);
	printf("| b-tagging eff = ( A1/A2 +2 )/( n1/n2 +2 ) = %5.3f|||\n", b_eff);

	printf("\n\n");
	printf("|Data Sample| n1 | n2 | n1Tight | n2Tight | \n");
	printf("|ttbar+jets | %.2e | %.2e | %.2e | %.2e |\n", 
			bTagging1JetEventsPassed_Proc[ChowderTTbar],
			bTagging2JetEventsPassed_Proc[ChowderTTbar],
			bTagging1JetEventsPassed_ZMassVeto_Proc[ChowderTTbar], 
			bTagging2JetEventsPassed_ZMassVeto_Proc[ChowderTTbar]);

	printf("|w+jets     | %.2e | %.2e | %.2e | %.2e |\n", 
			bTagging1JetEventsPassed_Proc[ChowderWJets],
			bTagging2JetEventsPassed_Proc[ChowderWJets],
			bTagging1JetEventsPassed_ZMassVeto_Proc[ChowderWJets], 
			bTagging2JetEventsPassed_ZMassVeto_Proc[ChowderWJets]);

	printf("|z+jets     | %.2e | %.2e | %.2e | %.2e |\n", 
			bTagging1JetEventsPassed_Proc[ChowderZJets],
			bTagging2JetEventsPassed_Proc[ChowderZJets],
			bTagging1JetEventsPassed_ZMassVeto_Proc[ChowderZJets], 
			bTagging2JetEventsPassed_ZMassVeto_Proc[ChowderZJets]);

	printf("| all       | %.2e | %.2e | %.2e | %.2e |\n", 
			bTagging1JetEventsPassed,
			bTagging2JetEventsPassed,
			bTagging1JetEventsPassed_ZMassVeto, 
			bTagging2JetEventsPassed_ZMassVeto);

	printf("\n\n");
	//PlotDiMuonInvMassWith1Band2B(A1, A2, bTagging1JetEventsPassed, bTagging2JetEventsPassed);
}
int TTbarAnalysis::GetObservedDistributions(Long64_t jentry, int process,  double EventWeight)
{
	if(_Debug)cout<<"observing..."<<endl;
	fChain->GetEntry(jentry);

	int PrimaryProcess = NonSoupProcess;
	if(_ifRunOnSoup) 
	{
		if (process < 0 ) cout<<"Unkonw Process =" << process << endl;
		else if(process <= 10 ) PrimaryProcess = ChowderWJets;
		else if(process <= 21 ) PrimaryProcess = ChowderZJets;
		else if(process <= 26 ) PrimaryProcess = ChowderTTbar;
		else if(process == 28 ) PrimaryProcess = GumboMinbias;
		else if(process <= 47 ) PrimaryProcess = GumboDiJets;
		else if(process <= 57 ) PrimaryProcess = GumboPhotonJets;
		else if(process == 58 ) PrimaryProcess = GumboHiggs;
		else if(process == 59 ) PrimaryProcess = GumboZprime;
		else if(process == 60 ) PrimaryProcess = StewBBbarToJPsi;
		else if(process == 61 ) PrimaryProcess = StewDiJets;
		else if(process <= 63 ) PrimaryProcess = StewBottomonium;
		else if(process <= 65 ) PrimaryProcess = StewCharmonium;
		else if(process <= 68 ) PrimaryProcess = StewBBE;
		else if(process == 69 ) PrimaryProcess = StewPPElectronX;
		else if(process == 70 ) PrimaryProcess = StewPPMuonX;
		else if(process == 99 ) PrimaryProcess = NonSoupProcess;
		else {cout<<"Warning : Unknown Process"<<endl;}

	}

	// - - -count how many b tags
	int nBTags = 0;
	if( (int)jetTagDiscriminator->size() != jetSize) cout<<"Warning: jettags != jetsize"<<endl;
	for(int i = 0; i<jetSize; i++)
	{
		Float_t jeteta = fabs((*jetEta)[i]);
		Float_t jetpt = (*jetPt)[i];
		if(jeteta>_MaxJetEta) continue;
		if(jetpt <_MinJetPt ) continue;
		if( (*jetTagDiscriminator)[i] > _BTagDiscriminator ) nBTags++; 
	}

	bool passDiMuonSelection  = false;
	bool passDiMuonMassCut = false;
	bool passZMassVeto = false;
	for(int m1 =0; m1< gmrSize; m1++ )
	{
		double glbpt = (*gmrPt)[m1];   		
		if(glbpt <_MinMuonPt) continue;
		double glbp = (*gmrP)[m1];
		double glbeta = (*gmrEta)[m1]; 
		if(fabs(glbeta) > _MaxMuonEta ) continue;
		double glbphi = (*gmrPhi)[m1];
		double glbcharge = (*gmrCharge)[m1];
		double glbdxy = (*gmrDXY)[m1];
		double glbiso03sumpt = (*gmrIso03sumPt)[m1];

		if( _ifApplyImpactParCut &&  fabs(glbdxy) > _MaxDxy ) continue;
		if( _ifApplyIsolationCut && glbiso03sumpt > _MaxIso03sumPt )	continue;

		// - - - - - - - - - -
		if(_SkipMuonsWhosePAboveThreshold && glbp> _MaxMuonP ) continue;

		h_obsvMuonP  ->Fill(glbp,EventWeight);
		h_obsvMuonEta->Fill(glbeta,EventWeight);
		h_obsvMuonPhi->Fill(glbphi,EventWeight);
		h_obsvMuonPt  ->Fill(glbpt,EventWeight);

		for(int m2 =m1+1; m2< gmrSize; m2++ )
		{
			double glbpt2 = (*gmrPt)[m2];   		
			if(glbpt2 <_MinMuonPt) continue;
			double glbp2 = (*gmrP)[m2];
			double glbeta2 = (*gmrEta)[m2]; 
			if(fabs(glbeta2) > _MaxMuonEta ) continue;
			double glbphi2 = (*gmrPhi)[m2];
			double glbcharge2 = (*gmrCharge)[m2];
			double glbdxy2 = (*gmrDXY)[m2];
			double glbiso03sumpt2 = (*gmrIso03sumPt)[m2];

			if( _ifApplyImpactParCut &&  fabs(glbdxy2) > _MaxDxy ) continue;
			if( _ifApplyIsolationCut && glbiso03sumpt2 > _MaxIso03sumPt )	continue;

			// - - - - - - - - - -
			if(_SkipMuonsWhosePAboveThreshold && glbp2> _MaxMuonP ) continue;

			if(_OppositeSign && glbcharge ==glbcharge2) continue; //skip it with same sign muons
			if(!_OppositeSign && glbcharge !=glbcharge2) continue; //

			//			if(_ifSkipTwoMuonsAssociatedWithTheSameJet && jetIndex[m1]==jetIndex[m2]) continue;

			TLorentzVector vmu1, vmu2, vmm; // LorentzVectors for mu+ mu- 
			vmu1.SetPtEtaPhiM(glbpt, glbeta, glbphi, MUMASS);
			vmu2.SetPtEtaPhiM(glbpt2, glbeta2, glbphi2, MUMASS);
			vmm = vmu1 + vmu2;		

			double mass = vmm.M();

			float tmpDR = GetDeltaR(glbeta, glbeta2, glbphi, glbphi2);
			if(tmpDR < _DeltaRCutForTwoMuons) continue;

			passDiMuonSelection = true;


			if( mass > _DiMuonMassCut )	passDiMuonMassCut = true;
			if( mass < _ZMassVetoLow || mass > _ZMassVetoHigh ) passZMassVeto = true;

			//if(mass<1)cout<<"mass before ="<<mass<<endl;
			h_obsvDeltaR->Fill(tmpDR,EventWeight);
			h_obsvInvMass->Fill(mass,EventWeight);
			h_InvariantMassByPrimaryProcess[PrimaryProcess]->Fill(mass,EventWeight);
			if(nBTags == 1) 
			{
				h_obsvInvMass1B->Fill(mass, EventWeight);
				h_InvMass1BByPrimaryProcess[PrimaryProcess]->Fill(mass, EventWeight);
			} 
			if(nBTags >= 2) 
			{
				h_obsvInvMass2B->Fill(mass, EventWeight);
				h_InvMass2BByPrimaryProcess[PrimaryProcess]->Fill(mass, EventWeight);
			}
			//if(mass<1)cout<<"mass after ="<<mass<<endl;

			h_obsvDeltaEta->Fill(fabs(glbeta-glbeta2),EventWeight);
			h_obsvDeltaPhi->Fill(DeltaPhi(glbphi,glbphi2),EventWeight);
		}
	}
	if( !passDiMuonSelection ) return 0;
	if (nBTags == 1 )
	{
		bTagging1JetEventsPassed += 1*EventWeight;
		bTagging1JetEventsPassed_Proc[PrimaryProcess] += 1*EventWeight;
		if(passZMassVeto)
		{
			bTagging1JetEventsPassed_ZMassVeto_Proc[PrimaryProcess] += 1*EventWeight;
			bTagging1JetEventsPassed_ZMassVeto += 1*EventWeight;
			h_numbers->Fill(NEventsWith1BTagged_ZMassVeto, EventWeight);
			if(PrimaryProcess == ChowderTTbar) h_numbers->Fill(NTTbarEventsWith1BTagged_ZMassVeto, EventWeight);
		}
		h_numbers->Fill(NEventsWith1BTagged, EventWeight);
		if(PrimaryProcess == ChowderTTbar) h_numbers->Fill(NTTbarEventsWith1BTagged, EventWeight);
	}
	if ( nBTags == 2 ) 
	{
		bTagging2JetEventsPassed += 1*EventWeight;
		bTagging2JetEventsPassed_Proc[PrimaryProcess] += 1*EventWeight;
		if(passZMassVeto)
		{
			bTagging2JetEventsPassed_ZMassVeto_Proc[PrimaryProcess] += 1*EventWeight;
			bTagging2JetEventsPassed_ZMassVeto += 1*EventWeight;
			h_numbers->Fill(NEventsWith2BTagged_ZMassVeto, EventWeight);
			if(PrimaryProcess == ChowderTTbar) h_numbers->Fill(NTTbarEventsWith2BTagged_ZMassVeto, EventWeight);
		}
		h_numbers->Fill(NEventsWith2BTagged, EventWeight);
		if(PrimaryProcess == ChowderTTbar) h_numbers->Fill(NTTbarEventsWith2BTagged, EventWeight);
	}
	if (nBTags >  2 ) 
	{
		bTaggingMJetEventsPassed += 1*EventWeight;
		bTaggingMJetEventsPassed_Proc[PrimaryProcess] += 1*EventWeight;
		if(passZMassVeto)
		{
			bTaggingMJetEventsPassed_ZMassVeto_Proc[PrimaryProcess] += 1*EventWeight;
			bTaggingMJetEventsPassed_ZMassVeto += 1*EventWeight;
			h_numbers->Fill(NEventsWithMBTagged_ZMassVeto, EventWeight);
			if(PrimaryProcess == ChowderTTbar) h_numbers->Fill(NTTbarEventsWithMBTagged_ZMassVeto, EventWeight);
		}
		h_numbers->Fill(NEventsWithMBTagged, EventWeight);
		if(PrimaryProcess == ChowderTTbar) h_numbers->Fill(NTTbarEventsWithMBTagged, EventWeight);
	}
	return nBTags;
}

TTbarAnalysis::~TTbarAnalysis() {
	outputFile_->cd();
	WriteObsvHists();
	outputFile_->Write() ;
	outputFile_->Close() ;
	delete outputFile_;
}
void TTbarAnalysis::GetJetProbDistribution(Long64_t jentry, double EventWeight)
{
	fChain->GetEntry(jentry);
	for(int i = 0; i<jetSize; i++)
	{
		if((*jetPt)[i]<_MinJetPt) continue;
		if((*jetEta)[i]>_MaxJetEta) continue;	
		int thePartonIndex  = FindPartonMatchedToThisJet(i, jetEta, jetPhi, pytSize, pytId, pytEta, pytPhi, pytPt, pytNDaughter, pytDaughterPdgId);	
		int id = 0;
		if( thePartonIndex < 0 ) id =0 ;
		else id = abs((*pytId)[thePartonIndex]);
		if (id == 21) id = 6;
		float prob = (*jetTagDiscriminator)[i] ;
		h_jetProb[id] -> Fill(prob, EventWeight);
		if( id!=5 ) h_jetProb[7] -> Fill(prob, EventWeight);
	}

}
void TTbarAnalysis::GetJetPerformance()
{
	int nbins = h_jetProb[0]->GetNbinsX();		
	double tot[8];
	float eff = 0;
	for(int i = 0; i<8 ; i++)
	{
		tot[i] = h_jetProb[i]->Integral(0, nbins+1);
		for(int bin = 1; bin <= nbins; bin++)	
		{
			eff = 1 - h_jetProb[i]->Integral(0, bin -1 ) / tot[i] ; 
			h_bTagEffVsProbCut[i]->SetBinContent(bin, eff);
		}
	}
}
void TTbarAnalysis::BookingObsvHists()
{
	h_obsvDeltaR = new TH1F("h_obsvDeltaR","h_obsvDeltaR",628,0,6.28);	
	h_obsvDeltaEta = new TH1F("h_obsvDeltaEta","h_obsvDeltaEta;|#{Delta}#{eta}|",628,0,6.28);	
	h_obsvDeltaPhi = new TH1F("h_obsvDeltaPhi","h_obsvDeltaPhi",800,-4,4);	
	h_obsvInvMass = new TH1F("h_obsvInvMass","h_obsvInvMass",7000,0,7000);
	//	h_obsvInvMassPassLooseCut = new TH1F("h_obsvInvMassPassLooseCut","h_obsvInvMassPassLooseCut",7000,0,7000);
	//	h_obsvInvMassPassLikelihoodRatioCut = new TH1F("h_obsvInvMassPassLikelihoodRatioCut","h_obsvInvMassPassLikelihoodRatioCut",7000,0,7000);
	h_obsvMuonP = new TH1F("h_obsvMuonP","h_obsvMuonP",7000,0,7000);
	h_obsvMuonPt = new TH1F("h_obsvMuonPt","h_obsvMuonPt",7000,0,7000);
	h_obsvMuonEta = new TH1F("h_obsvMuonEta","h_obsvMuonEta",800,-4,4);
	h_obsvMuonPhi = new TH1F("h_obsvMuonPhi","h_obsvMuonPhi",800,-4,4);
	h_obsvInvMass1B = new TH1F("h_obsvInvMass1B","h_obsvInvMass1B",7000,0,7000);
	h_obsvInvMass2B = new TH1F("h_obsvInvMass2B","h_obsvInvMass2B",7000,0,7000);

	h_obsvMuonP->Sumw2();
	h_obsvMuonPt->Sumw2();
	h_obsvMuonEta->Sumw2();
	h_obsvMuonPhi->Sumw2();
	h_obsvDeltaR->Sumw2();
	h_obsvDeltaEta->Sumw2();
	h_obsvDeltaPhi->Sumw2();
	h_obsvInvMass->Sumw2();
	h_obsvInvMass1B->Sumw2();
	h_obsvInvMass2B->Sumw2();
	//	h_obsvInvMassPassLooseCut->Sumw2();
	//	h_obsvInvMassPassLikelihoodRatioCut->Sumw2();

	TString s_partonname[8] = {"non", "d" , "u", "s", "c", "b" , "g", "nonb"};
	for(int i=0; i<8; i++)
	{
		TString s_title = "h_jetProb_"; s_title+=s_partonname[i];
		h_jetProb[i] = new TH1F(s_title,s_title, 120, 0, 1.2);
		h_jetProb[i]->Sumw2();
		s_title = "h_bTagEffVsProbCut_"; s_title+=s_partonname[i];
		h_bTagEffVsProbCut[i] = new TH1F(s_title, s_title, 120, 0, 1.2); // the x axis should be same as above hist
	}

	h_numbers = new TH1F("h_numbers","h_numbers", 20, 1, 21);
	h_numbers->Sumw2();

	outputFile_->mkdir("MassByProcess");
	outputFile_->cd("MassByProcess");
	for(int i=0; i<20; i++)
	{
		if(S_PrimaryProcess[i]=="")continue;
		TString s_title = "h_InvariantMassByPrimaryProcess"; s_title+=S_PrimaryProcess[i];
		h_InvariantMassByPrimaryProcess[i]= new TH1F(s_title,s_title,7000,0,7000);
		h_InvariantMassByPrimaryProcess[i]->Sumw2();
		s_title = "h_InvMass1BByPrimaryProcess"; s_title+=S_PrimaryProcess[i];
		h_InvMass1BByPrimaryProcess[i] = new TH1F(s_title, s_title, 7000, 0, 7000);
		h_InvMass1BByPrimaryProcess[i]->Sumw2();
		s_title = "h_InvMass2BByPrimaryProcess"; s_title+=S_PrimaryProcess[i];
		h_InvMass2BByPrimaryProcess[i] = new TH1F(s_title, s_title, 7000, 0, 7000);
		h_InvMass2BByPrimaryProcess[i]->Sumw2();
	}
}
void TTbarAnalysis::WriteObsvHists()
{
	h_obsvMuonP->Write();
	h_obsvMuonPt->Write();
	h_obsvMuonEta->Write();
	h_obsvMuonPhi->Write();
	h_obsvDeltaR->Write();
	h_obsvDeltaEta->Write();
	h_obsvDeltaPhi->Write();
	h_obsvInvMass->Write();
	h_obsvInvMass2B->Write();
	h_obsvInvMass1B->Write();
	//	h_obsvInvMass1BPassLooseCut->Write();
	//	h_obsvInvMassPassLikelihoodRatioCut->Write();
	for(int i=0; i<8; i++)
	{
		h_jetProb[i]->Write();
		h_bTagEffVsProbCut[i]->Write();
	}
}
int TTbarAnalysis::GetNumberOfBQuarkInAcceptance( Long64_t jentry, double _BQuarkAcceptanceEta)
{
	fChain->GetEntry(jentry);
	// - -  only count the  b  which comes from t
	int numb = 0;
	for(int i=0; i<pytSize; i++)
	{
		int b_id = abs((*pytId)[i]);
		int t_id = abs((*pytMotherPdgId)[i]);
		Float_t b_eta = fabs((*pytEta)[i]); 
		if( b_id==5 && t_id==6 && b_eta < _BQuarkAcceptanceEta ) numb++;
	}
	return numb;
}
void TTbarAnalysis::PlotDiMuonInvMassWith1Band2B(double A1, double A2, double n1, double n2)
{
	if(A2==0 || n2==0) 
	{
		cout<<"A2 =" <<A2 <<", n2 ="<<n2<<endl;
		return;
	}

	double b_eff = (A1/A2 +2 )/(n1/n2 +2 );
	//---------------
	//	N1 = N*A1;
	//	N2 = N*A2;
	//	n2 = N2*b_eff*b_eff;
	//	n1 = N1*b_eff+2*N2*b_eff*(1-b_eff);
	//  
	//  N = N2/A2; N2=n2/b_eff/b_eff; --> N = n2/b_eff/b_eff/A2;
	//  n1 = N*A1*b_eff + 2*N*A2*b_eff*(1-b_eff); --> N = n1/(A1*b_eff + 2*A2*b_eff*(1-b_eff));
	//---------------	


	double n2ScaleFactor = 	1/b_eff/b_eff/A2;
	double n1ScaleFactor =  1/(A1*b_eff + 2*A2*b_eff*(1-b_eff));

	int binwidth = 20;
	int RangeMin = 0;
	int RangeMax = 1000;

	TCanvas *can = new TCanvas("can","invariant mass");
	can->SetLogy(1);

	h_obsvInvMass->Rebin(binwidth);	
	h_obsvInvMass1B->Rebin(binwidth);	
	h_obsvInvMass2B->Rebin(binwidth);	

	h_obsvInvMass->GetXaxis()->SetRangeUser(RangeMin, RangeMax);
	h_obsvInvMass1B->GetXaxis()->SetRangeUser(RangeMin, RangeMax);
	h_obsvInvMass2B->GetXaxis()->SetRangeUser(RangeMin, RangeMax);

	h_obsvInvMass1B->Scale(n1ScaleFactor);
	h_obsvInvMass2B->Scale(n2ScaleFactor);

	char yTitle[100];
	sprintf(yTitle, "Entries / %.0f GeV/c", h_obsvInvMass->GetBinWidth(0));
	axis1F(h_obsvInvMass, "M_{#mu#mu} [GeV/c]",yTitle);
	Draw(h_obsvInvMass, kBlue, 1001, "HIST", kBlue, 2, kBlue);
	Draw(h_obsvInvMass1B, 0,1, "same", 2, 2, 2, 1.5, 22);
	Draw(h_obsvInvMass2B, 0,1, "same", 3, 2, 3, 1.5, 23);

	TLegend *leg = (TLegend*)SetLegend(0.58, 0.53, 0.94, 0.91);
	leg->AddEntry(h_obsvInvMass, " all ", "lfp");
	leg->AddEntry(h_obsvInvMass1B, " 1 b estimate", "lp");
	leg->AddEntry(h_obsvInvMass2B, " 2 b estimate", "lp");
	leg->Draw();
}


