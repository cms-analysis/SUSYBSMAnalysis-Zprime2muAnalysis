#define FakeMuonParameterizer_cxx
#include "FakeMuonParameterizer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH3.h>
#include "TSystem.h"
#include <TString.h>
#include <TF1.h>

#include "NtupleUtilities.C"
#include "LikelihoodUtilities.C"
#include "Configure.C"

TF1* fDxyS(0), *fDxyB(0), *fIsoSumPtS(0), *fIsoSumPtB(0);
double ResponseCut_Likelihood  = 0;
void FakeMuonParameterizer::Initialize(TString s_input)
{
	if(s_input=="") s_input = _InputFileToDoParameterizing;
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
}
void FakeMuonParameterizer::Loop(TString s_output, Long64_t entryToStart , Long64_t entriesToAnalyze)
{
	if (fChain == 0) return;

	if(s_output=="") s_output = _OutputFileParameterizedHistos;
	outputFile_  = new TFile( s_output, "RECREATE" ) ;
	outputFile_ -> cd();
	BookingHists();

	if(Use_Likelihood)
	{	
		cout<<"Start to initialize the Likelihood response..."<<endl;
		TFile *infileLLR = new TFile(_InputFileLikelihoodCalibrated,"READ");
		fDxyS = (TF1*)infileLLR->Get("fDxyS");
		fDxyB = (TF1*)infileLLR->Get("fDxyB");
		fIsoSumPtS = (TF1*)infileLLR->Get("fIsoSumPtS");
		fIsoSumPtB = (TF1*)infileLLR->Get("fIsoSumPtB");
		if(Use_Likelihood)ResponseCut_Likelihood = GetLRCutAtEff(_InputFileLikelihoodCalibrated, _KeepPromptMuonEfficiency);
	}


	Long64_t nentries = fChain->GetEntriesFast();
	cout<<"All Entries_"<<nentries<<endl;
	if(entriesToAnalyze < 0 ) entriesToAnalyze = nentries;
	cout<<"Starting from "<<entryToStart<<"th entry"<<endl;
	cout<<"Entries to be analyzed ="<<entriesToAnalyze<<endl;
	Long64_t nbytes = 0, nb = 0;
	entriesToAnalyze += entryToStart;
	if(entriesToAnalyze > nentries ) entriesToAnalyze = nentries;
	for (Long64_t jentry=entryToStart; jentry< entriesToAnalyze ;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry%10000==0) cout<<"Entry: "<<jentry<<endl;
		if(_Debug)cout<<"Enter entries loop"<<endl;

		// -  - - - -half event to analyze
		if((_TakeHalfEvents==2) && (jentry%2!=0)) continue;
		if((_TakeHalfEvents==1) && (jentry%2==0)) continue;

		// - - - do we want to take into accound evt weight ?
		if( _ifTakeIntoAccountWeight_Parameterizing == 0 )evtWeight=1;

		//- - - the dijets from soup ? 
		int process = GetProcess(pytSize, pytStatus, pytId, pytNMother);
		bool isDijetsEvent = (process==1);
		if( _ifRunOnGumboSoup && (isDijetsEvent == false)) continue;

		// count the jets above 20 GeV/c
		Int_t numOfJetsWithPtLargerThan20GeV =0;
		for(int i =0; i<jetSize; i++){
			double jetPt_ = (*jetPt)[i];
			if (jetPt_<_MinJetPt) continue;
			if (fabs((*jetEta)[i]) > _MaxJetEta ) continue;
			numOfJetsWithPtLargerThan20GeV++;
		}
		if (numOfJetsWithPtLargerThan20GeV< _MinNumberOfJetsAboveThreshold ) continue;


		for(int i =0; i<jetSize; i++){
			//-------------------------------------------
			if((*jetPt)[i]<_MinJetPt) continue;
			if (fabs((*jetEta)[i]) > _MaxJetEta ) continue;
			//-------------------------------------------
			h_allJetPEta	       ->Fill((*jetP)[i],(*jetEta)[i], evtWeight);
		}

		vector<Int_t> jetIndex;//matchted jet Index to order by  muon index
		jetIndex.clear();
		if(_Debug) cout<<"gmrSize_"<<gmrSize<<endl;
		if(_Debug) cout<<"gmrVectorSize_"<<gmrEta->size()<<endl;
		for(Int_t m =0; m<gmrSize; m++){
			jetIndex.push_back(-1);

			double glbpt = (*gmrPt)[m];   		
			if(glbpt <_MinMuonPt) continue;
		//	double glbp = (*gmrP)[m];
			double glbeta = (*gmrEta)[m]; 
			if(fabs(glbeta) > _MaxMuonEta ) continue;
		//	double glbphi = (*gmrPhi)[m];
		//	double glbcharge = (*gmrCharge)[m];
			double glbdxy = (*gmrDXY)[m];
			double glbiso03sumpt = (*gmrIso03sumPt)[m];

			if( _ifApplyImpactParCut &&  (fabs(glbdxy) > _MaxDxy) ) continue;
			if( _ifApplyIsolationCut && (glbiso03sumpt > _MaxIso03sumPt) )	continue;
			Bool_t passLikelihood1 = false;
			if( Use_Likelihood && _ifApplyLikelihoodRatioCut ) 
			{
				double llr1	= GetLikelihoodRatio( glbiso03sumpt, glbdxy, fDxyS, fDxyB, fIsoSumPtS, fIsoSumPtB );
				passLikelihood1 = (llr1 >= ResponseCut_Likelihood );
				if(!passLikelihood1) continue;
			}

			float minDeltaR = 999.;
			if(_Debug)cout<<"jetSize: "<<jetSize<<endl;
			for(Int_t i =0; i<jetSize; i++){

				//-------------------------------------------
				if((*jetPt)[i]<_MinJetPt) continue;
				if (fabs((*jetEta)[i]) > _MaxJetEta ) continue;
				//-------------------------------------------

				float deltaR = GetDeltaR((*gmrEta)[m],(*jetEta)[i], (*gmrPhi)[m],(*jetPhi)[i]);
				if (jetIndex[m] < 0 ){ //if the matched muon hasn't been filled
					jetIndex[m]=i;
					minDeltaR = deltaR;
				}
				else {// if it has filled a matched muon already.  then make a selection between previous muon and this muon
					Int_t j = jetIndex[m];
					//	  float deltaRj = GetDeltaR((*gmrEta)[m],(*jetEta)[j], (*gmrPhi)[m],(*jetPhi)[j]);

					if (minDeltaR>deltaR){
						minDeltaR=deltaR;
						//  if (deltaR <= deltaRj) {//matchPairs[irec].second=ijet;}
						jetIndex[m]=i;
					}
					else {//matchPairs[irec].second=j;
						jetIndex[m]=j;
					}
				}
			}//end jet size loop
			h_dRMuonJet ->Fill(minDeltaR, evtWeight);
			
			Int_t j = jetIndex[m];
			double jPt = (*jetPt)[j];
			double jP =  (*jetP)[j];
			double jEta =  (*jetEta)[j];
			double jPhi =  (*jetPhi)[j];

			h_muonPVsJetP->Fill(jP,(*gmrP)[m],evtWeight);
			h_muonEtaVsJetEta->Fill(jEta,(*gmrEta)[m],evtWeight);
			h_muonPhiVsJetPhi->Fill(jPhi,(*gmrPhi)[m],evtWeight);
			h_muonPtVsJetPt->Fill(jPt,(*gmrPt)[m],evtWeight);
			h_matchedJetPEta       ->Fill(jP,jEta,evtWeight);

		}//end gmr size loop
	}//end entries loop
}//end Loop() 
FakeMuonParameterizer::~FakeMuonParameterizer() {
	cout<<"~FakeMuonParameterizer()"<<endl;
	outputFile_->cd();
	WriteHists();
	outputFile_->Write() ;
	outputFile_->Close() ;
	delete outputFile_;
}
void FakeMuonParameterizer::BookingHists()
{
	h_allJetPEta = new TH2F("allJetPEta","allJetPEta;Jet P(GeV/c);Jet #eta",6000,0,6000,1100,-5.5,5.5);//pt vs eta
	h_matchedJetPEta = new TH2F("matchedJetPEta","matchedJetPEta;Jet P(GeV/c);Jet #eta",6000,0,6000,1100,-5.5,5.5);//pt vs eta
	h_muonPVsJetP = new TH2F("muonPVsJetP","muonPVsJetP",6000,0,6000,6000,0,6000);//xaxis is JetP, yaxis is MuonP
	h_muonEtaVsJetEta = new TH2F("muonEtaVsJetEta","muonEtaVsJetEta",1100,-5.5,5.5,600,-3.,3.);//xaxis is JetEta, yaxis is MuonEta
	h_muonPhiVsJetPhi = new TH2F("muonPhiVsJetPhi","muonPhiVsJetPhi",800,-4,4,800,-4,4);//xaxis is JetPhi, yaxis is MuonPhi
	h_muonPtVsJetPt = new TH2F("muonPtVsJetPt","muonPtVsJetPt",6000,0,6000,6000,0,6000);//xaxis is JetP, yaxis is MuonP
	h_dRMuonJet=new TH1F("dRMuonJet","dRMuonJet",10000,0,2);
	h_allJetPEta->Sumw2();
	h_matchedJetPEta->Sumw2();
	h_muonPVsJetP->Sumw2();
	h_muonEtaVsJetEta->Sumw2();
	h_muonPhiVsJetPhi->Sumw2();
	h_muonPtVsJetPt->Sumw2();
	h_dRMuonJet->Sumw2();

}
void FakeMuonParameterizer::WriteHists()
{
	h_allJetPEta->Write();
	h_matchedJetPEta->Write();
	h_muonPVsJetP->Write();
	h_muonEtaVsJetEta->Write();
	h_muonPhiVsJetPhi->Write();
	h_muonPtVsJetPt->Write();
	h_dRMuonJet->Write();
}
