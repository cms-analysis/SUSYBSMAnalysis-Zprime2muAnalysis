#define FakeMuonAnalyzer_cxx
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

#include "LikelihoodUtilities.C"
#include "FakeMuonAnalyzer.h"
#include "FakeMuonGenerator.h"
#include "NtupleUtilities.C"
#include "Configure.C"

const double MUMASS = 0.10566;
TF1* fDxyS(0), *fDxyB(0), *fIsoSumPtS(0), *fIsoSumPtB(0);
double ResponseCut_Likelihood  = 0;
void FakeMuonAnalyzer::Initialize(TString s_input)
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
}
//void FakeMuonAnalyzer::Loop()
void FakeMuonAnalyzer::Loop(TString s_output, Long64_t entryToStart , Long64_t entriesToAnalyze)
{
	if (fChain == 0) return;

	//-- Booking hists and creating file
	if(s_output=="") s_output = _OutputFileComparison;
	outputFile_  = new TFile( s_output, "RECREATE" ) ;
	//outputFile_  = new TFile( _OutputFileComparison, "RECREATE" ) ;
	outputFile_->mkdir("Comparison");
	outputFile_->cd("Comparison");
	BookingFakeHists();
	BookingObsvHists();

	//-- if apply likelihood cut, then intialize it
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

	//-- FakeMuonGenerator initilizing
	myGenerator = new FakeMuonGenerator(_FakeMuonGeneratorSeed);
	myGenerator->Initialize(_InputFileParameterizedHistos);

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

		// -  - - - -half event to analyze
		if((_TakeHalfEvents==1) && (jentry%2!=0)) continue;
		if((_TakeHalfEvents==2) && (jentry%2==0)) continue;

		// - - - do we want to take into accound evt weight ?
		if( _ifTakeIntoAccountWeight_Analyzing == 0 )evtWeight=1;

		//- - - the dijets from soup ? 
		int process = GetProcess(pytSize, pytStatus, pytId, pytNMother);
		bool isDijetsEvent = (process==1);
		if( _ifRunOnGumboSoup && isDijetsEvent == false) continue;

		// - - - run on W+jets sample ?
		bool isWjetsEvent = (process/10 == 2);

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

		// count the jets above 20 GeV/c
		Int_t numOfJetsWithPtLargerThan20GeV =0;
		for(int i =0; i<jetSize; i++){
			double jetPt_ = (*jetPt)[i];
			if (jetPt_<_MinJetPt) continue;
			if (fabs((*jetEta)[i]) > _MaxJetEta ) continue;
			numOfJetsWithPtLargerThan20GeV++;
		}
		if (numOfJetsWithPtLargerThan20GeV< _MinNumberOfJetsAboveThreshold) continue;

		Double_t EventWeight = evtWeight;
		GetPredictedDistributions(jentry, EventWeight);	
		GetObservedDistributions(jentry, EventWeight);

	}//end entries loop

	cout<<"Print weights:"<<endl;
	for(int i=0; i<(int)vweight.size();i++)
	{
		cout<<i <<" th weight= "<<vweight[i]<<endl;
	}
}//end Loop() 
void FakeMuonAnalyzer::GetPredictedDistributions( Long64_t jentry, double EventWeight)
{
	if(_Debug)cout<<"Producing fake muons...."<<endl;
	vector<FakeMuonGenerator::fakeMuon> blahs;
	blahs.clear();

	fChain->GetEntry(jentry);
	for(int i = 0; i<jetSize; i++)
	{
		float tmpJetP 	=(*jetP)[i];
		float tmpJetPt 	=(*jetPt)[i]; 
		float tmpJetEta =(*jetEta)[i];
		float tmpJetPhi =(*jetPhi)[i];

		if( tmpJetPt< _MinJetPt ) continue;
		if( fabs(tmpJetEta) > _MaxJetEta ) continue;

		if( tmpJetP   > MaxJetPInRefHist_ ) tmpJetP = MaxJetPInRefHist_;
		if( tmpJetEta > MaxJetEtaInRefHist_ ) tmpJetEta = MaxJetEtaInRefHist_;
		if( tmpJetEta < MinJetEtaInRefHist_ ) tmpJetEta = MinJetEtaInRefHist_;

		FakeMuonGenerator::fakeMuon blah;
		myGenerator->GenerateFakeMuon( tmpJetP,tmpJetPt, tmpJetEta, tmpJetPhi, blah );
		blahs.push_back(blah);
	}
	if(_Debug)cout<<blahs.size()<<" fake muons produced."<<endl;

	FillPredictedHists(blahs, EventWeight);
}
void FakeMuonAnalyzer::FillPredictedHists(vector<FakeMuonGenerator::fakeMuon> vfakemu, double EventWeight)
{
	for(int m1=0; m1<(int)vfakemu.size(); m1++)
	{
		float eta1=vfakemu[m1].eta;
		float phi1=vfakemu[m1].phi;
		float p1=vfakemu[m1].p;
		float pt1 = vfakemu[m1].pt;
		double ww = vfakemu[m1].weight * EventWeight ;
		h_fakeMuonP  ->Fill(p1, ww);
		h_fakeMuonEta->Fill(eta1,ww);
		h_fakeMuonPhi->Fill(phi1,ww);
		h_fakeMuonPt  ->Fill(pt1, ww);

		/*
		   float cphi1=TMath::Cos(phi1);
		   float sphi1=TMath::Sin(phi1);
		   float taz1=1./TMath::Tan(2.0*TMath::ATan(TMath::Exp(-eta1)));
		   float den1= TMath::Sqrt(cphi1*cphi1+sphi1*sphi1+taz1*taz1);

		   double Cpt1=0;
		   double Cp1=0;
		   if(den1!=0)
		   {
		   Cpt1=p1/den1;
		   Cp1 = pt1*den1;
		   double ww = vfakemu[m1].weight * EventWeight;
		   h_fakeMuonCPt  ->Fill(Cpt1,ww );
		   h_fakeMuonCP  ->Fill(Cp1,ww );
		   }
		   if(fakeNumW==0) fakeNumW=vfakemu[m1].weight;
		   else	fakeNumW*=vfakemu[m1].weight;
		 */			for(int m2=m1+1; m2<(int)vfakemu.size(); m2++)
		{
			if (_OppositeSign && vfakemu[m1].charge == vfakemu[m2].charge) continue;
			if (!_OppositeSign && vfakemu[m1].charge != vfakemu[m2].charge) continue;
			TLorentzVector vmu1, vmu2, vmm; // LorentzVectors for mu+ mu- 

			float eta2=vfakemu[m2].eta;
			float phi2=vfakemu[m2].phi;
			//float p2=vfakemu[m2].p;
			float pt2=vfakemu[m2].pt;
			/*
			   float cphi2=TMath::Cos(phi2);
			   float sphi2=TMath::Sin(phi2);
			   float taz2=1./TMath::Tan(2.0*TMath::ATan(TMath::Exp(-eta2)));
			   float den2= TMath::Sqrt(cphi2*cphi2+sphi2*sphi2+taz2*taz2);

			   double Cpt2=0;
			//	double Cp2=0;
			if(den2!=0)		Cpt2=p2/den2;
			 */
			vmu1.SetPtEtaPhiM(pt1,eta1,phi1,MUMASS);
			vmu2.SetPtEtaPhiM(pt2,eta2,phi2,MUMASS);
			vmm = vmu1 + vmu2;		
			//float w = (vfakemu[m2].weight * vfakemu[m1].weight);
			float w = (vfakemu[m2].weight * vfakemu[m1].weight) * EventWeight;

			vmu1.SetPtEtaPhiM(pt1,eta1,phi1,MUMASS);
			vmu2.SetPtEtaPhiM(pt2,eta2,phi2,MUMASS);
			vmm = vmu1 + vmu2;
			float tmpDR = GetDeltaR(eta1,eta2,phi1,phi2);
			if( tmpDR < _DeltaRCutForTwoMuons ) continue;
			h_fakeDeltaR->Fill(tmpDR,w);
			h_fakeDeltaEta->Fill(fabs(eta1-eta2),w);
			h_fakeDeltaPhi->Fill(DeltaPhi(phi1,phi2),w);
			h_fakeInvMass->Fill(vmm.M(), w);
			//	h_fakeInvMassC->Fill(vmm.M(), w);
		}

	}
}
void FakeMuonAnalyzer::GetObservedDistributions(Long64_t jentry, double EventWeight)
{
	if(_Debug)cout<<"observing..."<<endl;
	fChain->GetEntry(jentry);

	vector<Int_t> jetIndex;//matchted jet Index to order by  muon index
	jetIndex.clear();
	for(int m= 0; m<gmrSize; m++)
	{
		double glbpt = (*gmrPt)[m];   		
		double glbeta =(*gmrEta)[m]; 
		double glbphi =(*gmrPhi)[m]; 

		jetIndex.push_back(-1);

		//-------------------------------------------
		if(glbpt <_MinMuonPt) continue;
		if(fabs(glbeta) > _MaxMuonEta ) continue;
		//-------------------------------------------

		float minDeltaR = 999.;
		for(int i = 0; i<jetSize; i++)
		{
			double jetpt =  (*jetPt)[i];
			double jeteta = (*jetEta)[i];
			double jetphi = (*jetPhi)[i];

			//-------------------------------------------
			if(jetpt < _MinJetPt ) continue;
			if(fabs(jeteta) > _MaxJetEta ) continue;
			//-------------------------------------------

			float dR = GetDeltaR(glbeta,jeteta, glbphi, jetphi);
			if (jetIndex[m] < 0 ){ //if the matched muon hasn't been filled
				jetIndex[m]=i;
				minDeltaR = dR;
			}
			else {// if it has filled a matched muon already.  then make a selection between previous muon and this muon
				Int_t j = jetIndex[m];
				//	  float deltaRj = GetDeltaR((*gmrEta)[m],(*jetEta)[j], (*gmrPhi)[m],(*jetPhi)[j]);

				if (minDeltaR>dR){
					minDeltaR=dR;
					//  if (deltaR <= deltaRj) {//matchPairs[irec].second=ijet;}
					jetIndex[m]=i;
				}
				else {//matchPairs[irec].second=j;
					jetIndex[m]=j;
				}
			}
		}//end jet size loop
	}


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
		Bool_t passLikelihood1 = false;
		if( Use_Likelihood && _ifApplyLikelihoodRatioCut ) 
		{
			double llr1	= GetLikelihoodRatio( glbiso03sumpt, glbdxy, fDxyS, fDxyB, fIsoSumPtS, fIsoSumPtB );
			passLikelihood1 = (llr1 >= ResponseCut_Likelihood );
			if(!passLikelihood1) continue;
		}
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
			Bool_t passLikelihood2 = false;
			if( Use_Likelihood && _ifApplyLikelihoodRatioCut ) 
			{
				double llr2	= GetLikelihoodRatio( glbiso03sumpt2, glbdxy2, fDxyS, fDxyB, fIsoSumPtS, fIsoSumPtB );
				passLikelihood2 = (llr2 >= ResponseCut_Likelihood );
				if(!passLikelihood2) continue;
			}
			// - - - - - - - - - -
			if(_SkipMuonsWhosePAboveThreshold && glbp2> _MaxMuonP ) continue;

			if(_OppositeSign && glbcharge ==glbcharge2) continue; //skip it with same sign muons
			if(!_OppositeSign && glbcharge !=glbcharge2) continue; //

			if(_ifSkipTwoMuonsAssociatedWithTheSameJet && jetIndex[m1]==jetIndex[m2]) continue;

			TLorentzVector vmu1, vmu2, vmm; // LorentzVectors for mu+ mu- 
			vmu1.SetPtEtaPhiM(glbpt, glbeta, glbphi, MUMASS);
			vmu2.SetPtEtaPhiM(glbpt2, glbeta2, glbphi2, MUMASS);
			vmm = vmu1 + vmu2;		

			float tmpDR = GetDeltaR(glbeta, glbeta2, glbphi, glbphi2);
			if(tmpDR < _DeltaRCutForTwoMuons) continue;

			h_obsvDeltaR->Fill(tmpDR,EventWeight);
			h_obsvInvMass->Fill(vmm.M(),EventWeight);
			h_obsvDeltaEta->Fill(fabs(glbeta-glbeta2),EventWeight);
			h_obsvDeltaPhi->Fill(DeltaPhi(glbphi,glbphi2),EventWeight);
		}
	}
}

FakeMuonAnalyzer::~FakeMuonAnalyzer() {
	outputFile_->cd("Comparison");
	WriteFakeHists();
	WriteObsvHists();
	outputFile_->Write() ;
	outputFile_->Close() ;
	delete outputFile_;
	if(myGenerator) delete myGenerator;
}
void FakeMuonAnalyzer::BookingFakeHists()
{
	h_fakeDeltaR = new TH1F("h_fakeDeltaR","h_fakeDeltaR",628,0,6.28);	
	h_fakeDeltaEta = new TH1F("h_fakeDeltaEta","h_fakeDeltaEta;|#{Delta}#{eta}|",628,0,6.28);	
	h_fakeDeltaPhi = new TH1F("h_fakeDeltaPhi","h_fakeDeltaPhi",800,-4,4);	
	h_fakeInvMass = new TH1F("h_fakeInvMass","h_fakeInvMass",7000,0,7000);
	h_fakeMuonP = new TH1F("h_fakeMuonP","h_fakeMuonP",7000,0,7000);
	h_fakeMuonPt = new TH1F("h_fakeMuonPt","h_fakeMuonPt",7000,0,7000);
	h_fakeMuonEta = new TH1F("h_fakeMuonEta","h_fakeMuonEta",800,-4,4);
	h_fakeMuonPhi = new TH1F("h_fakeMuonPhi","h_fakeMuonPhi",800,-4,4);
	h_fakeMuonP->Sumw2();
	h_fakeMuonPt->Sumw2();
	h_fakeMuonEta->Sumw2();
	h_fakeMuonPhi->Sumw2();
	h_fakeDeltaR->Sumw2();
	h_fakeDeltaEta->Sumw2();
	h_fakeDeltaPhi->Sumw2();
	h_fakeInvMass->Sumw2();
}
void FakeMuonAnalyzer::BookingObsvHists()
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
	h_obsvMuonP->Sumw2();
	h_obsvMuonPt->Sumw2();
	h_obsvMuonEta->Sumw2();
	h_obsvMuonPhi->Sumw2();
	h_obsvDeltaR->Sumw2();
	h_obsvDeltaEta->Sumw2();
	h_obsvDeltaPhi->Sumw2();
	h_obsvInvMass->Sumw2();
	//	h_obsvInvMassPassLooseCut->Sumw2();
	//	h_obsvInvMassPassLikelihoodRatioCut->Sumw2();
}
void FakeMuonAnalyzer::WriteFakeHists()
{
	h_fakeMuonP->Write();
	h_fakeMuonPt->Write();
	h_fakeMuonEta->Write();
	h_fakeMuonPhi->Write();
	h_fakeDeltaR->Write();
	h_fakeDeltaEta->Write();
	h_fakeDeltaPhi->Write();
	h_fakeInvMass->Write();
}
void FakeMuonAnalyzer::WriteObsvHists()
{
	h_obsvMuonP->Write();
	h_obsvMuonPt->Write();
	h_obsvMuonEta->Write();
	h_obsvMuonPhi->Write();
	h_obsvDeltaR->Write();
	h_obsvDeltaEta->Write();
	h_obsvDeltaPhi->Write();
	h_obsvInvMass->Write();
	//	h_obsvInvMassPassLooseCut->Write();
	//	h_obsvInvMassPassLikelihoodRatioCut->Write();
}
