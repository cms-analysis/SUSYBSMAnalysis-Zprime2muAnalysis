/*
 * =====================================================================================
 *
 *       Filename:  LLR.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/16/2008 02:25:35 PM CEST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  IHEP, Beijing
 *
 * =====================================================================================
 */
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "LikelihoodUtilities.C"
#include "NtupleUtilities.C"
#include "MyFitters.h"
#include <vector>
void LikelihoodCalibration(TString inFileNameS, TString inFileNameB,  TString outputFileName)
{
	gSystem->Load("libvectorDict.so");
//	TStyle *tdrStyle = gROOT->GetStyle("tdrStyle");
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	// Loose Cut 
	double MinMuonPt_ = 20.;
	double MaxMuonDxy_ = 0.1;
	double MaxMuonIsoSumPt_ = 20.;

	TFile *infileS = new TFile(inFileNameS,"READ");
	TFile *infileB = new TFile(inFileNameB,"READ");

	TTree* theTreeS = (TTree*)infileS->Get("myTree");
	TTree* theTreeB = (TTree*)infileB->Get("myTree");
	Float_t evtWeightS =0, evtWeightB=0;
	Int_t gmrSizeS =0, gmrSizeB=0;
	vector<Float_t> * gmrPtS=0, *gmrPtB=0;
	vector<Float_t> * gmrDXYS=0, *gmrDXYB=0;
	vector<Float_t> * gmrIso03sumPtS=0, *gmrIso03sumPtB=0;

	Int_t pytSizeS =0;
	vector<int> *pytStatusS=0;
	vector<int> *pytIdS=0;
	vector<int> *pytNMotherS=0;	

	theTreeS->SetBranchAddress( "evtWeight", &evtWeightS);
	theTreeS->SetBranchAddress( "gmrSize", &gmrSizeS);
	theTreeS->SetBranchAddress( "gmrIso03sumPt", &gmrIso03sumPtS);
	theTreeS->SetBranchAddress( "gmrDXY", &gmrDXYS);
	theTreeS->SetBranchAddress( "gmrPt", &gmrPtS);
	theTreeB->SetBranchAddress( "evtWeight", &evtWeightB );
	theTreeB->SetBranchAddress( "gmrSize", &gmrSizeB );
	theTreeB->SetBranchAddress( "gmrIso03sumPt", &gmrIso03sumPtB );
	theTreeB->SetBranchAddress( "gmrDXY", &gmrDXYB );
	theTreeB->SetBranchAddress( "gmrPt", &gmrPtB );

	theTreeS->SetBranchAddress( "pytSize", &pytSizeS);
	theTreeS->SetBranchAddress( "pytStatus", &pytStatusS);
	theTreeS->SetBranchAddress( "pytId", &pytIdS);
	theTreeS->SetBranchAddress( "pytNMother", &pytNMotherS);



	TH1F *hDxyB = new TH1F("hDxyB","Dxy from Dijets; Dxy [cm]",1000,-0.1,0.1);
	TF1 *fDxyB =new TF1("fDxyB",fitBkgDxy,-0.1,0.1,11);
	TF1 *fDxyS =new TF1("fDxyS",fitSigDxy,-0.1,0.1,3);
	TH1F *hDxyS = new TH1F("hDxyS","hDxyS",1000,-0.1,0.1);
	TF1 *fIsoSumPtB =new TF1("fIsoSumPtB",fitBkgSumPt,0,20,3);
	TH1F *hIsoSumPtB = new TH1F("hIsoSumPtB","Iso sumPt from Dijets; Dxy [cm]",200, 0, 20);
	TF1 *fIsoSumPtS =new TF1("fIsoSumPtS",fitSigSumPt,0,20,5);
	TH1F *hIsoSumPtS = new TH1F("hIsoSumPtS","Iso sumPt from Signal; Dxy [cm]",200, 0, 20);




	Int_t entriesS = theTreeS->GetEntries();

	entriesS = 1000000;

	for(Long64_t ievt=0; ievt<entriesS; ievt++)
	{
		if(ievt%10000==0)cout<<ievt<<endl;

		theTreeS->GetEntry(ievt);

		int process = GetProcess(pytSizeS, pytStatusS, pytIdS, pytNMotherS);
		int p = process/10;
		if(p!=2) continue;// W events

		for(int i=0; i<gmrSizeS; i++)
		{
			if( (*gmrPtS)[i] < MinMuonPt_ ) continue;
			if( (*gmrDXYS)[i] > MaxMuonDxy_ ) continue;
			if( (*gmrIso03sumPtS)[i] > MaxMuonIsoSumPt_ ) continue;

			hDxyS->Fill((*gmrDXYS)[i]);
			hIsoSumPtS->Fill((*gmrIso03sumPtS)[i]);
		}
	}
	Int_t entriesB = theTreeB->GetEntries();
	for(Long64_t ievt=0; ievt<entriesB; ievt++)
	{
		if(ievt%10000==0)cout<<ievt<<endl;

		theTreeB->GetEntry(ievt);
		for(int i=0; i<gmrSizeB; i++)
		{
			if( (*gmrPtB)[i] < MinMuonPt_ ) continue;
			if( (*gmrDXYB)[i] > MaxMuonDxy_ ) continue;
			if( (*gmrIso03sumPtB)[i] > MaxMuonIsoSumPt_ ) continue;

			hDxyB->Fill((*gmrDXYB)[i]);
			hIsoSumPtB->Fill((*gmrIso03sumPtB)[i]);
		}
	}

	hDxyB->Scale(1./hDxyB->Integral(0,1001));
	fDxyB->SetParameters(1.65,-0.00028,0.0045, 0.065,-0.0004,0.018,0.3, 1, 0., 0.2,0.1);// initialize with resonable values
	//fDxyB->SetParNames("MainGauss Const","Mean","Sigma","TailGauss Const","Mean","Sigma","Frac","Pol Const");
	hDxyB->Fit("fDxyB","VR");

	hDxyS->Scale(1./hDxyS->Integral(0,1001));
	fDxyS->SetParameters(hDxyS->GetMaximum(),hDxyS->GetMean(),hDxyS->GetRMS());// initialize with resonable values
	hDxyS->Fit("fDxyS","VR");

	hIsoSumPtB->Scale(1./hIsoSumPtB->Integral(0,1001));
	fIsoSumPtB->SetParameters(0.5,0.05, 0.00001);// initialize with resonable values
	fIsoSumPtB->SetParNames("Constant At 0","Constant");
	hIsoSumPtB->Fit("fIsoSumPtB","VR");

	hIsoSumPtS->Scale(1./hIsoSumPtS->Integral(0,1001));
	fIsoSumPtS->SetParameters(0.5,1,2,2,0.00001);// initialize with resonable values
	//	fIsoSumPtS->SetParNames("Constant At 0","Constant");
	hIsoSumPtS->Fit("fIsoSumPtS","VR");


	TH1F *h_llrS = new TH1F("h_llrS","h_llrS",1000,0,1);
	TH1F *h_llrB = new TH1F("h_llrB","h_llrB",1000,0,1);

	TH1F *h_SigEffVsLLRCut = new TH1F("h_SigEffVsLLRCut","h_SigEffVsLLRCut",1000,0,1);
	TH1F *h_BkgRejVsLLRCut = new TH1F("h_BkgRejVsLLRCut","h_BkgRejVsLLRCut",1000,0,1);



	for(Long64_t ievt=0; ievt<entriesS; ievt++)
	{
		if(ievt%10000==0)cout<<ievt<<endl;

		theTreeS->GetEntry(ievt);
		int process = GetProcess(pytSizeS, pytStatusS, pytIdS, pytNMotherS);
		int p = process/10;
		if(p!=2) continue;// W events
		for(int i=0; i<gmrSizeS; i++)
		{
			if( (*gmrPtS)[i] < MinMuonPt_ ) continue;
			if( (*gmrDXYS)[i] > MaxMuonDxy_ ) continue;
			if( (*gmrIso03sumPtS)[i] > MaxMuonIsoSumPt_ ) continue;
			double llrS = GetLikelihoodRatio( (*gmrIso03sumPtS)[i], (*gmrDXYS)[i], fDxyS, fDxyB, fIsoSumPtS, fIsoSumPtB );
			h_llrS->Fill(llrS, evtWeightS);
		}

	}
	for(Long64_t ievt=0; ievt<entriesB; ievt++)
	{
		if(ievt%10000==0)cout<<ievt<<endl;

		theTreeB->GetEntry(ievt);
		for(int i=0; i<gmrSizeB; i++)
		{
			if( (*gmrPtB)[i] < MinMuonPt_ ) continue;
			if( (*gmrDXYB)[i] > MaxMuonDxy_ ) continue;
			if( (*gmrIso03sumPtB)[i] > MaxMuonIsoSumPt_ ) continue;

			double llrB = GetLikelihoodRatio( (*gmrIso03sumPtB)[i], (*gmrDXYB)[i], fDxyS, fDxyB, fIsoSumPtS, fIsoSumPtB );
			h_llrB->Fill(llrB,evtWeightB);
		}
	}

	double SEff[1000],BRej[1000];
	double totS = h_llrS->Integral(0,1001);
	double totB = h_llrB->Integral(0,1001);
	for(int i=1; i<=1000; i++)	
	{
		SEff[i-1]=1-h_llrS->Integral(0,i)/totS;	
		BRej[i-1]=h_llrB->Integral(0,i)/totB;	
		h_SigEffVsLLRCut->SetBinContent(i,SEff[i-1]);
		h_BkgRejVsLLRCut->SetBinContent(i,BRej[i-1]);
	}

	TCanvas *c_bkgsig =  new TCanvas("c_bkgsig","c_bkgsig");
	c_bkgsig->SetTopMargin(0.08);
	c_bkgsig->SetLeftMargin(0.18);
	c_bkgsig->SetRightMargin(0.05);
	c_bkgsig->SetBottomMargin(0.2);

	TGraph *g_BkgRejVsSigEff = new TGraph(1000, SEff, BRej);
	g_BkgRejVsSigEff->SetName("g_BkgRejVsSigEff");
	g_BkgRejVsSigEff->SetTitle("g_BkgRejVsSigEff;Sig Efficiency;Bkg Rejection");
	TAxis *xaxis = g_BkgRejVsSigEff->GetXaxis();
	TAxis *yaxis = g_BkgRejVsSigEff->GetYaxis();

  xaxis->SetLabelFont  (    42);
  xaxis->SetLabelOffset( 0.015);
  xaxis->SetNdivisions (   505);
  //xaxis->SetTitle      (xtitle);
  xaxis->SetTitleColor (kBlack);
  xaxis->SetTitleFont  (    42);
  xaxis->SetTitleOffset(  1.25);
  xaxis->SetTitleSize  (  0.05);

  yaxis->CenterTitle   ( kTRUE);
  yaxis->SetLabelFont  (    42);
  yaxis->SetLabelOffset( 0.020);
  yaxis->SetNdivisions (   505);
  //yaxis->SetTitle      (ytitle);
  yaxis->SetTitleColor (kBlack);  
  yaxis->SetTitleFont  (    42);
  yaxis->SetTitleOffset(  1.70);
  yaxis->SetTitleSize  (  0.05);
	g_BkgRejVsSigEff->Draw("AC*");

	c_bkgsig->Update();


	//	TFile *output = new TFile("MyLLRTrain.root","RECREATE");
	TFile *output = new TFile(outputFileName,"RECREATE");

	output->cd();
	h_llrB->Write();
	h_llrS->Write();
	g_BkgRejVsSigEff->Write();
	c_bkgsig->Write();
	h_SigEffVsLLRCut->Write();
	h_BkgRejVsLLRCut->Write();

	fIsoSumPtB -> Write();
	fIsoSumPtS -> Write();
	fDxyB		-> Write();
	fDxyS		-> Write();

	hDxyS->Write();
	hDxyB->Write();
	hIsoSumPtS->Write();
	hIsoSumPtB->Write();

	output->Write();
	output->Close();
}















