#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include <iostream>
#include "TGraphAsymmErrors.h"
//#include "PlotUtils/Utils.C"

bool _debug = true;
const int nVariables =  8;
enum enumHistComparison {enumPhi=0, enumEta, enumP, enumPt, enumMass, enumDeltaR, enumDeltaEta, enumDeltaPhi} ;
TString s_hist[nVariables] = {"MuonPhi", "MuonEta", "MuonP", "MuonPt", "InvMass", "DeltaR", "DeltaEta", "DeltaPhi"};
Int_t binWidth[nVariables] = {  10	,		10,			50,		50,			50, 	10,			10,			10	  };
TString s_xtitle[nVariables] = {"#phi_{#mu}", "#eta_{#mu}", "P_{#mu} [GeV/c]", "p_{T#mu} [GeV/c]", "M_{#mu#mu} [GeV/c^{2}]", "#deltaR_{#mu#mu}",
	"#delta#eta_{#mu#mu}", "#delta#phi_{#mu#mu}"};
TString s_histFake[nVariables], s_histOrin[nVariables];
TCanvas *canCompare[nVariables]; 
TCanvas *canRatio[nVariables];
TH1F *histFake[nVariables];
TH1F *histOrin[nVariables];

TFile *infile;
TGraphAsymmErrors *gr_asy; //= new TGraphAsymmErrors(n, x, y, ex, ex, ey, ey);
void ComposeTwoHists(TCanvas *can, TH1F *hist1, TH1F *hist2, bool SetNormFactorTo1=false){
	can->cd();
	if(SetNormFactorTo1 == true)
	{	
		double tot1 = hist1->GetSumOfWeights();
		double tot2 = hist2->GetSumOfWeights();
		if(_debug)cout<<hist1->GetTitle()<<" tot1="<<tot1<<" tot2="<<tot2<<endl;
		if(tot1 && tot2) {
			hist1 -> Scale(1./tot1);
			hist2 -> Scale(1./tot2);
		}
	}

	double max = TMath::Max(hist1->GetMaximum(), hist2->GetMaximum());
	hist1 -> SetMaximum(max*1.2);
	hist2 -> SetMaximum(max*1.2);
	double min = TMath::Min(hist1->GetMinimum(), hist2->GetMinimum());
//	hist1 -> SetMinimum(min*1.2);
//	hist2 -> SetMinimum(min*1.2);
}

void DrawRatio(TCanvas*can, TH1F *hist1, TH1F* hist2, TString s_pre, TString xTitle)//1-obsv; 2-pred;
{
	TString s_title = s_pre+"ratio";
	can->cd();
	can->SetLogy(1);
	TH1F *h_ratio = new TH1F(s_title,s_title, hist1->GetNbinsX(), hist1->GetBinLowEdge(1), hist1->GetBinLowEdge(hist1->GetNbinsX()+1));
	h_ratio = (TH1F*) hist1->Clone(s_title);
	h_ratio->Divide(hist2);
	h_ratio->SetBinContent(h_ratio->GetNbinsX(), 0);
	h_ratio->SetBinError(h_ratio->GetNbinsX(), 0);
	char yTitle[100];
	sprintf(yTitle, "ratio= obsv/pred");
	axis1F(h_ratio, xTitle.Data(),yTitle);
	Draw(h_ratio, 0,1, "e same", kBlue,2, kBlue );
	TF1 *f_pol = new TF1(s_title,"[0]");
	f_pol->SetParNames("Constant");
	//f_pol->SetRange(50,7000);
	f_pol->SetLineColor(kBlack);
	f_pol->SetParameter(0,1.);
	f_pol->SetParLimits(0,0.1,20);
	h_ratio->Fit(s_title,"","same");
	char s_fitstatsphi[100];
	sprintf(s_fitstatsphi, "Constant = %.2f +/- %.2f", f_pol->GetParameter(0), f_pol->GetParError(0));
	SetPave(s_fitstatsphi, 0.4, 0.8, 0.9, 0.9, kWhite, kBlue);
}

bool LoadHist (TH1F **hist, const char *hname, 
		TFile **file, const char *fname ){

	*file = new TFile( fname );

	TFile *temp = *file;

	*hist = (TH1F*) temp -> Get(hname);

	if(!*hist) 
	{
		cout<<"histogram "<<hname<<" doesn't exist"<<endl;
		return 0;
	}
	else return 1;
}
void Rebin(TH1F* hOrinInvMass, TH1F*& obsvInvMClone, TH1F*& h_variablebins , vector<int>& vmass)
{
	//--------------------
	//the idea is to use the TH1 to calculate and store the bin error automatically
	// 
	// if input hOrinInvMass, binwidth is 50 GeV,  obsvInvMClone will be drawn clonely, and not affect the hOrinInvMass 
	// the make variablebins:
	// ie, from 0 to 300 GeV, binwidth is 50 GeV, 
	// from 300 to 400 GeV, binwidth is 100 GeV: 
	// 			1. erase the bincontent and errors of below 300 and above 400. 
	//			2. rebin all bins into one bin 
	//			3. get the bincontent and error, assign them as the bin from 300 to 400
	// do the same for other bins as previous stetp
	//	vector<int> vmass;
	vmass.clear();
	vmass.push_back(300);
	vmass.push_back(400);
	vmass.push_back(1000);
	//vmass.push_back(1500);
	vmass.push_back(7000);

	int a0 = hOrinInvMass->FindBin(vmass[0]);	
	for (int i =1; i<a0; i++)
	{
		obsvInvMClone->SetBinContent(i,hOrinInvMass->GetBinContent(i) );
		obsvInvMClone->SetBinError(i,hOrinInvMass->GetBinError(i) );
		h_variablebins->SetBinContent(i,hOrinInvMass->GetBinContent(i) );
		h_variablebins->SetBinError(i,hOrinInvMass->GetBinError(i) );
	}


	TH1F *h_tmp;
	float ave, err;
	for(int Ncombined = 1; Ncombined<(int)vmass.size(); Ncombined++)
	{
		h_tmp = (TH1F*)hOrinInvMass->Clone("h_tmp");
		int a2 = hOrinInvMass->FindBin(vmass[Ncombined]);
		int a1 = hOrinInvMass->FindBin(vmass[Ncombined-1]);
		for(int i=1; i<a1; i++)
		{
			h_tmp->SetBinContent(i, 0);
			h_tmp->SetBinError(i, 0);
		}
		for(int i=a2; i<=hOrinInvMass->GetNbinsX(); i++)
		{
			h_tmp->SetBinContent(i, 0);
			h_tmp->SetBinError(i, 0);
		}
		//		cout<<"debug1"<<endl;
		h_tmp->Rebin(hOrinInvMass->GetNbinsX());
		h_variablebins->SetBinContent(a2,h_tmp->GetBinContent(1));
		h_variablebins->SetBinError(a2,h_tmp->GetBinError(1));

		ave = h_tmp->GetBinContent(1) / (float)(a2-a1);
		err = sqrt(pow(h_tmp->GetBinError(1),2)/(float)(a2-a1));
		for ( int i=a1; i<a2; i++)
		{
			obsvInvMClone->SetBinContent(i,ave);
			obsvInvMClone->SetBinError(i,err );
		}
		delete h_tmp;
	}

}

void DrawMassRatioAsymetry(TH1F *hFakeInvMass, TH1F *hOrinInvMass)
{
	TH1F *obsvInvMClone= new TH1F("ObservedInvariantMass","Invariant Mass; Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Entries Per 50GeV/c^{2}",140,0,7000);
	TH1F *fakeInvMClone= new TH1F("PredictedInvariantMass","Invariant Mass; Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Entries Per 50GeV/c^{2}",140,0,7000);
	obsvInvMCloneVariable= new TH1F("ObservedInvariantMassVariable","Invariant Mass; Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Entries Per 50GeV/c^{2}",140,0,7000);
	fakeInvMCloneVariable= new TH1F("PredictedInvariantMassVariable","Invariant Mass; Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Entries Per 50GeV/c^{2}",140,0,7000);
	vector<int> vmass;
	Rebin(hOrinInvMass, obsvInvMClone, obsvInvMCloneVariable, vmass);
	Rebin(hFakeInvMass, fakeInvMClone, fakeInvMCloneVariable, vmass);

	TCanvas *c1 = new TCanvas("asy","asy");
	c1->cd();
	c1->SetLogy(1);
	c1->SetLogx(1);
	obsvInvMCloneVariable->Divide(fakeInvMCloneVariable);
	float massbinwidth = hOrinInvMass->GetBinWidth(0);
	int nbinsinitial = (int)vmass[0]/massbinwidth;
	const int n = 100; 
	int trueN = (int)vmass.size() + nbinsinitial-1;
	Float_t x[n], y[n], ex[n], ey[n];
	for(int i=0; i<100; i++)
	{
		x[i]=0; y[i]=0; ex[i]=0; ey[i]=0;
	}
	for(int i=0; i<trueN; i++)
	{
		if(i < nbinsinitial) 
		{
			x[i]= massbinwidth * (i+1) - massbinwidth/2.;
			ex[i]=massbinwidth/2.;
			y[i]=obsvInvMCloneVariable->GetBinContent(i+1);
			ey[i]=obsvInvMCloneVariable->GetBinError(i+1);
		}
		else 
		{
			int massh = vmass[i-nbinsinitial+1];
			int massl = vmass[i-nbinsinitial];
			x[i] = massl + (massh -massl)/2.;
			ex[i]= (massh-massl)/2.;
			y[i] = obsvInvMCloneVariable->GetBinContent(obsvInvMCloneVariable->FindBin(massh));
			ey[i] = obsvInvMCloneVariable->GetBinError(obsvInvMCloneVariable->FindBin(massh));

		}

	}

	gr_asy= new TGraphAsymmErrors(n, x, y, ex, ex, ey, ey);
	gr_asy->SetTitle(";M_{#mu#mu} [GeV/c^{2}]; ratio = observed / predicted");
	gr_asy->Draw("AP");

	c1->Print("plots/massratio_asy.gif");
	c1->Print("plots/massratio_asy.eps");

}

/*

   mass pull
   TH1F *h_pull = new TH1F("hi","(ratio-1)/#sigma_{ratio};mass;(ratio-1)/#sigma_{ratio}",140,0,7000);
   Double_t par0 = f_pol->GetParameter(0); //value of 1st parameter
   cout << "Constant="<<par0<<endl;
   for(int i=1; i<=140;i++)
   {
//cout<<i<<"th bin error="<<h_ratio->GetBinError(i)<<endl;
if(h_ratio->GetBinError(i)==0) continue;
Double_t tmp = (h_ratio->GetBinContent(i)-1)/h_ratio->GetBinError(i);
h_pull->SetBinContent(i,tmp);
}
axis1F(h_pull, "m_{#mu#mu} [GeV/c^{2}]","(ratio-1)/#sigma_{ratio}");
Draw(h_pull, 0,1, "e", kRed,2, kRed );
 */

void SetPave(Char_t* title,
		Float_t x1,
		Float_t y1,
		Float_t x2,
		Float_t y2,
		Int_t   fColor,
		Int_t   tColor)
{
	TPaveText *pv = new TPaveText(x1,y1,x2,y2,"ndc");

	pv->SetBorderSize(     0);
	pv->SetFillColor (fColor);
	pv->SetTextAlign (    22);
	pv->SetTextFont  (    42);
	pv->SetTextSize  (  0.04);
	pv->SetTextColor (tColor);
	pv->AddText      ( title);

	pv->Draw();
}

//--Main Funcion
void CompareEtaPhiMomInv( TString fileName="", int pthat = 22, TString s_pre="")
{
	TString s_pthatbin = ""; s_pthatbin+=pthat;
	if(pthat>20) s_pthatbin="";

	for(int i=0; i<nVariables; i++)
	{
		s_histFake[i] = "h_fake"+s_hist[i]+s_pthatbin;
		s_histOrin[i] = "h_obsv"+s_hist[i]+s_pthatbin;
		//	canCompare[i] = 
	}
	for(int i=0; i<nVariables; i++)
	{
		canCompare[i]= new TCanvas(s_hist[i], s_hist[i]);
		canCompare[i]->cd();
		canCompare[i]->SetLogy(1);
		LoadHist( &histFake[i], "Comparison/"+s_histFake[i], &infile, fileName);
		LoadHist( &histOrin[i], "Comparison/"+s_histOrin[i], &infile, fileName);
		histFake[i]->Rebin(binWidth[i]);
		histOrin[i]->Rebin(binWidth[i]);
		ComposeTwoHists(canCompare[i], histFake[i], histOrin[i], 0);
		char yTitle[100];
		sprintf(yTitle, "Entries / %0.2f ", histFake[i]->GetBinWidth(0));
		axis1F(histFake[i], s_xtitle[i].Data(), yTitle);
		Draw(histFake[i], 0,kRed, "same", kRed,2, kRed );
		Draw(histOrin[i], kBlue,1007, "same", kBlue,2, kBlue );
		Draw(histFake[i], 0,kRed, "same", kRed,2, kRed );

		canCompare[i]->Print("plots/"+s_pre+s_hist[i]+s_pthatbin+".gif");
		canCompare[i]->Print("plots/"+s_pre+s_hist[i]+s_pthatbin+".eps");

		canRatio[i] = new TCanvas(s_hist[i]+"_ratio", s_hist[i]+"_ratio");
		if(i==enumMass) canRatio[i]->SetLogx(1);
		DrawRatio(canRatio[i], histOrin[i], histFake[i], s_hist[i], s_xtitle[i]);
		canRatio[i]->Print("plots/"+s_pre+s_hist[i]+s_pthatbin+"_ratio.gif");
		canRatio[i]->Print("plots/"+s_pre+s_hist[i]+s_pthatbin+"_ratio.eps");

		if(i==enumMass) DrawMassRatioAsymetry(histFake[i],histOrin[i]);


	}

}
//-----FIXME how to automatically rebin with desired different binwidths
