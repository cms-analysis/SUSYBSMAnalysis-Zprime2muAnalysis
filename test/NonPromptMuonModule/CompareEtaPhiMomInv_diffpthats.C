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
//#include "PlotUtils/Utils.C"
bool _UsingDerivedPt=false;//got P from distribution and pt from formula

TFile *file1;
TFile *file2;

TH1F  *hFakePhi, *hFakeEta, *hFakeP, *hFakePt, *hFakeInvMass;
TH1F  *hOrinPhi, *hOrinEta, *hOrinP, *hOrinPt, *hOrinInvMass;

//TH1F *h_eta_ratio;

TH1F *obsvInvMClone;
TH1F *fakeInvMClone;
TH1F *h_ratio;
TH1F *obsvDeltaR;
TH1F *fakeDeltaR;
TH1F *obsvDeltaEta;
TH1F *fakeDeltaEta;
TH1F *obsvDeltaPhi;
TH1F *fakeDeltaPhi;
TGraphAsymmErrors *gr_asy; //= new TGraphAsymmErrors(n, x, y, ex, ex, ey, ey);
void SetMaximum(TH1F* h1, TH1F* h2)
{
	double max = TMath::Max(h1->GetMaximum(), h2->GetMaximum());
	h1->SetMaximum(max*1.2);
	h2->SetMaximum(max*1.2);
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

void CompareEtaPhiMomInv( TString fileName="", int pthat = 20, TString s_pre=""){
	//	gInterpreter->ExecuteMacro("PlotUtils/FloridaStyle.C");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	TString s_pthatbin = ""; s_pthatbin+=pthat;
	if(pthat>20) s_pthatbin="";

	const int NCanvas = 19;
	TCanvas *c1[NCanvas];
	TString canTitle[NCanvas]={"phi","eta","p","pt","Mass","p_error_ratio","mass_ratio","ratioPull","deltaR", "deltaPhi", "deltaEta", "mass_graph_asymetricerror", "phi_ratio", "phi_error_ratio", "eta_ratio", "eta_error_ratio", "deltaEta_ratio","p_ratio", "deltaPhi_ratio"};
	for(int i =0; i<NCanvas; i++)
	{
		c1[i] = new TCanvas(canTitle[i],canTitle[i]);
	}

	c1[0]->cd();
	c1[0]->SetLogy(0);
	LoadHist( &hFakePhi, "Comparison/h_fakeMuonPhi"+s_pthatbin, &file1,fileName);
	LoadHist( &hOrinPhi, "Comparison/h_obsvMuonPhi"+s_pthatbin, &file2,fileName);
	hOrinPhi ->Rebin(10);
	hFakePhi ->Rebin(10);
	TH1F *h_OrinPhiClone = (TH1F*) hOrinPhi->Clone("h_OrinPhiClone");
	TH1F *h_FakePhiClone = (TH1F*) hFakePhi->Clone("h_FakePhiClone");
	SetMaximum(hFakePhi, hOrinPhi);
	char yTitle[100];
	sprintf(yTitle, "Entries / %0.2f ", hFakePhi->GetBinWidth(0));
	axis1F(hFakePhi, "#phi_{#mu}", yTitle);
	Draw(hFakePhi, 0,1, "", kRed,2, kRed );
	Draw(hOrinPhi, kBlue,1007, "e same", kBlue,2, kBlue );

	c1[12]->cd();
	c1[12]->SetLogy();
	TH1F *h_phi_ratio = new TH1F("h_phi_ratio","h_phi_ratio", h_OrinPhiClone->GetNbinsX(), h_OrinPhiClone->GetBinLowEdge(1), h_OrinPhiClone->GetBinLowEdge(h_OrinPhiClone->GetNbinsX()+1));
	h_phi_ratio = (TH1F*) h_OrinPhiClone->Clone("h_phi_ratio");
	h_phi_ratio->Divide(h_FakePhiClone);
	sprintf(yTitle, "ratio= #Phi_{obsv}/#Phi_{pred} ");
	axis1F(h_phi_ratio, "#phi_{#mu}",yTitle);
	Draw(h_phi_ratio, 0,1, "e same", kBlue,2, kBlue );
	TF1 *f_pol_phi = new TF1("f_pol_phi","[0]");
	f_pol_phi->SetParNames("Constant");
	//f_pol_phi->SetRange(50,7000);
	f_pol_phi->SetLineColor(kBlack);
	f_pol_phi->SetParameter(0,1.);
	f_pol_phi->SetParLimits(0,0.1,20);
	h_phi_ratio->Fit("f_pol_phi","","same");
	char s_fitstatsphi[100];
	sprintf(s_fitstatsphi, "Constant = %.2f +/- %.2f", f_pol_phi->GetParameter(0), f_pol_phi->GetParError(0));
	SetPave(s_fitstatsphi, 0.4, 0.8, 0.9, 0.9, kWhite, kBlue);

	c1[13]->cd();
	TH1F *h_phierror_ratio = new TH1F("h_phierror_ratio","h_phierror_ratio", h_OrinPhiClone->GetNbinsX(),h_OrinPhiClone->GetBinLowEdge(1), h_OrinPhiClone->GetBinLowEdge(h_OrinPhiClone->GetNbinsX()+1));
	for(int i=1; i<hOrinPhi->GetNbinsX(); i++)
	{
		float phi_error_orin = h_OrinPhiClone->GetBinError(i);
		float phi_error_fake = h_FakePhiClone->GetBinError(i);
		float ratio = 0;
		if(phi_error_fake!=0) ratio = phi_error_orin/phi_error_fake;
		h_phierror_ratio->SetBinContent(i, ratio);
	}
	sprintf(yTitle, "ratio= #sigma_{Phi}(obsv)/#sigma_{Phi}(pred) ");
	axis1F(h_phierror_ratio, "#phi_{#mu}",yTitle);
	Draw(h_phierror_ratio, kBlue,1007, "e same", kBlue,2, kBlue );

	c1[1]->cd();
	LoadHist( &hFakeEta, "Comparison/h_fakeMuonEta"+s_pthatbin, &file1, fileName);
	LoadHist( &hOrinEta, "Comparison/h_obsvMuonEta"+s_pthatbin, &file2, fileName);
	hOrinEta ->Rebin(10);
	hFakeEta ->Rebin(10);
	TH1F *h_OrinEtaClone = (TH1F*) hOrinEta->Clone("h_OrinEtaClone");
	TH1F *h_FakeEtaClone = (TH1F*) hFakeEta->Clone("h_FakeEtaClone");
	SetMaximum(hOrinEta, hFakeEta);
	sprintf(yTitle, "Entries / %0.2f ", hFakeEta->GetBinWidth(0));
	axis1F(hFakeEta, "#eta_{#mu}", yTitle);
	Draw(hFakeEta, 0,1, "", kRed,2, kRed );
	Draw(hOrinEta, kBlue,1007, "e same", kBlue,2, kBlue );


	c1[14]->cd();
	c1[14]->SetLogy();
	TH1F *h_eta_ratio = new TH1F("h_eta_ratio","h_eta_ratio", h_OrinEtaClone->GetNbinsX(), h_OrinEtaClone->GetBinLowEdge(1), h_OrinEtaClone->GetBinLowEdge(h_OrinEtaClone->GetNbinsX()+1));
	h_eta_ratio    = (TH1F*) h_OrinEtaClone->Clone("h_eta_ratio");
	h_eta_ratio -> Divide(h_FakeEtaClone); 
	sprintf(yTitle, "ratio= #eta_{obsv}/#eta_{pred} ");
	axis1F(h_eta_ratio, "#eta_{#mu}",yTitle);
	Draw(h_eta_ratio, 0, 1, "e same", kBlue,2, kBlue );
	TF1 *f_pol_eta = new TF1("f_pol_eta","[0]");
	f_pol_eta->SetParNames("Constant");
	f_pol_eta->SetLineColor(kBlack);
	f_pol_eta->SetParameter(0,1.);
	f_pol_eta->SetParLimits(0,0.1,20);
	h_eta_ratio->Fit("f_pol_eta","","same");
	char s_fitstats[100];
	sprintf(s_fitstats, "Constant = %.2f +/- %.2f", f_pol_eta->GetParameter(0), f_pol_eta->GetParError(0));
	SetPave(s_fitstats, 0.4, 0.8, 0.9, 0.9, kWhite, kBlue);



	c1[15]->cd();
	TH1F *h_etaerror_ratio = new TH1F("h_etaerror_ratio","h_etaerror_ratio", h_OrinEtaClone->GetNbinsX(),h_OrinEtaClone->GetBinLowEdge(1), h_OrinEtaClone->GetBinLowEdge(h_OrinEtaClone->GetNbinsX()+1));
	for(int i=1; i<hOrinEta->GetNbinsX(); i++)
	{
		float eta_error_orin = h_OrinEtaClone->GetBinError(i);
		float eta_error_fake = h_FakeEtaClone->GetBinError(i);
		float ratio = 0;
		if(eta_error_fake!=0) ratio = eta_error_orin/eta_error_fake;
		h_etaerror_ratio->SetBinContent(i, ratio);
	}
	sprintf(yTitle, "ratio= #sigma_{Eta}(obsv)/#sigma_{Eta}(pred) ");
	axis1F(h_etaerror_ratio, "#eta_{#mu}",yTitle);
	Draw(h_etaerror_ratio, kBlue,1007, "e same", kBlue,2, kBlue );


	c1[2]->cd();
	c1[2]->SetLogy(1);
	if(_UsingDerivedPt==true)
	{
		LoadHist( &hFakeP, "Comparison/h_fakeMuonP"+s_pthatbin, &file1, fileName);
	}
	else
	{
		LoadHist( &hFakeP, "Comparison/h_fakeMuonP"+s_pthatbin, &file1, fileName);
	}
	LoadHist( &hOrinP, "Comparison/h_obsvMuonP"+s_pthatbin, &file2, fileName);
	hFakeP ->Rebin(20);
	hOrinP ->Rebin(20);
	TH1F *h_OrinPClone = (TH1F*) hOrinP->Clone("h_OrinPClone");
	TH1F *h_FakePClone = (TH1F*) hFakeP->Clone("h_FakePClone");
	SetMaximum(hFakeP, hOrinP);
	sprintf(yTitle, "Entries / %0.f GeV/c", hFakeP->GetBinWidth(0));
	axis1F(hOrinP, "P_{#mu} [GeV/c]",yTitle);
	Draw(hOrinP, kBlue,1007, "e same", kBlue,2, kBlue );
	Draw(hFakeP, 0,1, "e same", kRed,2, kRed );

	c1[17]->cd();
	c1[17]->SetLogy(1);
	TH1F *h_p_ratio ;//= new TH1F("h_p_ratio","h_p_ratio", h_OrinPClone->GetNbinsX(),0,7000);
	h_p_ratio = (TH1F*)h_OrinPClone->Clone("h_p_ratio");
	h_p_ratio->Divide(h_FakePClone);
	sprintf(yTitle, "ratio= P_{obsv}/P_{pred} ");
	axis1F(h_p_ratio, "P_{#mu} [GeV/c]",yTitle);
	Draw(h_p_ratio, kBlue,1007, "e same", kBlue,2, kBlue );
	TF1 *f_pol_p = new TF1("f_pol_p","[0]");
	f_pol_p->SetParNames("Constant");
	f_pol_p->SetLineColor(kBlack);
	f_pol_p->SetParameter(0,1.);
	f_pol_p->SetParLimits(0,0.1,20);
	h_p_ratio->Fit("f_pol_p","","same");
	sprintf(s_fitstats, "Constant = %.2f +/- %.2f", f_pol_p->GetParameter(0), f_pol_p->GetParError(0));
	SetPave(s_fitstats, 0.4, 0.8, 0.9, 0.9, kWhite, kBlue);

	c1[5]->cd();
	TH1F *h_perror_ratio = new TH1F("h_perror_ratio","h_perror_ratio", h_OrinPClone->GetNbinsX(),0,7000);
	for(int i=1; i<hOrinP->GetNbinsX(); i++)
	{
		float p_error_orin = h_OrinPClone->GetBinError(i);
		float p_error_fake = h_FakePClone->GetBinError(i);
		float ratio = 0;
		if(p_error_fake!=0) ratio = p_error_orin/p_error_fake;
		h_perror_ratio->SetBinContent(i, ratio);
	}
	sprintf(yTitle, "ratio= #sigma_{P}(obsv)/#sigma_{P}(pred) ");
	axis1F(h_perror_ratio, "P_{#mu} [GeV/c]",yTitle);
	Draw(h_perror_ratio, kBlue,1007, "e same", kBlue,2, kBlue );



	c1[3]->cd();
	c1[3]->SetLogy(1);
	if(_UsingDerivedPt==true)
	{
		LoadHist( &hFakePt, "Comparison/h_fakeMuonPt"+s_pthatbin, &file1, fileName);
	}
	else
	{
		LoadHist( &hFakePt, "Comparison/h_fakeMuonPt"+s_pthatbin,  &file1,  fileName);
	}
	LoadHist( &hOrinPt, "Comparison/h_obsvMuonPt"+s_pthatbin, &file2, fileName);
	hOrinPt->Rebin(50);
	hFakePt->Rebin(50);
	SetMaximum(hOrinPt, hFakePt);
	sprintf(yTitle, "Entries / %0.f GeV/c", hFakePt->GetBinWidth(0));
	axis1F(hFakePt, "p_{T#mu} [GeV/c]",yTitle);
	Draw(hFakePt, 0,1, "", kRed,2, kRed );
	Draw(hOrinPt, kBlue,1007, "e same", kBlue,2, kBlue );


	c1[4]->cd();
	c1[4]->SetLogy(1);
	LoadHist( &hFakeInvMass, "Comparison/h_fakeInvMass"+s_pthatbin, &file1, fileName);
	LoadHist( &hOrinInvMass, "Comparison/h_obsvInvMass"+s_pthatbin, &file2, fileName);
	hFakeInvMass-> Rebin(50);
	hOrinInvMass-> Rebin(50);
	SetMaximum(hFakeInvMass, hOrinInvMass);
	Double_t max = TMath::Max(hFakeInvMass->GetMaximum(), hOrinInvMass->GetMaximum());
	hFakeInvMass -> SetMaximum(max*1.2);
	sprintf(yTitle, "Entries / %0.f GeV/c^{2}", hFakeInvMass->GetBinWidth(0));
	axis1F(hFakeInvMass, "m_{#mu#mu} [GeV/c^{2}]",yTitle);
	Draw(hFakeInvMass, 0,1, "", kRed,2, kRed );
	Draw(hOrinInvMass, kBlue,1007, "e same", kBlue,2, kBlue );



	obsvInvMClone= new TH1F("ObservedInvariantMass","Invariant Mass; Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Entries Per 50GeV/c^{2}",140,0,7000);
	fakeInvMClone= new TH1F("PredictedInvariantMass","Invariant Mass; Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Entries Per 50GeV/c^{2}",140,0,7000);
	obsvInvMCloneVariable= new TH1F("ObservedInvariantMassVariable","Invariant Mass; Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Entries Per 50GeV/c^{2}",140,0,7000);
	fakeInvMCloneVariable= new TH1F("PredictedInvariantMassVariable","Invariant Mass; Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Entries Per 50GeV/c^{2}",140,0,7000);
	vector<int> vmass;
	Rebin(hOrinInvMass, obsvInvMClone, obsvInvMCloneVariable, vmass);
	Rebin(hFakeInvMass, fakeInvMClone, fakeInvMCloneVariable, vmass);

	/*
	   c1[5]->cd();
	   sprintf(yTitle, "Entries / %0.f [GeV/c^{2}]", fakeInvMassClone->GetBinWidth(0));
	   axis1F(fakeInvMassClone, "m_{#mu#mu} [GeV/c^{2}]",yTitle);
	   Draw(fakeInvMassClone, 0,1, "", kRed,2, kRed );
	   Draw(obsvInvMassClone, kBlue,1007, "e same", kBlue,2, kBlue );
	 */

	c1[6]->cd();
	c1[6]->SetLogy(1);
	c1[6]->SetLogx(1);
	h_ratio = new TH1F("ObservedOverPredicted","ObservedOverPredicted;Invariant Mass of #mu^{+}#mu^{-} GeV/c^{2};Ratio",140,0,7000);// =obsv/fake
	h_ratio = (TH1F*)obsvInvMClone->Clone("ObservedOverPredicted");
	/*  for (int i=1; i<=140; i++)
		{
		h_ratio->SetBinContent( i, obsvInvMClone->GetBinContent(i)/fakeInvMClone->GetBinContent(i) );
		float eo = obsvInvMClone->GetBinError(i);
		float ef = fakeInvMClone->GetBinError(i);
		float ao = obsvInvMClone->GetBinContent(i);
		float af = fakeInvMClone->GetBinContent(i);

		float er = sqrt(eo/ao*eo/ao+ef/af*ef/af)*ao/af;
		h_ratio->SetBinError(i,er);
		}
	 */
	h_ratio -> Divide(fakeInvMClone);
	axis1F(h_ratio, "m_{#mu#mu} [GeV/c^{2}]","Ratio (Obs/Pred)");
	Draw(h_ratio, 0,1, "e same", kRed,2, kRed );
	TF1 *f_pol = new TF1("f_pol","[0]");
	f_pol->SetParNames("Constant");
	f_pol->SetRange(50,7000);
	f_pol->SetLineColor(kBlack);
	f_pol->SetParameter(0,1.);
	f_pol->SetParLimits(0,0.1,10);
	h_ratio->Fit("f_pol","VR","same");
	c1[6]->Update();
	TPaveStats *statRatio = (TPaveStats*)h_ratio->GetListOfFunctions()->FindObject("stats");
	statRatio->SetY1NDC(0.2); //new x start position
	statRatio->SetY2NDC(0.4); //new x end position
	statRatio->SetX1NDC(0.2); //new x start position
	statRatio->SetX2NDC(0.7); //new x end position
	statRatio->SetTextColor(kBlue) ;//text color for stats box
	h_ratio->Draw("sames");

	c1[11]->cd();
	c1[11]->SetLogy(1);
	c1[11]->SetLogx(1);
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



	c1[7]->cd();
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


	c1[8]->cd();
	c1[8]->SetLogy(1);
	LoadHist( &fakeDeltaR, "Comparison/h_fakeDeltaR"+s_pthatbin, &file1, fileName);
	LoadHist( &obsvDeltaR, "Comparison/h_obsvDeltaR"+s_pthatbin, &file2, fileName);
	obsvDeltaR ->Rebin(10);
	fakeDeltaR ->Rebin(10);
	obsvDeltaR -> SetMinimum(1e-4);
	SetMaximum(obsvDeltaR, fakeDeltaR);
	sprintf(yTitle, "Entries / %0.2f ", fakeDeltaR->GetBinWidth(0));
	axis1F(fakeDeltaR, "#deltaR_{#mu#mu}",yTitle);
	Draw(fakeDeltaR, 0,1, "", kRed,2, kRed );
	Draw(obsvDeltaR, kBlue,1007, "e same", kBlue,2, kBlue );


	c1[10]->cd();
	if ( LoadHist( &fakeDeltaEta, "Comparison/h_fakeDeltaEta"+s_pthatbin, &file1,  fileName)  && 
			LoadHist( &obsvDeltaEta, "Comparison/h_obsvDeltaEta"+s_pthatbin, &file2, fileName) )
	{
		obsvDeltaEta ->Rebin(10);
		fakeDeltaEta ->Rebin(10);
	TH1F *h_obsvDeltaEtaClone = (TH1F*) obsvDeltaEta->Clone("h_obsvDeltaEtaClone");
	TH1F *h_fakeDeltaEtaClone = (TH1F*) fakeDeltaEta->Clone("h_fakeDeltaEtaClone");
//		obsvDeltaEta -> SetMinimum(1e-4);

		SetMaximum(obsvDeltaEta, fakeDeltaEta);
		sprintf(yTitle, "Entries / %0.2f ", fakeDeltaEta->GetBinWidth(0));
		axis1F(fakeDeltaEta, "#delta#eta_{#mu#mu}",yTitle);
		Draw(fakeDeltaEta, 0,1, "", kRed,2, kRed );
		Draw(obsvDeltaEta, kBlue,1007, "e same", kBlue,2, kBlue );

	c1[16]->cd();
	c1[16]->SetLogy(1);
	TH1F *h_deltaEta_ratio = new TH1F("h_deltaEta_ratio","h_deltaEta_ratio", h_obsvDeltaEtaClone->GetNbinsX(),h_obsvDeltaEtaClone->GetBinLowEdge(1), h_obsvDeltaEtaClone->GetBinLowEdge(h_obsvDeltaEtaClone->GetNbinsX()+1));
	h_deltaEta_ratio = (TH1F*) h_obsvDeltaEtaClone->Clone("h_deltaEta_ratio");
	h_deltaEta_ratio->Divide(h_fakeDeltaEtaClone);
	cout<<h_deltaEta_ratio->GetMaximum()<<endl;
	int maxbin = h_deltaEta_ratio->FindBin(h_deltaEta_ratio->GetMaximum());
	cout<<maxbin<<endl;
	//h_deltaEta_ratio->SetBinContent(maxbin,0);//h_deltaEta_ratio->GetNbinsX(),0);
	//h_deltaEta_ratio->SetBinError(maxbin,0);//h_deltaEta_ratio->GetNbinsX(),0);
//	h_deltaEta_ratio->SetMaximum(100);
	sprintf(yTitle, "ratio= #delta#eta_{obsv}/#delta#eta_{pred} ");
	axis1F(h_deltaEta_ratio, "#delta#eta_{#mu}",yTitle);
	Draw(h_deltaEta_ratio, kBlue,1007, "e same", kBlue,2, kBlue );
	h_deltaEta_ratio->Print("all");
	cout<<h_deltaEta_ratio->GetMaximum()<<endl;

	}

	c1[9]->cd();
	if(  LoadHist( &fakeDeltaPhi, "Comparison/h_fakeDeltaPhi"+s_pthatbin, &file1, fileName) &&
			LoadHist( &obsvDeltaPhi, "Comparison/h_obsvDeltaPhi"+s_pthatbin, &file2, fileName))
	{
		obsvDeltaPhi ->Rebin(10);
		fakeDeltaPhi ->Rebin(10);
	TH1F *h_obsvDeltaPhiClone = (TH1F*) obsvDeltaPhi->Clone("h_obsvDeltaPhiClone");
	TH1F *h_fakeDeltaPhiClone = (TH1F*) fakeDeltaPhi->Clone("h_fakeDeltaPhiClone");
		obsvDeltaPhi -> SetMinimum(1e-4);

		SetMaximum(obsvDeltaPhi, fakeDeltaPhi);
		sprintf(yTitle, "Entries / %0.2f ", fakeDeltaPhi->GetBinWidth(0));
		axis1F(fakeDeltaPhi, "#delta#phi_{#mu#mu}",yTitle);
		Draw(fakeDeltaPhi, 0,1, "", kRed,2, kRed );
		Draw(obsvDeltaPhi, kBlue,1007, "e same", kBlue,2, kBlue );
	c1[18]->cd();
	c1[18]->SetLogy(1);
	TH1F *h_deltaPhi_ratio = new TH1F("h_deltaPhi_ratio","h_deltaPhi_ratio", h_obsvDeltaPhiClone->GetNbinsX(),h_obsvDeltaPhiClone->GetBinLowEdge(1), h_obsvDeltaPhiClone->GetBinLowEdge(h_obsvDeltaPhiClone->GetNbinsX()+1));
	h_deltaPhi_ratio = (TH1F*) h_obsvDeltaPhiClone->Clone("h_deltaPhi_ratio");
	h_deltaPhi_ratio->Divide(h_fakeDeltaPhiClone);
	sprintf(yTitle, "ratio= #delta#phi_{obsv}/#delta#phi_{pred} ");
	axis1F(h_deltaPhi_ratio, "#delta#phi_{#mu}",yTitle);
	Draw(h_deltaPhi_ratio, kBlue,1007, "e same", kBlue,2, kBlue );
	}

	for(int i =0; i<NCanvas; i++)
	{
		c1[i]->cd();
		c1[i]->Print("plots/"+s_pre+canTitle[i]+".eps");
		c1[i]->Print("plots/"+s_pre+canTitle[i]+".gif");
	}
}
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
//-----FIXME how to automatically rebin with desired different binwidths
