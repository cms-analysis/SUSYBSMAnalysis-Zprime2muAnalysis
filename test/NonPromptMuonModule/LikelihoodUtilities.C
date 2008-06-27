#ifndef LLR_C
#define LLR_C
#include "TF1.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include <iostream>
#include "TH1F.h"
Double_t GetLikelihoodRatio( Double_t sumpt, Double_t dxy, TF1 *fDxyS, TF1 *fDxyB, TF1 *fIsoSumPtS, TF1 *fIsoSumPtB )
{
		double PS_dxy = fDxyS->Eval(dxy);
		double PB_dxy = fDxyB->Eval(dxy);
		double PS_sumpt = fIsoSumPtS->Eval(sumpt);
		double PB_sumpt = fIsoSumPtB->Eval(sumpt);
		double llr = 0;
		double PSB = PS_dxy*PS_sumpt + PB_dxy*PB_sumpt;
		if(PSB!=0) llr = PS_dxy*PS_sumpt /  PSB;
		return llr ; //PS_dxy*PS_sumpt/(PS_dxy*PS_sumpt + PB_dxy*PB_sumpt);
}
double GetLRCutAtEff(TString trainedFileName, double effS, TString hSigVsLRCutName="h_SigEffVsLLRCut", TString hBkgVsLRCutName="h_BkgRejVsLLRCut")
{
	double ResponseCut_Likelihood = 0;
	TFile *TrainingResponseInput(0);
	if (!gSystem->AccessPathName( trainedFileName )) {
		cout << "--- Accessing Train OutPut file: " << trainedFileName << endl;
		TrainingResponseInput = TFile::Open(trainedFileName);
	} 
	if (!TrainingResponseInput) {
		cout << "ERROR: could not open trained file: " << trainedFileName << endl;
		exit(1);
	}
	TH1F *htmp = dynamic_cast<TH1F*>(TrainingResponseInput->Get(hSigVsLRCutName));
	Int_t binx;
	htmp->GetBinWithContent(effS,binx,0,0,0.001);
	ResponseCut_Likelihood = htmp->GetBinLowEdge(1) + htmp->GetBinWidth(1)*binx;
	if(htmp) delete htmp;
	htmp = dynamic_cast<TH1F*>(TrainingResponseInput->Get(hBkgVsLRCutName));
	Double_t bkgrej = htmp->GetBinContent(binx);
	cout<<"Keep Signal Efficiency="<<effS<<" and Reject bkg ="<<bkgrej<<"  llr cut at "<<ResponseCut_Likelihood<<endl;
	if(htmp) delete htmp;
	if(TrainingResponseInput) delete TrainingResponseInput;
	return ResponseCut_Likelihood;
}
#endif
