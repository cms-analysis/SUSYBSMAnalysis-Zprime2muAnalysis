/*
 * =====================================================================================
 *
 *       Filename:  plotJetPerformance.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/13/2008 02:40:44 PM CEST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  IHEP, Beijing
 *
 * =====================================================================================
 */
plotJetPerformance( TString s_infile = "a.root")
{
	TFile *inflie = new TFile(s_infile);	

	TH1F * h_jetProb[8];
	TH1F * h_bTagEffVsProbCut[8];
	TString s_partonname[8] = {"non", "d" , "u", "s", "c", "b" , "g", "nonb"};
	TCanvas *canJetProb = new TCanvas("canJetProb","canJetProb",  0,0,550,600);
	TCanvas *canBTagEff = new TCanvas("canBTagEff","canBTagEff",600,0,550,600);

	TLegend *leg = (TLegend*)SetLegend(0.78, 0.53, 0.94, 0.91);
	for(int i =0; i<8; i++)
	{
		canJetProb -> cd();
		TString s_title = "h_jetProb_" ; s_title+=s_partonname[i];
		h_jetProb[i] = (TH1F*)inflie->Get(s_title); 
		h_jetProb[i] -> Scale(1./h_jetProb[i]->Integral(0, h_jetProb[i]->GetNbinsX()+1));
		if (i == 0) axis1F(h_jetProb[i], "jet probability", "entries (normalized)");	
		Draw(h_jetProb[i], 0, 1, "same ", i+1, 2, i+1);

		canBTagEff -> cd();
		s_title = "h_bTagEffVsProbCut_"; s_title+=s_partonname[i];
		h_bTagEffVsProbCut[i]=(TH1F*)inflie->Get(s_title);
		if (i == 0) axis1F(h_bTagEffVsProbCut[i], "jet probability cut", "tag efficiency");	
		Draw(h_bTagEffVsProbCut[i], 0, 1, "same", i+1, 2, i+1);

		s_partonname[i] = " " + s_partonname[i];
		leg->AddEntry(h_jetProb[i], s_partonname[i], "lp");
	}
	float bJetProbabilityCut = 0.70; 
	float bTaggingEfficiency = h_bTagEffVsProbCut[5]->GetBinContent(h_bTagEffVsProbCut[5]->FindBin(bJetProbabilityCut));
	printf("|B-Tagging Efficiency (MC) = %4.2f, when set jet probability cut = %3.2f |\n",bTaggingEfficiency, bJetProbabilityCut);

	canJetProb->SetLogy(1);
	canJetProb->cd(); leg->Draw();
	canJetProb->Print("plots/"+s_infile+"JetProbabilityForDifferentFlavors.gif");

	canBTagEff->cd(); leg->Draw();
	canBTagEff->Print("plots/"+s_infile+"BTaggingEfficiencyForDifferentFlavors.gif");
}

