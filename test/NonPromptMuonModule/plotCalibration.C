plotCalibration(TString s_infile)
{
gStyle->SetPalette(1);
gStyle->SetPadRightMargin(0.15);

	TFile *infile = new TFile(s_infile);

    TCanvas *c1 = new TCanvas("c1","c1", 0, 0, 1350, 900);
  c1->Divide(3,2);
//c1.Divide(2,2);
c1->cd(1);
muonPhiVsJetPhi->Draw("colz");
muonPhiVsJetPhi->GetXaxis()->SetTitle("#phi_{jet}");
muonPhiVsJetPhi->GetYaxis()->SetTitle("#phi_{#mu}");
muonPhiVsJetPhi->Rebin2D(10,10);
//c1->Print("muonPhiVsJetP.gif");
c1->cd(2);

muonEtaVsJetEta->Draw("colz");
muonEtaVsJetEta->Rebin2D(10,10);
muonEtaVsJetEta->GetXaxis()->SetTitle("#eta_{jet}");
muonEtaVsJetEta->GetYaxis()->SetTitle("#eta_{#mu}");
//c1->Print("muonEtaVsJetEta.gif");

c1->cd(4);

dRMuonJet->Draw();
dRMuonJet->GetXaxis()->SetTitle("#deltaR between muon and jet");
dRMuonJet->GetYaxis()->SetTitle("Entries Per Bin");
dRMuonJet->Rebin(10);
dRMuonJet->GetYaxis()->SetTitleOffset(1.4);

c1->cd(5);

muonPtVsJetPt->Draw("colz");
muonPtVsJetPt->Rebin2D(50,50);
muonPtVsJetPt->GetXaxis()->SetTitle("Jet pT [GeV/c]");
muonPtVsJetPt->GetYaxis()->SetTitle("#mu pT [GeV/c]");
muonPtVsJetPt->GetYaxis()->SetTitleOffset(1.4);
//c1->Print("muonPtVsJetPt.gif");


c1->cd(3);
gStyle->SetPalette(1);
//gStyle->SetPadRightMargin(0.14);
//gStyle->SetPadGridY(1);

TH2F *allJet=infile->Get("allJetPEta");
TH2F *matchedJet=infile->Get("matchedJetPEta");

allJet->Rebin2D(100,20);
matchedJet->Rebin2D(100,20);

TH1F *h_allJetP = (TH1F*) allJet->ProjectionX("_allJetP");
TH1F *h_matchedJetP = (TH1F*) matchedJet->ProjectionX("_matchedJetP");

matchedJet->Divide(allJet);

matchedJet->SetTitle("Fake Rate Vs Jet P&Eta");
matchedJet->Draw()           ;

matchedJet->Draw("colz");
//matchedJet->SetMaximum(1.);
//c1->Print("FakeRateVsEtaP.gif");
//c1->Print("FakeRateVsEtaP.eps");

c1->cd(6);
h_matchedJetP->Divide(h_allJetP);
h_matchedJetP->Draw();

}
