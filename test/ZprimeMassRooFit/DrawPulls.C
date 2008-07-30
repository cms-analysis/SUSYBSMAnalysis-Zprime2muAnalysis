Int_t    rebin = 20;
Double_t xMin  = -4;
Double_t xMax  =  4;

void DrawPulls(Char_t fName[100] = "pulls.root")
{
  TFile* f = new TFile(fName, "read");

  TH1F* pull_N_ZP = (TH1F*)f->Get("pull_N_ZP");
  TH1F* pull_N_DY = (TH1F*)f->Get("pull_N_DY");

  pull_N_ZP->Rebin(rebin);
  pull_N_DY->Rebin(rebin);

  TF1* fitPull_ZP = new TF1("fitPull_ZP", "gaus", xMin, xMax);
  fitPull_ZP->SetLineWidth(3);

  TF1* fitPull_DY = new TF1("fitPull_DY", "gaus", xMin, xMax);
  fitPull_DY->SetLineWidth(3);

  Char_t yTitle_ZP[100];
  sprintf(yTitle_ZP, "entries / %.2f", pull_N_ZP->GetBinWidth(0));

  Char_t yTitle_DY[100];
  sprintf(yTitle_DY, "entries / %.2f", pull_N_DY->GetBinWidth(0));

  Char_t title[200] = "#left(N_{fit} - N_{Poisson}#right) / #sigma_{N_{fit}}";

  //----------------------------------------------------------------------------
  // Draw
  //----------------------------------------------------------------------------
  TCanvas* c1 = new TCanvas("c1","c1");

  axis1F(pull_N_ZP, title, yTitle_ZP);
  Draw(pull_N_ZP, kAzure-4, 1001, "hist", kAzure-4, 1);

  pull_N_ZP->Fit("fitPull_ZP","noqr");
  fitPull_ZP->Draw("same");

  TPad* pad_ZP = (TPad*)c1->GetPad(0);
  PRLabel(pad_ZP,"N_{Z'}",0,11,0.05);

  pad_ZP->GetFrame()->DrawClone();

  //----------------------------------------------------------------------------
  TCanvas* c2 = new TCanvas("c2","c2",575,0,550,600);

  axis1F(pull_N_DY, title, yTitle_DY);
  Draw(pull_N_DY, kAzure-4, 1001, "hist", kAzure-4, 1);

  pull_N_DY->Fit("fitPull_DY","noqr");
  fitPull_DY->Draw("same");

  TPad* pad_DY = (TPad*)c2->GetPad(0);
  PRLabel(pad_DY,"N_{DY}",0,11,0.05);

  pad_DY->GetFrame()->DrawClone();

  //----------------------------------------------------------------------------
  // Legends
  //----------------------------------------------------------------------------
  Char_t fitPrint_ZP[4][100];
  Char_t fitPrint_DY[4][100];

  TPaveText* pv_ZP = new TPaveText(0.63,0.64,0.92,0.88,"ndc");
  TPaveText* pv_DY = new TPaveText(0.63,0.64,0.92,0.88,"ndc");

  pv_ZP->SetBorderSize(   0);
  pv_ZP->SetFillColor (   0);
  pv_ZP->SetTextAlign (  13);
  pv_ZP->SetTextFont  (  42);
  pv_ZP->SetTextSize  (0.03);

  pv_DY->SetBorderSize(   0);
  pv_DY->SetFillColor (   0);
  pv_DY->SetTextAlign (  13);
  pv_DY->SetTextFont  (  42);
  pv_DY->SetTextSize  (0.03);

  sprintf(fitPrint_ZP[3], "Experiments = %d", pull_N_ZP->GetEntries());
  sprintf(fitPrint_DY[3], "Experiments = %d", pull_N_DY->GetEntries());

  pv_ZP->AddText(fitPrint_ZP[3]);
  pv_DY->AddText(fitPrint_DY[3]);

  for (Int_t i=0; i<3; i++) {

    sprintf(fitPrint_ZP[i],"%s = %.2f #pm %.2f",
	    fitPull_ZP->GetParName(i),
	    fitPull_ZP->GetParameter(i),
	    fitPull_ZP->GetParError(i));

    sprintf(fitPrint_DY[i],"%s = %.2f #pm %.2f",
	    fitPull_DY->GetParName(i),
	    fitPull_DY->GetParameter(i),
	    fitPull_DY->GetParError(i));

    pv_ZP->AddText(fitPrint_ZP[i]);
    pv_DY->AddText(fitPrint_DY[i]);
  }

  c1->cd();
  pv_ZP->Draw("same");
  c1->SaveAs("pull_ZP.gif");

  c2->cd();
  pv_DY->Draw("same");
  c2->SaveAs("pull_DY.gif");
}

