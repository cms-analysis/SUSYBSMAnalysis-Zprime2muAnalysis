const Double_t loZ      =  70.0;
const Double_t hiZ      = 110.0;
const Int_t    BinWidth =    10;
const Int_t    RangeMin =     0;
const Int_t    RangeMax =  1000;

Char_t   title[100];
Double_t lumi = 1e2;

//------------------------------------------------------------------------------
// plotTTbarMass
//------------------------------------------------------------------------------
void plotTTbarMass(Bool_t iso10 = true)
{
  if (!iso10) {
    cout << "\n plotTTbarMass_v2 not ready for iso10 = false\n" << endl;
    return;
  }

  enum ENUMCountingNumbers
    {
      NTTbarEvents=1, 
      NTTbarEventsWith1BPartonAccepted,
      NTTbarEventsWith2BPartonAccepted, 
      NEventsWith1BTagged, 
      NEventsWith2BTagged, 
      NEventsWithMBTagged,
      NEventsWith1BTagged_ZMassVeto, 
      NEventsWith2BTagged_ZMassVeto, 
      NEventsWithMBTagged_ZMassVeto,
      NTTbarEventsWith1BTagged,
      NTTbarEventsWith2BTagged, 
      NTTbarEventsWithMBTagged,
      NTTbarEventsWith1BTagged_ZMassVeto,
      NTTbarEventsWith2BTagged_ZMassVeto, 
      NTTbarEventsWithMBTagged_ZMassVeto,
    };
  
  //----------------------------------------------------------------------------
  // Get histograms
  //----------------------------------------------------------------------------
  TString s_mschen = (iso10 ? "_IsoCut" : "");
  TString s_path   = "/home/mschen/Analysis/TTbarAnalysis_test";
  TString s_infile = s_path + s_mschen + "/rootfiles/Chowder.root";

  TFile* infile = new TFile(s_infile);
  TFile* file0B = new TFile("/data/1a/piedra/analysis/rootfiles/all-dimuon-mass-process_iso10.root"  );
  TFile* file1B = new TFile("/data/1a/piedra/analysis/rootfiles/all-dimuon-mass-process1B_iso10.root");
  TFile* file2B = new TFile("/data/1a/piedra/analysis/rootfiles/all-dimuon-mass-process2B_iso10.root");

  TH1F* h_numbers       = (TH1F*)infile->Get("h_numbers");
  TH1F* h_invMass_ttbar = (TH1F*)infile->Get("MassByProcess/h_InvariantMassByPrimaryProcessChowderTTbar");
  TH1F* h_obsvInvMass   = (TH1F*)file0B->Get("hMassAll");
  TH1F* h_obsvInvMass1B = (TH1F*)file1B->Get("hMassAll");
  TH1F* h_obsvInvMass2B = (TH1F*)file2B->Get("hMassAll");

  //----------------------------------------------------------------------------
  // Remove 70 <= mass < 110 window
  //----------------------------------------------------------------------------
  for (Int_t j=0; j<h_obsvInvMass->GetNbinsX(); j++) {
    if (j >= h_obsvInvMass->FindBin(loZ) && j < h_obsvInvMass->FindBin(hiZ)) {
      h_obsvInvMass1B->SetBinContent(j,0.0);
      h_obsvInvMass2B->SetBinContent(j,0.0);

      h_obsvInvMass1B->SetBinError(j,0.0);
      h_obsvInvMass2B->SetBinError(j,0.0);
    }
  }

  //----------------------------------------------------------------------------
  // Compute n1ScaleFactor and n2ScaleFactor
  //----------------------------------------------------------------------------
  Double_t ttbarEvents1bAccepted       = h_numbers->GetBinContent(NTTbarEventsWith1BPartonAccepted);
  Double_t ttbarEvents2bAccepted       = h_numbers->GetBinContent(NTTbarEventsWith2BPartonAccepted);
  Double_t ttbarEventsTotal            = h_numbers->GetBinContent(NTTbarEvents);
  Double_t Error_ttbarEvents1bAccepted = h_numbers->GetBinError(NTTbarEventsWith1BPartonAccepted);
  Double_t Error_ttbarEvents2bAccepted = h_numbers->GetBinError(NTTbarEventsWith2BPartonAccepted);
  Double_t Error_ttbarEventsTotal      = h_numbers->GetBinError(NTTbarEvents);

  Double_t A1 = ttbarEvents1bAccepted / ttbarEventsTotal;
  Double_t A2 = ttbarEvents2bAccepted / ttbarEventsTotal;

  Double_t Error_A1 = 0;
  Error_A1 += pow(Error_ttbarEvents1bAccepted / ttbarEvents1bAccepted,2);
  Error_A1 += pow(Error_ttbarEventsTotal      / ttbarEventsTotal,     2);
  Error_A1  = A1 * sqrt(Error_A1);

  Double_t Error_A2 = 0;
  Error_A2 += pow(Error_ttbarEvents2bAccepted / ttbarEvents2bAccepted,2);
  Error_A2 += pow(Error_ttbarEventsTotal      / ttbarEventsTotal,     2);
  Error_A2  = A1 * sqrt(Error_A2);

  Double_t n1       = h_obsvInvMass1B->GetSumOfWeights();
  Double_t n2       = h_obsvInvMass2B->GetSumOfWeights();
  Double_t Error_n1 = sqrt(n1);
  Double_t Error_n2 = sqrt(n2);
	  
  Double_t b_eff       = (A1/A2 + 2) / (n1/n2 + 2);
  Double_t Error_b_eff = b_eff * sqrt(pow(Error_A1/A1,2) + pow(Error_A2/A2,2) +
				      pow(Error_n1/n1,2) + pow(Error_n2/n2,2));

  Double_t n2ScaleFactor = 1/b_eff/b_eff/A2;
  Double_t n1ScaleFactor = 1/(A1*b_eff + 2*A2*b_eff*(1-b_eff));

  //----------------------------------------------------------------------------
  // Draw
  //----------------------------------------------------------------------------
  TCanvas* can = new TCanvas("can","invariant mass");
  can->SetLogy(1);

  h_obsvInvMass  ->Rebin(BinWidth);	
  h_obsvInvMass1B->Rebin(BinWidth);	
  h_obsvInvMass2B->Rebin(BinWidth);	
  h_invMass_ttbar->Rebin(BinWidth);
  
  h_obsvInvMass  ->GetXaxis()->SetRangeUser(RangeMin,RangeMax);
  h_obsvInvMass1B->GetXaxis()->SetRangeUser(RangeMin,RangeMax);
  h_obsvInvMass2B->GetXaxis()->SetRangeUser(RangeMin,RangeMax);
  h_invMass_ttbar->GetXaxis()->SetRangeUser(RangeMin,RangeMax);

  h_obsvInvMass1B->Scale(n1ScaleFactor);
  h_obsvInvMass2B->Scale(n2ScaleFactor);

  char yTitle[100];
  sprintf(yTitle, "entries / %.0f GeV/c^{2}", h_obsvInvMass->GetBinWidth(0));
  axis1F(h_obsvInvMass, "m_{#mu#mu} [GeV/c^{2}]", yTitle);

  Double_t maxY = h_obsvInvMass->GetMaximum();
  if (h_obsvInvMass1B->GetMaximum() > maxY) maxY = h_obsvInvMass1B->GetMaximum();
  if (h_obsvInvMass2B->GetMaximum() > maxY) maxY = h_obsvInvMass2B->GetMaximum();
  if (h_invMass_ttbar->GetMaximum() > maxY) maxY = h_invMass_ttbar->GetMaximum();

  h_obsvInvMass->SetMaximum(10*maxY);

  Draw(h_obsvInvMass  ,kGreen -6,1001,"hist"     ,kGreen -6);
  Draw(h_invMass_ttbar,kOrange-3,1001,"hist,same",kOrange-3);
  Draw(h_obsvInvMass1B,        0,   1,"same"     ,kAzure -5,2.0,kAzure-5,1.2,kFullTriangleUp  );
  Draw(h_obsvInvMass2B,        0,   1,"same"     ,kRed   -3,2.0,kRed  -3,1.2,kFullTriangleDown);

  TLegend* leg = (TLegend*)SetLegend(0.73, 0.71, 0.89, 0.90);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_obsvInvMass  , " all "         , "f" );
  leg->AddEntry(h_invMass_ttbar, " t#bar{t}+jets", "f" );
  leg->AddEntry(h_obsvInvMass1B, " 1 b estimate" , "lp");
  leg->AddEntry(h_obsvInvMass2B, " 2 b estimate" , "lp");
  leg->Draw();

  //----------------------------------------------------------------------------
  // Cosmetics
  //----------------------------------------------------------------------------
  TPad* pad = (TPad*)can->GetPad(0);

  sprintf(title,"CMS Monte Carlo");
  PRLabel(pad,title,0,11,0.05);

  sprintf(title,"L #approx %.0f pb^{-1}",lumi);
  PRLabel(pad,title,1,31,0.05);
  
  Char_t effPrint[100];
  sprintf(effPrint,"#varepsilon_{b} = %.2f #pm %.2f",b_eff,Error_b_eff);
  SetPave(effPrint,0.03,0.03,0.30,0.09,50,kWhite);

  //----------------------------------------------------------------------------
  // Save
  //----------------------------------------------------------------------------
  can->GetFrame()->DrawClone();
  can->Print("tt_estimate_without_stew-bbjpsi-bbe-ppmux_v3.gif");	

  TFile* ttbarFile = new TFile("ttbar-estimates" + s_mschen + ".root",
			       "recreate");
  h_invMass_ttbar->Write("h_invMass_tt_truth");
  h_obsvInvMass1B->Write("h_invMass_tt_1bEst");
  h_obsvInvMass2B->Write("h_invMass_tt_2bEst");
  ttbarFile->Close();

  //----------------------------------------------------------------------------
  // Draw ratio
  //----------------------------------------------------------------------------
  TCanvas* canRatio = new TCanvas("canRatio","canRatio",600,0,550,600);
  canRatio->SetLogy(1);

  h_ratioInvMass1BoverGenuine = (TH1F*)h_obsvInvMass1B->Clone("h_ratioInvMass1BoverGenuine");
  h_ratioInvMass2BoverGenuine = (TH1F*)h_obsvInvMass2B->Clone("h_ratioInvMass2BoverGenuine");
  h_ratioInvMass1Bover2B      = (TH1F*)h_obsvInvMass1B->Clone("h_ratioInvMass1Bover2B");

  h_ratioInvMass1BoverGenuine->Divide(h_invMass_ttbar);
  h_ratioInvMass2BoverGenuine->Divide(h_invMass_ttbar);
  h_ratioInvMass1Bover2B     ->Divide(h_obsvInvMass2B);	

  axis1F(h_ratioInvMass1BoverGenuine,"m_{#mu#mu} [GeV/c^{2}]","ratio");

  Int_t col1 = h_obsvInvMass1B->GetLineColor();
  Int_t col2 = h_obsvInvMass2B->GetLineColor();
  Int_t col3 = h_invMass_ttbar->GetLineColor();

  Draw(h_ratioInvMass1BoverGenuine, 0, 1, "p"     , col1, 2, col1, 1, kFullTriangleUp  );
  Draw(h_ratioInvMass2BoverGenuine, 0, 1, "p,same", col2, 2, col2, 1, kFullTriangleDown);
  Draw(h_ratioInvMass1Bover2B,      0, 1, "p,same", col3, 2, col3, 1, kOpenCircle      );

  TLine* lineat1 = new TLine(RangeMin,1,RangeMax,1);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineStyle(     3);
  lineat1->SetLineWidth(     3);
  lineat1->Draw        ("same");

  TF1* f_pol1BoverGen = new TF1("f_pol1BoverGen", "[0]");
  TF1* f_pol2BoverGen = new TF1("f_pol2BoverGen", "[0]");
  TF1* f_pol1Bover2B  = new TF1("f_pol1Bover2B" , "[0]");

  h_ratioInvMass1BoverGenuine->Fit("f_pol1BoverGen","mnq","same");
  h_ratioInvMass2BoverGenuine->Fit("f_pol2BoverGen","mnq","same");
  h_ratioInvMass1Bover2B     ->Fit("f_pol1Bover2B" ,"mnq","same");

  char s_1Bover[100];
  sprintf(s_1Bover," 1b est / obs = %.2f #pm %.2f",
	  f_pol1BoverGen->GetParameter(0),f_pol1BoverGen->GetParError(0));

  char s_2Bover[100];
  sprintf(s_2Bover," 2b est / obs = %.2f #pm %.2f",
	  f_pol2BoverGen->GetParameter(0),f_pol2BoverGen->GetParError(0));

  char s_1Bover2B[100];
  sprintf(s_1Bover2B," 1b est / 2b est = %.2f #pm %.2f",
	  f_pol1Bover2B->GetParameter(0),f_pol1Bover2B->GetParError(0));

  TLegend* legRatio = (TLegend*)SetLegend(0.20, 0.77, 0.43, 0.92);
  legRatio->SetTextSize(0.03);
  legRatio->AddEntry(h_ratioInvMass1BoverGenuine, s_1Bover  , "lp");
  legRatio->AddEntry(h_ratioInvMass2BoverGenuine, s_2Bover  , "lp");
  legRatio->AddEntry(h_ratioInvMass1Bover2B,      s_1Bover2B, "lp");
  legRatio->Draw();

  //----------------------------------------------------------------------------
  // Cosmetics
  //----------------------------------------------------------------------------
  TPad* pad2 = (TPad*)canRatio->GetPad(0);

  sprintf(title,"CMS Monte Carlo");
  PRLabel(pad2,title,0,11,0.05);

  sprintf(title,"L #approx %.0f pb^{-1}",lumi);
  PRLabel(pad2,title,1,31,0.05);
  
  //----------------------------------------------------------------------------
  // Save
  //----------------------------------------------------------------------------
  canRatio->GetFrame()->DrawClone();
  canRatio->Print("DiMuonInvMass_Ratio.gif");
}

//------------------------------------------------------------------------------
// SetPave v2 (to be modified in PlotUtils)
//------------------------------------------------------------------------------
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
