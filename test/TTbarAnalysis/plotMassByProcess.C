const Bool_t   debug    = false;
const Bool_t   logScale =  true;

const Double_t boxY     =  0.04;
const Double_t evtMin   =  0.01;

const Int_t    nProcess =    16;
const Int_t    nSignal  =     4;
const Int_t    nTotal   = nProcess + nSignal;
const Int_t    BinWidth =    20;
const Int_t    RangeMin =     0;
const Int_t    RangeMax =  1000;

//------------------------------------------------------------------------------
// plotMassByProcess
//------------------------------------------------------------------------------
plotMassByProcess(Int_t  nBTags = 0,
		  Bool_t iso10  = false)
{
  if (nBTags < 0 || nBTags > 2) nBTags = 0;

  //----------------------------------------------------------------------------
  // Input file names
  //----------------------------------------------------------------------------
  TString s_piedra     = (iso10 ? "_iso10" : "");
  TString s_file_gumbo = "/data/1a/piedra/analysis/rootfiles/gumbo_skim" + s_piedra + ".root";

  TString s_kypreos    = (iso10 ? "Pt10" : "");
  TString s_file_tw    = "/data/0b/kypreos/ttbar/output" + s_kypreos + "/twOut.root";
  TString s_file_ww    = "/data/0b/kypreos/ttbar/output" + s_kypreos + "/wwOut.root";
  TString s_file_zz    = "/data/0b/kypreos/ttbar/output" + s_kypreos + "/zzOut.root";
  TString s_file_zw    = "/data/0b/kypreos/ttbar/output" + s_kypreos + "/zwOut.root";

  TString s_mschen     = (iso10 ? "_IsoCut" : "");
  TString s_path       = "/home/mschen/Analysis/TTbarAnalysis_test";
  TString s_file_ttbar = s_path + s_mschen + "/rootfiles/Chowder.root";
  TString s_file_wmunu = s_path + s_mschen + "/rootfiles/WmunuJets.root";
  TString s_file_zmumu = s_path + s_mschen + "/rootfiles/Zmumu.root";

  // pcufl1:~kkotov/CMSSW_2_0_5/src/stew.C

  TString s_kkotov     = (iso10 ? "iso2" : "new");
  TString s_file_stew  = "/data/1a/kkotov/stew_cuts_" + s_kkotov + ".root";

  //----------------------------------------------------------------------------
  // Define the process
  //----------------------------------------------------------------------------
  enum enumProcess
    {
      NonSoupProcess=0,
      GumboDiJets=1,
      ChowderWJets,
      ChowderZJets,
      ChowderTTbar,
      GumboPhotonJets,
      GumboMinbias,
      GumboHiggs,
      GumboZprime,
      StewBBbarToJPsi,
      StewDiJets,
      StewBottomonium,
      StewCharmonium,
      StewBBE,
      StewPPElectronX,
      StewPPMuonX,
      SignalTW,
      SignalWW,
      SignalZW,
      SignalZZ
    };

  TString S_PrimaryProcess[nProcess];

  S_PrimaryProcess[NonSoupProcess ] = "NonSoupProcess";
  S_PrimaryProcess[GumboDiJets    ] = "GumboDiJets";
  S_PrimaryProcess[GumboHiggs     ] = "GumboHiggs";
  S_PrimaryProcess[GumboZprime    ] = "GumboZprime";
  S_PrimaryProcess[GumboPhotonJets] = "GumboPhotonJets";
  S_PrimaryProcess[GumboMinbias   ] = "GumboMinbias";
  S_PrimaryProcess[ChowderTTbar   ] = "ChowderTTbar";
  S_PrimaryProcess[ChowderZJets   ] = "ChowderZJets";
  S_PrimaryProcess[ChowderWJets   ] = "ChowderWJets";
  S_PrimaryProcess[StewPPMuonX    ] = "StewPPMuonX";
  S_PrimaryProcess[StewPPElectronX] = "StewPPElectronX";
  S_PrimaryProcess[StewBBE        ] = "StewBBE";
  S_PrimaryProcess[StewBBbarToJPsi] = "StewBBbarToJPsi";
  S_PrimaryProcess[StewDiJets     ] = "StewDiJets";
  S_PrimaryProcess[StewBottomonium] = "StewBottomonium";
  S_PrimaryProcess[StewCharmonium ] = "StewCharmonium";

  TString label[nTotal];

  label[NonSoupProcess ] = " non-soup";
  label[StewPPMuonX    ] = " pp#mu+X";
  label[StewPPElectronX] = " ppe+X";
  label[StewBBE        ] = " bbe";
  label[StewCharmonium ] = " charmonium";
  label[StewBottomonium] = " bottomonium";
  label[StewDiJets     ] = " QCD 0<pt<15";
  label[StewBBbarToJPsi] = " b#bar{b} #rightarrow J/#psi";
  label[ChowderTTbar   ] = " t#bar{t}";
  label[ChowderZJets   ] = " Z #rightarrow #mu#mu";
  label[ChowderWJets   ] = " W #rightarrow #mu#nu";
  label[SignalTW       ] = " tW";
  label[SignalZW       ] = " ZW";
  label[SignalZZ       ] = " ZZ";
  label[SignalWW       ] = " WW";
  label[GumboDiJets    ] = " QCD";
  label[GumboPhotonJets] = " #gamma+jets";
  label[GumboHiggs     ] = " Higgs 150";
  label[GumboMinbias   ] = " min bias";
  label[GumboZprime    ] = " Z'";

  TString twlabel[nTotal];

  twlabel[NonSoupProcess ] = " non-soup";
  twlabel[StewPPMuonX    ] = " ppmu+X";
  twlabel[StewPPElectronX] = " ppe+X";
  twlabel[StewBBE        ] = " bbe";
  twlabel[StewCharmonium ] = " charmonium";
  twlabel[StewBottomonium] = " bottomonium";
  twlabel[StewDiJets     ] = " QCD 0&lt;pt&lt;15";
  twlabel[StewBBbarToJPsi] = " bbbar -> J/psi";
  twlabel[ChowderTTbar   ] = " ttbar+jets";
  twlabel[ChowderZJets   ] = " Z -> mumu";
  twlabel[ChowderWJets   ] = " W -> munu";
  twlabel[SignalTW       ] = " tW";
  twlabel[SignalZW       ] = " ZW";
  twlabel[SignalZZ       ] = " ZZ";
  twlabel[SignalWW       ] = " WW";
  twlabel[GumboDiJets    ] = " QCD";
  twlabel[GumboPhotonJets] = " photon+jets";
  twlabel[GumboHiggs     ] = " Higgs 150";
  twlabel[GumboMinbias   ] = " min bias";
  twlabel[GumboZprime    ] = " Z'";

  Color_t col[nTotal];
  
  col[NonSoupProcess ] = kBlack;
  col[SignalTW       ] = kAzure-6;
  col[SignalZW       ] = kAzure-7;
  col[SignalZZ       ] = kAzure-8; 
  col[SignalWW       ] = kAzure-9;
  col[GumboDiJets    ] = kRed+2;
  col[GumboPhotonJets] = kRed-3;
  col[GumboHiggs     ] = kRed-7;
  col[ChowderTTbar   ] = kGreen+2;
  col[ChowderZJets   ] = kOrange-3;
  col[ChowderWJets   ] = kOrange-2;
  col[StewPPMuonX    ] = kGreen+4;
  col[StewPPElectronX] = kGreen-1;
  col[StewBBE        ] = kGreen-5;
  col[StewCharmonium ] = kGreen-8;
  col[StewBottomonium] = kGreen-10;
  col[StewDiJets     ] = kRed+2;
  col[StewBBbarToJPsi] = kRed-3;

  Bool_t   notEmpty[nTotal];
  Int_t    procId  [nTotal];
  Double_t evtsId  [nTotal];

  for (int i=0; i<nTotal; i++) {
    notEmpty[i] = false;
    procId  [i] =  -999;
    evtsId  [i] =  -999;
  }

  TString bTagLabel    [3] = {"", "1B", "2B"};
  TString bTagLabelStew[3] = {"", "_1b", "_2b"};
  TString bTagDraw     [3] = {"no b-tagging", "1 b-tag", "2 b-tag"};

  Double_t legY = 0.0;

  TH1F* hMass[nTotal];

  THStack* hs = new THStack("hs","hs");

  //----------------------------------------------------------------------------
  // Get the Gumbo + Chowder histograms
  //----------------------------------------------------------------------------
  TFile* file_gumbo = new TFile(s_file_gumbo);		
  TFile* file_ttbar = new TFile(s_file_ttbar);		
  TFile* file_zmumu = new TFile(s_file_zmumu);		
  TFile* file_wmunu = new TFile(s_file_wmunu);		

  for (int i=0; i<nProcess; i++) {

    if (S_PrimaryProcess[i].Contains("Stew")) continue;

    TString s_title;

    if (i == ChowderZJets || i == ChowderWJets) {
      s_title = "h_obsvInvMass" + bTagLabel[nBTags];
    }
    else {
      s_title = "MassByProcess/h_Inv";

      if (nBTags == 0) s_title += "ariant";

      s_title += "Mass";
      s_title += bTagLabel[nBTags];
      s_title += "ByPrimaryProcess";
      s_title += S_PrimaryProcess[i];
    }

    if (s_title.Contains("Gumbo") || s_title.Contains("NonSoupProcess")) {
      hMass[i] = (TH1F*)file_gumbo->Get(s_title);

      int gumbo_events = 43341218;
      int ntupl_events = 36623868;

      double gumbo_scale = (double)ntupl_events / gumbo_events;

      hMass[i]->Scale(gumbo_scale);
    }

    if (i == ChowderTTbar) {
      hMass[i] = (TH1F*)file_ttbar->Get(s_title);
    }
    else if (i == ChowderZJets) {
      hMass[i] = (TH1F*)file_zmumu->Get(s_title);
      hMass[i]->Scale(8.293e-02);
    }
    else if (i == ChowderWJets) {
      hMass[i] = (TH1F*)file_wmunu->Get(s_title);
    }
  }

  //----------------------------------------------------------------------------
  // Get the Signal histograms
  //----------------------------------------------------------------------------
  TFile* file_tw = new TFile(s_file_tw);		
  TFile* file_ww = new TFile(s_file_ww);		
  TFile* file_zw = new TFile(s_file_zw);		
  TFile* file_zz = new TFile(s_file_zz);		

  //       weight    =  xsec * accept * lumi / Nevents
  Double_t weight_tw =  62.0 *    1.0 *  100 / 423791;
  Double_t weight_zw = 0.585 *    1.0 *  100 /  58897;
  Double_t weight_zz = 70.38 * 0.4177 *  100 /  78000;
  Double_t weight_ww = 114.3 *    1.0 *  100 / 805261;

  TString s_title = "h_obsvInvMass";

  s_title += bTagLabel[nBTags];

  hMass[SignalTW] = (TH1F*)file_tw->Get(s_title);
  hMass[SignalTW]->Scale(weight_tw);

  hMass[SignalWW] = (TH1F*)file_ww->Get(s_title);
  hMass[SignalWW]->Scale(weight_ww);

  hMass[SignalZW] = (TH1F*)file_zw->Get(s_title);
  hMass[SignalZW]->Scale(weight_zw);

  hMass[SignalZZ] = (TH1F*)file_zz->Get(s_title);
  hMass[SignalZZ]->Scale(weight_zz);

  //----------------------------------------------------------------------------
  // Get the Stew histograms
  //----------------------------------------------------------------------------
  TFile* file_stew = new TFile(s_file_stew);

  for (int i=0; i<nProcess; i++) {
    if (!S_PrimaryProcess[i].Contains("Stew")) continue;

    TString s_title = bTagLabelStew[nBTags];

    if (i == StewBBbarToJPsi) hMass[i] = (TH1F*)file_stew->Get("m60" + s_title);
    if (i == StewDiJets     ) hMass[i] = (TH1F*)file_stew->Get("m61" + s_title);
    if (i == StewBottomonium) {
      hMass[i] = (TH1F*)file_stew->Get("m62" + s_title);
      hMass[i]->Add((TH1F*)file_stew->Get("m63" + s_title));
    }
    if (i == StewCharmonium ) {
      hMass[i] = (TH1F*)file_stew->Get("m64" + s_title);
      hMass[i]->Add((TH1F*)file_stew->Get("m65" + s_title));
    }
    if (i == StewBBE        ) {
      hMass[i] = (TH1F*)file_stew->Get("m66" + s_title);
      hMass[i]->Add((TH1F*)file_stew->Get("m67" + s_title));
      hMass[i]->Add((TH1F*)file_stew->Get("m68" + s_title));
    }

    if (i == StewPPElectronX) hMass[i] = (TH1F*)file_stew->Get("m69" + s_title);
    if (i == StewPPMuonX    ) hMass[i] = (TH1F*)file_stew->Get("m70" + s_title);
  }

  //----------------------------------------------------------------------------
  // Write the sum into a file
  //----------------------------------------------------------------------------
  TH1F* hMassAll = new TH1F("hMassAll","hMassAll",7000,0,7000);
  hMassAll->Sumw2();

  for (int i=0; i<nTotal; i++) {

//    if (i == StewBBbarToJPsi) continue;
//    if (i == StewBBE        ) continue;
//    if (i == StewPPMuonX    ) continue;

    hMassAll->Add(hMass[i]);
  }

  TString s_outfile = "s-all-dimuon-mass-process";
  s_outfile += bTagLabel[nBTags];
  s_outfile += s_piedra;
  s_outfile += ".root";

  TFile* outfile = new TFile(s_outfile.Data(),"recreate");
  hMassAll->Write();
  outfile->Close();

  //----------------------------------------------------------------------------
  // Rebin
  //----------------------------------------------------------------------------
  for (int i=0; i<nTotal; i++) {

    hMass[i]->Rebin(BinWidth);	
    hMass[i]->GetXaxis()->SetRangeUser(RangeMin,RangeMax);

    if (hMass[i]->GetSumOfWeights() > evtMin) {
      legY += boxY;
      notEmpty[i] = true;
    }
  }

  //----------------------------------------------------------------------------
  // Order histograms
  //----------------------------------------------------------------------------
  Int_t iposition = 0;

  while (iposition < nTotal) {

    for (Int_t iprocess=0; iprocess<nTotal; iprocess++) {
      
      bool skip = false;
      
      for (int k=0; k<iposition; k++)
	if (iprocess == procId[k]) skip = true;
      
      if (skip) continue;

      Double_t nEvents = hMass[iprocess]->GetSumOfWeights();
      
      if (nEvents > evtsId[iposition]) {
	procId[iposition] = iprocess;
	evtsId[iposition] = nEvents;
      }
    }

    iposition++;
  }

  //----------------------------------------------------------------------------
  // Draw
  //----------------------------------------------------------------------------
  TCanvas* can = new TCanvas("can","invariant mass");

  if (logScale) can->SetLogy(1);

  TLegend *leg = (TLegend*)SetLegend(0.74, 0.91-legY, 0.87, 0.91);
  if (debug) leg->SetFillStyle(1001);

  leg->SetTextSize(0.03);

  for (Int_t iposition=nTotal-1; iposition>-1; iposition--) {

    if (!notEmpty[procId[iposition]]) continue;

    if (iposition == 0) {

      char yTitle[100];
      sprintf(yTitle, "entries / %.0f GeV/c",
	      hMass[procId[iposition]]->GetBinWidth(0));
      hMass[procId[iposition]]->SetTitle("");
      hMass[procId[iposition]]->GetYaxis()->SetTitle(yTitle);
      hMass[procId[iposition]]->GetYaxis()->SetTitleFont(42);
      hMass[procId[iposition]]->GetYaxis()->SetTitleOffset(  1.70);
      hMass[procId[iposition]]->GetYaxis()->SetTitleSize  (  0.05);
      hMass[procId[iposition]]->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");

      Draw(hMass[procId[iposition]],
	   col[procId[iposition]],1001,"hist",col[procId[iposition]]);
    }
    else
      Draw(hMass[procId[iposition]],
	   col[procId[iposition]],1001,"hist,same",col[procId[iposition]]);

    hs->Add(hMass[procId[iposition]]);
  }

  hs->Draw("hist,same");

  if (nBTags == 0 && logScale)
    hMass[procId[0]]->SetMaximum(1e3*hs->GetMaximum());

  //----------------------------------------------------------------------------
  // Labels
  //----------------------------------------------------------------------------
  printf("\n| *process*  |  *%s* |\n",
	 bTagDraw[nBTags].Data());

  for (Int_t iposition=0; iposition<nTotal; iposition++) {

    char printEvents[100];

    if (procId[iposition] == NonSoupProcess) continue;

    if (hMass[procId[iposition]]->GetSumOfWeights() > evtMin)
      sprintf(printEvents,"%8.1f",
	      hMass[procId[iposition]]->GetSumOfWeights());
    else
      sprintf(printEvents,"-");

    printf("|%-30s  |  %9s |\n",
	   twlabel[procId[iposition]].Data(),printEvents);

    if (!notEmpty[procId[iposition]]) continue;

    leg->AddEntry(hMass[procId[iposition]],label[procId[iposition]],"f");
  }

  printf("\n");

  leg->Draw();

  //----------------------------------------------------------------------------
  // Cosmetics
  //----------------------------------------------------------------------------
  TPad* pad = (TPad*)can->GetPad(0);

  char title[100];
  sprintf(title,"CMS Monte Carlo");
  PRLabel(pad,title,0,11,0.05);

  double lumi = 1e2;

  sprintf(title,"L #approx %.0f pb^{-1}",lumi);
  PRLabel(pad,title,1,31,0.05);
  
  if (nBTags > 0)
    SetPave(bTagDraw[nBTags].Data(),0.03,0.03,0.17,0.09,50,kWhite);

  if (iso10)
    SetPave("iso10",0.17,0.03,0.28,0.09,30,kWhite);

  //----------------------------------------------------------------------------
  // Save
  //----------------------------------------------------------------------------
  can->GetFrame()->DrawClone();

  sprintf(title,"all-dimuon-mass-process%s%s.gif",
	  bTagLabel[nBTags].Data(),s_piedra.Data());
  can->SaveAs(title);
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
