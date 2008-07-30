#include "RooConfigure.C"

Int_t    nFits     = 25;
Double_t residYMax = -1;
Double_t residYMin = -1;

//------------------------------------------------------------------------------
// ZprimeMassRooFit
//------------------------------------------------------------------------------
void ZprimeMassRooFit(Bool_t doTheFit = true,
		      Bool_t dyFit    = false,
		      Bool_t draw     = true,
		      Bool_t dyIn     = true,
		      Bool_t zpIn     = true,
		      Bool_t qqIn     = true,
		      Bool_t ttIn     = true)
{
  if (dyFit) {
    doTheFit = true;
    dyIn     = true;
    zpIn     = false;
    qqIn     = false;
    ttIn     = false;
  }

  for (Int_t iFit=0; iFit<nFits; iFit++) {

    Double_t zpEvents = rndm->Poisson(zp_pmean);
    Int_t    zpStart  = rndm->Uniform(zpAllData->numEntries());

    RooDataSet zpData("zpData", "zpData", RooArgSet(mass));
  
    Int_t iPeak =       0;
    Int_t iLoop = zpStart;

    while (iPeak < zpEvents) {
      
      Double_t zpMass = zpAllData->get(iLoop)->getRealValue("mass");

      iLoop++;

      if (iLoop == zpAllData->numEntries()) iLoop = 0;

      if (zpMass > zpPeak-zpWindow && zpMass < zpPeak+zpWindow) {

	mass.setVal(zpMass);
	zpData.add(RooArgSet(mass));
	
	iPeak++;
      }
    }

    Double_t dyEvents = rndm->Poisson(dy_pmean);
    Int_t dyStart     = rndm->Uniform(dyAllData->numEntries() - dyEvents);
    
    RooDataSet dyData("dyData", "dyData", RooArgSet(mass));
    
    for (Int_t i=dyStart; i<dyStart+dyEvents; i++) {
      mass.setVal(dyAllData->get(i)->getRealValue("mass"));
      dyData.add(RooArgSet(mass));
    }

    //----------------------------------------------------------------------------
    // QCD dataset
    //----------------------------------------------------------------------------
    N_QQ.setVal(rndm->Poisson(qq_pmean));

    RooDataSet* qqData = qqPdf.generate(mass,N_QQ.getVal());

    //----------------------------------------------------------------------------
    // ttbar dataset
    //----------------------------------------------------------------------------
    N_TT.setVal(rndm->Poisson(tt_pmean));

    RooDataSet* ttData = ttPdf.generate(mass,N_TT.getVal());

    //----------------------------------------------------------------------------
    // All datasets together
    //----------------------------------------------------------------------------
    RooDataSet* data = new RooDataSet("data","data",RooArgSet(mass));
  
    if (dyIn) data->append( dyData);
    if (zpIn) data->append( zpData);
    if (qqIn) data->append(*qqData);
    if (ttIn) data->append(*ttData);

    //--------------------------------------------------------------------------
    // QCD and ttbar histograms
    //--------------------------------------------------------------------------
    TH1F* qqMass = CreateHistogram("qqMass", mass, qqData, kCyan  +1);
    TH1F* ttMass = CreateHistogram("ttMass", mass, ttData, kOrange-3);

    //----------------------------------------------------------------------------
    // Total function
    //----------------------------------------------------------------------------
    if (!dyIn) {
      N_DY .setVal(0);
      N_DY .setConstant(kTRUE);
      dySetConstant(kTRUE);
    }

    if (!zpIn) {
      N_ZP.setVal(0);
      N_ZP.setConstant(kTRUE);
      mV  .setConstant(kTRUE);
      gV  .setConstant(kTRUE);
      sV  .setConstant(kTRUE);
    }

    if (!qqIn) {
      N_QQ  .setVal(0);
      N_QQ  .setConstant(kTRUE);
      qq_be1.setConstant(kTRUE);
      qq_ga1.setConstant(kTRUE);
      qq_be2.setConstant(kTRUE);
      qq_ga2.setConstant(kTRUE);
    }

    if (!ttIn) {
      N_TT  .setVal(0);
      N_TT  .setConstant(kTRUE);
      tt_alp.setConstant(kTRUE);
      tt_bet.setConstant(kTRUE);
      tt_gam.setConstant(kTRUE);
      tt_nor.setConstant(kTRUE);
      tt_mpv.setConstant(kTRUE);
      tt_sig.setConstant(kTRUE);
    }

    Double_t N_ALL = 0;

    N_ALL += N_DY.getVal();
    N_ALL += N_ZP.getVal();
    N_ALL += N_QQ.getVal();
    N_ALL += N_TT.getVal();

    if (N_ALL == 0) {
      printf(" zprime.C -- N_ALL = %f --> EXIT\n", N_ALL);
      return;
    }
    
    RooAddPdf totalPdf("totalPdf", "totalPdf",
		       RooArgList(zpPdf, dyPdf, qqPdf, ttPdf),
		       RooArgSet(N_ZP, N_DY, N_QQ, N_TT));
  
    //--------------------------------------------------------------------------
    // Fit
    //--------------------------------------------------------------------------
    if (dyFit) dySetConstant(kFALSE);
    else       dySetConstant(kTRUE );
    
    Int_t         fitTime   =    0;
    RooFitResult* fitResult = NULL;

    if (doTheFit) {

      fitTime = TTimeStamp().GetSec();

      RooNLLVar nll("nll","nll",totalPdf,*data,kTRUE);
      RooMinuit m(nll);

      m.setStrategy(0);
      m.migrad();

      fitResult = m.save();
      fitTime   = TTimeStamp().GetSec() - fitTime;
    }

    printf("\n");
    printf(" zprime.C -- the fit took %d seconds\n", fitTime);

    if (N_DY.getError() > 0)
      pull_N_DY->Fill((N_DY.getVal() - dy_pmean) / N_DY.getError());

    if (N_ZP.getError() > 0)
      pull_N_ZP->Fill((N_ZP.getVal() - zp_pmean) / N_ZP.getError());

    N_ALL = N_DY.getVal() + N_ZP.getVal() + N_QQ.getVal() + N_TT.getVal();

    Double_t f_ZP = N_ZP.getVal() / N_ALL;
    Double_t f_DY = N_DY.getVal() / N_ALL;

    if (!draw || iFit != 0) continue;

    //--------------------------------------------------------------------------
    // Draw - fit
    //--------------------------------------------------------------------------
    TCanvas* c1 = new TCanvas("c1", "c1");

    c1->SetLogy();

    RooPlot* frame = mass.frame();

    data->plotOn(frame,
		 Name("dat"),
		 MarkerColor(kBlack),
		 DrawOption("pz"));

    totalPdf.plotOn(frame,
		    Name("fit"),
		    LineStyle(kSolid),
		    LineWidth(3));

    SetRooPlot(frame,"m_{#mu#mu}","GeV");
  
    if (ttIn) ttMass->Draw("hist,same");
    if (qqIn) qqMass->Draw("hist,same");

    frame->SetMaximum(5e+5);
    frame->SetMinimum(1e-1);

    TLegend* leg = SetLegend(0.68,0.66,0.90,0.90);

    RooHist*  dataMass = (RooHist* )frame->getHist ("dat");
    RooCurve* fitCurve = (RooCurve*)frame->getCurve("fit");

    leg->AddEntry(dataMass, " MC data", "lp");
    leg->AddEntry(fitCurve, " fit" , "l" );

    if (ttIn) leg->AddEntry(ttMass, " t#bar{t}", "f");
    if (qqIn) leg->AddEntry(qqMass, " QCD"     , "f");

    leg->Draw();

    c1->GetFrame()->DrawClone();
    c1->SaveAs("dy_fit_50_1550.gif");

    //--------------------------------------------------------------------------
    // Draw - residuals
    //--------------------------------------------------------------------------
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 0, 550, 600);

    c2->SetLeftMargin(1.1*c2->GetLeftMargin());

    RooHist* resid = (RooHist*)dataMass->makeResidHist(*fitCurve);
    resid->SetMarkerColor(  50);
    resid->SetMarkerSize (0.75);
    resid->SetMarkerStyle(  20);

    RooPlot* yframe = mass.frame();
    yframe->addPlotable(resid,"pz");
    yframe->Draw();

    if (residYMax > 0) yframe->SetMaximum( residYMax);
    if (residYMin > 0) yframe->SetMinimum(-residYMin);

    SetRooPlot(yframe,"#mu#mu invariant mass","GeV");
    
    PRLabel((TPad*)c2->GetPad(0), "CMS Monte Carlo", 0, 11, 0.05);
    PRLabel((TPad*)c2->GetPad(0), "fit residuals"  , 1, 31, 0.05);
  }

  //----------------------------------------------------------------------------
  // Pulls
  //----------------------------------------------------------------------------
  TFile* out = new TFile("pulls.root", "recreate");
  pull_N_DY->Write();
  pull_N_ZP->Write();
  out->Close();

  //----------------------------------------------------------------------------
  // Print Poisson means for 100 pb-1
  //----------------------------------------------------------------------------
  printf("\n");
  if (zpIn) printf(" zp Poisson mean = %8.1f\n", zp_pmean);
  if (dyIn) printf(" dy Poisson mean = %8.1f\n", dy_pmean);
  if (ttIn) printf(" tt Poisson mean = %8.1f\n", tt_pmean);
  if (qqIn) printf(" qq Poisson mean = %8.1f\n", qq_pmean);
  printf("\n");
}

//------------------------------------------------------------------------------
// dySetConstant
//------------------------------------------------------------------------------
void dySetConstant(Bool_t flag)
{
  A    .setConstant(flag);
  B    .setConstant(flag);
  C    .setConstant(flag);
  Gamma.setConstant(flag);
  M    .setConstant(flag);
  kappa.setConstant(flag);
  theta.setConstant(flag);
}

//------------------------------------------------------------------------------
// CreateHistogram
//------------------------------------------------------------------------------
TH1F* CreateHistogram(Char_t*     name,
		      RooRealVar  x,
		      RooDataSet* dataset,
		      Color_t     color)
{
  TH1F* hist = x.createHistogram(name);
  dataset->fillHistogram(hist,RooArgList(x));

  hist->SetFillStyle( 1001);
  hist->SetFillColor(color);
  hist->SetLineColor(color);

  return hist;
}
