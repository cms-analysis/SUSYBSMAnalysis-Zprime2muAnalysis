#ifndef RooConfigure_C
#define RooConfigure_C

TRandom1* rndm = new TRandom1();

//------------------------------------------------------------------------------
// Pull histograms
//------------------------------------------------------------------------------
TH1F* pull_N_DY = new TH1F("pull_N_DY", "pull_N_DY", 1e3, -5, 5);
TH1F* pull_N_ZP = new TH1F("pull_N_ZP", "pull_N_ZP", 1e3, -5, 5);

//------------------------------------------------------------------------------
// Global variables
//------------------------------------------------------------------------------
Double_t kFactor      =  1.35;
Double_t intLumi      =   0.1;  // [fb]
Double_t zp_XSec      =   624;  // [fb]
Double_t dy_XSec      = 1.8e6;  // [fb]
Double_t zp_effFilter =   1.0;
Double_t dy_effFilter =  0.42;  // eta acceptance
Double_t zp_effReco   =   0.7;
Double_t dy_effReco   =   0.8;  // pt > 20
  
Double_t zp_pmean     = zp_XSec * kFactor * zp_effFilter * zp_effReco * intLumi;
Double_t dy_pmean     = dy_XSec * kFactor * dy_effFilter * dy_effReco * intLumi;

Double_t zpPeak       =  1000;
Double_t zpWindow     =   400;

Double_t xMin         =    50;
Double_t xMax         =  1550;

RooRealVar mass  ("mass"  , "mass"  , xMin, xMax);
RooRealVar weight("weight", "weight",    0,    1);

mass.setRange("fitRange", xMin, xMax);

//------------------------------------------------------------------------------
// Z' dataset
//------------------------------------------------------------------------------
RooDataSet* zpAllData =
RooDataSet::read("input/Zssm1000.txt", mass);

//------------------------------------------------------------------------------
// DY dataset
//------------------------------------------------------------------------------
RooDataSet* dyAllData =
RooDataSet::read("input/DrellYan_ll_40.txt", mass);

//------------------------------------------------------------------------------
// Get QCD and ttbar number of events
//------------------------------------------------------------------------------
TFile* bkgFile = new TFile("input/backgrounds.root", "read");

TH1F*  qqTrue = (TH1F*)bkgFile->Get("qqMass");

Double_t qq_pmean = qqTrue->Integral(qqTrue->FindBin(xMin),
				     qqTrue->FindBin(xMax));

TH1F*  ttTrue = (TH1F*)bkgFile->Get("ttMass");

Double_t tt_pmean = ttTrue->Integral(ttTrue->FindBin(xMin),
				     ttTrue->FindBin(xMax));

bkgFile->Close();

//------------------------------------------------------------------------------
// QCD function
//------------------------------------------------------------------------------
RooRealVar qq_be1("qq_be1", "qq_be1", 3.40194e+01);
RooRealVar qq_ga1("qq_ga1", "qq_ga1", 5.94078e+21);
RooRealVar qq_be2("qq_be2", "qq_be2", 3.46213e+01);
RooRealVar qq_ga2("qq_ga2", "qq_ga2", 9.99886e-02);

RooQCD qqPdf("qqPdf", "qqPdf",
	     mass, qq_be1, qq_ga1, qq_be2, qq_ga2);

//------------------------------------------------------------------------------
// ttbar function
//------------------------------------------------------------------------------
RooRealVar tt_alp("tt_alp", "tt_alp",  9.22795e-03);
RooRealVar tt_bet("tt_bet", "tt_bet", -3.85820e+01);
RooRealVar tt_gam("tt_gam", "tt_gam",  4.16641e+01);
RooRealVar tt_nor("tt_nor", "tt_nor",  6.28152e+01);
RooRealVar tt_mpv("tt_mpv", "tt_mpv",  7.78551e+01);
RooRealVar tt_sig("tt_sig", "tt_sig",  1.38399e+01);

RooTTbar ttPdf("ttPdf", "ttPdf",
	       mass, tt_alp, tt_bet, tt_gam, tt_nor, tt_mpv, tt_sig);

//------------------------------------------------------------------------------
// DY function
//------------------------------------------------------------------------------
RooRealVar A    ("A"    , "A"    ,  3.8745e+04, 0.0, 1e6);
RooRealVar B    ("B"    , "B"    ,  7.6023e+02, 0.0, 1e4);
RooRealVar C    ("C"    , "C"    ,  4.3701e+01, 0.0, 1e2);
RooRealVar Gamma("Gamma", "Gamma",  4.4435e+00);
RooRealVar M    ("M"    , "M"    ,  9.1359e+01);
RooRealVar kappa("kappa", "kappa", -3.7505e-04);
RooRealVar theta("theta", "theta",  1.4817e-02);

RooDrellYanC0 dy0("dy0", "dy0", mass, M,Gamma, theta);
RooDrellYanC1 dy1("dy1", "dy1", mass, M,Gamma, theta);
RooDrellYanC2 dy2("dy2", "dy2", mass, kappa);

RooAddPdf dyPdf("dyPdf", "dyPdf",
		RooArgList(dy0, dy1, dy2),
		RooArgList(A, B, C));

//------------------------------------------------------------------------------
// Z' function
//------------------------------------------------------------------------------
RooRealVar mV("mV", "mV", 1.0000e+03, 800, 1200);
RooRealVar gV("gV", "gV", 5.7193e+01,   0,  300);
RooRealVar sV("sV", "sV", 3.1278e-01,   0,  100);

RooVoigtian zpPdf("zpPdf", "zpPdf", mass, mV, gV, sV);

//------------------------------------------------------------------------------
// N
//------------------------------------------------------------------------------
RooRealVar N_DY("N_DY", "N_DY", 1e5, 0, 2e5);
RooRealVar N_ZP("N_ZP", "N_ZP", 1e2, 0, 200);
RooRealVar N_QQ("N_QQ", "N_QQ", 1e2);
RooRealVar N_TT("N_TT", "N_TT", 1e2);

#endif
