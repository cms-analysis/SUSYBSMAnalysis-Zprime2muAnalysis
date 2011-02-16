// File containing our user-defined functions which we use in root for
// histogram fitting. 
//
// Authors: Jason Mumford, Jordan Tucker, Slava Valuev, UCLA
//

#include <fstream>

#include "TF1.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Functions.h"

using namespace std;

//=============================================================================

double Lorentzian(double *x, double *par) {
  // Simple Lorentzian
  const bool debug = false;

  // par[0] = normalization constant
  // par[1] = FWHM
  // par[2] = mean resonance mass
  double f = par[0]*(par[1]/2.)/(TMath::Pi()*((x[0]-par[2])*
				      (x[0]-par[2])+(par[1]*par[1]/4.)));

  if (debug)
    LogDebug("Lorentzian") << " x " << x[0] << "\n"
			   << " p " << par[0] << " " << par[1] 
			   << " " << par[2] << " f " << f;

  return f;
}

//=============================================================================

double Relativistic_Lorentzian(double *x, double *par){
  // The formula is from "The Z Boson" by C. Caso and A. Gurtu, PDG (2002).
  // Special form of Lorentzian used to fit Z's at LEP.
  // Normalization seems to be OK, but it is not in the PDG, and I cannot
  // derive it analytically (SV, Oct 9, 2003).

  // par[0] = Normalization constant
  // par[1] = FWHM
  // par[2] = mean Z' mass

  double s = x[0]*x[0]; //square of dilepton mass
  double f = (s*par[0]*par[1]*par[1])/
   ((s-par[2]*par[2])*(s-par[2]*par[2]) + (s*s*par[1]*par[1]/(par[2]*par[2])));
  f *= 2.*s/(TMath::Pi()*par[1]*par[2]*par[2]);
  // Alternative expression, giving more symmetric Lorentzian
  // f *= 2./(TMath::Pi()*par[1]);

  return f;
}

//=============================================================================

double BkgndToLorentzian(double *x, double *par) {
  //Function for now one-parameter background, for now ~ 1/x^2
  // Put in factor of 1.e6 so coefficient is in units of TeV^2, which gives
  //  better numerics

  //if (x[0] != 0) return 1.e6*par[0]/(x[0]*x[0]);
  //if (x[0] != 0) return 4.e3*par[0]/(x[0]*x[0]); // (2000-4000 GeV)
  if (x[0] != 0) return 2.25e3*par[0]/(x[0]*x[0]); // (1500-4500 GeV)
  else return 0;
}

//=============================================================================

double expBckg(double *x, double *par) {
  // Attempt to describe the background (Drell-Yan) distributions with
  // exponential.
  // par[0] - normalization
  // par[1] - slope
  // par[2] = Int(exp(par[1]*x[0])) between mass_min and mass_max

  //double bckg = par[0]*exp(par[1]*x[0])/par[2];
  double bckg = par[0]*exp(par[1]*pow(x[0], 0.3))/par[2];

  // For tests of the background shape.
  // TF1 *f1 = new TF1("f1", "1.+0.1*0.5*(1+TMath::Erf((x-700)/50))",400,1600)
  //double coeff = 1. + 0.2*0.5*(1+TMath::Erf((x[0]-700.)/50.));
  //double coeff = 1. + 0.1*0.5*(1+TMath::Erf((x[0]-750.)/150.));
  //bckg = bckg*coeff;

  return bckg;
}

double expBckgNorm(double *x, double *par) {
  // Same as above, but self-normalized to unity.

  const bool debug  = false;
  const  int    npars = 3;
  static double parsave[npars] = {0.};
  static double bkgnorm = 0.;
  const  double mass_min = 1500.;
  const  double mass_max = 4500.;

  double bckg = 0.;

  bool renorm = false;
  /*
  const  double epsilon = 1.e-6;
  for (int ipar = 0; ipar < npars; ipar++) {
    if (fabs(par[ipar]) > epsilon && fabs(parsave[ipar]) > epsilon) {
      if (fabs(parsave[ipar]/par[ipar] - 1.) > epsilon) {
	renorm = true;  break;
      }
    }
    else { // par or parsave very close to 0
      if (fabs(parsave[ipar] - par[ipar]) > epsilon) {
	renorm = true;  break;
      }
    }
  } */
  for (int ipar = 0; ipar < npars; ipar++) {
    if (par[ipar] != parsave[ipar]) {
      renorm = true;  break;
    }
  }

  if (renorm) {
    TF1 *fbkg = new TF1("fbkg", expBckg, mass_min, mass_max, 3);
    fbkg->SetParameters(par[0], par[1], par[2]);
    bkgnorm = fbkg->Integral(mass_min, mass_max);
    delete fbkg;
    for (int ipar = 0; ipar < npars; ipar++)
      parsave[ipar] = par[ipar];

    if (debug)
      LogDebug("expBckgNorm")
	<< "new pars: " << par[0] << " " << par[1] << " " << par[2]
        << "\n bkgnorm = " << bkgnorm;
  }

  if (bkgnorm > 0.)
    bckg = expBckg(x, par)/bkgnorm;
  else
    edm::LogWarning("expBckgNorm") << "bkgnorm = " << bkgnorm << " <= 0";

  return bckg;
}

//=============================================================================

double lorentzianPlusBckg(double *x, double *par) {
  // First three parameters are passed to Lorentzian
  // 4th controls bkgrd. (So pass address of par[3] to BkgndToLorentzian -
  // see root manual pp. 76-77 for similar trick.)

  return Lorentzian(x, par) + BkgndToLorentzian(x, &par[3]);
}

double lorentzianPlusBckgNorm(double *x, double *par) {
  // An attempt to normalize the sum of the two components.

  double par_bckg[1];
  par_bckg[0] = 1. - par[0];
  return Lorentzian(x, par) + BkgndToLorentzian(x, par_bckg);
}

//=============================================================================

double lorentzianPlusExpbckg(double *x, double *par) {
  // First three parameters are passed to Lorentzian.
  // Next three parameters control exponential background.

  return Lorentzian(x, par) + expBckg(x, &par[3]);
}

double lorentzianPlusExpbckgNorm(double *x, double *par) {
  // An attempt to normalize the sum of the two components.
  double func = 0.;

  double par_bckg[3];
  par_bckg[0] = 1. - par[0];
  par_bckg[1] = par[3];
  par_bckg[2] = par[4];

  double sigfunc = Lorentzian(x, par);
  double bkgfunc = expBckg(x, par_bckg);
  func = sigfunc + bkgfunc;

#define RENORM
#ifdef RENORM
  // This numerical renormalization here is not really needed because both
  // Lorentzian and expBckg are analytically normalized to 1.  But I (Sl.)
  // think it may give a bit more accurate answer due to (-inf, mass_min)
  // and (mass_max; +inf) tails of Lorentzian.
  const  int    npars = 5;
  static double parsave[npars] = {0.};
  static double signorm = 0.0, bkgnorm = 0.0;
  const  double mass_min = 1500.;
  const  double mass_max = 4500.;
  if (par[1] != parsave[1] || par[2] != parsave[2]) {
    TF1 *flor = new TF1("flor", Lorentzian, mass_min, mass_max, 3);
    flor->SetParameters(1., par[1], par[2]);
    signorm = flor->Integral(mass_min, mass_max);
    parsave[1] = par[1];
    parsave[2] = par[2];
    delete flor;
  }

  if (par[3] != parsave[3] || par[4] != parsave[4]) {
    TF1 *fbkg = new TF1("fbkg", expBckg, mass_min, mass_max, 3);
    fbkg->SetParameters(1., par[3], par[4]);
    bkgnorm = fbkg->Integral(mass_min, mass_max);
    parsave[3] = par[3];
    parsave[4] = par[4];
    delete fbkg;
  }

  if (signorm > 0. && bkgnorm > 0.)
    func = sigfunc/signorm + bkgfunc/bkgnorm;
  else
    edm::LogWarning("lorentzianPlusExpbckgNorm")
      << "signorm = " << signorm << " bkgnorm = " << bkgnorm;
#endif

  return func;
}

double lorentzianPlusExpbckgN(double *x, double *par) {
  // lorentzianPlusExpbckgNorm times the normalization; used for plots.

  return par[0]*lorentzianPlusExpbckgNorm(x, &par[1]);
}

//=============================================================================

double lorengaufun_fixedsigma(double *x, double *par) {

  //Fit parameters:
  //par[0]=Normalization constant for Lorentzian
  //par[1]=FWHM for Lorentzian
  //par[2]=mean for Lorentzian

  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

  int bins = 100;      // number of convolution steps, bins
  double s = 5.0;      // convolution extends to +-s Gaussian sigmas

  double xx, xar[1];
  double sum = 0.0;

  // Mass resolution (sigma of a Gaussian) for a given mass.
  // double sigma = 48.; // GMR, for 1 TeV (sigma= 4.8%)
  double sigma = 39.; // Opt, for 1 TeV (sigma= 3.9%)
  // double sigma = 300.; // GMR, for 3 TeV (sigma=10.0%)
  // double sigma = 170.; // Opt, for 3 TeV (sigma= 5.7%)
  // double sigma = 700.; // GMR, for 5 TeV (sigma=14.0%)
  // double sigma = 425.; // Opt, for 5 TeV (sigma= 8.5%)

  // Range of convolution integral from -s*sigma to +s*sigma
  double xlow = x[0] - s*sigma;
  double xupp = x[0] + s*sigma;
  double width = (xupp-xlow)/bins; //width

  // Convolution integral of Lorentzian and Gaussian by sum (called
  // "extended trapezoidal rule" in English)
  for (int i = 1; i < bins; i++) {
    xx = xlow + (i*width);
    xar[0] = xx;
    sum += Lorentzian(xar, par) * TMath::Gaus(x[0], xx, sigma);
  }
  // Find the value of the first and last terms and add them on to sum.
  xar[0] = xlow;
  sum += 0.5*Lorentzian(xar, par)*TMath::Gaus(x[0], xlow, sigma);
  xar[0] = xupp;
  sum += 0.5*Lorentzian(xar, par)*TMath::Gaus(x[0], xupp, sigma);

  double Area = width*sum*invsq2pi/sigma;
  return Area;
}

double lorengaufun(double *x, double *par) {

  //Fit parameters:
  //par[0]=Normalization constant for Lorentzian
  //par[1]=FWHM for Lorentzian
  //par[2]=mean for Lorentzian

  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

  int bins = 100;      // number of convolution steps, bins
  double s = 5.0;      // convolution extends to +-s Gaussian sigmas

  double xx, xar[1];
  double sum = 0.0;

  // Range of convolution integral from -s*sigma to +s*sigma
  double sigma = mass_resolution(x, par);
  double xlow = x[0] - s*sigma;
  double xupp = x[0] + s*sigma;
  double width = (xupp-xlow)/bins; //width

  // Convolution integral of Lorentzian and Gaussian by sum (called
  // "extended trapezoidal rule" in English)
  for (int i = 1; i < bins; i++) {
    xx = xlow + (i*width);
    xar[0] = xx;
    sigma = mass_resolution(xar, par);
    sum += Lorentzian(xar, par) * TMath::Gaus(x[0], xx, sigma)/sigma;
  }
  // Find the value of the first and last terms and add them on to sum.
  xar[0] = xlow;
  sigma = mass_resolution(xar, par);
  sum += 0.5*Lorentzian(xar, par)*TMath::Gaus(x[0], xlow, sigma)/sigma;
  xar[0] = xupp;
  sigma = mass_resolution(xar, par);
  sum += 0.5*Lorentzian(xar, par)*TMath::Gaus(x[0], xupp, sigma)/sigma;

  double Area = width*sum*invsq2pi;
  return Area;
}

double mass_resolution(double *x, double *par) {
  // Returns mass resolution (sigma of a Gaussian, in GeV) for a given mass.
  double sigma;
  double mass = x[0];

  // JMTBAD this needs to be revisited
  // Parameterization below is the result of fitting the mean values of
  // (M_res - M_Z')/M_Z' at six points (0., 0.2, 0.4, 1, 3, and 5 TeV)
  // with two polynomials, one (2nd degree) at low masses, one (1st degree)
  // at high masses.
  //
  // Resolution values for Z's are calculated within (M_Z' +/- 0.1*M_Z').
  // Numbers for tracker and muon "long term" misalignment scenarios in
  // ORCA_8_13_2:
  // Mass value (GeV):   200,    400,    1000,   3000,   5000
  // Resolution:        0.021,  0.031,  0.042,  0.065,  0.090
  // Uncertainty:       0.001,  0.001,  0.001,  0.001,  0.002
  //
  // To obtain parameterizations, do the following:
  /* double x1[3] = {0., 200., 400.}
     double y1[3] = {0., 0.021, 0.032}
     double ex1[3] = {0., 0., 0.}
     double ey1[3] = {0.001, 0.001, 0.001}
     TGraph *gre1 = new TGraphErrors(3, x1, y1, ex1, ey1)
     gre1->Fit("pol2")
     double x2[4] = {400., 1000., 3000., 5000.}
     double y2[4] = {0.031, 0.042, 0.065, 0.090}
     double ex2[4] = {0., 0., 0., 0.}
     double ey2[4] = {0.001, 0.001, 0.001, 0.002}
     TGraph *gre2 = new TGraphErrors(4, x2, y2, ex2, ey2)
     gre2->Fit("pol1")
  */
  /* To plot them, do
     double x[6] = {0., 200., 400., 1000., 3000., 5000.}
     double y[6] = {0., 0.021, 0.031, 0.042, 0.065, 0.090}
     double ex[6] = {0., 0., 0., 0., 0., 0.}
     double ey[6] = {0., 0.001, 0.001, 0.001, 0.001, 0.002}
     TGraph *gre = new TGraphErrors(6,x,y,ex,ey)
     gre->Draw("A*")
     TF1 *f1 = new TF1("f1", "0.02751+1.258e-5*x", 447., 5000.)
     f1->Draw("same")
     TF1 *f2 = new TF1("f2", "0.00013*x-1.25e-7*x*x", 0., 447.)
     f2->Draw("same")
  */
  // ORCA_8_7_1
  //if (mass < 887.) sigma =          6.56e-5*mass - 3.07e-8*mass*mass;  // TMR
  //else             sigma = 0.0250 + 1.09e-5*mass - 8.12e-10*mass*mass;
  //if (mass < 567.) sigma =          7.e-5*mass   - 2.5e-8*mass*mass;  // GMR
  //else             sigma = 0.0244 + 1.33e-5*mass - 9.97e-10*mass*mass;
  //
  // ORCA_8_13_2
  // "First data" misalignment, low lumi pile-up.  Valid up to 3 TeV.
  // sigma =  0.03113 + 1.137e-4*mass - 1.563e-8*mass*mass;
  // "Long term" misalignment, both low-lumi and high-lumu pile-up.
  // if (mass < 447.) sigma =           1.3e-4*mass - 1.25e-7*mass*mass;
  // else             sigma = 0.02751 + 1.258e-5*mass;
  //
  // ORCA_6_3_1
  // Parameterization below is the result of fitting the mean values of
  // (M_rec - M_gen)/M_gen at four points ((M_gen +/- 3*FWHM) at 0.4, 1, 3
  // and 5 TeV).
  // Mass value (GeV):   400,   1000,   3000,   5000
  // TMR resolution:   0.019,  0.037,  0.063,  0.085
  // GMR resolution:   0.020,  0.043,  0.106,  0.130
  // sigma = 0.0131 + 2.16e-5*mass - 1.46e-9*mass*mass; // TMR
  // sigma = 0.00037 + 4.85e-5*mass - 4.52e-9*mass*mass; // GMR

  // Resolution in CMSSW_1_3_1 (ideal alignment) looks rather similar to
  // that in ORCA_6_3_1 (too bad!), so just use 6_3_1 parameterization for
  // now.  We should re-tune in the future.
  // sigma = 0.0131 + 2.16e-5*mass - 1.46e-9*mass*mass;

  // Numbers for 100pb-1 misalignment scenario in CMSSW_1_6_7:
  // Mass value (GeV):   300,    600,    1000,   1500,   2000,   2500,   3000
  // Resolution:        0.036   0.062,  0.081,  0.105,  0.122,  0.137,  0.165
  // Uncertainty:       0.001   0.001,  0.001,  0.001,  0.001,  0.002,  0.003
  /* To obtain parameterization:
     double x[7] = {300., 600., 1000., 1500., 2000., 2500., 3000.}
     double y[7] = {0.036, 0.062, 0.081, 0.105, 0.122, 0.137, 0.165}
     double ex[7] = {0., 0., 0., 0., 0., 0., 0.}
     double ey[7] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.003}
     TGraph *gre = new TGraphErrors(7, x, y, ex, ey)
     gre->Fit("pol2")

     To plot, fill x, y, ex, and ey arrays as above, and then:
     TGraph *gre = new TGraphErrors(7, x, y, ex, ey)
     gre->Draw("A*")
     TF1 *f1 = new TF1("f1", "0.01967+6.857e-5*x-8.074e-9*x*x", 200., 3500.)
     f1->Draw("same")
  */
  sigma = 0.0197 + 6.86e-5*mass - 8.07e-9*mass*mass;

  sigma *= mass;
  return sigma;
}

//=============================================================================

double lorengauPlusBckg(double *x, double *par) {
  // First three parameters are passed to lorengaufun, 4th controls background.

  return lorengaufun(x, par) + BkgndToLorentzian(x, &par[3]);
}

double lorengauPlusExpbckg(double *x, double *par) {
  // First three parameters are passed to Lorengaufun.
  // Next three parameters control exponential background.

  return lorengaufun(x, par) + expBckg(x, &par[3]);
}

// Accuracy of the Integral() (default is 10^{-6}).
const double EPS = 0.000001;

double lorengauPlusExpbckgNorm(double *x, double *par) {
  // An attempt to normalize the sum of the two components, the signal and
  // the background.
  static const bool debug = false;
#if defined(EML_1) || defined(EML_2)
  const  int npars = 6;
#else
  const  int npars = 5;
#endif
  static double parsave[npars] = {0.};
  static double signorm = 0.0, bkgnorm = 0.0;
#ifdef RENORMALL
  static double funorm = 0.0;
#endif
  bool renorm_sig = false, renorm_bkg = false;
  // const  double mass_min =  400.; // for 1 TeV
  // const  double mass_max = 1600.;
  const  double mass_min =  600.; // for 1.2 TeV
  const  double mass_max = 1600.;
  // const  double mass_min = 1000.; // for 1.5 TeV
  // const  double mass_min =  600.;
  // const  double mass_max = 2200.;
  // const  double mass_min = 1000.; // for 2 TeV
  // const  double mass_max = 3000.;
  // const  double mass_min = 1500.; // for 3 TeV
  // const  double mass_max = 4500.;
  // const  double mass_min = 3000.; // for 5 TeV
  // const  double mass_max = 7000.;

  double func = 0.;

  // Array with parameters for expBckg
  double par_bckg[3];
#ifndef EML_2
  par_bckg[0] = 1. - par[0];
#else
  par_bckg[0] = par[5];
#endif
  par_bckg[1] = par[3];
  par_bckg[2] = par[4];

  double sigfunc = lorengaufun(x, par);
  double bkgfunc = expBckg(x, par_bckg);

#ifdef RENORMALL
  // Integrate the sum of the two components.  This is OK since the
  // background has already been normalized to unity, hence
  // (Int[(1-N)*bckg + N*sig] = (1-N)*Int(bckg) + N*Int(sig) =
  //  1 - N + N*Int(sig) = 1 ==> Int(sig) = 1)
  for (int ipar = 0; ipar < npars; ipar++) {
    if (par[ipar] != parsave[ipar]) {
      renorm_sig = true;
      break;
    }
  }
  if (renorm_sig) {
    TF1 *flor = new TF1("flor", lorengaufun, mass_min, mass_max, 3);
    TF1 *fbkg = new TF1("fbkg", expBckg,     mass_min, mass_max, 3);
    flor->SetParameters(par[0], par[1], par[2]);
    fbkg->SetParameters(par_bckg[0], par_bckg[1], par_bckg[2]);
    signorm = flor->Integral(mass_min, mass_max);
    bkgnorm = fbkg->Integral(mass_min, mass_max);
    funorm  = signorm + bkgnorm;
    delete flor;
    delete fbkg;

    if (debug) {
      ostringstream out;
      out << "Renorm: " << endl;
      out << "old pars: ";
      for (int ipar = 0; ipar < npars; ipar++)
	out << parsave[ipar] << " ";
      out << endl;
      out << "new pars: ";
      for (int ipar = 0; ipar < npars; ipar++)
	out << par[ipar] << " ";
      out << endl;
      out << " signorm = " << signorm << " bkgnorm = " << bkgnorm
	   << " funorm = " << funorm;
      LogDebug("lorengauPlusExpbckgNorm") << out.str();
    }

    for (int ipar = 0; ipar < npars; ipar++)
      parsave[ipar] = par[ipar];
  }

  if (funorm > 0.)
    func = (sigfunc + bkgfunc)/funorm;
  else
    edm::LogWarning("lorengauPlusExpbckgNorm")
      << "funorm = " << funorm << " <= 0";
#else
  // Integrate and renormalize only the function parameters of which have
  // changed.  Since all background parameters (except for normalization
  // coefficient) are currently fixed, this is a faster way to proceed
  // (and we do not rely on the fact that the background has already been
  // appropriately normalized).

  for (int ipar = 1; ipar < npars; ipar++) { // exclude par[0] = norm
    if (par[ipar] != parsave[ipar]) {
      if      (ipar == 1 || ipar == 2) renorm_sig = true;
      else if (ipar == 3 || ipar == 4) renorm_bkg = true;
    }
  }

  ostringstream out;

  if (renorm_sig) {
    if (debug) out << "Renorm. signal: " << endl;
    TF1 *flor = new TF1("flor", lorengaufun, mass_min, mass_max, 3);
    double par_sig[3] = {1., par[1], par[2]};
    signorm = flor->Integral(mass_min, mass_max, par_sig, EPS);
    parsave[1] = par[1];
    parsave[2] = par[2];
    delete flor;
  }

  if (renorm_bkg) {
    if (debug) out << "Renorm. background: " << endl;
    TF1 *fbkg = new TF1("fbkg", expBckg, mass_min, mass_max, 3);
    double par_bkg[3] = {1., par[3], par[4]};
    bkgnorm = fbkg->Integral(mass_min, mass_max, par_bkg, EPS);
    parsave[3] = par[3];
    parsave[4] = par[4];
    delete fbkg;
  }
  if (debug && (renorm_sig || renorm_bkg)) {
    out << "new pars: ";
    for (int ipar = 0; ipar < npars; ipar++)
      out << par[ipar] << " ";
    out << endl;
    out << " signorm = " << signorm << " bkgnorm = " << bkgnorm;
  }
  
  if (debug) LogDebug("lorengauPlusExpbckgNorm") << out.str();

  if (signorm > 0. && bkgnorm > 0.)
    func = sigfunc/signorm + bkgfunc/bkgnorm;
  else
    edm::LogWarning("lorengauPlusExpbckgNorm") 
      << ": signorm = " << signorm << " bkgnorm = " << bkgnorm;
#endif

#ifdef EML_1
  func *= par[5]; // normalization
#endif

  return func;
}

double lorengauPlusExpbckgN(double *x, double *par) {
  // lorengauPlusExpbckgNorm times normalization; used for plots.

  return par[0]*lorengauPlusExpbckgNorm(x, &par[1]);
}

//=============================================================================

double expgaufun(double *x, double *par) {

  //returns convolution of exponential and gaussian
  //parameters:
  //par[0]=Extra Normalization constant for exponential, beyond 1/tau
  //par[1]=tau for exponential
  //par[2]=sigma for Gaussian at M=2000 GeV, dimensionless
  //par[3]=slope in sigma, Gev per 1000 GeV (typically ~ 0.025/1000

  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

  //!!! bins=10000 is very slow but gives nice precision for plot
  int bins = 10000;      // number of convolution steps, bins
  double s = 10.0;      // convolution extends to +-s Gaussian sigmas

  double xx, sig;
  double sum = 0.0;

  //Range of convolution integral from -s*sigma to +s*sigma
  double sigma = par[2];
  double slope = par[3];
  double xlow = x[0] - s*sigma*2000;
  double xupp = x[0] + s*sigma*2000;
  double width = (xupp-xlow)/bins; //width

  double tau = par[1];
  double tauinv = 1./tau;
  // Convolution integral of exponential and Gaussian by sum (called
  // "extended trapezoidal rule" in English)
  for (int i = 1; i < bins; i++) {
    xx = xlow + (i*width);
    //Use sigma at mass xx
    sig = sigma + (xx - 2000)*slope;
    if(sig<0.005) sig=0.005;
    sig = sig*xx; //convert from relative to absolute
    // LogDebug("expgaufun") << "sigma " << sigma << " slope " << slope
    //                       << " xx " << xx << " sig " << sig;
    if(xx>0) sum += tauinv*TMath::Exp(-xx/tau)*TMath::Gaus(x[0], xx, sig)/sig;
  }
  //find the value of the first and last terms and add them on to sum.
  xx = xlow;
  sig = sigma + (xx - 2000)*slope;
  if(sig<0.005) sig=0.005;
  sig = sig*xx; //convert from relative to absolute
  if(xx>0) sum += tauinv*TMath::Exp(-xx/tau) * TMath::Gaus(x[0], xx, sig)/sig;
  xx = xupp;
  sig = sigma + (xx - 2000)*slope;
  if(sig<0.005) sig=0.005;
  sig = sig*xx; //convert from relative to absolute
  if(xx>0) sum += tauinv*TMath::Exp(-xx/tau) * TMath::Gaus(x[0], xx, sig)/sig;

  double Area = par[0]*width*sum*invsq2pi;
  return Area;
}

//=============================================================================

TLorentzVector LorentzBoost(const TLorentzVector& boost_frame, const TLorentzVector& to_be_boosted) {
  TLorentzVector v = to_be_boosted;
  v.Boost(-boost_frame.BoostVector());
  return v;
}

double cos_angle(const TVector3& v1, const TVector3& v2) {
  return v1 * v2 / v1.Mag() / v2.Mag();
}

double cos_angle(const TLorentzVector& v1, const TLorentzVector& v2) {
  return cos_angle(v1.Vect(), v2.Vect());
}

//=============================================================================

double probks(const double alam) {
  // Kolmogorov-Smirnov probability function (D > observed).
  // From W. Press et al., "Numerical recipes in C", p.626.
  const double eps1 = 0.001;
  const double eps2 = 1.0e-8;
  double fac = 2.0, sum = 0.0, term, termbf = 0.0;

  double a2 = -2.0*alam*alam;
  for (int j = 1; j <= 100; j++) {
    term = fac*exp(a2*j*j);
    sum += term;
    if (fabs(term) <= eps1*termbf || fabs(term) <= eps2*sum)
      return sum;
    fac = -fac; // Alternating signs in sum.
    termbf = fabs(term);
  }
  return 1.0; // Get here only by failing to converge.
}
