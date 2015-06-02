// File containing our user-defined functions which we use in root for
// histogram fitting. 
//
// Authors: Jason Mumford, Jordan Tucker, Slava Valuev, UCLA
//

#include <fstream>

#include "TH2F.h"
#include "TFMultiD.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Functions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"

using namespace std;

AsymFitManager asymFitManager;

std::ostream& operator<< (std::ostream& out, const AsymFitData& data) {
  out << "pT = " << data.pT << ", rap = " << data.rapidity;
  out << ", phi = " << data.phi << ", mass = " << data.mass;
  out << ", cos_true " << data.cos_true << ", cos_cs = " << data.cos_cs;
  out << ", phi_cs = " << data.phi_cs << std::endl;
  out << "mistag_true = " << data.mistag_true 
      << ", mistag_cs = " << data.mistag_cs;
  return out;
}

//=============================================================================

double asym_2_PDF(double *x, double *par) {
  // Version of previous function, Asym_3_PDF, with b=1 fixed.
  // So b_func = 3/8 

  double func;
  const bool debug = false;

  func = par[0]*((3./8.) + par[1]*x[0] + (3./8.)*x[0]*x[0]);
  if (debug)
    LogDebug("asym_2_PDF")
      << "Parameters N: " << par[0] << "  A_fb: " << par[1]
      << "  x: " << x[0] << " f " << func;

  return func;
}

//=============================================================================

// renamed from asym_3_PDF for flexibility in choosing the angular 
// fit function
double asym_3_PDF(double *x, double *par) {
  // Returns P(x) proportional to 1 + ax + bx^2, for use in fitting
  //   angular distributions where x = cos(theta), -1 <= x <= 1.
  // In this case, it is useful to rewrite P as
  // K (B + A*x + B*b*x^2), where the term in parenthesis is normalized
  // to unit integral and 
  // K = par[0] is the normalization constant,
  // A = par[1] is the forward backward asymmetry,
  // b = par[2] controls difference in magnitude between the
  //            x^0 term and the x^2 term. (b=1 for typical S.M. distributions)
  // Input x = x[0]. (It's in an array for compatibility with fits)
  // For b fixed to 1, see following function Asym_2_PDF

  double func, b_func;
  const bool debug = false;

  b_func = 3./(2.*(3.+par[2]));
  func = par[0]*(b_func + par[1]*x[0] + b_func*par[2]*x[0]*x[0]);
  if (debug)
    LogDebug("asym_3_PDF")
      << "Parameters N: " << par[0] << "  b: " << par[2]
      << "  x: " << x[0] << "  A_fb: " << par[1]
      << "  f: " << func;

  if(func<asymFitManager.epsilon()) func=asymFitManager.epsilon();
  return func;
}

//=============================================================================

double asym_mistag_PDF(double *x, double *pars) {
  double rap = pars[3];
  double cos = x[0];
  double neg_cos = -cos;
  double w = 0.;
  if (asymFitManager.correct_mistags()) { w = mistagProb(rap, cos, asymFitManager.peak_mass()); }
  double func = (1. - w)*asym_3_PDF(&cos, pars) + 
    w*asym_3_PDF(&neg_cos, pars);

  return func;
}

//=============================================================================

double GravitonCos_2_PDF(double *x, double *par) {
  // Version of previous function,  GravitonCos_3_PDF, with a=0 fixed.
  // So only 1+b*x^4, and a=0 in C_norm.
  // So now par[1] is coefficient of x^4 rather than x^2.

  double func, C_norm;
  const bool debug = false;

  C_norm = 0.5/(1. + par[1]/5.);
  double xsq = x[0]*x[0];
  func = par[0]*C_norm*(1. + par[1]*xsq*xsq);
  if (debug)
    LogDebug("GravitonCos_2_PDF")
      << "Parameters N: " << par[0] << "  b: " << par[1]
      << "  x: " << x[0] 
      << "  f: " << func;

  if(func<asymFitManager.epsilon()) func=asymFitManager.epsilon();
  return func;
}

//=============================================================================

double GravitonCos_th_PDF(double *x, double *par) {
  // Returns P(x) =  K*C_norm*((s+1) - (3s)*x^2 + (4s-1)x^4),
  //   for use in fitting angular distributions where
  //   x = cos(theta), -1 <= x <= 1.
  // C_norm is normalization for integral from -1 to 1:
  //    C_norm = 5/8 = 0.625
  // K = par[0] is the additional input normalization constant
  // s = par[1] is the parameter (ratio of qqbar to gg origin)
  // Input x = x[0]. (It's in an array for compatibility with fits)
  double func, C_norm;
  const bool debug = false;

  C_norm = 0.625;
  double xsq = x[0]*x[0];
  func = par[0]*C_norm*(par[1]+1. - 3.*par[1]*xsq + (4.*par[1]-1.)*xsq*xsq);
  if (debug)
    LogDebug("GravitonCos_th_PDF")
      << "Parameters N: " << par[0] << "  s: " << par[1]
      << "  x: " << x[0] << "  f: " << func << endl;

  if(func<asymFitManager.epsilon()) func=asymFitManager.epsilon();
  return func;
}

//=============================================================================

double GravitonCos_3_PDF(double *x, double *par) {
  // Returns P(x) =  K*C_norm*(1 + ax^2 + bx^4), for use in fitting
  //   angular distributions where x = cos(theta), -1 <= x <= 1.
  // C_norm is normalization for integral from -1 to 1:
  //    C_norm = 0.5*(1 + a/3 + b/5)
  // K = par[0] is the additional input normalization constant
  // A = par[1] is the ratio of x^2 term to x^0 term
  // b = par[2] is the ratio of x^4 term to x^0 term
  // Input x = x[0]. (It's in an array for compatibility with fits)
  // For a fixed to 0, see following function GravitonCos_2_PDF

  double func, C_norm;
  const bool debug = false;

  C_norm = 0.5/(1. + par[1]/3. + par[2]/5.);
  double xsq = x[0]*x[0];
  func = par[0]*C_norm*(1. + par[1]*xsq + par[2]*xsq*xsq);
  if (debug)
    LogDebug("GravitonCos_3_PDF")
      << "Parameters N: " << par[0] << "  a: " << par[1]
      << "  b " << par[2] << "  x: " << x[0] 
      << "  f: " << func;

  if(func<asymFitManager.epsilon()) func=asymFitManager.epsilon();
  return func;
}

//=============================================================================
double asym2D(double *x, double *par) {
  const bool debug = false;
  bool paramsChanged = false;
  cout<<paramsChanged<<endl;//raffa
  static bool first_event = true;
  static double parsave[3];
  if (first_event) {
    for (int i = 0; i < 3; i++) parsave[i] = 0;
    first_event = false;
  }

  // Enter into this loop if parameters have changed or this is first time
  // in loop.
  if ( par[0]!=parsave[0] || par[1]!=parsave[1] || par[2]!=parsave[2] ) {
    paramsChanged = true;
    for (int i = 0; i < 3; i++) parsave[i] = par[i];
    if (asymFitManager.debug())
      edm::LogVerbatim("asym2D") 
	<< "Norm = " << par[0] << ", A_FB = " << par[1] << ", b = " << par[2];
  }

  double cos_cs = x[0];
  double neg_cos_cs = -x[0];
  double rap = x[1];
  //double mass = x[2]; // passed in just for mistagProb; will not be fit to.
  // JMTBAD fake the mass for the 2D fit for now
  double mass = asymFitManager.peak_mass();
  //double pL     = x[2];
  //double qpL    = x[3];
  //double pT     = x[4];

  static double etalim = asymFitManager.eta_lim();

  // See if dimuon is within detector acceptance.
  if (fabs(rap) > fabs(rapMaxAccept(&cos_cs, &etalim))) return asymFitManager.epsilon();

  double w = 0.;
  // A Few options for calculating mistag probability:
  // 1. Default version using rapidity and cos_cs
  if (asymFitManager.correct_mistags()) { w = mistagProb(rap, cos_cs, mass); }
  // 2. Most precise version using quark and dilepton pL
  //if (asymFitManager.correct_mistags()) { w = mistagProbVsPL(qpL, pL); }
  // 3. Version using dilepton pT and rap
  //if (asymFitManager.correct_mistags()) { w = mistagProbVsPtRap(pT, rap); }
  // 4. Version using dilepton pL
  //if (asymFitManager.correct_mistags()) { w = mistagProbVsPL(pL); }
  // 5. Version using only dilepton rapidity
  //if (asymFitManager.correct_mistags()) { w = mistagProbVsRap(rap); }
 

  double func = (1.-w)*asym_3_PDF(&cos_cs, par) + 
    w*asym_3_PDF(&neg_cos_cs, par);

  func *= yDistRTC(&rap, par);

  // Some debug statements
  if (debug) {
    ostringstream out;
    out << "x = ";
    for (int i = 0; i < 2; i++) out << x[i] << " ";
    out << "  ypdf = " << yDistRTC(&rap, par) << " w = " << w 
	<< " asymPDF = " << asym_3_PDF(&cos_cs, par) << endl;
    out << "  N = " << par[0] << ", b = " << par[2] << ", A_fb = " 
	<< par[1] << ", f = " << func << endl;
    LogDebug("asym2D") << out.str();
  }

  return func;
}
  

//=============================================================================

double execAsym2D(double *x, double *par) {
  // Same as asym2D, except that this function can be used by the 
  // unbinned fitter, since it uses a monte-carlo integration method to 
  // normalize the function.
  const bool debug = false;
  static bool first_event = true;
  static int nchange = 0;
  static unsigned int nNormPoints = 0;
  static double lo_lim[2] = {-1., -asymFitManager.max_rap()};
  static double hi_lim[2] = { 1.,  asymFitManager.max_rap()};

  if (debug)
    LogDebug("execAsym2D") << "cos, rap = " << x[0] << ", " << x[1];

  static double anorm, parsave[3];
  // Start with asym2D function.
  double func = asym2D(x, par);

  // This flag is used to keep track of when the parameters have changed.
  // We only need to recalculate the normalization when the parameters
  // have changed.
  bool paramsChanged = false;

  // Initialize parameters on first call
  if (first_event) {
    for (int i = 0; i < 3; i++)
      parsave[i] = 0;
    nNormPoints = (unsigned int) par[3];

    if (asymFitManager.debug()) {
      ostringstream out;
      out << "initial parameters: ";
      for (int i = 0; i < 4; i++)
	out << par[i] << " ";
      out << endl;
      out << "  lo_lim = " << lo_lim[0] << ", " << lo_lim[1] << endl;
      out << "  hi_lim = " << hi_lim[0] << ", " << hi_lim[1];
      edm::LogInfo("execAsym2D") << out.str();
    }

    first_event = false;
  }

  // Enter into this loop if parameters have changed or this is first time
  // in loop.
  if (par[0]!=parsave[0] || par[1]!=parsave[1] || par[2]!=parsave[2]) {
    paramsChanged = true;
    for (int i = 0; i < 3; i++) parsave[i] = par[i];
    anorm = 0.;

    // If in 2D mode, normalization is done over all 2 variables using
    // nIntPoints in crude Monte-Carlo method.  It is important that the
    // same number of points be used each time this integral is performed.
    double relerr = 0.0;
    double epsilon = .10;
    unsigned int status = 999;

    // Can set this to true if you want to debug integration
    const bool debug_int = false;

    TFMultiD *f_2d = new TFMultiD("f_2d", asym2D, 2, 3);
    f_2d->SetParameters(par[0], par[1], par[2]);      
    anorm = f_2d->crudeMCIntegral(2, lo_lim, hi_lim, epsilon, relerr, 
				  status, debug_int, nNormPoints);

    // Dump if integration does not have enough points for specified accuracy
    if (status != 1)
      edm::LogWarning("execAsym2D")
	<< "INTEGRATION NOT TAKEN TO SPECIFIED ACCURACY\n"
	<< "  epsilon = " << epsilon << ", relerr = " << relerr 
	<< ", status = " << status;

    // Shows relative error of integration, and status indicated if this
    // error is smaller than what was set for epsilon.
    if (asymFitManager.debug()) 
      edm::LogVerbatim("execAsym2D")
	<< "anorm = " << anorm << ", relerr =  " << relerr 
	<< ", status = " << status << endl;

    delete f_2d;
  }

  // Dump values every time parameters change.
  if (paramsChanged && asymFitManager.debug())
    edm::LogVerbatim("execAsym2D")
      << "nchange " << nchange++ << ": func/anorm = "  << func
      << "/" << anorm << " = " << func/anorm << ", N = " << par[0] 
      << ", b = " << par[2] << ", A_fb = " << par[1]; 

  // Prevent division by zero
  if (anorm < asymFitManager.epsilon()) anorm = asymFitManager.epsilon();
  func /= anorm;
  
  // Return epsilon if function is too small (not necessary?)
  //if (func < asymFitManager.epsilon()) return asymFitManager.epsilon();
  return func;
}

//=============================================================================

double asymResSmear2D(double *x, double *par) {
  // This function is the same as asym2D described above, but also
  // has a resolution smear over all measured quantities.  For now, this
  // smear is a very simple gaussian smear.  At some point, can be improved
  // to have smear as 2D resolution plots.  'par' array contains 2 extra
  // variables which are used to set the X0 values in the Gaussian functions.
  // Example:  TMath::Gaus(par[i+3], x[i], SIGMA[i])
  const bool debug = false;
  double gaus[2];
  double gaus_tot = 1.;
  static double sigma[2] = {
    asymFitManager.rec_sigma(0), asymFitManager.rec_sigma(1)
  };

  // Could use either asym2D or normalized version (execAsym2D).
  // Since this function will be normalized again, not necessary to use
  // execAsmCSSmear (which takes longer).
  double asym = asym2D(x, par);
  //double asym = execAsym2D(x, par);

  // Do smear over all variables.
  for (int i = 0; i < 2; i++) {

    // Remember: par[i+3] is used to set X0 value in Gaussian. 
    gaus[i] = TMath::Gaus(par[i+3], x[i], sigma[i], true);
    gaus_tot *= gaus[i];
  }

  // Multiply gaussian-smear total and multiply by asym function.
  double result = asym*gaus_tot;

  // Some dumps.
  if (debug) {
    ostringstream out;
    out << "x0[] = ";
    for (int i = 0; i < 2; i++) out << par[i+3] << " ";
    out << endl << "  x[] = ";
    for (int i = 0; i < 2; i++) out << x[i] << " ";
    out << endl << "  dif[] = ";
    for (int i = 0; i < 2; i++) out << fabs(x[i] - par[i+3]) << " "; 
    out << endl << "  sigma[] = ";
    for (int i = 0; i < 2; i++) out << sigma[i] << " ";
    out << endl << "  n SIGMA[] = ";
    for (int i = 0; i < 2; i++) out << fabs(x[i] - par[i+3])/sigma[i] << " ";
    out << endl << "  gaus[] = ";
    for (int i = 0; i < 2; i++) out << gaus[i] << " ";
    out << endl << "  asym*gaus_tot = " << asym 
	<< "*" << gaus_tot << " = " << result;
    LogDebug("asymResSmear2D") << out.str();
  }

  return result;
}

//=============================================================================

double asymResSmearNorm2D(double *x, double *par) {
  // Routine which is used in conjunction with asymResSmear2D routine to make
  // normalized smeared PDF. 
  const bool debug = false;
  static double sigma[2] = {
    asymFitManager.rec_sigma(0), asymFitManager.rec_sigma(1)
  };

  // a_lim and b_lim are taken from PythiaSamples.h
  static double lo_limits[2] = {
    asymFitManager.a_lim(0), asymFitManager.a_lim(1)
  };
  static double hi_limits[2] = {
    asymFitManager.b_lim(0), asymFitManager.b_lim(1)
  };

  // Not necessary to use execAsym2D (see explanation in asym2D)
  double asym = asym2D(x, par);
  //double asym = execAsym2D(x, par);

  double lim1, lim2;
  double erf = 1.;

  // Flag to keep track of variables I need to normalize, because they
  // have non-cyclical physical boundaries (phi has a cyclical physical
  // boundary).

  // Loop over all dimensions
  for (int i = 0; i < 2; i++) {
    // This next calculation can be derived by starting with a general
    // smeared function (over limits from a_lim to b_lim), 
    // int(f(x)*Gaus(x-x0, SIGMA)*dx), and integrating
    // over all possible values of x0.  After a change a variables, you get 
    // .5*int(f(x)*erf(x-a_lim/(sqrt(2)*sigma))-
    //             erf(x-b_lim/(sqrt(2)*sigma))*dx).
    lim1 = ((x[i] - lo_limits[i])/(TMath::Sqrt(2.)*sigma[i]));
    lim2 = ((x[i] - hi_limits[i])/(TMath::Sqrt(2.)*sigma[i]));
    double tmp_erf = .5*(TMath::Erf(lim1) - TMath::Erf(lim2));
    erf *= tmp_erf;
    if (debug)
      LogDebug("asymResSmearNorm2D") << "  lim1 = " << lim1
				     << ", lim2 = " << lim2
				     << ", tmp_erf = " << tmp_erf << endl;
  }

  double result = asym*erf;
  if (debug) {
    ostringstream out;
    out << "  cos rap = " << x[0] << " " << x[1] << endl;
    out << "  lo_limits hi_limits " << lo_limits[0] << " " << lo_limits[1]
	<< " " << hi_limits[0] << " " << hi_limits[1] << endl;
    out << "  lim1 lim2 sigma1 sigma2 = " << lim1 << " " << lim2 
	<< " " << sigma[0] << " " << sigma[1] << endl;
    out << "  asym*erf = " << asym << "*" << erf << " = " << result;
    LogDebug("asymResSmearNorm2D") << out.str();
  }

  return result;
}

//=============================================================================

double slowAsymResSmearNorm2D(double *x, double *par) {
  // Slow routine (more accurate) which is used in conjunction with 
  // asymResSmear2D routine to make normalized smeared PDF. 
  const bool debug = false;
  static bool first_event = true;
  static int icount = 0;
  static double sigma[2] = {
    asymFitManager.rec_sigma(0), asymFitManager.rec_sigma(1)
  };
  static double nSigma = 0;
  static unsigned int nSmearPoints = 0;
  double x_low[2], x_upp[2];
  if (first_event) {
    nSmearPoints = (unsigned int) par[4];
    nSigma = (unsigned int) par[5];

    if (asymFitManager.debug()) {
      ostringstream out;
      out << "initial parameters: ";
      for (int i = 0; i < 5; i++) out << par[i] << " ";
      edm::LogInfo("slowAsymResSmearNorm2D") << out.str();
    }

    first_event = false;
  }

  // Set limits of area to smear under
  for (int i = 0; i < 2; i++) {
    x_low[i] = x[i] - nSigma*sigma[i]; 
    x_upp[i] = x[i] + nSigma*sigma[i]; 
  }
  // Limits of cos_cs must not exceed -1.  1.
  if (x_low[0] < -1.) { x_low[0] = -1.; }
  if (x_upp[0] >  1.) { x_upp[0] =  1.; }

  TFMultiD *f_sARS = new TFMultiD("f_sARS", asymResSmear2D, 2, 5);
  f_sARS->SetParameters(par[0], par[1], par[2]);
  for (int i = 0; i < 2; i++) { f_sARS->FixParameter(i+3, x[i]); }
  f_sARS->SetParNames("Norm", "A_fb", "b", "cos_0", "rap_0");

  // Some initialization of flags used in monte-carlo integration.
  double relerr = 0.0;
  double res_err = 0.0;
  unsigned int status = 999;
  const bool debug_int = false;
  double epsilon = .5;

  // Use Monte-Carlo to do smear.
  double res_smear = f_sARS->crudeMCIntegral(2, x_low, x_upp, epsilon, 
					       relerr, status, debug_int, 
					       nSmearPoints);
  res_err = res_smear*relerr;

  if (debug) {
    ostringstream out;
    out << "icount " << icount++ << endl;
    const string s_var[2] = {"cos", "rap"};
    for (int i = 0; i < 2; i++) {
      out << " " <<  s_var[i] <<" = " << x[i] << ", " << s_var[i] 
	  << "_sigma = " << sigma[i]
	  << ", " << s_var[i] << "_low = " << x_low[i] << ", " 
	  << s_var[i] << "_upp = " << x_upp[i] << endl;
    }
    out << " res_smear " << res_smear << "+/-" << res_err 
	<< ", relerr " << relerr << ", status " << status;
    LogDebug("slowAsymResSmearNorm2D") << out.str();
  }

  delete f_sARS;
  return res_smear;
}
  
//=============================================================================

double recAsym2D(double *x, double *par) {
  // Function which uses asymResSmear2D and asymResSmearNorm2D to make a 
  // normalized smeared version of asym2D (or execAsym2D).
  const bool debug = false;

  // Flag used for doing occasional dumps
  //int count_dump = 9999999;
  static int icount  = 0;
  static bool first_event = true;
  static double anorm, anorm_err, res_err, parsave[3];
  static double sigma[2] = {
    asymFitManager.rec_sigma(0), asymFitManager.rec_sigma(1)
  };
  static double lo_limits[2] = {
    asymFitManager.a_lim(0), asymFitManager.a_lim(1)
  };
  static double hi_limits[2] = {
    asymFitManager.b_lim(0), asymFitManager.b_lim(1)
  };
  static double nSigma = 0;
  static unsigned int nSmearPoints = 0;
  static unsigned int nNormPoints = 0;

  // Only need to recalculate normalization when parameters have changed
  bool paramsChanged = false;

  double epsilonSmear = .5;
  double epsilonNorm = .5;
  double x_low[2], x_upp[2];
  double res_smear = 0.;

  if (debug)
    LogDebug("recAsym2D") << "  dil rap, cos_cs = " << x[0] << " " << x[1];

  // Initialize saved parameters to zero on first call.
  if (first_event) { 
    nNormPoints = (unsigned int)par[3];
    nSmearPoints = (unsigned int)par[4];
    nSigma = (unsigned int)par[5];
    for (int i = 0; i < 3; i++) parsave[i] = 0;

    if (asymFitManager.debug()) {
      ostringstream out;
      out << "initial parameters: ";
      for (int i = 0; i < 6; i++) out << par[i] << " ";
      edm::LogInfo("recAsym2D") << out.str();
    }
  }

  // Save parameters from fit every time they change, so function can know
  // when they change.
  if (par[0]!=parsave[0] || par[1]!=parsave[1] || par[2]!=parsave[2]) {
    paramsChanged = true;
    anorm = 0.;
    for (int i = 0; i < 3; i++) parsave[i] = par[i];
  }

  // Set limits of area to smear under
  for (int i = 0; i < 2; i++) {
    x_low[i] = x[i] - nSigma*sigma[i];
    x_upp[i] = x[i] + nSigma*sigma[i];
  }
  // Limits of cos_cs must not exceed -1.  1.
  if (x_low[0] < -1.) { x_low[0] = -1.; }
  if (x_upp[0] >  1.) { x_upp[0] =  1.; }

  if (debug) {
    ostringstream out;
    for (int i = 0; i < 2; i++) {
      out << "x = " << x[i] << ", SIGMA = " << sigma[i] << endl;
      out << "x_low = " << x_low[i] << ", x_upp = " << x_upp[i] << endl;
    }
    LogDebug("recAsym2D") << out.str();
  }

  if (first_event) {
    if (asymFitManager.debug()) {
      ostringstream out;
      out << "Doing recAsym2D fit over Mass range "
	  << asymFitManager.fit_win(0) << "-"
	  << asymFitManager.fit_win(1) << endl
	  << "  Smearing over " << nSigma << " sigma...";
      if (nSmearPoints == 0)
	out << "  MC integration for smear taken to " << epsilonSmear
	    << " precision" << endl;
      else
	out << "  " << nSmearPoints
	    << " points used in MC integration for smear" << endl;

      if (nNormPoints == 0)
	out << "  MC integration for normalization taken to " << epsilonNorm 
	    << " precision" << endl;
      else
	out << "  " << nNormPoints
	    << " points used in MC integration for normalization " << endl;
      out << "  lo_limits = " << lo_limits[0] << " " << lo_limits[1] << endl;
      out << "  hi_limits = " << hi_limits[0] << " " << hi_limits[1] << endl;
      out << "  x_low     = " << x_low[0] << " " << x_low[1] << endl;
      out << "  x_upp     = " << x_upp[0] << " " << x_upp[1];
      edm::LogInfo("recAsym2D") << out.str();
    }

    first_event = false;
  }

  // Smear part of PDF.  Smear is over all generated values of variables
  // The reconstructed values are stored as extra parameters in asymResSmear2D
  // function (see comments in this function for more details).
  TFMultiD *f_aRS = new TFMultiD("f_aRS", asymResSmear2D, 2, 5);
  f_aRS->SetParameters(par[0], par[1], par[2], x[0], x[1]);
  for (int i = 0; i < 2; i++) { f_aRS->FixParameter(i+3, x[i]); }
  f_aRS->SetParNames("Norm", "A_fb", "b", "cos_0", "rap_0");

  // Some initialization of flags used in monte-carlo integration.
  double relerr = 0.0;
  unsigned int status = 999;
  const bool debug_int = false;
  // Use Monte-Carlo to do smear.
  res_smear = f_aRS->crudeMCIntegral(2, x_low, x_upp, epsilonSmear, relerr, 
  			     status, debug_int, nSmearPoints);
  if (debug)
    LogDebug("recAsym2D") << icount << " res_smear "
			  << res_smear << ", relerr " << relerr 
			  << ", status " << status << ", x[1] ";

  res_err = res_smear*relerr;
  
  // Dump if integration does not have enough points for specified accuracy
  if (status != 1 && debug) {
    ostringstream out;
    out << "res_smear integration not up to specified accuracy!" << endl;
    out << "  res_smear = " << res_smear << endl;
    out << "  epsilonSmear = " << epsilonSmear << ", relerr = " << relerr 
	<< ", status = " << status << endl;
    out << "  icount = " << icount << ", x[] = ";
    for (int i = 0; i < 2; i++) out << x[i] << " ";
    out << endl << "  x_low[] = ";
    for (int i = 0; i < 2; i++) out << x_low[i] << " ";
    out << endl << "  x_upp[] = ";
    for (int i = 0; i < 2; i++) out << x_upp[i] << " ";
    out << endl << "  SIGMA[] = ";
    for (int i = 0; i < 2; i++) out << sigma[i] << " ";
    edm::LogWarning("recAsym2D") << out.str();
  }

  // Must calculate normalization part every time parameters change
  if (paramsChanged) {
    TFMultiD *f_aRSN = new TFMultiD("f_aRSN", asymResSmearNorm2D, 2, 3);
    f_aRSN->SetParameters(par[0], par[1], par[2]);
    //TFMultiD *f_aRSN = new TFMultiD("f_aRSN", slowAsymResSmearNorm2D, 2, 5);
    //f_aRSN->SetParameters(par[0], par[1], par[2]);
    //f_aRSN->FixParameter(3, nSigma);
    //f_aRSN->FixParameter(4, nSmearPoints);
    //f_aRSN->SetParNames("Norm", "A_fb", "b", "nSigma", "nSmearPoints");

    // Different initialization for this integration
    relerr = 0.0;
    status = 999;
    anorm = f_aRSN->crudeMCIntegral(2, lo_limits, hi_limits, epsilonNorm, 
				    relerr, status, debug_int, nNormPoints);

    if (asymFitManager.debug())
      edm::LogVerbatim("recAsym2D") << "anorm = " << anorm
				    << ", relerr =  " << relerr 
				    << ", status = " << status;

    if (status != 1 && debug)
      edm::LogWarning("recAsym2D")
	<< "  NORMALIZATION INTEGRATION NOT TAKEN TO SPECIFIED ACCURACY" 
	<< endl << "  epsilonNorm = " << epsilonNorm
	<< ", relerr = " << relerr << ", status = " << status;

    anorm_err = anorm*relerr;
    delete f_aRSN;
  }

  if ((paramsChanged && asymFitManager.debug()) || debug)
    edm::LogVerbatim("recAsym2D")
     << "icount " << icount++
     << "\n  res_smear = " << res_smear << "+/-" << res_err
     << ", anorm = " << anorm << "+/-" << anorm_err  
     << "\n  res_smear/anorm = " << res_smear/anorm
     << "\n  Norm = " << par[0] << ", b = " << par[2] << ", A_fb = " << par[1];

  // Prevent division by zero
  if (anorm < asymFitManager.epsilon()) anorm = asymFitManager.epsilon();
  //if (res_smear/anorm < asymFitManager.epsilon()) return asymFitManager.epsilon();

  // divide by anorm to get final answer
  res_smear /= anorm;

  delete f_aRS; 
  return res_smear;
}  

//=============================================================================

double asym6D(double *x, double *par) {
  const bool debug = false;

  TLorentzVector v_dil, v_mum_star, v_mup_star;

  AsymFitData data;
  data.cos_cs = x[0];
  data.rapidity = x[1];
  data.pT = x[2];
  data.phi = x[3];
  data.mass = x[4];
  data.phi_cs = x[5];

  // Calculate 4-vectors of dilepton in lab frame, mu+, mu- in dilepton
  // rest frame
  calc4Vectors(data, v_dil, v_mum_star, v_mup_star, debug);
  data.cut_status = diRapAccept(v_dil, v_mum_star, v_mup_star);

  if (debug)
    LogDebug("asym6D") << data;

  // See if dimuon is within detector acceptance.
  if (data.cut_status != NOTCUT) return asymFitManager.epsilon();

  // Multiply function by PDF's of the following variables.  This is only
  // necessary for multid mode.
  double mass   = data.mass;
  double pt_sqr = data.pT*data.pT;
  double rap    = data.rapidity;
  double phi_cs = data.phi_cs;
  double cos_cs = data.cos_cs;

  double w = 0.;
  if (asymFitManager.correct_mistags()) { w = mistagProb(fabs(rap), cos_cs, mass); }
  //if (asymFitManager.correct_mistags()) { w = mistagProbVsRap(fabs(rap)); }

  double neg_cos_cs = -data.cos_cs;
  double func = (1.-w)*asym_3_PDF(&cos_cs, par) + 
    w*asym_3_PDF(&neg_cos_cs, par);

  // Use PDF of pT^2, but then convert to PDF for pT. 
  func *= 2.*data.pT*ptSqrDist(&pt_sqr, par);
  func *= massDist(&mass, par);
  //LogDebug("asym6D") << "mass " << mass 
  //                   << " massDist " << massDist(&mass, par);
    
  // Revised-Thermalized model for now.
  func *= yDistRTC(&rap, par); 
  
  // PDF of Phi Collins-Soper
  func *= phiCSDist(&phi_cs, par);

  // Some debug statements
  if (debug) {
    ostringstream out;
    out << "  pTPdf = "<< 2*data.pT*ptSqrDist(&pt_sqr, par)
	<< ", massPdf = " << massDist(&mass, par)
	<< ", YPdf = " << yDistRTC(&rap, par) 
	<< ", phiPdf = " << phiCSDist(&phi_cs, par) << endl;
    out << "x = ";
    for (int i = 0; i < 6; i++) out << x[i] << " ";
    out << "  N = " << par[0] << ", b = " << par[2] << ", A_fb = " 
	<< par[1] << ", f = " << func;
    LogDebug("asym6D") << out.str();
  }

  return func;
}

//=============================================================================

double execAsym6D(double *x, double *par) {
  // Same as asym6D, except that this function can be used by the 
  // unbinned fitter, since it uses a monte-carlo integration method to 
  // normalize the function.
  const bool debug = false;
  static bool first_event = true;
  static int nchange = 0;
  static unsigned int nNormPoints = 0;

  if (debug)
    LogDebug("execAsym6D") 
      << "  dil (rap, pt, phi, m) = (" << x[1] << ", " << x[2] << ", "
      << x[3] << ", " << x[4] << ")\n"
      << "  cs_cos(theta), cs_phi = " << x[0] << ", " << x[5]; 

  static double anorm, parsave[3];

  // Start with asym6D function.
  double func = asym6D(x, par);

  // This flag is used to keep track of when the parameters have changed.
  // We only need to recalculate the normalization when the parameters
  // have changed.
  bool paramsChanged = false;

  // Initialize parameters on first call
  if (first_event) {
    // Number of points used in normalization is taken from 4th parameter
    // sent to this function
    nNormPoints = (unsigned int) par[3];
    for (int i = 0; i < 3; i++) parsave[i] = 0;

    if (asymFitManager.debug()) {
      ostringstream out;
      out << "initial parameters" << endl;
      for (int i = 0; i < 4; i++) out << par[i] << " ";
      out << endl;
      out << "  a_lim = {" << asymFitManager.a_lim(0);
      for (int i = 0; i < 5; i++)
	out << ", " << asymFitManager.a_lim(i+1);
      out << "}" << endl << "  b_lim = {" << asymFitManager.b_lim(0);
      for (int i = 0; i < 5; i++)
	out << ", " << asymFitManager.b_lim(i+1);
      out << "}";
      edm::LogInfo("execAsym6D") << out.str();
    }

    first_event = false;
  }

  // Enter into this loop if parameters have changed or this is first time
  // in loop.
  if (par[0]!=parsave[0] || par[1]!=parsave[1] || par[2]!=parsave[2]) {
    paramsChanged = true;
    anorm = 0.;
    for (int i = 0; i < 3; i++) parsave[i] = par[i];

    // If in 6D mode, normalization is done over all 6 variables using
    // nIntPoints in crude Monte-Carlo method.  It is important that the
    // same number of points be used each time this integral is performed.
    double relerr = 0.0;
    double epsilon = .10;
    unsigned int status = 999;

    // Can set this to true if you want to debug integration
    const bool debug_int = false;

    // Limits of integration (a_lim, b_lim) are found in PythiaSamples.h
    // and AsymFunctions.h.  These will be different for every sample.
    TFMultiD *f_6d = new TFMultiD("f_6d", asym6D, 6, 3);
    f_6d->SetParameters(par[0], par[1], par[2]);      

    anorm = f_6d->crudeMCIntegral(6, asymFitManager.a_lim_arr(),
				  asymFitManager.b_lim_arr(), epsilon, relerr, 
				  status, debug_int, nNormPoints);

    // Dump if integration does not have enough points for specified accuracy
    if (status != 1)
      edm::LogWarning("execAsym6D")
	<< "INTEGRATION NOT TAKEN TO SPECIFIED ACCURACY\n"
	<< "  epsilon = " << epsilon << ", relerr = " << relerr 
	<< ", status = " << status;

    // Can also try Gaussian-Quadrature method (fit wouldn't converge for me)
    /* anorm = f_6d->Integral(6, asymFitManager.a_lim_arr(),
	                         asymFitManager.b_lim_arr(), 
                                 epsilon, relerr, status); */

    // Shows relative error of integration, and status indicated it this
    // error is smaller than what was set for epsilon.
    if (asymFitManager.debug()) 
      edm::LogVerbatim("execAsym6D") 
	<< "anorm = " << anorm << ", relerr =  " << relerr 
	<< ", status = " << status;

    delete f_6d;
  }


  // Dump values every time parameters change.
  if ((paramsChanged && asymFitManager.debug()) || debug)
    edm::LogVerbatim("execAsym6D")
      << "nchange " << nchange++ << ": func/anorm = " << func 
      << "/" << anorm << " = " << func/anorm << ", N = " << par[0] 
      << ", b = " << par[2] << ", A_fb = " << par[1];

  // Prevent division by zero
  if (anorm < asymFitManager.epsilon()) anorm = asymFitManager.epsilon();
  func /= anorm;
  
  // Return epsilon if function is too small (not necessary?)
  //if (func < asymFitManager.epsilon()) return asymFitManager.epsilon();
  return func;
}

//=============================================================================

double asymResSmear6D(double *x, double *par) {
  // This function is the same as asym6D described above, but also
  // has a resolution smear over all measured quantities.  For now, this
  // smear is a very simple gaussian smear.  At some point, can be improved
  // to have smear as 2D resolution plots.  'par' array contains 6 extra
  // variables which are used to set the X0 values in the Gaussian functions.
  // Example:  TMath::Gaus(par[i+3], x[i], SIGMA[i])
  const bool debug = false;
  double gaus[6];
  double gaus_tot = 1.;
  double result = 0.;

  // I use these variables to dump contents of routine every iteration of 
  // count_dump.  It's not so informative to see it dumped on every iteration.
  static int count = 0;
  static int count_used = 0;
  //int count_dump = 9999999;
  static double gaus_tot_max;
  if (count%1000000==0) {
    if (debug) LogDebug("asymResSmear6D") << "zeroing gaus_tot_max";
    gaus_tot_max=0.;
    count_used=0;
  }
  count++;

  // Do smear over all variables.
  for (int i = 0; i < 6; i++) {
    // Remember: par[i+3] is used to set X0 value in Gaussian.
    // rec_sigma is read in from the config file
    // (and will be different for every sample).
    //gaus[i] = TMath::Gaus(par[i+3], x[i], asymFitManager.rec_sigma(i), true);
    double arg = (par[i+3]-x[i])/asymFitManager.rec_sigma(i);
    double res = exp(-0.5*arg*arg);
    gaus[i] = res/(2.50662827463100024*asymFitManager.rec_sigma(i));
    gaus_tot *= gaus[i];
  }

  if (gaus_tot < gaus_tot_max*1.e-6) {return result;}
  count_used++;
  if (gaus_tot > gaus_tot_max) {
    gaus_tot_max=gaus_tot;
    if (debug)
      LogDebug("asymResSmear6D") << "new gaus_tot_max " << gaus_tot_max
				 << " count " << count
				 << " used " << count_used;
  }

  // Could use either asym6D or normalized version (execAsym6D).
  // Since this function will be normalized again, not necessary to use
  // execAsmCSSmear (which takes longer).
  double asym = asym6D(x, par);
  //double asym = execAsym6D(x, par);

  // Multiply gaussian-smear total and multiply by asym function.
  result = asym*gaus_tot;

  // Some dumps.
  if (debug) {
    ostringstream out;
    out << "x0[] = ";
    for (int i = 0; i < 6; i++) out << par[i+3] << " ";
    out << endl << "  x[] = ";
    for (int i = 0; i < 6; i++) out << x[i] << " ";
    out << endl << "  dif[] = ";
    for (int i = 0; i < 6; i++) out << fabs(x[i] - par[i+3]) << " ";
    out << endl << "  SIGMA[] = ";
    for (int i = 0; i < 6; i++) out << asymFitManager.rec_sigma(i) << " ";
    out << endl << "  n SIGMA[] = ";
    for (int i = 0; i < 6; i++)
      out << fabs(x[i] - par[i+3])/asymFitManager.rec_sigma(i) << " ";
    out << endl << "  gaus[] = ";
    for (int i = 0; i < 6; i++) out << gaus[i] << " ";
    out << endl << "  asym*gaus_tot = " << asym << "*" << gaus_tot << " = "
	<< result;
    LogDebug("asymResSmear6D") << out.str();
  }

  return result;
}

//=============================================================================

double asymResSmearNorm6D(double *x, double *par) {
  // Routine which is used in conjunction with asymResSmear6D routine to make
  // normalized smeared PDF.  This routine allows me to avoid a 12 dimensional
  // integral, which would naively be expected since I have to smear over 
  // 6 dimensions, and then integrate over all possible reconstructed values
  // for those 6 dimensions.  Because of a calculation trick, I can reduce
  // the problem to a separate 6D integral (6D is the upper limit; I really
  // only have to integrate over the variables which have non-cyclical physical
  // cut-offs, such as pT).
  const bool debug = false;

  // Not necessary to use execAsym6D (see explanation in asym6D)
  double asym = asym6D(x, par);
  //double asym = execAsym6D(x, par);

  double lim1, lim2;
  double erf = 1.;

  // Flag to keep track of variables I need to normalize, because they
  // have non-cyclical physical boundaries (phi has a cyclical physical
  // boundary).
  static bool do_norm_of_var[6] = { true, false, true,
				    false, false, false };
 
  // Loop over all dimensions
  for (int i = 0; i < 6; i++) {
    // Is this variable one that needs normalization done
    if (do_norm_of_var[i]) {
      // This next calculation can be derived by starting with a general
      // smeared function (over limits from a_lim to b_lim), 
      // int(f(x)*Gaus(x-x0, SIGMA)*dx), and integrating
      // over all possible values of x0.  After a change a variables, you get 
      // .5*int(f(x)*erf(x-a_lim/(sqrt(2)*SIGMA))-
      //             erf(x-b_lim/(sqrt(2)*SIGMA))*dx).
      // Again, a_lim and b_lim are taken from PythiaSamples.h
      lim1 = ((x[i] - asymFitManager.a_lim(i))/(TMath::Sqrt(2.)*asymFitManager.rec_sigma(i)));
      lim2 = ((x[i] - asymFitManager.b_lim(i))/(TMath::Sqrt(2.)*asymFitManager.rec_sigma(i)));
      double tmp_erf = .5*(TMath::Erf(lim1) - TMath::Erf(lim2));
      erf *= tmp_erf;
      if (debug)
	LogDebug("asymResSmearNorm6D") << "  lim1 = " << lim1 
				  << ", lim2 = " << lim2
				  << ", tmp_erf = " << tmp_erf;
    }
  }

  double result = asym*erf;

  if (debug) {
    ostringstream out;
    out << "  x[] = ";
    for (int i = 0; i < 6; i++) out << x[i] << " ";
    out << endl << "  asym*erf = " << asym << "*" << erf << " = " << result;
    LogDebug("asymResSmearNorm6D") << out.str();
  }

  return result;
}
  
//=============================================================================

double recAsym6D(double *x, double *par) {
  // Function which uses asymResSmear6D and asymResSmearNorm6D to make a 
  // normalized smeared version of asym6D (or execAsym6D).
  const bool debug = false;

  // Flag used for doing occasional dumps
  //int count_dump = 9999999;
  static int icount  = 0;
  static bool first_event = true;
  static double nSigma = 0;
  static unsigned int nSmearPoints = 0;
  static unsigned int nNormPoints = 0;

  // When function is smeared, some variables have physical boundaries which
  // must never be passed (example, pT can never be less than 0).  Keep track
  // of those with boundaries.
  //                        cos_cs  rap    pT    phi     M    phi_cs
  bool set_max_limits[6] = {true,  false, true, false, false, false};

  static double anorm, anorm_err, res_err, parsave[3];

  // Only need to recalculate normalization when parameters have changed
  bool paramsChanged = false;

  double epsilonSmear = .5;
  double epsilonNorm = .5;
  double x_low[6], x_upp[6];
  double res_smear = 0.;

  if (debug)
    LogDebug("recAsym6D")
      << "  dil (rap, pt, phi, m) = (" << x[1] << ", " << x[2]
      << ", " << x[3] << ", " << x[4] << ")\n"
      << "  cs_cos(theta), cs_phi = " << x[0] << ", " << x[5]; 

  // Initialize saved parameters to zero on first call.
  if (first_event) { 
    for (int i = 0; i < 3; i++) parsave[i] = 0;

    // These values are stored in parameters sent to function
    nNormPoints = (unsigned int) par[3];
    nSmearPoints = (unsigned int) par[4];
    nSigma = par[5];
  }

  // Save parameters from fit every time they change, so function can know
  // when they change.
  if (par[0]!=parsave[0] || par[1]!=parsave[1] || par[2]!=parsave[2]) {
    paramsChanged = true;
    anorm = 0.;
    for(int i = 0; i < 3; i++) parsave[i] = par[i];
  }

  for (int i = 0; i < 6; i++) {
    // Set limits of area to smear under
    x_low[i] = x[i] - nSigma*asymFitManager.rec_sigma(i); 
    x_upp[i] = x[i] + nSigma*asymFitManager.rec_sigma(i); 
    
    // Set boundaries for those variables which have physical, non-cyclical
    // boundaries.
    if (set_max_limits[i]) {
      if (x_low[i] < asymFitManager.limit_low(i))
	x_low[i] = asymFitManager.limit_low(i);
      if (x_upp[i] > asymFitManager.limit_upp(i))
      x_upp[i] = asymFitManager.limit_upp(i);
    }
    if (debug)
      LogDebug("recAsym6D")
	<< "x = " << x[i] << ", SIGMA = " << asymFitManager.rec_sigma(i)
	<< "\nx_low = " << x_low[i] << ", x_upp = " << x_upp[i];
  }

  if (first_event) {
    if (asymFitManager.debug()) {
      ostringstream out;
      out << "Doing recAsym6D fit over Mass range "
	  << asymFitManager.fit_win(0) << "-"
	  << asymFitManager.fit_win(1) << endl;
      out << "  Smearing over " << nSigma << " sigma..." << endl;
      if (nSmearPoints == 0)
	out << "  MC integration for smear taken to " << epsilonSmear 
	    << " precision" << endl;
      else
	out << "  " << nSmearPoints 
	    << " points used in MC integration for smear" << endl;
      if (nSmearPoints == 0)
	out << "  MC integration for normalization taken to " << epsilonNorm 
	    << " precision" << endl;
      else
	out << "  " << nNormPoints
	    << " points used in MC integration for normalization " << endl;
      out << "  a_lim = {" << asymFitManager.a_lim(0);
      for (int i = 0; i < 5; i++) out << ", " << asymFitManager.a_lim(i+1);
      out << "}" << endl << "  b_lim = {" << asymFitManager.b_lim(0);
      for (int i = 0; i < 5; i++) out << ", " << asymFitManager.b_lim(i+1);
      out << "}" << endl << "  x_low = {" << x_low[0];
      for (int i = 0; i < 5; i++) out << ", " << x_low[i+1];
      out << "}" << endl << "  x_upp = {" << x_upp[0];
      for (int i = 0; i < 5; i++) out << ", " << x_upp[i+1];
      out << "}";
      edm::LogInfo("recAsym6D") << out.str();
    }

    first_event = false;
  }

  // Smear part of PDF.  Smear is over all generated values of variables
  // The reconstructed values are stored as extra parameters in asymResSmear6D
  // function (see comments in this function for more details).
  TFMultiD *f_aRS = new TFMultiD("f_aRS", asymResSmear6D, 6, 9);
  f_aRS->SetParameters(par[0], par[1], par[2]);
  for (int i = 0; i < 6; i++) f_aRS->FixParameter(i+3, x[i]);
  f_aRS->SetParNames("Norm", "A_fb", "b", "cos_cs_0", "rap_0", "pT_0", 
		     "phi_0", "mass_0", "phi_cs_0");

  // Some initialization of flags used in monte-carlo integration.
  double relerr = 0.0;
  unsigned int status = 999;
  const bool debug_int = false;
  // Use Monte-Carlo to do smear.
  res_smear = f_aRS->crudeMCIntegral(6, x_low, x_upp, epsilonSmear, 
				     relerr, status, debug_int, nSmearPoints);
  if (debug)
    LogDebug("recAsym6D") << "res_smear " << res_smear
			  << ", relerr " << relerr
			  << ", status " << status;

  res_err = res_smear*relerr;
  
  // Dump if integration does not have enough points for specified accuracy
  if (status != 1 && debug) {
    ostringstream out;
    out << "recAsym6D:";
    out << "  res_smear integration not up to specified accuracy!" << endl;
    out << "  res_smear = " << res_smear << endl;
    out << "  epsilonSmear = " << epsilonSmear << ", relerr = " << relerr 
	<< ", status = " << status << endl;
    out << "  icount = " << icount << ", x[] = ";
    for (int i = 0; i < 6; i++) out << x[i] << " ";
    out << endl << "  x_low[] = ";
    for (int i = 0; i < 6; i++) out << x_low[i] << " ";
    out << endl << "  x_upp[] = ";
    for (int i = 0; i < 6; i++) out << x_upp[i] << " ";
    out << endl << "  SIGMA[] = ";
    for (int i = 0; i < 6; i++) out << asymFitManager.rec_sigma(i) << " ";
    edm::LogWarning("recAsym6D") << out.str();
  }

  // Must calculate normalization part every time parameters change
  if (paramsChanged) {
    TFMultiD *f_aRSN = new TFMultiD("f_aRSN", asymResSmearNorm6D, 6, 3);
    f_aRSN->SetParameters(par[0], par[1], par[2]);

    // Different initialization for this integration
    relerr = 0.0;
    status = 999;
    anorm = f_aRSN->crudeMCIntegral(6, asymFitManager.a_lim_arr(),
				    asymFitManager.b_lim_arr(), epsilonNorm, 
				    relerr, status, debug_int, nSmearPoints);

    if (asymFitManager.debug())
      edm::LogVerbatim("recAsym6D") << "anorm = " << anorm
				    << ", relerr =  " << relerr 
				    << ", status = " << status;

    if (status != 1)
      edm::LogWarning("recAsym6D")
	<< "  NORMALIZATION INTEGRATION NOT TAKEN TO SPECIFIED ACCURACY\n" 
	<< "  epsilonNorm = " << epsilonNorm << ", relerr = " << relerr 
	<< ", status = " << status;
    
    anorm_err = anorm*relerr;
    delete f_aRSN;
  }

  if ((paramsChanged && asymFitManager.debug()) || debug)
    edm::LogVerbatim("recAsym6D") 
      << "icount " << icount++
      << "\n  res_smear = " << res_smear << "+/-" << res_err
      << ", anorm = " << anorm << "+/-" << anorm_err  
      << "\n  res_smear/anorm = " << res_smear/anorm
      << "\n  Norm = " << par[0] << ", b = " << par[2] << ", A_fb = " << par[1]
      << "\n  nSigma = " << nSigma << ", nSmearPoints = " << nSmearPoints 
      << ", nNormPoints = " << nNormPoints;

  // Prevent division by zero
  if (anorm < asymFitManager.epsilon()) anorm = asymFitManager.epsilon();
  //if (res_smear/anorm < asymFitManager.epsilon()) return asymFitManager.epsilon();

  // divide by anorm to get final answer
  res_smear /= anorm;

  delete f_aRS; 
  return res_smear;
}  

//=============================================================================

double rapMaxAccept(double *x, double *par) {
  // Calculate maximum y of z-prime given maximum lab_pseudorap of muons
  // variable x[0] is cos theta* for which we want acceptance of muons
  // at plus and minus cos theta*.
  // theta* = angle between mu- and the quark in the Z' CMS frame
  // ("true" theta*).
  const bool debug = false;
  
  double cos = fabs(x[0]);
  double a = fabs(par[0]);
  
  // Use eta = -ln(tan(theta/2)) to derive 
  // cos(theta) = (exp(2*eta) - 1)/(exp(2*eta) + 1)
  
  //return 0 if cos(theta) too big even for y=0;
  if (cos > (exp(2*a) - 1.)/(exp(2*a) + 1.)) {
    //LogDebug("rapMaxAccept") << "returning 0 for a = " << a;
    return 0.;
  }
  
  // To derive, start with -ln(tan(theta/2)) = -0.5*ln(tan(theta/2)^2).
  // Then use following trig identities.
  // sin(theta/2)^2 = 1 - cos(theta)
  // cos(theta/2)^2 = 1 + cos(theta) to obtain
  // eta = .5*ln((1+cos(theta))/(1-cos(theta)))
  //
  // Since rapidity adds when Lorentz boosted (and rapidity and pseudorapidity
  // are very close) 
  // Y_mu_lab = Y_mu_cms + Y_zprime_lab  (_cms is Z' center of mass)
  // Y_zprime_lab_max = Y_mu_lab_max - Y_mu_cms
  double func = a - 0.5*log((1+cos)/(1-cos));
  if (debug)
    LogDebug("rapMaxAccept")
      << "a " << a << " cos " << cos << " func " << func;

  return func;
}

//=============================================================================

CUTSTATUS diRapAccept(TLorentzVector v_dil, TLorentzVector v_mum, 
		      TLorentzVector v_mup) {
  // Function which returns true or false depending on if an event is 
  // within the detector acceptance. It takes as input the 4-vectors of the
  // dilepton in the lab frame, and the mu+ and mu- in the dilepton rest frame.
  // The muon 4-vector for each muon is then boosted into the lab frame to 
  // see if the eta is less than the maximum lab eta.
  const bool debug = false;

  CUTSTATUS returnval = NOTCUT;
  if (debug)
    LogDebug("diRapAccept")
      << "(pT, eta, phi, mass) in Dilepton Rest Frame:"
      << "\n    mu- = (" << v_mum.Pt() << ", " << v_mum.Eta() << ", " 
      << v_mum.Phi() << ", " << v_mum.M() << ")"
      << "\n    mu+ = (" << v_mup.Pt() << ", " << v_mup.Eta() << ", " 
      << v_mup.Phi() << ", " << v_mup.M() << ")";

  // Boost muon 4-vectors into lab frame.
  v_mum.Boost(v_dil.BoostVector());
  v_mup.Boost(v_dil.BoostVector());

  // check if muons are within eta acceptance
  double eta_mum = v_mum.Eta();
  double eta_mup = v_mup.Eta();
  if (eta_mum < asymFitManager.mum_eta_lim_lo() || eta_mum > asymFitManager.mum_eta_lim_hi() ||
      eta_mup < asymFitManager.mup_eta_lim_lo() || eta_mup > asymFitManager.mup_eta_lim_hi())
    returnval = ETACUT;

  // Simple pT cut (currently disabled by setting to 0).
  if (v_mum.Pt() < asymFitManager.mum_pt_min() || v_mup.Pt() < asymFitManager.mup_pt_min())
    returnval = CUTSTATUS(returnval | PTCUT);

  if (debug) {
    ostringstream out;
    out << "(pT, eta, phi, mass) in Lab Frame:" << endl;
    out << "    mu- = (" << v_mum.Pt() << ", " << v_mum.Eta() << ", " << v_mum.Phi() << ", " << v_mum.M() << ")" << endl;
    out << "    mu+ = (" << v_mup.Pt() << ", " << v_mup.Eta() << ", " << v_mup.Phi() << ", " << v_mup.M() << ")" << endl;
    out << "  MUM_ETA_LIM = (" << asymFitManager.mum_eta_lim_lo() << ", " << asymFitManager.mum_eta_lim_hi() << ")" << endl;
    out << "  MUP_ETA_LIM = (" << asymFitManager.mup_eta_lim_lo() << ", " << asymFitManager.mup_eta_lim_hi() << ")" << endl;
    out << "  MUP_PT_LIM = " << asymFitManager.mup_pt_min() << endl;
    out << "  MUM_PT_LIM = " << asymFitManager.mum_pt_min() << endl;
    out << "  returnval = " << returnval;
    LogDebug("diRapAccept") << out.str();
  }

  return returnval;
}

//=============================================================================

void calc4Vectors(AsymFitData& x,  TLorentzVector& v_dil, 
		  TLorentzVector& v_mum_prime, 
		  TLorentzVector& v_mup_prime, bool debug) {
  // Routine which calculates the 4-vectors, of the dilepton in the lab frame,
  // and the mu- and mu+ in the dilepton rest frame.  It uses as input the
  // pt, eta, phi, mass of the dilepton, and the theta* and phi* of the 
  // mu- in the dilepton rest frame.
  //debug = false;
  TVector3 v3_mum_star;
  TLorentzVector v_mum_star;
  double pt_dil       = x.pT;
  double rap_dil      = x.rapidity;
  double phi_dil      = x.phi;
  double mass_dil     = x.mass;
  double cos_theta_cs = x.cos_cs;
  double phi_cs       = x.phi_cs;
  double p_mu         = mass_dil/2.;
  double eta_dil = calcEta(pt_dil, rap_dil, mass_dil, debug);

  // Set 4-vector for dilepton
  v_dil.SetPtEtaPhiM(pt_dil, eta_dil, phi_dil, mass_dil);

  // Boost beam and target into dilepton CM frame
  TLorentzVector v_pp1(0., 0., -asymFitManager.beam_energy(), asymFitManager.beam_energy());
  TLorentzVector v_pp2(0., 0.,  asymFitManager.beam_energy(), asymFitManager.beam_energy());
  TLorentzVector v_beam, v_targ;
  v_pp1.Boost(-v_dil.BoostVector());
  v_pp2.Boost(-v_dil.BoostVector());
  if (v_dil.Pz() > 0.) {
    v_targ = v_pp1;
    v_beam = v_pp2;
  }
  else {
    v_targ = v_pp2;
    v_beam = v_pp1;
  }    
  TVector3 v3_beam  = v_beam.Vect().Unit();
  TVector3 v3_targ  = v_targ.Vect().Unit();

  // Define z* axis to be line which bisects boosted beam and target lines
  TVector3 v3_zstar = (v3_beam-v3_targ).Unit();

  // Define y* axis to be perpendicular to plane formed by boosted beam 
  // and target lines
  TVector3 v3_ystar = (v3_targ.Cross(v3_beam)).Unit();

  // x* axis is now defined by y* and z* axis
  TVector3 v3_xstar = (v3_ystar.Cross(v3_zstar)).Unit();

  // get pt, theta, phi of mu- with respect to these axis
  double theta_cs = acos(cos_theta_cs);
  double pt_mu = fabs(p_mu*sin(theta_cs));
  if (phi_cs < 0.) phi_cs += 2*TMath::Pi();

  // Set 3-vector for mu- in star frame according to star axis
  v3_mum_star.SetPtThetaPhi(pt_mu, theta_cs, phi_cs);  

  // Rotate normal unboosted axis into new rotation
  TRotation a;
  a.RotateAxes(v3_xstar, v3_ystar, v3_zstar);

  // Transform boosted componenets of mu- into rotated axis
  TVector3 v3_mum_prime = v3_mum_star;
  v3_mum_prime.Transform(a);

  // Set 4-vectors of mu+, mu- (which are opposite to each other) in dilepton
  // CM frame.
  v_mum_prime.SetVectM(v3_mum_prime, asymFitManager.lepton_mass());
  v_mup_prime.SetVectM(-v3_mum_prime, asymFitManager.lepton_mass());

  if (debug) {
    ostringstream out;
    out << "  dil(pt,rap,phi,mass) = (" << pt_dil << ", " << rap_dil << ", "
	<< phi_dil << ", " << mass_dil << ")" << endl;
    out << "  cos_theta_cs = " << cos_theta_cs << ", theta_cs = " << theta_cs 
	<< ", phi_cs = " << phi_cs << ", p_mu = " << p_mu << endl;
    out << "  v3_beam = (" << v3_beam.X() << ", " << v3_beam.Y() << " "
	<< v3_beam.Z() << ")" << endl;
    out << "  v3_targ = (" << v3_targ.X() << " " 
	<< v3_targ.Y() << ", " << v3_targ.Z() << ")" << endl;
    out << "  v3_xstar = (" << v3_xstar.X() << ", " << v3_xstar.Y() << " "
	<< v3_xstar.Z() << ")" << endl;
    out << "  v3_ystar = (" << v3_ystar.X() << ", " << v3_ystar.Y() << " "
	<< v3_ystar.Z() << ")" << endl;
    out << "  v3_zstar = (" << v3_zstar.X() << ", " << v3_zstar.Y() << " "
	<< v3_zstar.Z() << ")" << endl;
    out << "  v3_mum_star = (" << v3_mum_star.X() << ", " << v3_mum_star.Y() 
	<< " " << v3_mum_star.Z() << ")" << endl;
    out << "  v3_mum_prime = (" << v3_mum_prime.X() << ", " 
	<< v3_mum_prime.Y() << " " << v3_mum_prime.Z() << ")" << endl;
    out << "  mu- (px*,py*,pz*,E*) = (" << v_mum_prime.Px() << ", "
	<< v_mum_prime.Py() << ", " << v_mum_prime.Pz() << ", " 
	<< v_mum_prime.E() << ")" << endl;
    out << "  mu+ (px*,py*,pz*,E*) = (" << v_mup_prime.Px() << ", "
	<< v_mup_prime.Py() << ", " << v_mup_prime.Pz() << ", " 
	<< v_mup_prime.E() << ")";
    LogDebug("calc4Vectors") << out.str();
  }
}

//=============================================================================

double calcCosThetaTrue(TLorentzVector v_quark, TLorentzVector v_mum, 
			TLorentzVector v_dil, bool debug) {
  v_quark.Boost(-v_dil.BoostVector());
  v_mum.Boost(-v_dil.BoostVector());
  double num = v_quark.Vect()*v_mum.Vect();
  double den = v_quark.Vect().Mag()*v_mum.Vect().Mag();
  double cos = num/den;
  if (debug) {
    ostringstream out;
    out << "mu-   = (" << v_mum.Px() << ", " << v_mum.Py() << ", " 
	 << v_mum.Pz() << ", " << v_mum.E() << ")" << endl;
    out << "quark = (" << v_quark.Px() << ", " << v_quark.Py() << ", " 
	 << v_quark.Pz() << ", " << v_quark.E() << ")" << endl;
    out << "dil   = (" << v_dil.Px() << ", " << v_dil.Py() << ", " 
	 << v_dil.Pz() << ", " << v_dil.E() << ")" << endl;
    out << "cos_true = " << num << "/" << den << " = " << cos;
    LogDebug("calcCosThetaTrue") << out.str();
  }
  return cos;
}

//=============================================================================

double calcCosThetaCSAnal(TLorentzVector v_dil, TLorentzVector v_mum, 
			  TLorentzVector v_mup, bool debug) {
  // Function to return value of cos(theta*) in Collins-Soper frame
  // takes as input 4-vector of dilepton in lab frame, and 4-vectors of mu+
  // and mu- in dilepton CM frame.
  //debug = false;

  // Get pz and E components of mu+ and mu- in lab frame.
  double pz_mum = v_mum.Pz();
  double e_mum  = v_mum.E();
  double pz_mup = v_mup.Pz();
  double e_mup  = v_mup.E();

  // Get mass and pt of dilepton in lab frame
  double pt_dil   = v_dil.Pt();
  double pl_dil   = v_dil.Pz();
  double mass_dil = v_dil.M();

  // Calculate cos_theta and return
  double cos_theta_cs = calcCosThetaCSAnal(pz_mum, e_mum, pz_mup, e_mup, 
					     pt_dil, pl_dil, mass_dil, debug);
  return cos_theta_cs;
}

//=============================================================================

double calcCosThetaCSAnal(double pz_mum, double e_mum, double pz_mup, 
			  double e_mup, double pt_dil, double pl_dil,
			  double mass_dil, bool debug) {
  // Analytical calculation of Collins-Soper cos(theta).  Uses pz, e of mu+
  // and mu-, and pt, pl, and mass of dilepton in lab frame.
  //debug = false;

  double mum_minus = (1./sqrt(2.))*(e_mum - pz_mum);
  double mum_plus  = (1./sqrt(2.))*(e_mum + pz_mum);
  double mup_minus = (1./sqrt(2.))*(e_mup - pz_mup);
  double mup_plus  = (1./sqrt(2.))*(e_mup + pz_mup);
  double dil_term  = 2./(mass_dil*sqrt((mass_dil*mass_dil) + 
					 (pt_dil*pt_dil)));
  double mu_term   = (mum_plus*mup_minus) - (mum_minus*mup_plus);
  double cos_cs    = dil_term*mu_term;

  // The above calculation assumed dilepton pL > 0. Flip the sign of
  // cos_cs if this isn't true.
  if (pl_dil < 0.)
    cos_cs *= -1.;

  if (debug) {
    ostringstream out;
    out << "  mu- (pz, E) = (" << pz_mum << ", " << e_mum << "), mu+ = ("
	<< pz_mup << ", " << e_mup << ")" << endl;
    out << "  pt_dil = " << pt_dil << ", pl_dil = " << pl_dil 
	<< ", mass_dil = " << mass_dil << endl;
    out << "  mum_minus = " << mum_minus << ", mum_plus = " 
	<< mum_plus << endl;
    out << "  mup_minus = " << mup_minus << ", mup_plus = " 
	<< mup_plus << endl;
    out << "  dil_term = " << dil_term << ", mu_term = " << mu_term 
	<< ", cos_cs = " << cos_cs;
    LogDebug("calcCosThetaCSAnal") << out.str();
  }
  return cos_cs;
}

//=============================================================================

double calcPhiCSAnal(double px_mum, double py_mum, double px_mup, 
		     double py_mup, double pt_dil, double eta_dil,
		     double phi_dil, double mass_dil, bool debug) {
  // Analytical calculation of collins-soper phi using 4-vector of dilepton 
  // and the px,py components of mu+, and mu- in the lab frame (can also use 
  // dilepton CM frame).
  //debug = true;

  TLorentzVector v_dil;
  TVector3 v3_beam, v3_R_T;

  // Create 4-vector out of Z' compononets.
  v_dil.SetPtEtaPhiM(pt_dil, eta_dil, phi_dil, mass_dil);

  // Approximate the longitudinal momentum of the quark to be in the same 
  // direction as that of the dilepton (here the beam is defined as the 
  // direction of the quark).
  if (v_dil.Pz() > 0.) v3_beam.SetXYZ(0., 0., asymFitManager.beam_energy());
  else  v3_beam.SetXYZ(0., 0., -asymFitManager.beam_energy());

  // Make a transverse unit vector in the direction of beam x dilepton
  v3_R_T = (v3_beam.Cross(v_dil.Vect())).Unit();

  // Store transverse components of vectors in appropriate containers.
  TVector2 v2_delta_T(px_mum-px_mup, py_mum-py_mup);
  TVector2 v2_Q_T(v_dil.X(), v_dil.Y());
  v2_Q_T = v2_Q_T.Unit();
  TVector2 v2_R_T(v3_R_T.X(), v3_R_T.Y());

  double Q_term = (sqrt((mass_dil*mass_dil + (pt_dil*pt_dil))))/mass_dil;
  double delta_R_term = v2_delta_T*v2_R_T;
  double delta_Q_term = v2_delta_T*v2_Q_T;
  double phi_cs = atan2(Q_term*delta_R_term, delta_Q_term);
  if (phi_cs < 0.) phi_cs += 2*TMath::Pi();

  if (debug) {
    ostringstream out;
    out << "  mu- (px,py) = (" << px_mum << ", " << py_mum << ")  mu+ = (" 
	<< px_mup << ", " << py_mup << ")" << endl;
    out << "  Z'(pT,eta,phi,m) = (" << pt_dil << ", " << eta_dil << ", "
	<< phi_dil << ", " << mass_dil << ")" << endl;
    out << "  v3_R_T = (" << v3_R_T.X() << ", " << v3_R_T.Y() << ", " 
	<< v3_R_T.Z() << ")" << endl;
    out << "  v2_delta_T = (" << v2_delta_T.X() << ", " << v2_delta_T.Y() 
	<< ")" << endl;
    out << "  v2_Q_T = (" << v2_Q_T.X() << ", " << v2_Q_T.Y() 
	<< "), v2_R_T = (" << v2_R_T.X() << ", " << v2_R_T.Y() << ")" << endl;
    out << "  Q_term = " << Q_term << ", delta_R_term = " << delta_R_term 
	<< ", delta_Q_term = " << delta_Q_term  << endl;
    out << "  phi_cs = " << phi_cs << endl;
    LogDebug("calcPhiCSAnal") << out.str();
  }

  return phi_cs;
}

//=============================================================================

void calcCSQuantities(TLorentzVector v_dil, TLorentzVector v_mum, 
		      double &cos_theta_cs, double &phi_cs, bool debug) {
  // Calculate Collins-Soper cos(theta) and phi using 4-vector of 
  // dilepton and mu- in lab frame. This is the calculation using the 
  // definition of Collins-Soper.
  //debug = false;

  // Boost mu- into dilepton CM frame.
  v_mum.Boost(-v_dil.BoostVector());
  TVector3 v3_mum = v_mum.Vect();

  // Boost beam and targe 4-vectors into dilepton CM frame
  TLorentzVector v_pp1(0., 0., -asymFitManager.beam_energy(), asymFitManager.beam_energy());
  TLorentzVector v_pp2(0., 0.,  asymFitManager.beam_energy(), asymFitManager.beam_energy());
  TLorentzVector v_beam, v_targ;
  v_pp1.Boost(-v_dil.BoostVector());
  v_pp2.Boost(-v_dil.BoostVector());

  // Tag beam so that it is in same direction as longitudinal momentum of
  // dilepton
  if (v_dil.Pz() > 0.) {
    v_targ = v_pp1;
    v_beam = v_pp2;
  }
  else {
    v_targ = v_pp2;
    v_beam = v_pp1;
  }  

  // Make boosted beam and target vectors unit vectors.
  TVector3 v3_beam = v_beam.Vect().Unit();
  TVector3 v3_targ = v_targ.Vect().Unit();

  // Find Collins-Soper z-axis and q_t axis (which I call x axis since phi* 
  // will be measured in plane transverse to z-axis and starting from the 
  // q_t axis).
  TVector3 v3_zstar = (v3_beam-v3_targ).Unit();
  TVector3 v3_ystar = (v3_targ.Cross(v3_beam)).Unit();
  TVector3 v3_xstar = (v3_ystar.Cross(v3_zstar)).Unit();

  // Calculate Collins-Soper cos-theta
  cos_theta_cs = v3_mum.Dot(v3_zstar)/(v3_mum.Mag());
  double sin_theta_cs = sqrt(1. - cos_theta_cs*cos_theta_cs);

  // Derive and return Collins-Soper phi
  double cos_phi_cs = v3_mum.Dot(v3_xstar)/(v3_mum.Mag()*sin_theta_cs);
  phi_cs = acos(cos_phi_cs);
  if (v3_ystar.Dot(v3_mum) < 0.) phi_cs *= -1.;
  if (phi_cs < 0.) phi_cs += 2*TMath::Pi();
  if (debug) {
    ostringstream out;
    out << "  Z' (pT,eta,phi,M) = (" << v_dil.Pt() << ", " << v_dil.Eta() 
	<< ", " << v_dil.Phi() << ", " << v_dil.M() << ")" << endl;
    out << "  mu- (px*,py*,pz*,E) = (" << v_mum.Px() << ", " << v_mum.Py()
	<< ", " << v_mum.Pz() << ", " << v_mum.E() << endl;
    out << "  v3_beam = (" << v3_beam.X() << ", " << v3_beam.Y() << ", "
	<< v3_beam.Z() << "), mag = " << v3_beam.Mag() << endl;
    out << "  v3_targ = (" << v3_targ.X() << ", " << v3_targ.Y() << ", "
	<< v3_targ.Z() << "), mag = " << v3_targ.Mag() << endl;
    out << "  v3_zstar = (" << v3_zstar.X() << ", " << v3_zstar.Y() << ", "
	<< v3_zstar.Z() << "), mag = " << v3_zstar.Mag() << endl;
    out << "  v3_xstar = (" << v3_xstar.X() << ", " << v3_xstar.Y() << ", "
	<< v3_xstar.Z() << "), mag = " << v3_xstar.Mag() << endl;
    out << "  v3_mum = (" << v3_mum.X() << ", " << v3_mum.Y() << ", "
	<< v3_mum.Z() << "), mag = " << v3_mum.Mag() << endl;
    out << "  cos_theta_cs = " << cos_theta_cs << ", sin_theta_cs = " 
	<< sin_theta_cs << endl;
    out << "  cos_phi_cs = " << cos_phi_cs << ", phi_cs = " << phi_cs;
    LogDebug("calcCSQuantities") << out.str();
  }
}

//=============================================================================

double calcEta(double pt, double rap, double mass, bool debug) {
  // Calculate eta given pt, rapidity, mass
  double e_trans = sqrt((pt*pt) + (mass*mass));
  double e_tot = e_trans*cosh(rap);
  double p = sqrt((e_tot*e_tot) - (mass*mass));
  double beta = p/e_tot;
  double eta = atanh((1./beta)*tanh(rap));

  if (debug) {
    ostringstream out;
    out << "calcEta:" << endl
	<< "  pt = " << pt << ", rap = " << rap << ", mass = " << mass << endl;
    out << "  e_trans = " << e_trans << ", e_tot = " << e_tot << ", p = " 
	<< p << ", beta = " << beta << endl;
    out << "  eta = " << eta;
    LogDebug("calcEta") << out.str();
  }

  return eta;
}

//=============================================================================

double mistagProb(double rap, double cos, double mass) {
  //Option to use either parameterizations or 2D histogram
  static const bool debug = false;
  double prob;

  if (asymFitManager.use_mistag_hist()) {
    if (fabs(rap) > 3.5) return 0.;
    static int nCosBins = asymFitManager.h2_mistagProb->GetNbinsX();
    static int nRapBins = asymFitManager.h2_mistagProb->GetNbinsY();
    int cos_bin = int(cos*nCosBins) + 1;
    int rap_bin = int(fabs(rap)/3.5*nRapBins) + 1;
    prob = asymFitManager.h2_mistagProb->GetBinContent(cos_bin, rap_bin);

    // Alternative way of calculating using separate cos and rap 1-D
    // histograms. JMTBAD should probably go back to using #ifdefs
    // instead of that if (asymFitManager.use_mistag_hist()) above,
    // and have one for MISTAG_HIST_2D vs _1D ...
    /*
    static int nCosBins = h_cos_true_mistag_prob->GetNbinsX();
    static int nRapBins = h_rap_mistag_fraction->GetNbinsX();
    int cos_bin = int(-cos*nCosBins) + 1;
    int rap_bin = int(fabs(rap)/7.*nRapBins) + 1;
    double cos_prob = h_cos_true_mistag_prob->GetBinContent(cos_bin);
    double rap_prob = h_rap_mistag_fraction->GetBinContent(rap_bin);
    prob = cos_prob + rap_prob - 2.*cos_prob*rap_prob;
    */

    if (debug) {
      LogDebug("mistagProb") 
	<< "cos rap " << cos << "(" << cos_bin << ") " << rap
	<< "(" << rap_bin << ") prob " << prob;
      //  << "\n   cos_prob = " << cos_prob << " rap_prob = " << rap_prob;
    }
  }
  else {
    double prob_v_rap = mistagProbVsRap(rap, mass);
    double prob_v_cos = mistagProbVsCos(cos);
    prob = prob_v_rap + prob_v_cos - 2.*prob_v_rap*prob_v_cos;
    if (debug)
      LogDebug("mistagProb")
	<< "cos rap prob " << cos << " " << rap << " " << prob;
  }

  return prob;
}

double mistagProbVsRap(double rap, double mass) {
  if (asymFitManager.calculate_mistag() && mass > 0) 
    return asymFitManager.mistag_calc()->omega(rap, mass);

  // Compute mistag probability as function of y
  double absy = fabs(rap);
  double prob;
  const bool debug = false;

  if (asymFitManager.use_mistag_hist()) {
    static int nRapBins = asymFitManager.h_rap_mistag_prob->GetNbinsX();
    int rap_bin = int(fabs(rap)/3.5*nRapBins) + 1;
    prob = asymFitManager.h_rap_mistag_prob->GetBinContent(rap_bin);
    if (debug)
      LogDebug("mistagProbVsRap") 
	<< "rap rap_bin prob " << rap << " "<< rap_bin << " " << prob;
    return prob;
  }
  
   //quadratic fit
  prob = 0.5 + asymFitManager.mistag_pars[0]*absy + asymFitManager.mistag_pars[1]*absy*absy;
  if (prob > .5) prob = .5;
  if (absy > asymFitManager.mistag_pars[2] || prob < 0.) prob = 0.;

  return prob;
}

double mistagProbVsCos(double cos) {
  double abs_cos = fabs(cos);
  double prob = asymFitManager.mistag_pars[3] + asymFitManager.mistag_pars[4]*exp(-asymFitManager.mistag_pars[5]*abs_cos);
  return prob;
}

double mistagProbVsPL(double pL) {
  //Calculate mistag probability as of function of longitudinal momentum
  const bool debug = false;
  if (fabs(pL) > 5000.) return 0.;
  static int nBins = asymFitManager.h_pL_mistag_prob->GetNbinsX();
  int pL_bin = int(pL/5000.*nBins) + 1;
  double prob = asymFitManager.h_pL_mistag_prob->GetBinContent(pL_bin);
  if (debug)
    LogDebug("mistagProbVsPL") << "pL " << pL << " pL_bin " << pL_bin
			       << " prob " << prob;
  return prob;
}

double mistagProbVsPL(double qpL, double dilpL) {
  // Calculate mistag probability as function of quark pL and dilepton pL
  // This is the most correct version of mistag
  const bool debug = false;
  if (fabs(qpL) > 500. || fabs(dilpL) > 5000.) return 0.;
  static int qplBins   = asymFitManager.h2_pL_mistag_prob->GetNbinsX();
  static int dilplBins = asymFitManager.h2_pL_mistag_prob->GetNbinsY();
  int qpL_bin = int(fabs(qpL)/500.*qplBins) + 1;
  int dilpL_bin = int(fabs(dilpL)/5000.*dilplBins) + 1;
  double prob = asymFitManager.h2_pL_mistag_prob->GetBinContent(qpL_bin, dilpL_bin);
  if (debug)
    LogDebug("mistagProbVsPL") 
      << "qpL dilpL, qpL_bin dilpL_bin prob " << qpL << " " << dilpL << " "
      << qpL_bin << " " << dilpL_bin << " " << prob;
  return prob;
}

double mistagProbVsPtRap(double pT, double rap) {
  // Calculate mistag probability as function of dilepton pT and rapidity
  const bool debug = false;
  if (fabs(pT*pT) > 700. || fabs(rap) > 3.5) return 0.;
  static int ptBins = asymFitManager.h2_pTrap_mistag_prob->GetNbinsX();
  static int rapBins = asymFitManager.h2_pTrap_mistag_prob->GetNbinsY();
  int pT_bin = int(fabs(pT*pT)/700.*ptBins) + 1;
  int rap_bin = int(fabs(rap)/3.5*rapBins) + 1;
  double prob = asymFitManager.h2_pTrap_mistag_prob->GetBinContent(pT_bin, rap_bin);
  if (debug)
    LogDebug("mistagProbVsPtRap")
      << "pT rap, pT_bin rap_bin prob " << pT << " " << rap << " "
      << pT_bin << " " << rap_bin << " " << prob;
  return prob;
}
  
//=============================================================================

double massDist(double *x, double *par) {
  double mass = x[0];

  if (mass < 0.)
    throw cms::Exception("massDist") << "mass = " << mass << " is negative\n";

  const int mass_type = asymFitManager.mass_type();

  if (mass_type < 0 || mass_type > 3)
    throw cms::Exception("massDist") << "mass_type = " << mass_type
				     << " is unknown\n";
  double f = -999.;
  if (mass_type == AsymFitManager::MASS_EXP)
    f = expBckg(&mass, asymFitManager.mass_pars);
  else if (mass_type == AsymFitManager::MASS_LOR)
    f = Lorentzian(&mass, asymFitManager.mass_pars);
  else if (mass_type == AsymFitManager::MASS_LOREXP)
    f = lorentzianPlusExpbckg(&mass, asymFitManager.mass_pars);

  // Don't normalize, otherwise you suffer from roundoff problems in 6D fit
  //f /= asymFitManager.mass_pars[nPars];
  //LogDebug("massDist") << "f/norm = " << f << "/" << asymFitManager.mass_pars[nPars];
  return f;
}

//=============================================================================

double ptSqrDist(double *x, double *par) {
  double pt_sq = x[0];
  if (pt_sq < 0.)
    throw cms::Exception("ptSqrDist") << "pt_sq = " << pt_sq
				      << " is negative\n";

  double f = (asymFitManager.pt_pars[0]*exp(-asymFitManager.pt_pars[1]*sqrt(pt_sq))) +
    (asymFitManager.pt_pars[2]*exp(-asymFitManager.pt_pars[3]*sqrt(pt_sq)));
  f /= asymFitManager.pt_pars[4];
  return f;
}

//=============================================================================

double phiCSDist(double *x, double *par) {
  static bool first_event = true;
  static bool ibOn = false;

  if (first_event) {
    first_event = false;

    // This has a non-flat distribution when internal-Brem is turned on.  
    // So switch off if parametrizations have not been set (or if fit did
    // not converge properly).
    if (asymFitManager.phi_cs_pars[0] > 0. && asymFitManager.phi_cs_pars[1] > 0. && 
	asymFitManager.phi_cs_pars[2] > 0. && asymFitManager.phi_cs_pars[3] > 0. && asymFitManager.phi_cs_pars[4] > 0.)
      ibOn = true;
  }

  if (!ibOn) return 1.;

  double phi_cs = x[0];
  // Clamp in [0, pi)
  if (phi_cs < 0.) phi_cs += 2*TMath::Pi();
  if (phi_cs > 2*TMath::Pi()) phi_cs -= 2*TMath::Pi();
  if (phi_cs > TMath::Pi()) phi_cs -= TMath::Pi();
  //if (phi_cs > TMath::Pi()) phi_cs = 2*TMath::Pi() - phi_cs;
  if (phi_cs < 0. || phi_cs > TMath::Pi())
    throw cms::Exception("phiCSDist") << "phi_cs = " << phi_cs
				      << "is negative\n";

  double f = 1.;
  int y = (int) asymFitManager.phi_cs_pars[3];
  for (int i = 0; i < y; i++)
    f *= (phi_cs - asymFitManager.phi_cs_pars[2]);
  f *= asymFitManager.phi_cs_pars[1];
  f += asymFitManager.phi_cs_pars[0];
  f /= asymFitManager.phi_cs_pars[4];

  return f;
}

//=============================================================================

double yDistRTC(double *x, double *par) {
  // Return y production distribution of Revised Thermalized cylinder model.
  double a = fabs(x[0]);
  double f = asymFitManager.rap_pars[0]*(tanh((asymFitManager.rap_pars[1]*a)+asymFitManager.rap_pars[2]+asymFitManager.rap_pars[3])-
			  tanh((asymFitManager.rap_pars[1]*a)-asymFitManager.rap_pars[2]+asymFitManager.rap_pars[3])+
			  tanh((asymFitManager.rap_pars[1]*a)+asymFitManager.rap_pars[2]-asymFitManager.rap_pars[3])-
			  tanh((asymFitManager.rap_pars[1]*a)-asymFitManager.rap_pars[2]-asymFitManager.rap_pars[3]));
  f /= asymFitManager.rap_pars[4];
  return f;
}


//=============================================================================

double cosTrueVsCS(double *x, double *par) {
  // This function can be used for doing convolution between true cos-theta
  // and Collins-Soper cos-theta (Although this has been determined to be
  // unnecessary for the purpose of our fits.
  const bool debug = false;
  static int nCSBins   = asymFitManager.h2_cos_cs_vs_true->GetNbinsX();
  static int nTrueBins = asymFitManager.h2_cos_cs_vs_true->GetNbinsY();
  bool paramsChanged = false;
  static bool firstCall = true;
  static double parsave[3];

  // Enter into this loop if parameters have changed or this is first time
  // in loop.
  if(par[0]!=parsave[0] || par[1]!=parsave[1] || par[2]!=parsave[2] 
     || firstCall) {
    firstCall = false;
    paramsChanged = true;
    for (int i = 0; i < 3; i++) { parsave[i] = par[i]; }
  }

  int csBin = int(((x[0]+1.)/2.)*nCSBins) + 1;
  double func = 0.;
  double norm = asymFitManager.h2_cos_cs_vs_true->Integral(csBin, csBin, 1, nTrueBins);

  for (int i = 1; i < nTrueBins+1; i++) {
    double tempProb = asymFitManager.h2_cos_cs_vs_true->GetBinContent(csBin, i);
    double tempTrue = asymFitManager.h2_cos_cs_vs_true->GetYaxis()->GetBinLowEdge(i);
    func += asym_3_PDF(&tempTrue, par)*tempProb;
    if (debug)
      LogTrace("AsymFunctions") << "tempProb " << tempProb
				<< " tempTrue " << tempTrue
				<< " func " << func;
  }

  if (debug || paramsChanged)
    LogTrace("AsymFunctions") << "cos_cs csBin func/norm = " << x[0]
			      << " " << csBin << " "  << func 
			      << "/" << norm << " = " << (func/norm);

  //if (norm < asymFitManager.epsilon()) return asymFitManager.epsilon();
  return func/norm;
}

double testSmearAsym(double *x, double *par) {
  // Function used for smearing 1D asym_3_PDF (used for testing smears).
  const bool debug = false;

  // Mean value and sigma used in gaussian
  double x0 = par[3];
  double sigma = par[4];

  // Calculate probability of x given x[0] in a gaussian, multiply by PDF
  // to smear
  double asym = asym_3_PDF(x, par);
  double gaus = TMath::Gaus(x[0], x0, sigma, true);
  double func = asym*gaus;

  if (debug) {
    ostringstream out;
    out << "testSmearAsym pars = ";
    for (int i = 0; i < 5; i++) out << par[i] << " ";
    out << endl << "  x x0 asym*gaus = " << x[0] << " " << x0 << " " << asym 
	<< "*" << gaus << " = " << func;
    LogTrace("AsymFunctions") << out.str();
  }

  return func;
}

double testSmearAsymInt(double *x, double *par) {
  // sigma, number of sigma and number of smeared points are determined from
  // parameters sent to this function.
  static bool first_event = true;
  const bool debug = false;
  static double par_save[3];
  static int icount = 0;

  // Number of sigma which will be used in normalization of testSmearAsym
  // Normally would normalize over entire range of -1 to 1, but reducing
  // the number of sigma can potentially increase the accuracy in the 
  // calculation of the area under the gaussian
  double nSigma = par[3];

  // This is the sigma used in testSmearAsym
  double sigma  = par[4];

  // The number of points to use when integrating over the gaussian integral
  // and the seed used as the starting point in the MC integration
  unsigned int nSmearPoints  = (unsigned int) par[5];
  int seed = (unsigned int) par[7];

  // Store first 3 parameters so that if they change, we can dump information.
  if (first_event) { 
    for (int i = 0; i < 3; i++) par_save[i] = 0; 
    first_event = false;
  }

  // Initialize smear function so that the area can be calculated with 
  // MC integration
  TFMultiD* f_testSmear = new TFMultiD("f_testSmear", testSmearAsym, 1, 5);
  f_testSmear->SetParameters(par[0], par[1], par[2], x[0], sigma);
  f_testSmear->FixParameter(0, par[0]);
  f_testSmear->FixParameter(3, x[0]);
  f_testSmear->FixParameter(4, sigma);
  f_testSmear->SetParNames("Norm", "A_fb", "b", "cos_0", "sigma");

  double epsilon = .5;
  double relerr = 0.;
  unsigned int status = 999;

  // Limits of the integration are determined by the number of sigma away
  // from the mean.  Otherwise comment out for conditions B and C.
  double low = x[0] - nSigma*sigma;
  double up  = x[0] + nSigma*sigma;
  //low = -1.; up = 1.;

  // For data generated with boundary condition A, uncomment next 2 lines.
  // (see testSmear routine)
  if (low < -1.) { low = -1.; }
  if (up  >  1.) { up  =  1.; }


  // Calculation of smeared function with MC integral
  double f = f_testSmear->crudeMCIntegral(1, &low, &up, epsilon, relerr,
					    status, false, nSmearPoints, seed);

  // Some dumps when one of the first 3 parameters has changed.
  if (par_save[0] != par[0] || par_save[1] != par[1] || par_save[2] != par[2]
      || debug) {
    for (int i = 0; i < 3; i++) par_save[i] = par[i];

    ostringstream out;
    out << icount << "testSmearAsymInt pars = ";
    for (int i = 0; i < 8; i++) out << par[i] << " ";
    out << endl;
    out << "  x = " << x[0] << " low up status f relerr " << low << " " << up 
	<< " " << status << " " << f << " " << relerr;
    LogTrace("AsymFunctions") << out.str();
  }
  
  delete f_testSmear;
  return f;
}

double testSmearAsymNorm(double *x, double *par) {
  // This function is used to quickly normalize testSmearAsymInt.
  // Normally it would be necessary to integrate testSmearAsymInt over
  // all possible values of x.  Since testSmearAsymInt already has a 1D
  // integral, normalizing it would require a 2D integral.  Thanks to
  // an algebra trick, it is only necessary to integrate this function
  // over all possible values of x and then divide testSmearAsymInt by the 
  // answer
  static double sigma  = par[3];

  double asym = asym_3_PDF(x, par);
  double lim1 = ((x[0] + 1.)/(TMath::Sqrt(2.)*sigma));
  double lim2 = ((x[0] - 1.)/(TMath::Sqrt(2.)*sigma));
  double erf = .5*(TMath::Erf(lim1) - TMath::Erf(lim2));
  double f = asym*erf;

  //LogTrace("AsymFunctions") << " sigma = " << sigma << ", lim1 = "
  //<< lim1 << ", lim2 = " << lim2 << ", asym*erf = " << asym << "*"
  //<< erf << " = " << f << endl;
  return f;
}

double recTestSlow(double *x, double *par) {
  // This is the naive algorithm to normalize testSmearAsymInt: integrate
  // testSmearAsymInt over all possible values of x.  Hence you must calculate
  // a double integral.  If the number of variables were large (such as 6
  // for execAsym6D), then to normalize would require a 12 deep integral.  This
  // would be a time-prohibitive calculation.
  static bool firstEvent = true;
  static double parsave[3];

  // Number of sigma and value of sigma to apply gaussian smear over
  double nSigma = par[3]; 
  double sigma = par[4];

  // Number of MC points used in calculation of gaussian area
  unsigned int nSmearPoints = (unsigned int) par[5];

  // Number of MC points used to normalize testSmearAsymInt
  unsigned int nNormPoints = (unsigned int) par[6];

  // Seed used for starting MC integration for smear and normalization
  int seed = (int) par[7];

  // Store first 3 parameters so we can do dumps if they change.
  if (firstEvent) {
    firstEvent = false;
    for (int i = 0; i < 3; i++) parsave[i] = 0.;

    ostringstream out;
    out << "recTestSlow initial parameters" << endl;
    for (int i = 0; i < 8; i++) out << par[i] << " ";
    LogTrace("AsymFunctions") << out.str();
  }

  static double den = 1.;;

  // Calculate smeared function.
  double num = testSmearAsymInt(x, par);

  // if free parameters have changed, recalculate integral of testSmearAsymInt
  if (parsave[0] != par[0] || parsave[1] != par[1] || parsave[2] != par[2]) {
    TFMultiD* f_testSmearNorm = new TFMultiD("f_testSmearNorm", 
					     testSmearAsymInt, 1, 8);
    f_testSmearNorm->SetParameters(par[0], par[1], par[2], nSigma, sigma,
				   nSmearPoints, 0., seed);  

    // sigma, number of sigma, number of smear and integration points are
    // determined from parameters sent from routine which calls this function.
    f_testSmearNorm->FixParameter(0, par[0]);
    f_testSmearNorm->FixParameter(3, nSigma);
    f_testSmearNorm->FixParameter(4, sigma);
    f_testSmearNorm->FixParameter(5, nSmearPoints);
    f_testSmearNorm->FixParameter(6, 0.);
    f_testSmearNorm->FixParameter(7, seed);
    f_testSmearNorm->SetParNames("Norm", "A_fb", "b", "nSigma", 
				 "sigma", "nSmearPoints", "not used", "seed");

    double epsilon = .5;
    double relerr = 0.;
    unsigned int status = 999;

    // Limits of integration are all possible values of x0
    double low = -1.;
    double up  = 1.;
    den = f_testSmearNorm->crudeMCIntegral(1, &low, &up, epsilon, relerr, 
					   status, false, nNormPoints);

    for (int i = 0; i < 3; i++) parsave[i] = par[i];

    // Some occasional dumps
    ostringstream out;
    out << "recTestSlow Parameters = " ;
    for (int i = 0; i < 8; i++) out << par[i] << " ";
    out << endl;
    out << "  x low up status relerr " << x[0] << " " << low << " " << up 
	<< " " << status << " " << relerr << endl;
    out << "  num/den = " << num << "/" << den << " = " << num/den;
    LogTrace("AsymFunctions") << out.str();

    delete f_testSmearNorm;
  }

  double f = num/den;
  return f; 
}

double recTestFast(double *x, double *par) {
  // As opposed to recTestSlow, we have created this function which
  // does not use a double integral to normalize testSmearAsymInt, but
  // rather uses testSmearAsymNorm, so that 2 separate integrals are calculated
  // instead of a double integral.  So for 6 variables, we will only need
  // to calculate 2 separate 6 deep integrals instead of one 12 deep integral.
  static bool firstEvent = true;
  static double parsave[3];

  // Number of sigma and value of sigma used in gaussian smear.
  double nSigma = par[3]; 
  double sigma = par[4];

  // Number of points used for calcluation of gaussian smear.
  //unsigned int nSmearPoints = (unsigned int) par[5];

  // Number of poitns used for normalization of gaussian smear.  Note: 
  // experience has shown that nSmearPoints and nNormPoints should be the
  // same number.
  unsigned int nNormPoints = (unsigned int) par[6];

  // Starting seed to be used for MC integrations.  Note: same seed should
  // be used for smear and for normalization.
  int seed = (int) par[7];

  // Store free parameters so we can dump information if they change
  if (firstEvent) {
    firstEvent = false;
    for (int i = 0; i < 3; i++) parsave[i] = 0.;

    ostringstream out;
    out << "recTestFast initial parameters" << endl;
    for (int i = 0; i < 8; i++) out << par[i] << " ";
    LogTrace("AsymFunctions") << out.str();
  }

  // Smeared function
  double num = testSmearAsymInt(x, par);

  static double den;
  // If free parameters change, must recalculate denominator used for 
  // normalization
  if (parsave[0] != par[0] || parsave[1] != par[1] || parsave[2] != par[2]) {

    // Function used for normalization
    TFMultiD* f_testSmearNorm = new TFMultiD("f_testSmearNorm", 
					     testSmearAsymNorm, 1, 4);
    f_testSmearNorm->SetParameters(par[0], par[1], par[2], sigma);

    // sigma determined from parameters sent from routine which calls 
    // this function.
    f_testSmearNorm->FixParameter(0, par[0]);
    f_testSmearNorm->FixParameter(3, sigma);
    f_testSmearNorm->SetParNames("Norm", "A_fb", "b", "sigma");

    double epsilon = .5;
    double relerr = 0.;
    unsigned int status = 999;

    // Limits of integration should be applied over all possible values of x
    double low = -1. - (nSigma*sigma);
    double up  = 1. + (nSigma*sigma);
    //low = -1.; up = 1.;

    // For data generated with boundary condition A, uncomment next 2 lines.
    // (see testSmear routine).  Otherwise comment out for conditions B and C.
    if (low < -1.) low = -1.;
    if (up  >  1.)  up =  1.;

    den = f_testSmearNorm->crudeMCIntegral(1, &low, &up, epsilon, relerr, 
					   status, false, nNormPoints, seed);

    for (int i = 0; i < 3; i++) parsave[i] = par[i];

    // Some dumps
    ostringstream out;
    out << "recTestFast pars = " ;
    for (int i = 0; i < 8; i++) out << par[i] << " ";
    out << endl;
    out << "  x low up status relerr " << x[0] << " " << low << " " << up 
	<< " " << status << " " << relerr << endl;
    out << "  num/den = " << num << "/" << den << " = " << num/den << endl;
    LogTrace("AsymFunctions") << out.str();

    delete f_testSmearNorm;
  }

  double f = num/den;
  return f; 
}

double recTest4BinFitFast(double* x, double* par) {
  // Multiply recTestFast routine by a constant number to be used in binned
  // fits (this number is the number of binned events multiplied by the 
  // width of the bin).
  static bool first_event = true;
  double f = recTestFast(x, par); 
  double norm = par[8];
  double func = f*norm;

  if (first_event) {
    first_event = false;
    
    ostringstream out;
    out << "recTest4BinFitFast pars:  ";
    for (int i = 0; i < 9; i++) out << par[i] << " ";
    out << endl << "f*norm = " << f << "*" << norm << " = " << func;
    LogTrace("AsymFunctions") << out.str();
  }

  return func;
}

double recTest4BinFitSlow(double* x, double* par) {
  // Multiply recTestSlow routine by a constant number to be used in binned
  // fits (this number is the number of binned events multiplied by the 
  // width of the bin).
  static bool first_event = true;
  double f = recTestSlow(x, par); 
  double norm = par[8];
  double func = f*norm;

  if (first_event) {
    first_event = false;

    ostringstream out;
    out << "recTest4BinFitSlow pars:  ";
    for (int i = 0; i < 9; i++) out << par[i] << " ";
    out << endl << "f*norm = " << f << "*" << norm << " = " << func;
    LogTrace("AsymFunctions") << out.str();
  }

  return func;
}

