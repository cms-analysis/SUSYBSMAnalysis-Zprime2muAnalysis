//Class to do unbinned maximum likelihood fit.
//Inspired by TTreePlayer::UnbinnedFit in Root's TTreePlayer.cxx
//But this version gets passed arrays of data rather then getting it from tree
//
// Bob Cousins 8/2002

#include <fstream>

#include "TF1.h"
#include "TFMultiD.h"
#include "TROOT.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "UnbinnedFitter.h"

using namespace std;

// for now, these must be global so that unbinnedFitLikelihoodFCN can
// get at them; with ROOT 5.16, said function can be a class method
// and this can be cleaned up
TF1 *fitfunc = 0;
TFMultiD *fitfuncMultiD = 0;
int  nfitentries = 0;
double *fitdata1 = 0;
double *fitdata2 = 0;
double *fitdata3 = 0;
double *fitweight = 0;
double sum_fitweight = 0.;
TVirtualFitter *tFitter = 0;
vector<double*> fitdata;
int fcnmode = 0;

UnbinnedFitter::UnbinnedFitter() {
  fitfunc = 0;
  fitfuncMultiD = 0;
  nfitentries = 0;
  fitdata1 = 0;
  fitdata2 = 0;
  fitdata3 = 0;
  fitweight = 0;
  sum_fitweight = 0.;
  tFitter = 0;
  fitdata.clear();
}

// even though we couldn't make it private and we could get at this
// by using an extern declaration, have an accessor for niceness
TVirtualFitter* UnbinnedFitter::getFitter() {
  return tFitter;
}

void unbinnedFitLikelihoodFCN(int &npar, double *gin, double &r,
			      double *par, int flag)
{
  // The Minuit FCN function used by the unbinned likelihood fit.
  // Inspired by TreeUnbinnedFitLikelihood
  // Arguments are Minuit standard, function is returned in r
  // This function uses TF1 or TFMultiD depending on fcnmode variable.

  static bool debug = false;
  static double logEpsilon = -230.; // protect against negative probabilities
  double logL = 0.0, prob;

  double x[MAXDIM] = {0.};

  //LogTrace("UnbinnedFitter") << "fitdata1 " << fitdata[0][0] 
  //                               << " fitdata2 " << fitdata[0][1]
  //				   << "\n nfitentries " << nfitentries;

  for (int i = 0; i < nfitentries; i++) {
    if (fcnmode == 0) {
      x[0] = fitdata1[i];
      if (fitdata2) x[1] = fitdata2[i];
      if (fitdata3) x[2] = fitdata3[i];
      prob = fitfunc->EvalPar(x,par); 
    }
    else if (fcnmode == 1) {
      for (unsigned int j = 0; j < fitdata.size(); j++) {
	x[j] = fitdata[j][i];
      }
      prob = fitfuncMultiD->EvalPar(x,par); 
      //LogTrace("UnbinnedFitter") << prob << " : ";
    }
    else
      throw cms::Exception("unbinnedFitLikelihoodFCN")
	<< "fcnmode = " << fcnmode << " not valid";
      
    //LogTrace("UnbinnedFitter") 
    // << "fitdata1[" << i << "] = " << fitdata1[i] << ", pars = [" 
    // << par[0] << ", " << par[1] << "], prob = " << prob;

    if (prob > 0.) {
      prob = TMath::Log(prob);
      if (fitweight) prob *= fitweight[i];
      logL += prob;
    }
    else {
      ostringstream out;
      out << "UnbinnedFitLikelihoodFCN: prob = " << prob
	  << " <= 0 ; i " << i << " x[0] " << x[0];

      if (fcnmode == 0) {
	if (fitdata2) out << " x[1] " << x[1];  
	if (fitdata3) out << " x[2] " << x[2];  
      }
      else if (fcnmode == 1)
	for (unsigned int j = 1; j < fitdata.size(); j++)
	  out << " x[" << j << "] " << x[j];

      if (fitweight) out << " w = "    << fitweight[i];
      out << endl;
      for (int ipar = 0; ipar < fitfunc->GetNpar(); ipar++)
	out << " par[" << ipar << "] = " << par[ipar];
      edm::LogWarning("UnbinnedFitter") << out.str();

      //int j; cin >> j;
      logL += logEpsilon;
    }
  }

  // Const; irrelevant?
  // JMT yes? http://root.cern.ch/root/roottalk/roottalk04/2604.html
  //if (fitweight) logL -= (TMath::Log(sum_fitweight))*nfitentries;

  // Extra term for extended maximum likelihood
#if defined(EML_1) && !defined(EML_2)
  if (fitweight) logL -= par[5]*nfitentries/sum_fitweight;
  else           logL -= par[5];
#endif
#if !defined(EML_1) && defined(EML_2)
  if (fitweight) logL -= (par[0] + par[5])*nfitentries/sum_fitweight;
  else           logL -= (par[0] + par[5]);
#endif

  r = -2.*logL;
  //LogTrace("UnbinnedFitter") << "r: " << setprecision(9) << r
  //                               << ", flag = " << flag;

  if (debug) {
    static int ncall=0;
    if (ncall++ <= 20) {
      ostringstream out;
      out << " ncall " << ncall << " nfitentries " << nfitentries;
      if (fcnmode == 0) {
	if (fitdata1) out << " data1 " << fitdata1[0];
	if (fitdata2) out << " data2 " << fitdata2[0];
	if (fitdata3) out << " data3 " << fitdata3[0];
      }
      else if (fcnmode == 1) {
	out << " ndim " << fitdata.size();
	for (unsigned int j = 0; j < fitdata.size(); j++)
	  out << " fitdata[" << j << "][0] " << fitdata[j][0];
      }
      if (fitweight) out << " weight " << fitweight[0];      
      LogTrace("UnbinnedFitter") << out.str();
    }
  }
}

//_____________________________________________________________________________

int UnbinnedFitter::unbinnedFitExec(const char *funcname, 
					const Option_t *option, 
					const int nentries, 
					vector<double*> data, 
					double *weight, double& log_ML)
{
//*-*-*-*-*-*Unbinned fit of one or more variable(s) from arrays passed
//  based on TTreePlayer::UnbinnedFit
//*-*        ===================================================
//
//  funcname is a TF1 function.
//
//   The list of fit options is given in parameter option.
//      option = "Q" Quiet mode (minimum printing)
//             = "V" Verbose mode (default is between Q and V)
//             = "E" Perform better Errors estimation using Minos technique
//             = "M" More. Improve fit results
//             = "C" Just call the function once and return its return
//   For the fit to be meaningful, the function must be self-normalized.
//   1, 2 and 3 Dimensional fits are supported if the pointers are non-zero.
//   weights are used if the pointers are non-zero.

  int i,nvpar,nparx;
  double par, we, al, bl;
  double eplus,eminus,eparab,globcc,fmin,edm,errdef,werr;
  double arglist[10];
  log_ML = 0;
  fcnmode = 1;

  // Set the global fit function so that UnbinnedFitLikelihoodFCN can find it.
  fitfuncMultiD = (TFMultiD*)gROOT->GetFunction(funcname);
  if (!fitfuncMultiD) {
    edm::LogWarning("UnbinnedFitter") << "unknown function " << funcname;
    return 0;
  }
  int npar = fitfuncMultiD->GetNpar();
  if (npar <= 0) {
    edm::LogWarning("UnbinnedFitter") << "illegal number of parameters "
					  << "npar = " << npar; 
    return 0; 
  }

  //Store pointers to data
  fitdata = data;
  fitweight = weight;

  TString opt = option;
  opt.ToLower();

  bool verbose = opt.Contains("v");
  ostringstream out;

  //calculate dimension of input vector x of function
  int ndim = data.size();
  if (verbose)
    out << "In unbinnedFitExec, # of parameters = " << npar << ", # dim = " 
	<< ndim << ", # var = " << fitfuncMultiD->GetNdim() << endl;

  nfitentries = nentries;
  //Compute total sum of weights to set the normalization factor
  sum_fitweight = 0;
  if (weight)
    for (i = 0; i < nentries; i++) sum_fitweight += weight[i];
  else {
    if (verbose)
      out << "Unit weights used" << endl;
  }
  if (verbose)
    out << "Total events in the fit region = " << nfitentries
	<< "; sum of the weights = " << sum_fitweight << endl;

  // Create and set up the fitter.  Argument is maximum number of params
  tFitter = TVirtualFitter::Fitter(0,5);
  tFitter->Clear();
  tFitter->SetFCN(unbinnedFitLikelihoodFCN);

  tFitter->SetObjectFit(fitfuncMultiD); //Not needed? 
  //tFitter->SetPrecision(.000000001);

  // Some initialisations
  if (!verbose) {
    arglist[0] = -1;
    tFitter->ExecuteCommand("SET PRINT", arglist,1);
    arglist[0] = 0;
    tFitter->ExecuteCommand("SET NOW",   arglist,0);
  }

  // Setup the parameters (#, name, start, step, min, max)
  double min, max;
  for (i = 0; i < npar; i++) {
    fitfuncMultiD->GetParLimits(i, min, max);
    if (verbose)
      out << "par[" << i << "] = " << fitfuncMultiD->GetParName(i) 
	  << ", limits = [" << min << ", " << max <<"]" << endl;
    double start = fitfuncMultiD->GetParameter(i);
    double step  = fabs(start)/100.;
    if (start == 0.) step = 0.001;   
    if (min < max) {
      tFitter->SetParameter(i, fitfuncMultiD->GetParName(i), start, step, min, max);
    }
    else {
      tFitter->SetParameter(i, fitfuncMultiD->GetParName(i), start, step, 0., 0.);
    }

    // Check for a fixed parameter
    if (max <= min && fabs(min) > 0.0) {
      tFitter->FixParameter(i);
    }
  }  // end for loop through parameters

  if (verbose)
    edm::LogInfo("UnbinnedFitter") << out.str();

  // Reset Print level
  if (verbose) {
    arglist[0] = 0;
    tFitter->ExecuteCommand("SET PRINT", arglist,1);
  }

  if (opt.Contains("c")) {
    arglist[0] = 1;
    tFitter->ExecuteCommand("CALL", arglist, 1);
  }
  else {
    // Now ready for minimization step
    arglist[0] = TVirtualFitter::GetMaxIterations();
    arglist[1] = 1;
    // arglist[1] = 0.1; 
    tFitter->ExecuteCommand("MIGRAD", arglist, 2);
    // tFitter->ExecuteCommand("MINIMIZE", arglist, 2);
    if (opt.Contains("m")) {
      tFitter->ExecuteCommand("IMPROVE",arglist,0);
    }
    if (opt.Contains("e")) {
      tFitter->ExecuteCommand("HESSE",arglist,0);
      tFitter->ExecuteCommand("MINOS",arglist,0);
    }
    fitfuncMultiD->SetNDF(fitfuncMultiD->GetNumberFitPoints()-npar);
  }
    
  // Get return status into function
  char parName[50];
  for (i=0;i<npar;i++) {
    tFitter->GetParameter(i,parName, par,we,al,bl);
    if (opt.Contains("e")) werr = we;
    else {
      tFitter->GetErrors(i,eplus,eminus,eparab,globcc);
      if (eplus > 0 && eminus < 0) werr = 0.5*(eplus-eminus);
      else                         werr = we;
    }
    fitfuncMultiD->SetParameter(i,par);
    fitfuncMultiD->SetParError(i,werr);
  }
  
  //cov_status = a status integer indicating how good is the covariance matrix:
  //                    0= not calculated at all
  //                    1= approximation only, not accurate
  //                    2= full matrix, but forced positive-definite
  //                    3= full accurate covariance matrix
  // fmin is the best function value ever found.
  int cov_status = 0;
  cov_status = tFitter->GetStats(fmin,edm,errdef,nvpar,nparx);
  
  // Print final values of parameters.
  if (!opt.Contains("q") && !verbose) {
    edm = 0;
    tFitter->PrintResults(1, edm);
  }
  
  // Logarithm of maximum likelihood
  log_ML = -fmin/2.;
  
  return cov_status;
}

//_____________________________________________________________________________

int UnbinnedFitter::unbinnedFitExec(const char *funcname, 
					const Option_t *option, 
					const int nentries, double *data1, 
					double *data2, double *data3, 
					double *weight, double& log_ML){
//*-*-*-*-*-*Unbinned fit of one or more variable(s) from arrays passed
//  based on TTreePlayer::UnbinnedFit
//*-*        ===================================================
//
//  funcname is a TF1 function.
//
//   The list of fit options is given in parameter option.
//      option = "Q" Quiet mode (minimum printing)
//             = "V" Verbose mode (default is between Q and V)
//             = "E" Perform better Errors estimation using Minos technique
//             = "M" More. Improve fit results
//
//   For the fit to be meaningful, the function must be self-normalized.
//   1, 2 and 3 Dimensional fits are supported if the pointers are non-zero.
//   weights are used if the pointers are non-zero.

  int i,nvpar,nparx;
  double par, we, al, bl;
  double eplus,eminus,eparab,globcc,fmin,edm,errdef,werr;
  double arglist[10];
  log_ML = 0;
  fcnmode = 0;

  // Set the global fit function so that UnbinnedFitLikelihoodFCN can find it.
  fitfunc = (TF1*)gROOT->GetFunction(funcname);
  if (!fitfunc) {
    edm::LogWarning("UnbinnedFitter") << "unknown function " << funcname;
    return 0;
  }
  int npar = fitfunc->GetNpar();
  if (npar <= 0) {
    edm::LogWarning("UnbinnedFitter") << "illegal number of parameters "
					  << "npar = " << npar; 
    return 0;
  }

  TString opt = option;
  opt.ToLower();

  bool verbose = opt.Contains("v");
  ostringstream out;

  //Store pointers to data
  fitdata1 = data1;
  fitdata2 = data2;
  fitdata3 = data3;
  fitweight = weight;

  int ndim = 0; //calculate dimension of input vector x of function
  if(data1) ndim++;
  if(data2) ndim++;
  if(data3) ndim++;
  if (verbose)
    out << "In UnbinnedFitExec, # of parameters = " << npar << ", # dim = " 
	<< ndim << ", # var = " << fitfunc->GetNdim() << endl;

  // Check that function has same dimension as number of variables
  if (ndim != fitfunc->GetNdim()) {
    edm::LogWarning("UnbinnedFitter")
      << "function dimension " << fitfunc->GetNdim() 
      << " not equal to expression dimension " << ndim;
     return 0;
  }

  nfitentries = nentries;
  //Compute total sum of weights to set the normalization factor
  sum_fitweight = 0;
  if (weight) {
    for (i = 0; i < nentries; i++) sum_fitweight += weight[i];
  }
  else {
    if (verbose) out << "Unit weights used" << endl;
  }
  if (verbose)
    out << "Total events in the fit region = " << nfitentries
	<< "; sum of the weights = " << sum_fitweight << endl;

  // Create and set up the fitter.  Argument is maximum number of params
  tFitter = TVirtualFitter::Fitter(0, npar);
  tFitter->Clear();
  tFitter->SetFCN(unbinnedFitLikelihoodFCN);
  tFitter->SetObjectFit(fitfunc); //Not needed? 

  // Some initialisations
  if (!verbose) {
    arglist[0] = -1;
    tFitter->ExecuteCommand("SET PRINT", arglist,1);
    arglist[0] = 0;
    tFitter->ExecuteCommand("SET NOW",   arglist,0);
  }

  // Setup the parameters (#, name, start, step, min, max)
  double min, max;
  for (i = 0; i < npar; i++) {
    fitfunc->GetParLimits(i, min, max);
    if (verbose)
      out << "par[" << i << "] = " << fitfunc->GetParName(i) 
	  << ", limits = [" << min << ", " << max <<"]" << endl;
    double start = fitfunc->GetParameter(i);
    double step  = fabs(start)/100.;
    if (start == 0.) step = 0.001;
    if (min < max) {
      tFitter->SetParameter(i, fitfunc->GetParName(i), start, step, min, max);
    }
    else {
      tFitter->SetParameter(i, fitfunc->GetParName(i), start, step, 0., 0.);
    }

    // Check for a fixed parameter
    if (max <= min && fabs(min) > 0.0) {
      tFitter->FixParameter(i);
    }
  }  // end for loop through parameters

  if (verbose)
    edm::LogInfo("UnbinnedFitter") << out.str();

  // Reset Print level
  if (verbose) {
    arglist[0] = 0;
    tFitter->ExecuteCommand("SET PRINT", arglist,1);
  }

  // Now ready for minimization step
  arglist[0] = TVirtualFitter::GetMaxIterations();
  arglist[1] = 1;
  // arglist[1] = 0.1;
  int errstat = tFitter->ExecuteCommand("MIGRAD", arglist, 2);
  // tFitter->ExecuteCommand("MINIMIZE", arglist, 2);
  // Switch to SIMPLEX if MIGRAD fails to converge
  if (errstat != 0)
    tFitter->ExecuteCommand("MINIMIZE", arglist, 2);
  if (opt.Contains("m")) {
    tFitter->ExecuteCommand("IMPROVE",arglist,0);
  }
  if (opt.Contains("e")) {
    tFitter->ExecuteCommand("HESSE",arglist,0);
    tFitter->ExecuteCommand("MINOS",arglist,0);
  }
  fitfunc->SetNDF(fitfunc->GetNumberFitPoints()-npar);

  // Get return status into function
  char parName[50];
  for (i=0; i<npar; i++) {
    tFitter->GetParameter(i,parName, par,we,al,bl);
    if (opt.Contains("e")) werr = we;
    else {
      tFitter->GetErrors(i,eplus,eminus,eparab,globcc);
      if (eplus > 0 && eminus < 0) werr = 0.5*(eplus-eminus);
      else                         werr = we;
    }
    fitfunc->SetParameter(i,par);
    fitfunc->SetParError(i,werr);
  }

  //cov_status = a status integer indicating how good is the covariance matrix:
  //                    0= not calculated at all
  //                    1= approximation only, not accurate
  //                    2= full matrix, but forced positive-definite
  //                    3= full accurate covariance matrix
  // fmin is the best function value ever found.
  int cov_status = 0;
  cov_status = tFitter->GetStats(fmin,edm,errdef,nvpar,nparx);

  // Print final values of parameters.
  if (!opt.Contains("q") && !verbose) {
    edm = 0;
    tFitter->PrintResults(1, edm);
  }

  // Logarithm of maximum likelihood
  log_ML = -fmin/2.;

  return cov_status;
}
