#ifndef ASYM_FUNCTIONS_h
#define ASYM_FUNCTIONS_h

#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

#include "AsymFitData.h"
#include "AsymFitManager.h"

extern bool mistagInData;
extern bool doingGravFit;
extern const bool useMistagHist;
extern bool asymDebug;
extern double mistag_pars[6];
extern double mass_pars[7];
extern double rap_pars[5];
extern double pt_pars[5];
extern double phi_cs_pars[5];
extern int nMistagBins;
extern TH2F* h2_mistagProb;
extern TH1F* h_rap_mistag_prob;
extern TH1F* h_cos_true_mistag_prob;
extern TH1F* h_pL_mistag_prob;
extern TH2F* h2_pL_mistag_prob;
extern TH2F* h2_cos_cs_vs_true;
extern TH2F* h2_pTrap_mistag_prob;

extern AsymFitManager asymFitManager;

// small number to return instead of zero
const double EPSILON = 1.e-6;

// roughly model the acceptance with a cut in eta at 2.4
const double ETALIM  = 2.4;
//const double ETALIM  = 9999.;
// able to separately cut on negative and positive muons
const double MUP_ETA_LIM[2] = {-ETALIM, ETALIM};
const double MUM_ETA_LIM[2] = {-ETALIM, ETALIM};

// parameters for a simple pT cut for graviton studies
const double PTMIN = 20; // GeV/c 
const double MUP_PT_MIN = PTMIN;
const double MUM_PT_MIN = PTMIN;

const double MUMASS = 0.10566; // GeV/c^2

enum MASS_TYPE { MASS_EXP=1, MASS_LOR, MASS_LOREXP };

double asym_2_PDF(double *x, double *par);
// asym_3_PDF is now a function pointer, initialized in AsymFunctions.C
// we should give it a more generic name but that could introduce bugs
extern double (*asym_3_PDF)(double *x, double *par);
// renamed asym_3_PDF to asym_3_PDF_real
double asym_3_PDF_real(double *x, double *par);
double asym_mistag_PDF(double *x, double *pars);
double GravitonCos_2_PDF(double *x, double *par);
double GravitonCos_th_PDF(double *x, double *par);
double GravitonCos_3_PDF(double *x, double *par);

double asym2D(double *x, double *par);
double execAsym2D(double *x, double *par);
double asymResSmear2D(double *x, double *par);
double asymResSmearNorm2D(double *x, double *par);
double slowAsymResSmearNorm2D(double *x, double *par);
double recAsym2D(double *x, double *par);
double asym6D(double *x, double *par);
double execAsym6D(double *x, double *par);
double asymResSmear6D(double *x, double *par);
double asymResSmearNorm6D(double *x, double *par);
double recAsym6D(double *x, double *par);
CUTSTATUS diRapAccept(TLorentzVector v_dil, TLorentzVector v_mum, 
		      TLorentzVector v_mup);
double rapMaxAccept(double *x, double *par);
void calc4Vectors(AsymFitData& x,  TLorentzVector& v_dil, 
		  TLorentzVector& v_mum_prime, 
		  TLorentzVector& v_mup_prime, bool debug);
double calcEta(double pt, double rap, double mass, bool debug);
double calcCosThetaTrue(TLorentzVector v_quark, TLorentzVector v_mum, 
			  TLorentzVector v_dil, bool debug);
double calcCosThetaCSAnal(TLorentzVector v_dil, TLorentzVector v_mum, 
			    TLorentzVector v_mup, bool debug);
double calcCosThetaCSAnal(double pz_mum, double e_mum, double pz_mup, 
			    double e_mup, double pt_dil, double pl_dil,
			    double mass_dil, bool debug);
double calcPhiCSAnal(double px_mum, double py_mum, double px_mup, 
		       double py_mup, double pt_dil, double eta_dil,
		       double phi_dil, double mass_dil, bool debug);
void calcCSQuantities(TLorentzVector v_dil, TLorentzVector v_mum, 
		      double &cos_theta_cs, double &phi_cs, bool debug);
double mistagProb(double rap, double cos);
double mistagProbVsRap(double rap);
double mistagProbVsPtRap(double pT, double rap);
double mistagProbVsCos(double cos);
double mistagProbVsPL(double pL);
double mistagProbVsPL(double qpL, double dilqL);
double massDist(double *x, double *par);
double ptSqrDist(double *x, double *par);
double phiCSDist(double *x, double *par);
double yDistRTC(double *x, double *par);
double cosTrueVsCS(double *x, double *par);
double testSmearAsym(double *x, double *par);
double testSmearAsymInt(double *x, double *par);
double testSmearAsymNorm(double *x, double *par);
double recTestFast(double *x, double *par);
double recTestSlow(double *x, double *par);
double recTest4BinFitFast(double* x, double* par);
double recTest4BinFitSlow(double* x, double* par);

#endif
