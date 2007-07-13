#ifndef FUNCTIONS_h
#define FUNCTIONS_h

#include "TLorentzVector.h"

double Lorentzian(double *x, double *par);
double Relativistic_Lorentzian(double *x, double *par);
double expBckg(double *x, double *par);
double expBckgNorm(double *x, double *par);
double lorentzianPlusBckg(double *x, double *par);
double lorentzianPlusBckgNorm(double *x, double *par);
double lorentzianPlusExpbckg(double *x, double *par);
double lorentzianPlusExpbckgNorm(double *x, double *par);
double lorentzianPlusExpbckgN(double *x, double *par);

double lorengaufun(double *x, double *par);
double lorengauPlusBckg(double *x, double *par);
double lorengauPlusExpbckg(double *x, double *par);
double lorengauPlusExpbckgNorm(double *x, double *par);
double lorengauPlusExpbckgN(double *x, double *par);
double mass_resolution(double *x, double *par);

TLorentzVector LorentzBoost(TLorentzVector boost_frame,
			    TLorentzVector rest_frame);
double Momentum(double pT, double eta);
double probks(const double alam);

#endif
