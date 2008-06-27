#ifndef MyFitters_H
#define MyFitters_H
/*
 * =====================================================================================
 *
 *       Filename:  MyFitFunctions.C
 *
 *    Description: 	Fitter for Iso03sumPt and Dxy of muons from  Zmumu(sig) and Dijets(bkg)
 *
 *        Version:  1.0
 *        Created:  05/16/2008 02:44:02 PM CEST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  IHEP, Beijing
 *
 * =====================================================================================
 */
#include "TMath.h"

// - - - single gaus model
double fitSigDxy(double *x, double par[]) // --- for prompt muons  : single gaus
{
	double dxy = x[0], c=par[0], m=par[1], s=par[2]; //main gaussian 
	double g=0;
	if(s!=0)g = c*TMath::Exp(-0.5*((dxy-m)/s)*((dxy-m)/s));
	return g;
}

// - - - triple gauses model
double fitBkgDxy(double *x, double par[]) // - --  for jet muons
{
	double dxy = x[0], c=par[0], m=par[1], s=par[2]; //main gaussian 
	double c2=par[3], m2=par[4], s2=par[5]; //tail gaussian
	double frac = par[6];
	double c3=par[7], m3=par[8], s3=par[9]; //tail constant 
	double frac2 = par[10];
	double g=0;
	double g2=0;
	double g3=0;
	if(s!=0)g = c*exp(-0.5*((dxy-m)/s)*((dxy-m)/s));
	if(s2!=0)g2 = c2*exp(-0.5*((dxy-m2)/s2)*((dxy-m2)/s2));
	if(s3!=0)g3 = c3*exp(-0.5*((dxy-m3)/s3)*((dxy-m3)/s3));
	//	return g+c2*g2+c3;
	return g+c2*g2+c3*g3;
}

// - - - 3 regions
double fitSigSumPt(double *val, double par[])
{

	Bool_t norm = 0;
	Double_t x = val[0];
	Double_t c0 = par[0]; //when sumpt =0 , return it
	Double_t con = par[1];
	Double_t mpv = par[2];
	Double_t sigma = par[3];

	Double_t c1 = par[4]; // the constants between (0.1, 0.5)

	if(x<0.1) return c0;
	else if(x<0.5) return c1;
	else return con*TMath::Landau(x, mpv, sigma, norm);
}

// - - - 3 regions
double fitBkgSumPt(double *x, double par[])
{
	double result = 0;
	double pt = x[0];
	double c0 = par[0]; //when sumpt =0 , return it
	double c1 = par[1]; // when sumpt > 0.6 ,  pol0	
	double c2 = par[2]; // when sumpt  in (0.1, 0.6)
	if(pt<0.1) result = c0;
	else if(pt<0.6) result =c2;
	else result = c1;
	return result;	
}
#endif
