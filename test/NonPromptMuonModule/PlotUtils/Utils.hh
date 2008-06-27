#ifndef _UTILS_HH
#define _UTILS_HH (1)

#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

//------------------------------------------------------------------------------
// Colors
//------------------------------------------------------------------------------
const    Color_t signalColor =   2; // Red
const    Color_t otherColor  =   7; // Light green
const    Color_t bg1Color    =  38; // Light blue
const    Color_t bg2Color    =  11; // Light brown
const    Color_t bg3Color    =   5; // Yellow
const    Color_t bg4Color    =   6; // Magenta

const    Color_t pionColor   = 188; // Very light blue 
const    Color_t kaonColor   =  50; // My red
const    Color_t protColor   =   5; // Yellow

//------------------------------------------------------------------------------
// Methods
//------------------------------------------------------------------------------
void     SetAxis(TAxis *&xaxis, TAxis *&yaxis, char *xtitle, char *ytitle);

void     axisTF1(TF1  *f, char *xtitle, char *ytitle);

void     axis1F (TH1F *h, char *xtitle, char *ytitle);

void     axis1D (TH1D *h, char *xtitle, char *ytitle);

void     axis2F (TH2F *h, char *xtitle, char *ytitle);

void     Draw(TH1F *histo,        int    fColor = 0,   int fStyle = 1,
              char *drawOpt = "", int    lColor = 1,   int lWidth = 1,
	      int   mColor  = 1,  double mSize  = 0.5, int mStyle = 20);

void     Draw1D(TH1D    *histo,        Color_t fColor = 0,   int fStyle = 1,
		char    *drawOpt = "", Color_t lColor = 1,   int lWidth = 1,
		Color_t  mColor  = 1,  double  mSize  = 0.5, int mStyle = 20);

void     Draw(TH1F *hist, TF1* func, float lo=-999, float hi=-999);

TLegend* SetLegend(float x0, float y0, float x1, float y1);

void     SetPave(char *title, float x, float y);

void     cdfPave();

void     SetTF1(TF1 *f,
		int fStyle,
		int fColor,
		int lStyle,
		int lColor,
		int lWidth,
		int nPoints = 1000);

double   chi2RebinTest(TH1F   *hist1,
                       TH1F   *hist2,
                       double &chi2,
		       int    &NDF,
                       double  low       = 0.0,
		       double  high      = 0.0,
		       int     verbose   =   0,
		       int     threshold =  20,
		       int     fitParams =   0);

double   chi2RebinTest(TH1F   *hist,
                       TF1    *fun,
		       double &chi2,
		       int    &NDF,
                       double  low       = 0.0,
		       double  high      = 0.0,
		       int     verbose   =   0,
		       int     threshold =  20);

void     PrintFitProbability(TCanvas *canvas,
                             TH1F    *hist1,
			     TH1F    *hist2,
                             double   low       = 0.0,
			     double   high      = 0.0,
			     int      verbose   =   0,
			     int      threshold =  20,
			     int      padNum    =   0,
			     int      fitParams =   0);


void     PrintFitProbability(TCanvas *canvas,
                             TH1F    *hist,
			     TF1     *fun,
                             double   low       = 0.0,
			     double   high      = 0.0,
			     int      verbose   =   0,
			     int      threshold =  20,
			     int      padNum    =   0);

void     PrintFitProbability(TPad   *pad,
                             TH1F   *hist1,
			     TH1F   *hist2,
                             double  low       = 0.0,
			     double  high      = 0.0,
			     int     verbose   =   0,
			     int     threshold =  20,
			     int     fitParams =   0);

void     PrintFitProbability(TPad   *pad,
                             TH1F   *hist,
			     TF1    *fun,
                             double  low       = 0.0,
			     double  high      = 0.0,
			     int     verbose   =   0,
			     int     threshold =  20);

void     PrintIt(TPad *pad, char *title);

void     PRLabel(TPad   *pad,
		 char   *title,
		 double  delta,
		 int     tAlign,
		 double  tSize = 0.05);

void     ComparisonLegend(TH1F *h1,     TH1F *h2,
			  char *style1, char *style2,
			  char *title1, char *title2);

void     SetMultiGraph(TMultiGraph *gr,
		       TAxis *&x, TAxis *&y,
		       char  *xt, char  *yt);

void     SetTGraphErrors(TGraphErrors *gr,
			 int           lWidth =    2,
			 int           mColor =   50,
			 double        mSize  = 1.25,
			 int           mStyle =   20);

void     SetTGraphAsymmErrors(TGraphAsymmErrors *g,
			      Color_t            mColor = kBlue,
			      int                mStyle = 20,
			      double             mSize  = 1.0,
			      Color_t            lColor = kBlack,
			      int                lStyle = 1,
			      Width_t            lWidth = 1);

double   eff   (double rs, double ws, double nt);
double   dil   (double rs, double ws);
double   eD2   (double rs, double ws, double nt);

double   effErr(double rs, double rsErr,
		double ws, double wsErr,
		double nt, double ntErr);

double   dilErr(double rs, double rsErr,
		double ws, double wsErr);

double   eD2Err(double rs, double rsErr,
		double ws, double wsErr,
		double nt, double ntErr);


void SetMaximum(TH1F* h1, TH1F* h2);

#endif
