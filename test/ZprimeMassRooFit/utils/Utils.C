#include "Utils.hh"

#include "TLatex.h"
#include "TF1.h"
#include "TPaveText.h"

//------------------------------------------------------------------------------
// Axis settings
//------------------------------------------------------------------------------
void SetAxis(TAxis *&xaxis, TAxis *&yaxis, char *xtitle, char *ytitle)
{
  xaxis->SetLabelFont  (    42);
  xaxis->SetLabelOffset( 0.015);
  xaxis->SetNdivisions (   505);
  xaxis->SetTitle      (xtitle);
  xaxis->SetTitleColor (kBlack);
  xaxis->SetTitleFont  (    42);
  xaxis->SetTitleOffset(  1.25);
  xaxis->SetTitleSize  (  0.05);

  yaxis->CenterTitle   ( kTRUE);
  yaxis->SetLabelFont  (    42);
  yaxis->SetLabelOffset( 0.020);
  yaxis->SetNdivisions (   505);
  yaxis->SetTitle      (ytitle);
  yaxis->SetTitleColor (kBlack);  
  yaxis->SetTitleFont  (    42);
  yaxis->SetTitleOffset(  1.70);
  yaxis->SetTitleSize  (  0.05);
}

//------------------------------------------------------------------------------
// TF1 axis settings
//------------------------------------------------------------------------------
void axisTF1(TF1 *func, char *xtitle, char *ytitle)
{
  func->SetMarkerSize(0.5);
  func->SetMarkerStyle(kFullCircle);

  func->SetTitle("");
  
  TAxis *xaxis = (TAxis*)func->GetXaxis();
  TAxis *yaxis = (TAxis*)func->GetYaxis();
  
  SetAxis(xaxis,yaxis,xtitle,ytitle);
}

//------------------------------------------------------------------------------
// TH1F axis settings
//------------------------------------------------------------------------------
void axis1F(TH1F *histo, char *xtitle, char *ytitle)
{
  histo->SetMarkerSize(0.5);
  histo->SetMarkerStyle(kFullCircle);

  histo->SetTitle("");
  
  TAxis *xaxis = (TAxis*)histo->GetXaxis();
  TAxis *yaxis = (TAxis*)histo->GetYaxis();

  SetAxis(xaxis,yaxis,xtitle,ytitle);
}

//------------------------------------------------------------------------------
// TH1D axis settings
//------------------------------------------------------------------------------
void axis1D(TH1D *histo, char *xtitle, char *ytitle)
{
  histo->SetMarkerSize(0.5);
  histo->SetMarkerStyle(kFullCircle);

  histo->SetTitle("");
  
  TAxis *xaxis = (TAxis*)histo->GetXaxis();
  TAxis *yaxis = (TAxis*)histo->GetYaxis();

  SetAxis(xaxis,yaxis,xtitle,ytitle);
}

//------------------------------------------------------------------------------
// TH2F axis settings
//------------------------------------------------------------------------------
void axis2F(TH2F *histo, char *xtitle, char *ytitle)
{
  histo->SetMarkerSize(0.5);
  histo->SetMarkerStyle(kFullCircle);

  histo->SetTitle("");

  TAxis *xaxis = (TAxis*)histo->GetXaxis();
  TAxis *yaxis = (TAxis*)histo->GetYaxis();

  SetAxis(xaxis,yaxis,xtitle,ytitle);
}

//------------------------------------------------------------------------------
// Draw TH1F
// Draw(histo,fColor,fStyle,drawOpt,lColor,lWidth,mColor,mSize,mStyle)
//------------------------------------------------------------------------------
void Draw(TH1F *histo,   int    fColor, int fStyle,
          char *drawOpt, int    lColor, int lWidth,
	  int   mColor,  double mSize,  int mStyle)
{
  histo->SetDirectory(0     );
  histo->SetFillColor(fColor);
  histo->SetFillStyle(fStyle);

  histo->SetLineColor  (lColor);
  histo->SetLineWidth  (lWidth);
  histo->SetMarkerColor(mColor);
  histo->SetMarkerSize (mSize );
  histo->SetMarkerStyle(mStyle);

  if (histo != NULL) histo->Draw(drawOpt);
}

//------------------------------------------------------------------------------
// Draw TH1D
//------------------------------------------------------------------------------
void Draw1D(TH1D    *histo,   Color_t fColor, int fStyle,
	    char    *drawOpt, Color_t lColor, int lWidth,
	    Color_t  mColor,  double  mSize,  int mStyle)
{
  histo->SetDirectory(0     );
  histo->SetFillColor(fColor);
  histo->SetFillStyle(fStyle);

  histo->SetLineColor  (lColor);
  histo->SetLineWidth  (lWidth);
  histo->SetMarkerColor(mColor);
  histo->SetMarkerSize (mSize );
  histo->SetMarkerStyle(mStyle);

  if (histo != NULL) histo->Draw(drawOpt);
}

//------------------------------------------------------------------------------
// Draw TH1F and TF1
//------------------------------------------------------------------------------
void Draw(TH1F *hist, TF1* func, float lo, float hi)
{
  char yTitle[200];

  sprintf(yTitle,"Candidates per %.0f MeV/c^{2}",
          double(hist->GetBinWidth(1))*1e3);

  if (lo != -999)
    hist->GetXaxis()->SetRange(hist->FindBin(lo),hist->FindBin(hi));

  hist->SetDirectory (   0);
  hist->SetFillColor ( 188);
  hist->SetFillStyle (1001);
  hist->SetLineColor (  10);
  hist->SetLineWidth (   0);
  hist->SetMarkerSize(   0);
  
  hist->Draw("hist");
  func->SetLineWidth(3);
  func->Draw("same");

  axis1F(hist,"Candidate Mass [GeV/c^{2}]",yTitle);
}

//------------------------------------------------------------------------------
// Set TLegend
//------------------------------------------------------------------------------
TLegend *SetLegend(float x0, float y0, float x1, float y1)
{
  TLegend *leg = new TLegend(x0,y0,x1,y1);
  leg->SetFillStyle (0);            
  leg->SetBorderSize(0);
  leg->SetTextAlign(12);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  return leg;
}

//------------------------------------------------------------------------------
// Set TPaveText
//------------------------------------------------------------------------------
void SetPave(char *title, float x, float y)
{
  TPaveText *pv = new TPaveText(x,y,x+0.2,y+0.03,"ndc");
  pv->SetBorderSize(0);
  pv->SetFillColor(0);
  pv->SetTextAlign(33);
  pv->SetTextFont(42);
  pv->SetTextSize(0.035);
  pv->AddText(title);
  pv->Draw();
}

//------------------------------------------------------------------------------
// CDF Run II Monte Carlo
//------------------------------------------------------------------------------
void cdfPave()
{
  TPaveText *tPave = new TPaveText(0.35,0.95,0.55,0.98,"ndc");
  tPave->SetBorderSize(0);
  tPave->SetFillColor(0);
  tPave->SetTextAlign(23);
  tPave->SetTextFont(42);
  tPave->SetTextSize(0.045);
  tPave->AddText("CDF Run II Monte Carlo");
  tPave->Draw();
}

//------------------------------------------------------------------------------
// Set TF1
//------------------------------------------------------------------------------
void SetTF1(TF1 *func,
	    int fStyle,
	    int fColor,
	    int lStyle,
	    int lColor,
	    int lWidth,
	    int nPoints)
{
  func->SetFillStyle(fStyle);
  func->SetFillColor(fColor);
  func->SetLineStyle(lStyle);
  func->SetLineColor(lColor);
  func->SetLineWidth(lWidth);
  func->SetNpx(nPoints);
}

//------------------------------------------------------------------------------
// Chi-squared of two histograms
//------------------------------------------------------------------------------
double chi2RebinTest(TH1F   *hist1,
                     TH1F   *hist2,
                     double &chi2,
		     int    &NDF,
                     double  low,
		     double  high,
		     int     verbose,
		     int     threshold,
		     int     fitParams)
{
  double chi2Sum = 0.0;

  double nEvts1     = 0;
  double nEvts2     = 0;
  int    nBins      = 0;
  int    binStarted = 1;

  double evtErr1 = 0.0;
  double evtErr2 = 0.0;
  double binChi2 = 0.0;

  int minBin, maxBin; 
  
  if (low == 0 && high == 0) {
    minBin = 1;
    maxBin = hist1->GetNbinsX(); 
  }
  else {
    minBin = hist1->FindBin(low); 
    maxBin = hist1->FindBin(high); 
  }
  
  for (int i=minBin; i<=maxBin; i++) {
    nEvts1 += (double)hist1->GetBinContent(i);
    nEvts2 += (double)hist2->GetBinContent(i);
    
    evtErr1 += (double)hist1->GetBinError(i)*hist1->GetBinError(i);
    evtErr2 += (double)hist2->GetBinError(i)*hist2->GetBinError(i);
    
    if (nEvts1 < threshold && i !=  hist1->GetNbinsX()) continue;
    
    binChi2 = (nEvts1 - nEvts2) * (nEvts1 - nEvts2);
    
    if (evtErr1 != 0) binChi2 /= evtErr1;
    
    chi2Sum += binChi2;

    if (verbose >= 1)
      printf("%3d)[%3d,%3d] data: %7.2f th: %7.2f  chi2: %7.2f total: %8.2f\n",
             nBins,binStarted,i,nEvts1,nEvts2,binChi2,chi2Sum); 
 
    binStarted = i+1;
    nEvts1     = 0.0;
    nEvts2     = 0.0;
    evtErr1    = 0.0;
    evtErr2    = 0.0;
    nBins++;
  }
  
  chi2 = chi2Sum;
  NDF = nBins - 1 - fitParams;
  
  return TMath::Prob(chi2Sum,NDF);
}

//------------------------------------------------------------------------------
// Chi-squared of a histogram and a function
//------------------------------------------------------------------------------
double chi2RebinTest(TH1F   *hist,
                     TF1    *fun,
		     double &chi2,
		     int    &NDF,
                     double  low,
		     double  high,
		     int     verbose,
		     int     threshold)
{
  double chi2Sum = 0.0;

  double nEvts      = 0;
  double evtErr     = 0;
  int    nBins      = 0;
  int    binStarted = 1;

  double theory, lowEdge, highEdge, binChi2;
  
  double binWidth = hist->GetXaxis()->GetBinWidth(1);

  int minBin, maxBin; 
  
  if (low == 0 && high == 0) {
    minBin = 1;
    maxBin = hist->GetNbinsX();
  }
  else {
    minBin = hist->FindBin(low); 
    maxBin = hist->FindBin(high); 
  }

  if (verbose >= 1)
    printf("\n >>> Bins chi-squared information <<<\n\n");

  for (int i=minBin; i<maxBin; i++) {
    lowEdge  = hist->GetXaxis()->GetBinLowEdge(binStarted);
    highEdge = hist->GetXaxis()->GetBinUpEdge(i);

    if (lowEdge  < low ) lowEdge  = low;
    if (highEdge > high) highEdge = high;
    
    nEvts  += hist->GetBinContent(i);    
    evtErr += pow(hist->GetBinError(i),2);

    if (nEvts < threshold && i != hist->GetNbinsX()) continue;
    
    theory = fun->Integral(lowEdge,highEdge) / binWidth;
    
    if (evtErr > 1e-100) {
      binChi2 = (nEvts - theory) * (nEvts - theory) / evtErr;
      nBins++;
    }  
    else {
      binChi2 = 0.0;
      if (verbose > 0 && nEvts != 0)
        printf("\n ThesisStyle::chi2RebinTest - WARNING: error is 0");
    }

    chi2Sum += binChi2;

    if (verbose >= 1) {
      printf("%3d)[%3d,%3d] data:%6.2f th:%6.2f",
    	     nBins,binStarted,i,nEvts,theory);
      printf(" err^2:%7.4f  chi2:%6.2f total:%7.2f\n",
    	     evtErr,binChi2,chi2Sum);
    }

    binStarted = i+1;
    nEvts      = 0;
    evtErr     = 0;
  }

  int npar  = fun->GetNpar();
  int fixed = 0;
  double parmin, parmax;

  for (int i=0; i<npar; i++) {
    fun->GetParLimits(i,parmin,parmax);
    if ((parmin*parmax != 0) && (parmin >= parmax))
      fixed++;
  }

  if (verbose >= 0)
    printf("\n Function has %d parameters, %d fixed\n\n",npar,fixed);

  int NDFSum = nBins - npar + fixed;
  
  chi2 = chi2Sum;
  NDF  = NDFSum - 1;

  double prob = TMath::Prob(chi2Sum,NDFSum);
  return prob;
}

//------------------------------------------------------------------------------
// Print the probability between two histograms in a TCanvas
//------------------------------------------------------------------------------
void PrintFitProbability(TCanvas *canvas,
                         TH1F    *hist1,
			 TH1F    *hist2,
                         double   low,
			 double   high,
			 int      verbose,
			 int      threshold,
			 int      padNum,
			 int      fitParams)
{
  TPad *pad = (TPad*)canvas->GetPad(padNum);
  PrintFitProbability(pad,hist1,hist2,low,high,verbose,threshold,fitParams);
}

//------------------------------------------------------------------------------
// Print the probability between a histogram and a function in a TCanvas
//------------------------------------------------------------------------------
void PrintFitProbability(TCanvas *canvas,
                         TH1F    *hist,
			 TF1     *fun,
                         double   low,
			 double   high,
			 int      verbose,
			 int      threshold,
			 int      padNum)
{
  TPad *pad = (TPad*)canvas->GetPad(padNum);
  PrintFitProbability(pad,hist,fun,low,high,verbose,threshold);
}

//------------------------------------------------------------------------------
// Print the probability between two histograms in a TPad
//------------------------------------------------------------------------------
void PrintFitProbability(TPad   *pad,
                         TH1F   *hist1,
			 TH1F   *hist2,
                         double  low,
			 double  high,
			 int     verbose,
			 int     threshold,
			 int     fitParams)
{
  char   title[100];
  double chi2;
  int    ndf;
  double prob  = 1e2;
  double probK = 1e2;

  prob *=
    chi2RebinTest(hist1,hist2,chi2,ndf,low,high,verbose,threshold,fitParams);

  probK *= hist1->KolmogorovTest(hist2);

  sprintf(title,"#chi^{2} / NDF = %.2f / %d, Prob = %.2f%%",chi2,ndf,prob);
  sprintf(title,"%s, K-Prob = %.2f%%",title,probK);
  PrintIt(pad,title);
}

//------------------------------------------------------------------------------
// Print the probability between a histogram and a function in a TPad
//------------------------------------------------------------------------------
void PrintFitProbability(TPad   *pad,
                         TH1F   *hist,
			 TF1    *fun,
                         double  low,
			 double  high,
			 int     verbose,
			 int     threshold)
{
  char   title[100];
  double chi2;
  int    ndf;
  double prob = 1e2;

  prob *= chi2RebinTest(hist,fun,chi2,ndf,low,high,verbose,threshold);

  sprintf(title,"#chi^{2} / NDF = %.2f / %d, Prob = %.2f%%",chi2,ndf,prob);
  PrintIt(pad,title);
}

//------------------------------------------------------------------------------
// Real printing called by PrintFitProbability
//------------------------------------------------------------------------------
void PrintIt(TPad *pad, char *title)
{
  pad->Update();

  TLatex *latex = new TLatex();
  latex->SetTextFont(  42);  
  latex->SetTextSize(0.04);

  double xmin = pad->GetUxmin();
  double xmax = pad->GetUxmax();
  double ymin = pad->GetUymin();
  double ymax = pad->GetUymax();		     

  double xpos = xmin + 0.50 * (xmax - xmin);
  double ypos = ymax + 0.06 * (ymax - ymin);

  if (pad->GetLogy())
    ypos = pow(10,ypos);

  latex->SetTextAlign(22);
  latex->DrawLatex(xpos,ypos,title);
}

//------------------------------------------------------------------------------
// Draw a TLegend of two histograms TH1F
//------------------------------------------------------------------------------
void ComparisonLegend(TH1F *h1,     TH1F *h2,
		      char *style1, char *style2,
		      char *title1, char *title2)
{
  TLegend *tl = SetLegend(0.456044,
			  0.743007,
			  0.734432,
			  0.882867);

  tl->AddEntry(h1,title1,style1);
  tl->AddEntry(h2,title2,style2);
  tl->Draw();
}

//------------------------------------------------------------------------------
// Label for 'CDF Run II Preliminary' and 'L = 355 pb-1'
//------------------------------------------------------------------------------
void PRLabel(TPad   *pad, 
	     char   *title,
	     double  delta,
	     int     tAlign,
	     double  tSize)
{
  pad->Update();

  TLatex *latex = new TLatex();
  latex->SetTextFont(   42);  
  latex->SetTextSize(tSize);

  double xmin = pad->GetUxmin();
  double xmax = pad->GetUxmax();
  double ymin = pad->GetUymin();
  double ymax = pad->GetUymax();		     

  double xpos = xmin + delta * (xmax - xmin);
  double ypos = ymax +  0.03 * (ymax - ymin);

  if (pad->GetLogy())
    ypos = pow(10,ypos);

  latex->SetTextAlign(tAlign);
  latex->DrawLatex(xpos,ypos,title);
}

//------------------------------------------------------------------------------
// TMultiGraph settings
//------------------------------------------------------------------------------
void SetMultiGraph(TMultiGraph *gr,
		   TAxis       *&x, TAxis *&y,
		   char        *xt, char  *yt)
{
  gr->SetTitle("");
  gr->Draw("apz");

  x = gr->GetXaxis();
  y = gr->GetYaxis();

  x->SetLabelFont(42);
  y->SetLabelFont(42);
  x->SetLabelOffset(0.02);
  y->SetLabelOffset(0.02);
  x->SetNdivisions(507);
  y->SetNdivisions(507);
  x->SetTitle(xt);
  y->SetTitle(yt);
  x->SetTitleFont(42);
  y->SetTitleFont(42);
  x->SetTitleOffset(1.40);
  y->SetTitleOffset(1.85);
  x->SetTitleSize(0.05);
  y->SetTitleSize(0.05);
}

//------------------------------------------------------------------------------
// TGraphErrors settings
//------------------------------------------------------------------------------
void SetTGraphErrors(TGraphErrors *gr,
		     int           lWidth,
		     int           mColor,
		     double        mSize,
		     int           mStyle)
{
  gr->SetLineWidth  (lWidth);
  gr->SetMarkerColor(mColor);
  gr->SetMarkerSize ( mSize);
  gr->SetMarkerStyle(mStyle);
}

//------------------------------------------------------------------------------
// Set TGraphAsymmErrors
//------------------------------------------------------------------------------
void SetTGraphAsymmErrors(TGraphAsymmErrors *g,
			  Color_t            mColor,
			  int                mStyle,
			  double             mSize,
			  Color_t            lColor,
			  int                lStyle,
			  Width_t            lWidth)
{
  g->SetMarkerColor(mColor);
  g->SetMarkerStyle(mStyle);
  g->SetMarkerSize (mSize );
  g->SetLineColor  (lColor);
  g->SetLineStyle  (lStyle);
  g->SetLineWidth  (lWidth);
}

//------------------------------------------------------------------------------
// Flavor tagging efficiency
//------------------------------------------------------------------------------
double eff(double rs, double ws, double nt)
{
  double a = fabs(rs);
  double b = fabs(ws);
  double c = fabs(nt);

  if (a+b+c == 0)
    return 0.0;

  double feff = (a+b) / (a+b+c);
  return feff;
}

//------------------------------------------------------------------------------
// Flavor tagging dilution
//------------------------------------------------------------------------------
double dil(double rs, double ws)
{  
  double a = fabs(rs);
  double b = fabs(ws);

  if (a+b == 0)
    return 0.0;

  double fasym = (a-b) / (a+b);
  return fasym;
}

//------------------------------------------------------------------------------
// Flavor tagging eD2
//------------------------------------------------------------------------------

double eD2(double rs, double ws, double nt)
{
  double a = fabs(rs);
  double b = fabs(ws);
  double c = fabs(nt);

  double feD2 = eff(a,b,c) * dil(a,b) * dil(a,b);
  return feD2;
}

//------------------------------------------------------------------------------
// Flavor tagging efficiency error
//------------------------------------------------------------------------------
double effErr(double rs, double rsErr,
	      double ws, double wsErr,
	      double nt, double ntErr)
{  
  double a    = fabs(rs);
  double b    = fabs(ws);
  double c    = fabs(nt);
  double aErr = fabs(rsErr);
  double bErr = fabs(wsErr);
  double cErr = fabs(ntErr);

  if (a+b+c == 0)
    return 0.0;

  double Err = sqrt(pow(c,2)*(pow(aErr,2) + pow(bErr,2)) + pow((a+b)*cErr,2));
  Err /= ((a+b+c) * (a+b+c));
  return Err;
}

//------------------------------------------------------------------------------
// Flavor tagging dilution error
//------------------------------------------------------------------------------
double dilErr(double rs, double rsErr,
	      double ws, double wsErr)
{  
  double a    = fabs(rs);
  double b    = fabs(ws);
  double aErr = fabs(rsErr);
  double bErr = fabs(wsErr);

  if (a+b == 0)
    return 0.0;

  double Err = (2.0 / pow(a+b,2)) * sqrt(pow(a*bErr,2) + pow(b*aErr,2));
  return Err;
}

//------------------------------------------------------------------------------
// Flavor tagging eD2 error
//------------------------------------------------------------------------------
double eD2Err(double rs, double rsErr,
	      double ws, double wsErr,
	      double nt, double ntErr)
{
  double a    = fabs(rs);
  double b    = fabs(ws);
  double c    = fabs(nt);
  double aErr = fabs(rsErr);
  double bErr = fabs(wsErr);
  double cErr = fabs(ntErr);

  if (a+b+c==0 || a+b==0)
    return 0.0;

  double factor1 = dil(a,b) / pow(a+b+c,2);
  double factor2 = pow(4.0*(a+b+c)/(a+b),2) * (pow(a*bErr,2) + pow(b*aErr,2));
  factor2 += pow((a-b)*cErr,2);

  double Err = factor1 * sqrt(factor2);
  return Err;
}
