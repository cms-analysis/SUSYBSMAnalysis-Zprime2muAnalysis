#include <iostream>
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>
#include "TMinuit.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TStopwatch.h"

#include "TGraphErrors.h"
TF1* tempfit;
//bool const _useData = false;

bool const _useCaching = true;
bool const _isTest = true;
double const _lowerBound = 120;
double const _upperBound = 3000;
int const _numPseudo = 400;

double const _startingExpectedSignal = 0;
double const _maxExpectedSignal = 14.0;
double const _scanningStep = 0.1;


int _theZPrimeMass = 1750;
double _jitterLambdaB = 1.0;
double _jitterZPrimeWidth= 1.0;


TRandom3* rand3;

static double const drellYan(double* x, double* par);
static double const voigt(double* x, double* par);

TF1* fdrellYan, *fSignal, *fmix, *fBackground;
TF1* fGenSignal, *fGenBackground;
//TF1* fdrellYan, *fSignal, *fmix, *fBackground;
//double const avgBackground = 23.67;	//60 pb-1
//double avgBackground = 19.5241; //50 pb-1

TH1F* hbkg, *hsig, *hsnb;

double const GetMass(TF1* func, TH1F* reverseTable);

TH1F* h;
//TH1F* hBackgroundTable, *hSignalTable;

std::vector<double> GetToyMC(double const lambdaS, double lambdaB); 


TGraphErrors* grePvalueVsLambdaS;
TGraphErrors* grePvalueVsLambdaSprot;
TGraphErrors* grePvalueVsLambdaSPL;

TFile* _outFile;

TMinuit* minuit = 0;
//
// minimization function for Minuit
//
std::vector<double> _MASSES;
void minuitPdf(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag); 
double const EvalPdf(double const lambdaS, double const lambdaB); 
double const CalculateMedian(std::vector<double> vec);
double const CalculateMean(std::vector<double> vec);
double const CalculateRMS(std::vector<double> vec, double const mean);
double const CalculatePValue(std::vector<double> const& vec, double const t0);

std::pair<int,double> GetMinuitAnswer(std::vector<double> const toymc, double const lambdaS, double const lambdaB, bool const fixLambdaS);

//
// z pole fit
//
double ZPole(Double_t* x, Double_t* par);
//
// tail fit
//
double shitTail(double* x, double* par);
//
// mix the two
//
double MixFunc(double* x, double* par);

std::vector<double> const PutDataIntoVec(TTree* tree);

double Signal500(double* x, double* par); // this was named when I was adding templates...  

double ExpoConvGaus(double* x, double* par);

double BckgFunc(double* x, double* par);

TF1* fTtbar;

void const SetSignalParameters(int const mass);

double const GetLimitFromExpoFit(TGraphErrors* gre);
double const CalcLimitFomExpFitVal(TF1* thefit);

void const SetSignalParametersFixed(int const mass_);
double GetZssmGamma(int const _mass);
double GetResolution(int const _mass);
std::pair<double,double> GetBand(TGraphErrors* gre, double const pval);

std::vector<std::pair<double,double> > FitCacheVec;
double const EvalPdfCached(double const lambdaS, double const lambdaB); 
double const GetPValError(double const pval, int const num);



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void limits(int inZpmass = 1750, double jitterLambdaB = 0, double jitterZPrimeWidth = 0) {
	_theZPrimeMass = inZpmass;
	rand3 = new TRandom3(0);
	_jitterLambdaB 		+= jitterLambdaB;
	_jitterZPrimeWidth 	+= jitterZPrimeWidth;

	TGraphErrors* greLLRVsLambdaS= new TGraphErrors();

/*
//TF1* fmix;
   1  tail A       2.44298e+01   7.15120e-03  -4.75412e-07   1.70094e-04
   2  tail #alpha   6.99626e+00   2.47564e-02  -1.25056e-08  -3.19569e-02
   3  tail #kappa   2.03121e-01   3.95772e-04   4.13567e-09  -1.39444e-01
   4  pole A       9.45741e+06   1.91097e+05   1.52909e-10   5.22108e-03
   5  pole #sigma   2.63045e+00   1.41963e-02  -1.60154e-10  -3.02380e-05
   6  pole #theta   1.80690e-02   1.92189e-04   1.16035e-10  -2.45090e-02
   7  pole B       5.12014e+01   1.05762e+00   2.68504e-07   9.44424e-09
   8  pole C       0.00000e+00     fixed    
   9  pole #kappa   0.00000e+00     fixed    
  10  turnon width   2.19753e-02   5.53804e-04  -2.34669e-10   7.58831e-05
  11  Norm         1.00000e+00     fixed    

   1  tail A       2.39254e+01   8.60747e-02   1.17352e-04   2.07462e-04
   2  tail #alpha   6.68434e+00   5.22540e-02   1.29660e-06   4.61947e-02
   3  tail #kappa   2.07129e-01   7.02808e-04   7.35396e-07   1.68286e-01
   4  pole A       9.13142e+06   2.02697e+05   1.52567e-06   1.33285e-01
   5  pole #sigma   2.37072e+00   6.97259e-03   2.09848e-04   6.17631e-04
   6  pole #theta   1.77966e-02   2.11498e-04   1.27014e-06  -2.07846e-02
   7  pole B       5.58682e+01   1.00144e+00   1.86666e-02  -1.47247e-05
   8  pole C       0.00000e+00     fixed    
   9  pole #kappa   0.00000e+00     fixed    
  10  turnon width   2.24108e-02   5.47565e-04   2.41905e-05   2.06586e-02
  11  Norm         1.00000e+00     fixed    

*/

	
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// DY parts
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
	fmix = new TF1("fmix",MixFunc, _lowerBound, _upperBound, 11);
	fmix->SetRange(_lowerBound,_upperBound);
	fmix->FixParameter(0 	,	2.39254e+01	);
	fmix->FixParameter(1 	, 	6.68434e+00	); 
	fmix->FixParameter(2 	, 	2.07129e-01	); 
	fmix->FixParameter(3 	,	9.13142e+06	);  
	fmix->FixParameter(4 	, 	2.37072e+00	); 
	fmix->FixParameter(5 	, 	1.77966e-02	); 
	fmix->FixParameter(6 	,	5.58682e+01	);  
	fmix->FixParameter(7 	,	0.00000e+00	);  
	fmix->FixParameter(8 	, 	0.00000e+00	); 
	fmix->FixParameter(9	,  	2.24108e-02	);
	fmix->SetParameter(10,1);

	fmix->SetNpx(10000);

	fmix->SetParameter(10, 1./fmix->Integral(_lowerBound, _upperBound));
//	fmix->Draw();
//
//	hBackgroundTable->Draw();

///////////////////////////////////////////////////////////////////////////////
//
// ttbar 
//
///////////////////////////////////////////////////////////////////////////////
	fTtbar = new TF1("ttbar", ExpoConvGaus, _lowerBound, _upperBound, 4);  
	fTtbar->SetRange(_lowerBound,_upperBound);
	fTtbar->SetNpx(10000);
	fTtbar->SetParameters(1, 2.50938e+01, 1.73358e-02, 	4.80605e+01);

	fTtbar->SetParameter(0, 1./fTtbar->Integral(_lowerBound, _upperBound));



//	fTtbar->Draw();
//	return;


//	RooConstVar ttbarSigma	("ttbarSigma"	, "ttbarGamma"	,2.50938e+01);//,20, 200);
//	RooConstVar ttbarGamma	("ttbarGamma"	, "ttbarSigma"	,1.73358e-02);//,20, 200);
//	RooConstVar ttbarMean	("ttbarMean"	, "ttbarMean"	,4.80605e+01);//,20, 200);

//double ExpoConvGaus(double* x, double* par) {


	fdrellYan 	= fmix;
//	fBackground	= fmix;


//	double BckgFunc(double* x, double* par){

	fBackground = new TF1("fBackground",BckgFunc, _lowerBound, _upperBound,1); 
	fBackground->FixParameter(0,1.0);

//	fBackground->Draw();

///////////////////////////////////////////////////////////////////////////////
//
// get the signal 
//
///////////////////////////////////////////////////////////////////////////////

/*
	fSignal = new TF1("fSignal", Signal500, _lowerBound, _upperBound,8);
		fSignal->SetRange(_lowerBound,_upperBound);
		fSignal->SetParNames("A", "mean", "sigma", "lg", "alpha", "B", "sigma2", "Norm");

	SetSignalParameters(inZpmass);
*/
	fSignal = new TF1("fSignal",voigt,_lowerBound,_upperBound,4);
		fSignal->SetRange(_lowerBound,_upperBound);
		fSignal->SetParNames("A", "mean", "sigma", "lg");
		SetSignalParametersFixed(inZpmass);


	fSignal->SetParameter("Norm",1./fSignal->Integral(_lowerBound,_upperBound));
	fSignal->SetNpx(10000);
//	fSignal->Draw();

//
// clone the templates for generators
//
	fGenSignal = (TF1*)fSignal->Clone("fGenSignal");
	fGenBackground = (TF1*)fBackground->Clone("fGenBackground");
	fGenSignal->SetParameter("Norm",1);	
	fGenSignal->SetParameter("lg",fGenSignal->GetParameter("lg")*_jitterZPrimeWidth); 
	fGenSignal->SetParameter("Norm",1./fSignal->Integral(_lowerBound,_upperBound));

///////////////////////////////////////////////////////////////////////////////
//
//	


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//
// do the limit calculation
//
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



	int numPseudo = 100; // hold-over from when there was no data

//
// start up minuit one time
//
	minuit	= new TMinuit(2);
	minuit	-> SetPrintLevel(-1);
	minuit	-> SetFCN(minuitPdf);

	double lambdaS = 0;
	double lambdaB = 77.;
//
// calculate the reference likelihoods
//
	std::map<double,std::vector<double> > refTable;

//
// use toy MC to get an estimate
// actually... it's set to get data.  I hacked out the use of toy MC
//
	if (true){
		TString sDataName ="/Users/kypreos/physics/templates/limitSetting/data123.root"; 
		sDataName ="/Users/kypreos/physics/templates/limitSetting/merged.root"; 
		sDataName ="/Users/kypreos/physics/templates/limitSetting/mergedNewShit.root"; 
		sDataName = "dimuons_255nb.root";
		sDataName = "data/dimuons_merged_2900InvNb_nocosmics_new.root";
		sDataName = "data/muons1090.root";
//		sDataName = "dimuons_20100825.root";

		TFile inData(sDataName, "open");
		std::cout<<"Getting data file from: "<<sDataName<<std::endl;	

//	
// the the tree with masses.  roofit does it the same way. snap.
// 
//
	
//		TTree* inTree = (TTree*)inData.Get("tree");
		TTree* inTree = (TTree*)inData.Get("muons");
		std::vector<double> dataMasses; dataMasses.clear();
		dataMasses =  PutDataIntoVec(inTree); 
		lambdaB = dataMasses.size(); 
		numPseudo = 1;
		std::vector<double> const& toymc= dataMasses;//GetToyMC(0, lambdaB, hBackgroundTable, hSignalTable);



		for (int jexp = 0; jexp < numPseudo; ++jexp) {
//			std::vector<double> toymc=GetToyMC(0, lambdaB, hBackgroundTable, hSignalTable);
			lambdaS = 0;	
			double llfit = 0, llfix =0;
			std::pair<int,double> answerPair = GetMinuitAnswer(toymc, lambdaS, lambdaB,false);
			if (answerPair.first ==0) llfit = answerPair.second;
		
//double const _startingExpectedSingal = 1.0;
//double const _scanning step = 0.1;
//			for (lambdaS = 0.3; lambdaS < 5.5; lambdaS += 0.1){
			for (lambdaS = _startingExpectedSignal; lambdaS < _maxExpectedSignal; lambdaS += _scanningStep){
	
				std::pair<int,double> answerPairfix = GetMinuitAnswer(toymc, lambdaS, lambdaB,true);
				if (answerPairfix.first ==0) llfix = answerPairfix.second;
				double llr = llfit - llfix;
				
				printf("lambdaS = %3.2f\tfit = %5f\tfix = %5f\tllr = %5f\n",
					lambdaS, llfit,llfix,llr);
					std::vector<double> tmp; tmp.clear(); tmp.push_back(llr);
					refTable.insert(std::pair<double,std::vector<double> > (lambdaS, tmp));	
					greLLRVsLambdaS->SetPoint(greLLRVsLambdaS->GetN(),lambdaS,llr);
					
//				if (lambdaS < 2) lambdaS+=1;
			}
	
		}	
	} 
	
	


	
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//
// OK. here's where I really do the limit
//
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



	grePvalueVsLambdaS = new TGraphErrors();
	grePvalueVsLambdaSprot = new TGraphErrors();

	grePvalueVsLambdaSPL = new TGraphErrors();
//
// I decided to do it where I calculate the reference values and then iterate on the values that were calculated. 
// I could have calculated a reference for each point, 
// but that would have put more shit in shit and I wasn't in the mood to debug that.
//

	TString sOutname = TString::Format("out%d/fne_%04d_v1.root",_numPseudo,inZpmass);
	if (_isTest) sOutname = "test.root";
	_outFile = new TFile(sOutname,"recreate");
	_outFile->cd();
	fSignal->Write();
	fBackground->Write();
	_outFile->mkdir("pseudoExp");
	_outFile->GetDirectory("pseudoExp")->cd();

	TStopwatch timer;
	timer.Start();


	double refValPL	= refTable.begin()->second.front();
	refValPL = greLLRVsLambdaS->GetY()[0];
	std::cout<<"ref val for PL: "<<refValPL<<std::endl;
	for (std::map<double,std::vector<double> >::const_iterator tableit = refTable.begin(), 
		tableitend = refTable.end(); tableit!=tableitend; ++tableit) {

	

		rand3->SetSeed(0);	
		lambdaS = tableit->first;
		std::vector<double> vec = tableit->second;		
		TString dirname = TString::Format("lambdaS_%02d_%02d",(int)lambdaS, int(lambdaS*100)-100*int(lambdaS)  ); 
//(lambdaS-((int)lambdaS))*1000);
//		int assmunch = int(lambdaS*100);
//		std::cout<<assmunch<<std::endl;;	
		_outFile->cd();
		_outFile->GetDirectory("pseudoExp")->mkdir(dirname);
		_outFile->GetDirectory("pseudoExp")->GetDirectory(dirname)->cd();
		TTree* _llrTree = new TTree("llrTree","llrTree");
		float _storeLlr = 0;
		float _storeLlrPL = 0;
		_llrTree->Branch("llr", &_storeLlr,"llr/F"); 
		_llrTree->Branch("llrPL", &_storeLlrPL,"llr/F"); 

		hbkg = new TH1F("hbkg", "DY background generation", 1000,_lowerBound,_upperBound );
		hsig = new TH1F("hsig", "Zprime signal generation", 1000,_lowerBound,_upperBound );
		

		double refVal(0);

		std::cout<<"lambdaS = "<<tableit->first
			<<"\t test-val = "<<vec.front();

//		numPseudo = 1000;
		refVal = vec.front();

		std::vector<double> llrVals; llrVals.clear();
		std::vector<double> llrValsPL; llrValsPL.clear();

		for (int jexp = 0; jexp < _numPseudo; ++jexp) {
			std::vector<double> toymc=GetToyMC(lambdaS, lambdaB);//, hBackgroundTable, hSignalTable);
			std::pair<int,double> answerFit = GetMinuitAnswer(toymc, lambdaS, lambdaB,false);
			if (answerFit.first != 0) continue;
			std::pair<int,double> answerFix = GetMinuitAnswer(toymc, lambdaS, lambdaB,true);
			if (answerFix.first != 0) continue;
			std::pair<int,double> answerFix0 = GetMinuitAnswer(toymc, 0, lambdaB,true);
			if (answerFix0.first != 0) continue;
				
			double llr = answerFit.second - answerFix.second;
			llrVals.push_back(llr);
			double llr0 = answerFit.second - answerFix0.second;
			llrValsPL.push_back(llr0);


			_storeLlr = llr;
			_storeLlrPL = llr0;
			_llrTree->Fill();
		}
	
		sort(llrVals.begin(),llrVals.end());
		sort(llrValsPL.begin(),llrValsPL.end());
	

		double meanLlr 		= CalculateMean(llrVals); 
		double medianLlr 	= CalculateMedian(llrVals); 
		double rmsLlr		= CalculateRMS(llrVals,meanLlr);		

		std::cout
			<<"\t median = "<<medianLlr
			<<"\t mean = "<<meanLlr
			<<"\t rms = "<<rmsLlr;		


		double delta = 1.2;
		delta = rmsLlr;

		double pval		= CalculatePValue(llrVals,refVal);  //FC
		double pvalProt	= CalculatePValue(llrVals,refVal+rmsLlr); //protected FC (often suggested)
//		std::cout<<std::endl<<refValPL<<std::endl;

		double pvalPL	= 1.-CalculatePValue(llrValsPL,greLLRVsLambdaS->GetY()[0]); //profile likelihood


		TH1F* ht0 = new TH1F("ht0", "FC test statistic", 1, 0, 1);
		TH1F* ht0_prot = new TH1F("ht0_prot", "test statistic protected", 1, 0, 1);
		TH1F* ht1 = new TH1F("ht1", "Profile likelihood test statistic", 1, 0, 1);
		ht0->SetBinContent(1,refVal);
		ht0_prot->SetBinContent(1,refVal+rmsLlr);
		ht1->SetBinContent(1,refValPL);

		grePvalueVsLambdaS->SetPoint(grePvalueVsLambdaS->GetN(), lambdaS, pval);
		grePvalueVsLambdaS->SetPointError(grePvalueVsLambdaS->GetN()-1, 0.5*_scanningStep, GetPValError(pval,_numPseudo));
		
		grePvalueVsLambdaSprot->SetPoint(grePvalueVsLambdaSprot->GetN()			, lambdaS			, pvalProt);
		grePvalueVsLambdaSprot->SetPointError(grePvalueVsLambdaSprot->GetN()-1	, 0.5*_scanningStep	, GetPValError(pvalProt,_numPseudo));

		grePvalueVsLambdaSPL->SetPoint(grePvalueVsLambdaSPL->GetN(), lambdaS, pvalPL);
		grePvalueVsLambdaSPL->SetPointError(grePvalueVsLambdaSPL->GetN()-1, 0.5*_scanningStep, GetPValError(pvalPL,_numPseudo));


		std::cout<<"\tp-value (FC) = "<<pval;
		std::cout<<"\tp-value (PL) = "<<pvalPL;
		std::cout<<std::endl;
		_llrTree->Write();
		hsig->Write();
		hbkg->Write();
		ht0->Write();	
		ht0_prot->Write();	
		ht1->Write();	
		_outFile->cd();
//		fSignal->Write("fSignal");
//		fBackground->Write("fBackground");
	
		_llrTree->Delete();	
		hsig->Delete();
		hbkg->Delete();
		ht0->Delete();	
		ht0_prot->Delete();	
		ht1->Delete();	
	
		if (pval < 0.001 && pvalProt < 0.001 && lambdaS > 5) break;
	
	} // end of the loop on test-statistics

//
// clean up and go home to a nice glass of Glenfiddich -- at least the 15-year. That's tasty. 
//
	timer.Stop();
	timer.Print();
	_outFile->cd();

	std::pair<double,double> band1 = GetBand(grePvalueVsLambdaS		,0.05);
	std::pair<double,double> band2 = GetBand(grePvalueVsLambdaSprot	,0.05);
	std::pair<double,double> band3 = GetBand(grePvalueVsLambdaSPL	,0.05);

	printf("%10s: [%5.2f, %5.2f]\n","FC", band1.first,band1.second) ;
	printf("%10s: [%5.2f, %5.2f]\n","profile L", band3.first,band3.second); 
	printf("%10s: [%5.2f, %5.2f]\n","protected", band2.first,band2.second); 

	TH1F* ht0 		= new TH1F("band0", "FC conf band", 2, 0, 2);
	TH1F* ht0_prot 	= new TH1F("band0_prot", "conf band", 2, 0, 2);
	TH1F* ht1 		= new TH1F("band1", "profile likelihood conf band", 2, 0, 2);
	ht0->SetBinContent(1,band1.first);
	ht0->SetBinContent(2,band1.second);
	ht0_prot->SetBinContent(1,band2.first);
	ht0_prot->SetBinContent(2,band2.second);
	ht1->SetBinContent(1,band3.first);
	ht1->SetBinContent(2,band3.second);
	

	grePvalueVsLambdaS->Write("pvalueVsLambdaS");
	grePvalueVsLambdaSprot->Write("pvalueVsLambdaSprot");
	grePvalueVsLambdaSPL->Write("pvalueVsLambdaSPL");
	greLLRVsLambdaS->Write("LLRVsPvalue");
	

	ht0->Write();
	ht0_prot->Write();

	TH1F* hinfo = new TH1F("hinfo", "running info", 5,0,5);
		hinfo->GetXaxis()->SetBinLabel(1,"running time");
		hinfo->SetBinContent(1,timer.RealTime());

		hinfo->GetXaxis()->SetBinLabel(2,"z-prime mass");
		hinfo->SetBinContent(2,inZpmass);
	hinfo->Write();
	grePvalueVsLambdaS->Draw("ape");	

//	GetLimitFromExpoFit(grePvalueVsLambdaS);


	return;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
double const drellYan(double* x, double* par){

	double norm = par[0];
	double alpha = par[1];
	double k = par[2];
	
	double retVal = exp(-alpha*TMath::Power(x[0],k));

	retVal *= norm;
	return retVal;
}

double const voigt(double* x, double* par) {

	double norm 	= par[0];
	double mean 	= par[1];
	double sigma	= par[2];
	double gamma	= par[3];

	double retval = TMath::Voigt(x[0]-mean,sigma,gamma);
	retval *= norm;
	return  retval;
}

double const GetMass(TF1* func, TH1F* reverseTable){

	Double_t min, max;
	func->GetRange(min,max);

	double p 	= rand3->Rndm();
	int bin 	= reverseTable->FindBin(p);
	double val 	= reverseTable->GetBinContent(bin);
	double delta = (max - min)/(1.*reverseTable->GetNbinsX());		
	double valnew = val-delta;
	if (valnew < min) valnew = min;
	
	double mod = 2.*rand3->Rndm()*delta;
	double mass = valnew+mod;
	return mass;
}
//
// gives you a vector of toy MC
//
std::vector<double> GetToyMC(double const lambdaS, double lambdaB){
	std::vector<double> retvec; retvec.clear(); 


	lambdaB *= _jitterLambdaB;
	
	int numSig = (int)rand3->Poisson(lambdaS);
	int numBkg = (int)rand3->Poisson(lambdaB);
	if (lambdaS ==0) numSig =0;

	FitCacheVec.clear();
	
//
// this was set up to use look-up tables and stuff
// I got lazy. this works fine
// the tables were yelling at me because of precesion. it still worked fine but I hate warnings
// i think it was running into a float-precesion issue in filling the table so the binning got mad
// totally stupid, right?
//
//fGenSignal, *fGenBackground;
	for (int i = 0; i < numSig; ++i) {
//		double val = fSignal->GetRandom();
		double val = fGenSignal->GetRandom();
		retvec.push_back(val);
		hsig->Fill(retvec.back());
	}
	for (int i = 0; i < numBkg; ++i) {
//		double val = fBackground->GetRandom();
		double val = fGenBackground->GetRandom();
		retvec.push_back(val);
		hbkg->Fill(val);
	}
	sort(retvec.begin(), retvec.end()); // sort the masses.  this is generally good to do for precesion calculations
//	std::vector<std::pair<double,double> > FitCacheVec;

	for (std::vector<double>::const_iterator it = retvec.begin(); it !=retvec.end(); ++it) {
		double a = fSignal->Eval(*it);
		double b = fBackground->Eval(*it);
		FitCacheVec.push_back(std::pair<double,double>(a,b) );
	}

	return retvec;
} 
//
// a median calculator
//
double const CalculateMedian(std::vector<double> vec){
	sort(vec.begin(), vec.end());
	return vec[vec.size()/2+1];

}
//
// This does all that good "fitting stuff"
// it works by giving it the toy mc points, the expected signal and background
// and an option on whether or not to "fix" the expected signal
// returns an <int,double> pair with the fit flag (in case it fails) and the log-likelihood ratio
//
std::pair<int,double> GetMinuitAnswer(std::vector<double> const toymc, double const lambdaS, double const lambdaB, bool const fixLambdaS){

	Double_t arglist[2];
	Int_t ierflg = 0;    //error flag that gets spit out after minimization 

	arglist[0] = 1;      //this likes to have an initial error set otherwise it breaks itself
	minuit->mnexcm("SET ERR", arglist ,1,ierflg); //more not breaking.
//	minuit->mnparm(0,"stepLambdaB", lambdaB, 24, 0.,	100	,ierflg);
	minuit->mnparm(0,"stepLambdaB", 1.*rand3->Poisson(lambdaB), 24, 0.,	150	,ierflg);

	if (fixLambdaS) minuit->mnparm(1,"stepLambdaS", lambdaS, 0		, 0.,	50	,ierflg);
	else minuit->mnparm(1,"stepLambdaS", 1.*rand3->Poisson(lambdaS), 5		, 0.,	50	,ierflg);
//	else minuit->mnparm(1,"stepLambdaS", lambdaS, 5		, 0.,	50	,ierflg);

	arglist[0] = 5000;		// This is the number of function calls that minuit will run before it commits suicide.
	arglist[1] = 100;		// I think this is iterations over the second derivative but i'm really not at all sure
	_MASSES = toymc;
	minuit->mnexcm("MIGRAD", arglist ,2,ierflg); 
	minuit->mnexcm("MIGRAD", arglist ,2,ierflg); 

	double valLambdaB, errLambdaB;
	double valLambdaS, errLambdaS;
	double loglikelihood =0;//, llr =0;
	minuit->GetParameter(0, valLambdaB, errLambdaB);
	minuit->GetParameter(1, valLambdaS, errLambdaS);
	loglikelihood = EvalPdf(valLambdaS, valLambdaB);
//	printf("lambdaB = %7.4f\n", valLambdaB);
//	printf("lambdaS = %7.4f\n", valLambdaS);
//	printf("-log(likelihood) = %10.5f\n", loglikelihood);


	return std::pair<int,double>(ierflg,loglikelihood);
}

double const CalculateMean(std::vector<double> vec){
	double retval = 0;
	for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it) retval+= *it;
	retval /= vec.size();
	return retval;
}
double const CalculateRMS(std::vector<double> vec, double const mean){
	double retval = 0;
	for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it) retval+= (*it-mean)*(*it-mean);
	retval /=vec.size();
	return sqrt(retval);	
}

//
// z pole fit
//
double ZPole(Double_t* x, Double_t* par){

double const PDG_MASS_Z 	= 91.1876; 	//GeV
double const PDG_WIDTH_Z	= 2.4952;	//GeV

	double A	= par[0];
	double mu 	= PDG_MASS_Z;
	double gamma= PDG_WIDTH_Z;

	double sigma=par[1];
	double mass = x[0];
	double theta = par[2];
	
//
// first part
// 
	double poleval = A*TMath::Voigt(mass-mu,sigma,gamma);//gamma,sigma);
	double pdfTerm = TMath::Exp(-1.*theta*mass);

//
// interference
//	
	double B = par[3];
	double mass2 = mass*mass;
	double mu2	= mu*mu;
	double deltam2 = mass2-mu2;

	double interference = B*mu*deltam2*deltam2/ ( deltam2*deltam2 + mass2*mass2*gamma*gamma/mu2   ); 
//(mass2-mu2)*(mass2-mu2)/((mass2-mu2)*(mass2-mu2)+mass2*mass2*gamma*gamma/mu2);

//
// final expo
//
	double C = par[4];
	double kappa = par[5];

	double outexpo = C/(mass*mass)*TMath::Exp(-1.*kappa*mass);
	return poleval*pdfTerm + interference*pdfTerm + outexpo;
}
//
// tail fit
//
double shitTail(double* x, double* par) {
	double mass = x[0];
	double A = par[0];
	double alpha = par[1];
	double kappa = par[2];

	return TMath::Exp(A+-1.*alpha*TMath::Power(mass,kappa));
}
//
// mix the two
//
double MixFunc(double* x, double* par) {

	double scale = par[10];

	double turnonWidth = par[9];

	double mass = x[0];
	double tailA = par[0];
	double tailAlpha = par[1];
	double tailKappa = par[2];

	double sigA = par[3];
	double sigSigma= par[4];
	double sigTheta = par[5];
	double sigB	= par[6];
	double sigC	= par[7];
	double sigKappa= par[8];


	double erfval = 0.5+0.5*TMath::Erf(turnonWidth*(mass-120));//   120+0.3*x[0]);

	double valSig[] = {sigA,sigSigma,sigTheta,sigB,sigC,sigKappa};
	double val1a = ZPole(x,valSig);

	double valTail[] = {tailA,tailAlpha,tailKappa};

//	double val1b = shitTail(x,par);
	double val1b = shitTail(x,valTail);

	double	val1 = val1a * (1.-erfval);
	double 	val2	= val1b*erfval;

	
	return scale*(val1+val2);

}
//
// put the tree into a vector. I like vectors. It comes back sorted in terms of mass
//
std::vector<double> const PutDataIntoVec(TTree* tree){
	float mass;
	std::vector<double> retVec; retVec.clear();

//
// oh wow. with this setup, you can easily just use the 2 4-vectors to read mass... just a thought...
// i strongly suggest that you put in a pruned tree before this 
// step else you'll be hacking this until the end of time. 
//
//

	tree->SetBranchAddress("mass", &mass);	
	int numEntries = tree->GetEntries();
	
	for (int jEntry = 0; jEntry < numEntries; ++jEntry){
		tree->GetEntry(jEntry);
		if (mass < 1*_lowerBound || mass > _upperBound) continue;
		retVec.push_back(mass);
	}

	std::cout<<"number of events in data: "<<retVec.size()<<std::endl;	
	std::sort(retVec.begin(),retVec.end());
	return retVec;
}

double Signal500(double* x, double* par) {
	
//	fit = new TF1("fit", "[0]*TMath::Voigt(x-[1],[2],[3])*TMath::Exp(-[4]*x)+[5]*TMath::Gaus(x,[1],[6])");
//	fit ->SetParNames("A", "mean", "sigma", "lg", "alpha", "B", "sigma2");

	double mass = x[0];

	double A		= par[0];
	double mean		= par[1];
	double sigma	= par[2];
	double lg		= par[3];
	double alpha	= par[4];
	double B		= par[5];
	double sigma2	= par[6];
	double norm		= par[7];	
	
	double retVal = A*TMath::Voigt(mass-mean,sigma,lg)*TMath::Exp(-1.*alpha*mass)+B*TMath::Gaus(mass,mean,sigma2);

	return norm*retVal;
}

//
// pdf calculation for minuit
//
void minuitPdf(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
//
// we're only fitting for 2 parameters. There's no reason to see anything else...
//
	double lambdaB = par[0];
	double lambdaS = par[1];
	f = EvalPdfCached(lambdaS,lambdaB);
}
//
// pdf evaluation with lambdaS and lambdaB
//
// this is totally not optimized.  
// the TF1 evaluations could be cached
// this happened in one incarnation that was way to advanced
// the one-file script is so much easier
//
double const EvalPdf(double const lambdaS, double const lambdaB){ 

	double likelihood = 0;
	double sum = lambdaB+lambdaS;
	if (_MASSES.empty()) return -1.*sum;
	for (std::vector<double>::const_iterator mass = _MASSES.begin(); mass != _MASSES.end(); ++mass){
		double a = lambdaS*fSignal->Eval(*mass);
		double b = lambdaB*fBackground->Eval(*mass);
		double val = a+b;
		likelihood += log(val);
	}
	likelihood -= sum;
	return -2.*likelihood;	
}
//
// general form for an exponential convoluted with a gaussian
// i.e. the form for the t-tbar background
//
double ExpoConvGaus(double* x, double* par) {

	double mass = x[0];
	double norm	 = par[0];
	double sigma	= par[1];
	double gamma	= par[2];
	double mean		= par[3];

	double expo = 0.5*(gamma*gamma*sigma*sigma) + gamma*(mean-mass);
		expo = TMath::Exp(expo);
	double erf = TMath::Erf( (gamma*sigma*sigma + mean-mass) /(sigma*sqrt(2.)));	
	return norm*expo*(1-erf);///(2*gamma);
}
//
// if you want to add multiple backgrounds, this works wonders
// need to know the fraction that is ttbar. 
//
double BckgFunc(double* x, double* par){
	double f1 = par[0];
	double fn = 1.-f1;
	double mass = x[0];
	
	double val1 = fdrellYan->Eval(mass);
	double val2 = fTtbar->Eval(mass);

	return val1*f1+fn*val2;
}
//
// set the signal parameters for different z-prime models 
//
void const SetSignalParameters(int const mass_) {

//	double A		= par[0];
//	double mean		= par[1];
//	double sigma	= par[2];
//	double lg		= par[3];
//	double alpha	= par[4];
//	double B		= par[5];
//	double sigma2	= par[6];
//	double norm		= par[7];	

//
// for default parameters, use the info from the z
//
		fSignal->SetParameters(
		9.45741e+06,
		(1.*mass_),
		2.63045,
		2.4952,		// use the z width
		3.91833e-03,
		6.42994e+02,
		3.03020e+01,
		1.);
	if (mass_ == 500) {
		fSignal->SetParameters(
		1.05021e+05,
		4.93904e+02,
		7.90580e+01,
		4.48431e+01,
		3.91833e-03,
		6.42994e+02,
		3.03020e+01,
		1.);
	}
	if (mass_ == 750) {	fSignal->SetParameters(
		9.45534e+04,
		7.39962e+02,
		1.31342e+02,
		9.05346e+01,
		2.72171e-03,
		3.08867e+02,
		5.75050e+01,
		1.);}
	if (mass_ == 1000) {	fSignal->SetParameters(

			8.00731e+04,			
			9.81335e+02	,		
			1.89187e-01	,		
			9.99006e-02	,		
			1.03518e-01	,		
			2.46877e+02	,		
			1.00291e+02 ,
			1);	}

	if (mass_ == 1250) {	fSignal->SetParameters(
			 2.00000e+02,			
			 1.20649e+03,			
			 2.00000e-01,			
			 1.00000e-01,			
			 2.04984e+00,			
			 1.60941e+02,			
			 1.47620e+02,			
			1);	}
	if (mass_ == 1500) {	fSignal->SetParameters(
			2.00000e+02,				
			1.42002e+03	,			
			2.00000e-01	,			
			1.00000e-01	,			
			8.53613e+00	,			
			1.28112e+02	,			
			2.04584e+02	,			
			1);	}
	if (mass_ == 1750) {	fSignal->SetParameters(
			 2.00000e+02,			
			 1.65944e+03,			
			 2.00000e-01,			
			 1.00000e-01,			
			 1.45729e+00,			
			 9.12729e+01,			
			 2.56125e+02,			
			1);	}



}

void const SetSignalParametersFixed(int const mass_){

		double width = GetZssmGamma(mass_);
		double sigma = GetResolution(mass_);

		fSignal->SetParameters(1, 1.*mass_,width,sigma); 
}

//
// get the limit estimate
// 

double const GetLimitFromExpoFit(TGraphErrors* gre){

	tempfit = new TF1("tempfit", "expo");
//	tempfit->SetRange(_startingExpectedSignal, 6.);
	tempfit->SetLineColor(kRed);
	tempfit->SetParameters(-1,-1);
	tempfit->SetParLimits(0,0,-1000);
	tempfit->SetParLimits(1,0,-1000);
//	gre->Fit("tempfit","rb");
	gre->Fit("tempfit","b");

	gre->Draw("ape");

	double limitval = CalcLimitFomExpFitVal(tempfit);
//	double limitval = (log(0.05)-tempfit->GetParameter(0))/tempfit->GetParameter(1);

//	std::cout<<"95\% limit from expo fit: "<<limitval<<std::endl;
	return limitval;
}
// 
// 
// 
double const CalcLimitFomExpFitVal(TF1* thefit){
	double limitval = (log(0.05)-thefit->GetParameter(0))/thefit->GetParameter(1);
	std::cout<<"95\% limit from expo fit: "<<limitval<<std::endl;
	return limitval;	

}

double GetZssmGamma(int const _mass){

	double p0 =	-3.98166    ; 
	double p1 =	0.0417735    ;
	double p2 =	0.000148462  ;
	double val = p0+_mass*p1+_mass*_mass*p2;
	return val;

}
double GetResolution(int const _mass){
	double p0 = 0.009031;
	double p1 = 9.89e-05;
	double p2 = -2.579e-08;
	double sigma = p0+_mass*p1+_mass*_mass*p2;
	return sigma;
}
//
// calculate the p-value of a test-statistic
//
double const CalculatePValue(std::vector<double> const& vec, double const t0){
	double retVal = 0;
	for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it){
		if (*it < t0) retVal+=1;
	}
	return (retVal/(1.*vec.size()));
}

std::pair<double,double> GetBand(TGraphErrors* gre, double const pval){

	int numPoints = gre->GetN();
//#	std::pair<double,double> retPair(0,999);

	std::vector<double> vals; vals.clear();
	
	double* yvals = gre->GetY();
	double* xvals = gre->GetX();
	
	for (int i = 0; i < numPoints; ++i){	
		if (yvals[i] > pval) vals.push_back(xvals[i]); 
	}
	return std::pair<double,double>(vals.front(),vals.back());

}

double const EvalPdfCached(double const lambdaS, double const lambdaB){ 

	if (FitCacheVec.empty() || !_useCaching) return EvalPdf(lambdaS,lambdaB);	


	double likelihood = 0;
	double sum = lambdaB+lambdaS;
	if (_MASSES.empty()) return -1.*sum;
	for (std::vector<std::pair<double,double> >::const_iterator mass = FitCacheVec.begin(); mass != FitCacheVec.end(); ++mass){

		double a = lambdaS*(mass->first);
		double b = lambdaB*(mass->second);
		double val = a+b;
		likelihood += log(val);
	}
	likelihood -= sum;
	return -2.*likelihood;	
}
double const GetPValError(double const pval, int const num){
	return sqrt(pval*(1.-pval)/double(num));
}

