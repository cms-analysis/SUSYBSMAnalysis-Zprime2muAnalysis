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

#include "TGraphErrors.h"

//bool const _useData = false;

double const _lowerBound = 50;
double const _upperBound = 3000;
int const _numPseudo = 1000;


int _theZPrimeMass = 1750;
TRandom3* rand3;

static double const drellYan(double* x, double* par);
static double const voigt(double* x, double* par);

TF1* fdrellYan, *fSignal, *fmix, *fBackground;
TF1* fTtbar;
//double const avgBackground = 23.67;	//60 pb-1
//double avgBackground = 19.5241; //50 pb-1

TH1F* hbkg, *hsig, *hsnb;

double const GetMass(TF1* func, TH1F* reverseTable);
TH1F* GenerateIntegralTableNew (TF1* func, TString const sname, double const res);

TH1F* h;
//TH1F* hBackgroundTable, *hSignalTable;

//std::vector<double> GetToyMC(double const lambdaS, double const lambdaB, TH1F* hBkgTable, TH1F* hSigTable);// old version 
std::vector<double> GetToyMC(double const lambdaS, double const lambdaB); 
//int DoFit(bool const fixLambdaS, double const lambdaB, double const lambdaS);


TGraphErrors* grePvalueVsLambdaS;

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
double Signal1000(double* x, double* par); // there used to be a different template

double ExpoConvGaus(double* x, double* par);

double BckgFunc(double* x, double* par);



void limits(int inZpmass = 1750) {
	_theZPrimeMass = inZpmass;
	rand3 = new TRandom3(0);
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
*/

	
	fmix = new TF1("fmix",MixFunc, _lowerBound, _upperBound, 11);
	fmix->SetRange(_lowerBound,_upperBound);
	fmix->FixParameter(0 	,	2.44298e+01	);
	fmix->FixParameter(1 	, 	6.99626e+00	); 
	fmix->FixParameter(2 	, 	2.03121e-01	); 
	fmix->FixParameter(3 	,	9.45741e+06	);  
	fmix->FixParameter(4 	, 	2.63045e+00	); 
	fmix->FixParameter(5 	, 	1.80690e-02	); 
	fmix->FixParameter(6 	,	5.12014e+01	);  
	fmix->FixParameter(7 	,	0.00000e+00	);  
	fmix->FixParameter(8 	, 	0.00000e+00	); 
	fmix->FixParameter(9	,  	2.19753e-02	);
	fmix->SetParameter(10,1);

	fmix->SetNpx(10000);

	fmix->SetParameter(10, 1./fmix->Integral(_lowerBound, _upperBound));
	fmix->Draw();
//
//	hBackgroundTable->Draw();

	fTtbar = new TF1("ttbar", ExpoConvGaus, _lowerBound, _upperBound, 4);  
	fTtbar->SetRange(_lowerBound,_upperBound);
	fTtbar->SetNpx(10000);
	fTtbar->SetParameters(1, 2.50938e+01, 1.73358e-02, 	4.80605e+01);


	fTtbar->SetParameter(0, 1./fTtbar->Integral(_lowerBound, _upperBound));
	fTtbar->Draw();
//	return;


//	RooConstVar ttbarSigma	("ttbarSigma"	, "ttbarGamma"	,2.50938e+01);//,20, 200);
//	RooConstVar ttbarGamma	("ttbarGamma"	, "ttbarSigma"	,1.73358e-02);//,20, 200);
//	RooConstVar ttbarMean	("ttbarMean"	, "ttbarMean"	,4.80605e+01);//,20, 200);

//double ExpoConvGaus(double* x, double* par) {


	fdrellYan 	= fmix;
//	fBackground	= fmix;


//	double BckgFunc(double* x, double* par){

	fBackground = new TF1("fBackground",BckgFunc, _lowerBound, _upperBound,1); 
	fBackground->FixParameter(0,0.94);

	fBackground->Draw();
	
//	return;
//
// get the signal tables
//
	fSignal = new TF1("fSignal", Signal500, _lowerBound, _upperBound,8);
		fSignal->SetRange(_lowerBound,_upperBound);
		fSignal->SetParNames("A", "mean", "sigma", "lg", "alpha", "B", "sigma2", "Norm");
//
// this is the default setting-- i.e. for a mass of 500 GeV
//
		fSignal->SetParameters(
		1.05021e+05,
		4.93904e+02,
		7.90580e+01,
		4.48431e+01,
		3.91833e-03,
		6.42994e+02,
		3.03020e+01,
		1.);

//
// lots of cases
//
	if (_theZPrimeMass == 750) {	fSignal->SetParameters(
		9.45534e+04,
		7.39962e+02,
		1.31342e+02,
		9.05346e+01,
		2.72171e-03,
		3.08867e+02,
		5.75050e+01,
		1.);}
	if (_theZPrimeMass == 1000) {	fSignal->SetParameters(

			8.00731e+04,			
			9.81335e+02	,		
			1.89187e-01	,		
			9.99006e-02	,		
			1.03518e-01	,		
			2.46877e+02	,		
			1.00291e+02 ,
			1);	}

	if (_theZPrimeMass == 1250) {	fSignal->SetParameters(
			 2.00000e+02,			
			 1.20649e+03,			
			 2.00000e-01,			
			 1.00000e-01,			
			 2.04984e+00,			
			 1.60941e+02,			
			 1.47620e+02,			
			1);	}
	if (_theZPrimeMass == 1500) {	fSignal->SetParameters(
			2.00000e+02,				
			1.42002e+03	,			
			2.00000e-01	,			
			1.00000e-01	,			
			8.53613e+00	,			
			1.28112e+02	,			
			2.04584e+02	,			
			1);	}
	if (_theZPrimeMass == 1750) {	fSignal->SetParameters(
			 2.00000e+02,			
			 1.65944e+03,			
			 2.00000e-01,			
			 1.00000e-01,			
			 1.45729e+00,			
			 9.12729e+01,			
			 2.56125e+02,			
			1);	}


	fSignal->SetParameter("Norm",1./fSignal->Integral(_lowerBound,_upperBound));
	fSignal->SetNpx(10000);
	fSignal->Draw();

//	return;
//	hSignalTable		= GenerateIntegralTableNew(fSignal, "hSignalTable", 2.);
//	hSignalTable->Draw();
//
//	
	hbkg = new TH1F("hbkg", "DY background generation", 1000,_lowerBound,_upperBound );
	hsig = new TH1F("hsig", "Zprime signal generation", 1000,_lowerBound,_upperBound );


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
//		sDataName = "dimuons_20100825.root";

		TFile inData(sDataName, "open");
		std::cout<<"Getting data file from: "<<sDataName<<std::endl;	

//	
// the the tree with masses.  roofit does it the same way. snap.
// 
//
	
		TTree* inTree = (TTree*)inData.Get("tree");
		std::vector<double> dataMasses; dataMasses.clear();
// vectors are better... and the tree needs to be pruned. 
// the ntupler should be fixed to not keep like-sign stuff  
// we don't need to worry about the stats now
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
		
			for (lambdaS = 2.3; lambdaS < 5.5; lambdaS += 0.1){
	
				std::pair<int,double> answerPairfix = GetMinuitAnswer(toymc, lambdaS, lambdaB,true);
				if (answerPairfix.first ==0) llfix = answerPairfix.second;
				double llr = llfit - llfix;
				
				printf("lambdaS = %3.2f\tfit = %5f\tfix = %5f\tllr = %5f\n",
					lambdaS, llfit,llfix,llr);
					std::vector<double> tmp; tmp.clear(); tmp.push_back(llr);
					refTable.insert(std::pair<double,std::vector<double> > (lambdaS, tmp));	
					

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

//
// I decided to do it where I calculate the reference values and then iterate on the values that were calculated. 
// I could have calculated a reference for each point, 
// but that would have put more shit in shit and I wasn't in the mood to debug that.
//
	for (std::map<double,std::vector<double> >::const_iterator tableit = refTable.begin(), 
		tableitend = refTable.end(); tableit!=tableitend; ++tableit) {
	

	
		lambdaS = tableit->first;
		std::vector<double> vec = tableit->second;		

		double refVal(0);

		std::cout<<"lambdaS = "<<tableit->first
			<<"\t test-val = "<<vec.front();

//		numPseudo = 1000;
		refVal = vec.front();

		std::vector<double> llrVals; llrVals.clear();

		for (int jexp = 0; jexp < _numPseudo; ++jexp) {
			std::vector<double> toymc=GetToyMC(lambdaS, lambdaB);//, hBackgroundTable, hSignalTable);
			std::pair<int,double> answerFit = GetMinuitAnswer(toymc, lambdaS, lambdaB,false);
			if (answerFit.first != 0) continue;
			std::pair<int,double> answerFix = GetMinuitAnswer(toymc, lambdaS, lambdaB,true);
			if (answerFix.first != 0) continue;
				
			double llr = answerFit.second - answerFix.second;
			llrVals.push_back(llr);
		}
	
		sort(llrVals.begin(),llrVals.end());
	

//		double meanLlr 		= CalculateMean(llrVals); 
//		double medianLlr 	= CalculateMedian(llrVals); 
//		double rmsLlr		= CalculateRMS(llrVals,meanLlr);		

//		std::cout
//			<<"\t median = "<<medianLlr
//			<<"\t mean = "<<meanLlr
//			<<"\t rms = "<<rmsLlr;		

		double pval= 0;

		for (std::vector<double>::const_iterator it = llrVals.begin(); it != llrVals.end(); ++it){

			if (*it< refVal) pval+=1;

		}

//		std::cout<<"\tp-value = "<<pval;
		pval /= 1.*_numPseudo;
		std::cout<<"\tp-value = "<<pval;
		std::cout<<std::endl;
	
		double pvalerr = sqrt(pval*(1.-pval)/numPseudo);
		grePvalueVsLambdaS->SetPoint(grePvalueVsLambdaS->GetN(), lambdaS, pval);
		grePvalueVsLambdaS->SetPointError(grePvalueVsLambdaS->GetN()-1, 0.1, pvalerr);

//		break;
	}


	_outFile = new TFile(TString::Format("out_%d_1p1.root",_theZPrimeMass), "recreate");
	_outFile->cd();
	grePvalueVsLambdaS->Write("pvalueVsLambdaS");
	hbkg->Write();
	hsig->Write();


	return;
}

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
TH1F* GenerateIntegralTableNew (TF1* func, TString const sname, double const res){

	double min, max;
	func->GetRange(min,max);

	typedef std::pair<double,double> valInt;

	std::vector<valInt> vecd;			
	vecd.push_back(valInt(0,min));

	double intval = 0;

	double val= min;
	while (intval<0.9999 && val< max){
		intval = func->Integral(min, val);
//		intval = std::max(intval,vecd.back().second);
		vecd.push_back(valInt(intval,val));
		val += res;	
	}
	
	vecd.push_back(valInt(1,max));
	double bkgbins[vecd.size()];

	for (unsigned int i = 0; i < vecd.size(); ++i) bkgbins[i] = vecd[i].first;

	TH1F* massTable= new TH1F(sname, "masstable", vecd.size()-1, bkgbins);

	for (unsigned int i = 0; i < vecd.size(); ++i) massTable->SetBinContent(i, vecd[i].second);
	return massTable;
}
//
// gives you a vector of toy MC
// deprecated
//
/*
std::vector<double> GetToyMC(double const lambdaS, double const lambdaB, TH1F* hBkgTable, TH1F* hSigTable){
	std::vector<double> retvec; retvec.clear(); 
	
	int numSig = (int)rand3->PoissonD(lambdaS);
	int numBkg = (int)rand3->PoissonD(lambdaB);
	if (lambdaS ==0) numSig =0;

//
// this was set up to use look-up tables and stuff
// I got lazy. this works fine
// the tables were yelling at me because of precesion. it still worked fine but I hate warnings
// i think it was running into a float-precesion issue in filling the table so the binning got mad
// totally stupid, right?
//

	for (int i = 0; i < numSig; ++i) {
		double val = fSignal->GetRandom();
		retvec.push_back(val);
	}
	for (int i = 0; i < numBkg; ++i) {
		double val = fBackground->GetRandom();
		retvec.push_back(val);
	}
	sort(retvec.begin(), retvec.end());
	return retvec;
} 
*/
//
// gives you a vector of toy MC
//
std::vector<double> GetToyMC(double const lambdaS, double const lambdaB){
	std::vector<double> retvec; retvec.clear(); 
	
	int numSig = (int)rand3->PoissonD(lambdaS);
	int numBkg = (int)rand3->PoissonD(lambdaB);
	if (lambdaS ==0) numSig =0;

//
// this was set up to use look-up tables and stuff
// I got lazy. this works fine
// the tables were yelling at me because of precesion. it still worked fine but I hate warnings
// i think it was running into a float-precesion issue in filling the table so the binning got mad
// totally stupid, right?
//
	for (int i = 0; i < numSig; ++i) {
		double val = fSignal->GetRandom();
		retvec.push_back(val);
//		hsig->Fill(retvec.back());
	}
	for (int i = 0; i < numBkg; ++i) {
		double val = fBackground->GetRandom();
		retvec.push_back(val);
//		hbkg->Fill(val);
	}
	sort(retvec.begin(), retvec.end()); // sort the masses.  this is generally good to do for precesion calculations
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
	double poleval = A*TMath::Voigt(mass-mu,gamma,sigma);
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
//
//
std::vector<double> const PutDataIntoVec(TTree* tree){
	bool const _checkCharge = true;
//
// fast way to read in...
// 
	typedef struct {
		int charge;
		float pt;
		float eta;
		float phi;
	} _TrackInfo;

	float mass;
	_TrackInfo mu1, mu2;

	std::vector<double> retVec; retVec.clear();



//
// oh wow. with this setup, you can easily just use the 2 4-vectors to read mass... just a thought
//
	tree->SetBranchAddress("recoCandMass", &mass);	
	tree->SetBranchAddress("reco1", &mu1);	
	tree->SetBranchAddress("reco2", &mu2);	
	int numEntries = tree->GetEntries();
	
	for (int jEntry = 0; jEntry < numEntries; ++jEntry){
		tree->GetEntry(jEntry);
		if (_checkCharge && mu1.charge == mu2.charge) continue;

		// cut on the masses we don't want
		if (mass < _lowerBound || mass > _upperBound) continue;
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

double Signal1000(double* x, double* par) {

	double mass = x[0];
		
	double norm = par[0];
	double turnonMass = par[1];
	double turnonWidth = par[2];

	double erfval = 0.5+0.5*TMath::Erf(turnonWidth*(mass-turnonMass));//   120+0.3*x[0]);


//	fit = new TF1("fit", "[0]*TMath::Voigt(x-[1],[2],[3])*TMath::Exp(-[4]*x)+[5]*TMath::Gaus(x,[1],[6])");

	double A 	= par[3];
	double mean 	= par[4];
	double sigma1 	= par[5];
	double lz	 	= par[6];
	double alpha	= par[7];
	double B		= par[8];
	double sigma2	= par[9];


	double sideRight = A*TMath::Voigt(mass-mean,sigma1,lz)*TMath::Exp(-1.*alpha*mass)+B*TMath::Gaus(mass,mean,sigma2);
	double right1 = (erfval)*sideRight;

	double C = par[10];
	double d1 = par[11];
	double e1 = par[12];
	
	double sideLeft	= TMath::Exp(C+d1*mass+e1*mass*mass); 
	double left1 = (1.-erfval)*sideLeft;
//	fit2 = new TF1("fit2", "TMath::Exp([0]+[1]*x+[2]*x*x)");

	return norm*(right1+left1);

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
	f = EvalPdf(lambdaS,lambdaB);
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

