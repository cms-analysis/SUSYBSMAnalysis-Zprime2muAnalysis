#include "FakeMuonGenerator.h"

FakeMuonGenerator::FakeMuonGenerator( int seed ):
	_fakeRateHist (0),
	_deltaEtaHist (0),
	_deltaPhiHist (0),
	_deltaPHist   (0),
	_deltaPtHist   (0),
	_deltaEtaArray(0),
	_deltaPhiArray(0),
	_deltaPArray  (0),
	_deltaPtArray(0)
{
	_gRandom.SetSeed ( seed );
}

FakeMuonGenerator::~FakeMuonGenerator(){

	_inFile -> Close();

	_fakeRateHist = 0;
	_deltaEtaHist = 0;
	_deltaPhiHist = 0;
	_deltaPHist   = 0;
	_deltaPtHist   = 0;
	if (_deltaEtaHist ){
		delete[] _deltaEtaHist;
		_deltaEtaHist = 0;
	}

	if (_deltaPhiHist){
		delete[] _deltaPhiHist;
		_deltaPhiHist = 0;
	}

	if (_deltaPHist){
		delete[] _deltaPHist;
		_deltaPHist = 0;
	}
	if (_deltaPtHist){
		delete[] _deltaPtHist;
		_deltaPtHist = 0;
	}
	delete _inFile;
}

void FakeMuonGenerator::SetupArray( float **ptr, TH2F *hist, int *nEntry, float &min, float &width){

	long int dim = hist -> GetXaxis() -> GetNbins() * hist -> GetYaxis() -> GetNbins();

	printf("temporary array dimension: %ld floats\n", dim);

	*ptr = new float [ dim ];

	float *array = *ptr;

	nEntry[0] = hist -> GetYaxis() -> GetNbins();
	nEntry[1] = hist -> GetXaxis() -> GetNbins();
	min       = hist -> GetYaxis() -> GetBinLowEdge(1);
	width     = hist -> GetYaxis() -> GetBinWidth(1);

	//  printf("nEntry = %d\n", nEntry);

	int ncol = nEntry[0];

	float norm = 0;

	for (int i = 0; i < hist -> GetXaxis() -> GetNbins() ; i++){

		norm = 0;

		for (int j = 0; j < hist -> GetYaxis() -> GetNbins(); j++){

			norm += hist -> GetBinContent(i+1, j+1);
			array[i*ncol + j] = norm;

		}

		if (norm == 0.0)
			norm = 1.0;
		else
			norm = 1.0 / norm;

		for (int j = 0; j < hist -> GetYaxis() -> GetNbins(); j++){

			array[i*ncol + j] *= norm;

		}

	}

}

float FakeMuonGenerator::ThrowFromArray( float *array, int bin, int nPerBin, double min, double width ){

	//  printf("bin: %d nPerBin: %d, min: %lf, width: %lf\n", bin, nPerBin, min, width);

	int firstBin = (bin-1)*nPerBin;

	double rand = _gRandom . Rndm();

	int ibin;

	for ( ibin = 0; (rand >= array[firstBin+ibin]) && (ibin < nPerBin); ibin++ );

	if( ibin == nPerBin )
		ibin = nPerBin-1;

	rand = _gRandom . Rndm() * width;

	double smear = min + ibin*width + rand;

	return smear;
}

void FakeMuonGenerator::Initialize( const char *fileName ){

	_inFile = new TFile(fileName);
//FIXME
//-------those hists names need to be configurable 
//add check  also for file name
	TH2F *tempHist = (TH2F*) _inFile -> Get("allJetPEta");
	tempHist->Rebin2D(20,20); 
	_fakeRateHist = (TH2F*) _inFile -> Get("matchedJetPEta");
	_fakeRateHist -> Rebin2D(20,20);
	_fakeRateHist -> Divide( tempHist );

	_deltaEtaHist = (TH2F*) _inFile -> Get("muonEtaVsJetEta");
	_deltaPhiHist = (TH2F*) _inFile -> Get("muonPhiVsJetPhi");
	_deltaPHist   = (TH2F*) _inFile -> Get("muonPVsJetP");
	_deltaPtHist   = (TH2F*) _inFile -> Get("muonPtVsJetPt");

	_deltaPHist -> Rebin2D(10,10);
	_deltaPtHist -> Rebin2D(10,10);	

	printf("setting up delta Eta array...\n");

	SetupArray( &_deltaEtaArray, _deltaEtaHist, _nEtaPerRow, _etaMin, _etaWidth );

	printf("setting up delta Phi array...\n");

	SetupArray( &_deltaPhiArray, _deltaPhiHist, _nPhiPerRow, _phiMin, _phiWidth );

	printf("setting up delta P array...\n");

	SetupArray( &_deltaPArray,   _deltaPHist,   _nDelPPerRow, _pMin, _pWidth);

	printf("setting up delta Pt array...\n");

	SetupArray( &_deltaPtArray,   _deltaPtHist,   _nDelPtPerRow, _ptMin, _ptWidth);
}

int FakeMuonGenerator::GenerateFakeMuon ( double p, double pt, double eta, double phi, FakeMuonGenerator::fakeMuon  &muon ){

	//  printf("p: %f\n", p);
	//  printf("eta: %f\n", eta);
	//  printf("phi: %f\n", phi);

	double borderP = p-1e-6;

	int retVal = 1;

	double fakeRate=0;


	double random;


	int ibin;

	// first, for a given p and eta, calculate the fake rate..

	ibin = _fakeRateHist -> FindBin( p, eta );

	fakeRate = _fakeRateHist -> GetBinContent(ibin);




	// next, throw a momentum for our new muon..
	ibin = _deltaPHist -> GetXaxis() -> FindBin( borderP );

	double fakeP=0;
	double fakeEta=0;
	double fakePhi=0;
	double fakePt=0;
	int    charge = 1;


	double borderEta = eta-1e-6;
	double borderPhi = phi-1e-6;
	double borderPt = pt - 1e-6;
	int  ietabin = _deltaEtaHist->GetXaxis()->FindBin(borderEta);
	int  iphibin = _deltaPhiHist->GetXaxis()->FindBin(borderPhi);
	int  iptbin = _deltaPtHist->GetXaxis()->FindBin(borderPt);
	//  printf("ibin: %d\n", ibin );

	fakeP = -1;

	//at some slice, the projection only has muons with P<2.5, then it will stuck here .. loop  

	int n=0;
	while (fakeP < 2.5 && n<10000 ){

		fakeP = ThrowFromArray( _deltaPArray, ibin, _nDelPPerRow[0], 
				_pMin, _pWidth);
		n++;
	}

	fakePt = ThrowFromArray( _deltaPtArray, iptbin, _nDelPtPerRow[0],
			_ptMin, _ptWidth);


	// throw an eta for our new muon..
	fakeEta = ThrowFromArray( _deltaEtaArray, ietabin, _nEtaPerRow[0], 
			_etaMin, _etaWidth);

	// throw a phi for our new muon..
	fakePhi = ThrowFromArray( _deltaPhiArray, iphibin, _nPhiPerRow[0], 
			_phiMin, _phiWidth);

	// throw a charge for our new muon...
	random = _gRandom . Rndm();
	if (random > 0.5) charge = -charge;

	// finally, fill the muon...

	muon.weight = fakeRate; 
	muon.charge = charge;
	muon.p = fakeP;
	muon.eta = fakeEta;
	muon.phi = fakePhi;
	muon.pt=fakePt;

	//  printf("generated charge, weight, p, eta, phi: %d %f %f %f %f\n",
	// charge, fakeRate, fakeP, fakeEta, fakePhi);
	//	}
	//	printf("after gen: muons.size= %d\n",muons.size());
	return retVal;
}
