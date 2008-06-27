#ifndef __FAKEMUONGENERATOR_H
#define __FAKEMUONGENERATOR_H (1)
/*
	Author: Ivan K. Furic
*/
#include "TRandom.h"
#include "TFile.h"
#include "TH2F.h"
#include <vector>
class FakeMuonGenerator{

	public:

		// structure of fake muon
		struct fakeMuon{
			double weight;
			int charge;
			double p;
			double eta;
			double phi;
			double pt;
		};



		FakeMuonGenerator( int seed = 13 );
		~FakeMuonGenerator();

		// Initialize , get histograms and extract data arrays 
		void Initialize( const char *fileName = 0 );

		// Main function to produce a fake muon
		int  GenerateFakeMuon ( double p, double pt, double eta, double phi, fakeMuon& muon );

		// Throw muon properties according to the calibrations
		float ThrowFromArray( float *array, int bin, int nPerBin, 
				double min, double width );

	protected:

		void SetupArray( float **ptr, TH2F *hist, int *nEntry,
				float &min, float &width);

		TRandom _gRandom;

		TFile *_inFile;

		TH2F  *_fakeRateHist;
		TH2F  *_deltaEtaHist;
		TH2F  *_deltaPhiHist;
		TH2F  *_deltaPHist;
		TH2F  *_deltaPtHist;

		float   *_deltaEtaArray;
		float   *_deltaPhiArray;
		float   *_deltaPArray;
		float   *_deltaPtArray;

		int   _nEtaPerRow[2];
		int   _nPhiPerRow[2];
		int   _nDelPPerRow[2];
		int   _nDelPtPerRow[2];

		float _etaWidth;
		float _phiWidth;
		float _pWidth;
		float _ptWidth;

		float _etaMin;
		float _phiMin;
		float _pMin;
		float _ptMin;

};

#endif
