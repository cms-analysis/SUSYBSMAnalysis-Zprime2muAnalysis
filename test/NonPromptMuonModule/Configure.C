#ifndef Configure_C
#define Configure_C
//---common parameters for parameterizing and analyzing
Bool_t Use_Likelihood = 0; // if use likehood ratio to discriminate
Double_t _MinJetPt = 20; // - - minimum jet pt to produce a fake muon
Double_t _MaxJetEta= 2.5; // - - mamimum jet eta to produce a fake muon , havn't used it yet
Double_t _MinMuonPt = 20; // - - consistence cut with TeV Muon 2007 effort
Double_t _MaxMuonEta = 2.4 ; //--- in muon system's acceptance
Int_t _MinNumberOfJetsAboveThreshold = 2;// - - How many jets required to selected events
Bool_t _ifRunOnGumboSoup = true; // extract dijets
//Bool_t _ifTakeIntoAccountWeight = false;// evt weight
Int_t _TakeHalfEvents = 0; //0-all; 1- odd to para, even to analyze; 2- odd to analyze, even to para 
Bool_t _ifApplyLikelihoodRatioCut = false; //
Double_t _KeepPromptMuonEfficiency = 0.99; //
Bool_t _ifApplyIsolationCut = false; // for muon, traker based sumpt 
Double_t _MaxIso03sumPt = 20;
Bool_t _ifApplyImpactParCut = false; // for muon, impact parameter
Double_t _MaxDxy = 0.1;
Bool_t _Debug = false; // control the print out
TString _InputFileLikelihoodCalibrated = "MyLLRTrain.root"; // the LLR calibration from W+jets and dijets 

Double_t MaxJetPInRefHist_  = 6000. ; 
Double_t MaxJetEtaInRefHist_ = 5.5 ;
Double_t MinJetEtaInRefHist_ = -5.5 ;

Double_t _MinPtHatToSkip = 30;// not used yet
Double_t _MaxPtHatToSkip = 50;// not used yet

// - - - - parameterizing
TString _InputFileToDoParameterizing = "/data/0b/mschen/CSA07/FullSim100pb/Gumbo_Test.root";
TString _OutputFileParameterizedHistos = "DijetsParameterizing.root";
Bool_t _ifTakeIntoAccountWeight_Parameterizing = false;// evt weight

// - - - -analyzing
Bool_t _OppositeSign = true; // 
Bool_t _ifSkipTwoMuonsAssociatedWithTheSameJet = true;
Double_t _DeltaRCutForTwoMuons = 0.0;
Bool_t _SkipMuonsWhosePAboveThreshold = true;
Double_t _MaxMuonP = 6000.;
TString _InputFileToDoAnalyzing = "/data/0b/mschen/CSA07/FullSim100pb/Gumbo_Test.root";
TString _InputFileParameterizedHistos = "DijetsParameterizing.root";
Int_t _FakeMuonGeneratorSeed = 1000;
TString _OutputFileComparison = "DijetsPredObsv.root";
Bool_t _ifTakeIntoAccountWeight_Analyzing = false;// evt weight


//-- CSA07 QCD Dijets weights ~100/pb
// - - - QCD Dijets pt-hats can be a discriminant
//- - /data/1b/mschen/Gumbo_Test_Dijets.root

	double DiJetsWeight[21] = {
		378571 , 
		87687.7,
		23706.5,
		6415.43,
		862.138,
		8.7069 ,
		260.847,
		39.52  ,
		2.04167,
		0.14958,
		0.528814, 
		0.0351   ,
		0.0555285,
		0.0408   ,
		1.08e-06 ,
		4.29e-05 ,
		0.000238 ,
		0.00533333,
		0.000725  ,
		8.44e-06  ,
		0.0363333 
	};
	TString s_pthat[21] = {
		"0~15", 
		"15~20", 
		"20~30", 
		"30~50", 
		"50~80", 
		"170~230",
		"80~120",
		"120~170",
		"230~300",
		"380~470",
		"300~380", 
		"800~1000", 
		"470~600",  
		"600~800", 
		"3500~inf", 
		"2600~3000",
		"2200~2600",
		"1400~1800",
		"1800~2200",
		"3000~3500",
		"1000~1400" 
	};

#endif 

