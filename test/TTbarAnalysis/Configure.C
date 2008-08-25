#ifndef Configure_C
#define Configure_C
Float_t _MinJetPt = 20; // - - minimum jet pt to be studied
Float_t _MaxJetEta= 2.4; // - - mamimum jet eta to be studied 
Float_t _MinMuonPt = 20; // - - consistence cut with TeV Muon 2007 effort
Float_t _MaxMuonEta = 2.4 ; //--- in muon system's acceptance
Int_t _MinNumberOfJetsAboveThreshold = 0;// - - How many jets required to select events
Bool_t _ifRunOnSoup = true; // extract dijets
Bool_t _ifApplyIsolationCut = true; // for muon, traker based sumpt 
Float_t _MaxIso03sumPt = 10.;
Bool_t _ifApplyImpactParCut = false; // for muon, impact parameter
Float_t _MaxDxy = 0.1;
Bool_t _Debug = false; // control the print out


// - - - -analyzing
Bool_t _OppositeSign = true; // 
Bool_t _ifSkipTwoMuonsAssociatedWithTheSameJet = true;
Float_t _DeltaRCutForTwoMuons = 0.0;
Bool_t _SkipMuonsWhosePAboveThreshold = false;
Float_t _MaxMuonP = 6000.;
TString _InputFileToDoAnalyzing = "/data/1b/mschen/Chowder8_TxAtLeastOneMuonPtGT20GeV.root";
TString _OutputFile = "test.root";
Bool_t _ifTakeIntoAccountWeight_Analyzing = true;// evt weight

Float_t  _BTagDiscriminator = 0.7;
Float_t _BQuarkAcceptanceEta = 2.4;

Float_t _DiMuonMassCut = 0;

Float_t _ZMassVetoLow = 70.;
Float_t _ZMassVetoHigh = 110.;

#endif 

