//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun  2 23:38:08 2008 by ROOT version 5.14/00f
// from TTree myTree/myTree
// found on file: /data/0b/mschen/CSA07/FullSim100pb/Chowder_Test_TxAtLeastOneMuonPtGT20GeV.root
//////////////////////////////////////////////////////////

#ifndef NtupleBase_h
#define NtupleBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

class NtupleBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Float_t         evtWeight;
   Float_t         genEventWeight;
   Float_t         genEventScale;
   Int_t           genEventProcID;
   Float_t         evtProcessId;
   Int_t           evtRun;
   Int_t           evtNum;
   vector<bool>    *hltResults;
   vector<bool>    *l1Results;
   Int_t           gmrSize;
   vector<float>   *gmrPt;
   vector<float>   *gmrEta;
   vector<float>   *gmrPhi;
   vector<float>   *gmrP;
   vector<float>   *gmrPx;
   vector<float>   *gmrPy;
   vector<float>   *gmrPz;
   vector<float>   *gmrTheta;
   vector<float>   *gmrVx;
   vector<float>   *gmrVy;
   vector<float>   *gmrVz;
   vector<float>   *gmrCharge;
   vector<float>   *gmrNDoF;
   vector<float>   *gmrChi2;
   vector<float>   *gmrChi2Norm;
   vector<float>   *gmrDXY;
   vector<float>   *gmrDTheta;
   vector<float>   *gmrDPt;
   vector<float>   *gmrDEta;
   vector<float>   *gmrDPhi;
   vector<float>   *gmrDDXY;
   vector<float>   *gmrIso03nTracks;
   vector<float>   *gmrIso03sumPt;
   vector<float>   *gmrDz;
   vector<float>   *gmrD0;
   vector<float>   *gmrDsz;
   vector<float>   *gmrDDz;
   vector<float>   *gmrDD0;
   vector<float>   *gmrDDsz;
   vector<float>   *gmrInnerX;
   vector<float>   *gmrInnerY;
   vector<float>   *gmrInnerZ;
   vector<float>   *gmrNumberOfValidHits;
   vector<float>   *gmrNumberOfLostHits;
   vector<float>   *gmrTrackerHits;
   Int_t           pytSize;
   vector<float>   *pytEnergy;
   vector<float>   *pytPt;
   vector<float>   *pytEta;
   vector<float>   *pytPhi;
   vector<float>   *pytP;
   vector<float>   *pytPx;
   vector<float>   *pytPy;
   vector<float>   *pytPz;
   vector<float>   *pytTheta;
   vector<float>   *pytVx;
   vector<float>   *pytVy;
   vector<float>   *pytVz;
   vector<float>   *pytCharge;
   vector<float>   *pytMass;
   vector<int>     *pytIndex;
   vector<int>     *pytNMother;
   vector<int>     *pytMotherIndex;
   vector<int>     *pytMotherPdgId;
   vector<int>     *pytId;
   vector<int>     *pytStatus;
   vector<int>     *pytNDaughter;
   vector<int>     *pytDaughterIndex;
   vector<int>     *pytDaughterPdgId;
   Int_t           jetSize;
   vector<float>   *jetPt;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetP;
   vector<float>   *jetPx;
   vector<float>   *jetPy;
   vector<float>   *jetPz;
   vector<float>   *jetTheta;
   vector<float>   *jetEnergy;
   vector<float>   *jetEt;
   vector<float>   *metEta;
   vector<float>   *metPhi;
   vector<float>   *metPx;
   vector<float>   *metPy;
   vector<float>   *metPz;
   vector<float>   *metSumEt;
   vector<float>   *metEt;

   // List of branches
   TBranch        *b_evtWeight;   //!
   TBranch        *b_genEventWeight;   //!
   TBranch        *b_genEventScale;   //!
   TBranch        *b_genEventProcID;   //!
   TBranch        *b_evtProcessId;   //!
   TBranch        *b_evtRun;   //!
   TBranch        *b_evtNum;   //!
   TBranch        *b_hltResults;   //!
   TBranch        *b_l1Results;   //!
   TBranch        *b_gmrSize;   //!
   TBranch        *b_gmrPt;   //!
   TBranch        *b_gmrEta;   //!
   TBranch        *b_gmrPhi;   //!
   TBranch        *b_gmrP;   //!
   TBranch        *b_gmrPx;   //!
   TBranch        *b_gmrPy;   //!
   TBranch        *b_gmrPz;   //!
   TBranch        *b_gmrTheta;   //!
   TBranch        *b_gmrVx;   //!
   TBranch        *b_gmrVy;   //!
   TBranch        *b_gmrVz;   //!
   TBranch        *b_gmrCharge;   //!
   TBranch        *b_gmrNDoF;   //!
   TBranch        *b_gmrChi2;   //!
   TBranch        *b_gmrChi2Norm;   //!
   TBranch        *b_gmrDXY;   //!
   TBranch        *b_gmrDTheta;   //!
   TBranch        *b_gmrDPt;   //!
   TBranch        *b_gmrDEta;   //!
   TBranch        *b_gmrDPhi;   //!
   TBranch        *b_gmrDDXY;   //!
   TBranch        *b_gmrIso03nTracks;   //!
   TBranch        *b_gmrIso03sumPt;   //!
   TBranch        *b_gmrDz;   //!
   TBranch        *b_gmrD0;   //!
   TBranch        *b_gmrDsz;   //!
   TBranch        *b_gmrDDz;   //!
   TBranch        *b_gmrDD0;   //!
   TBranch        *b_gmrDDsz;   //!
   TBranch        *b_gmrInnerX;   //!
   TBranch        *b_gmrInnerY;   //!
   TBranch        *b_gmrInnerZ;   //!
   TBranch        *b_gmrNumberOfValidHits;   //!
   TBranch        *b_gmrNumberOfLostHits;   //!
   TBranch        *b_gmrTrackerHits;   //!
   TBranch        *b_pytSize;   //!
   TBranch        *b_pytEnergy;   //!
   TBranch        *b_pytPt;   //!
   TBranch        *b_pytEta;   //!
   TBranch        *b_pytPhi;   //!
   TBranch        *b_pytP;   //!
   TBranch        *b_pytPx;   //!
   TBranch        *b_pytPy;   //!
   TBranch        *b_pytPz;   //!
   TBranch        *b_pytTheta;   //!
   TBranch        *b_pytVx;   //!
   TBranch        *b_pytVy;   //!
   TBranch        *b_pytVz;   //!
   TBranch        *b_pytCharge;   //!
   TBranch        *b_pytMass;   //!
   TBranch        *b_pytIndex;   //!
   TBranch        *b_pytNMother;   //!
   TBranch        *b_pytMotherIndex;   //!
   TBranch        *b_pytMotherPdgId;   //!
   TBranch        *b_pytId;   //!
   TBranch        *b_pytStatus;   //!
   TBranch        *b_pytNDaughter;   //!
   TBranch        *b_pytDaughterIndex;   //!
   TBranch        *b_pytDaughterPdgId;   //!
   TBranch        *b_jetSize;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetP;   //!
   TBranch        *b_jetPx;   //!
   TBranch        *b_jetPy;   //!
   TBranch        *b_jetPz;   //!
   TBranch        *b_jetTheta;   //!
   TBranch        *b_jetEnergy;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_metEta;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_metPz;   //!
   TBranch        *b_metSumEt;   //!
   TBranch        *b_metEt;   //!

   NtupleBase(TTree *tree=0);
//	NtupleBase() {};
   virtual ~NtupleBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(){};
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupleBase_cxx
NtupleBase::NtupleBase(TTree *tree)
{

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/0b/mschen/CSA07/FullSim100pb/Chowder_Test_TxAtLeastOneMuonPtGT20GeV.root");
      if (!f) {
         f = new TFile("/data/0b/mschen/CSA07/FullSim100pb/Chowder_Test_TxAtLeastOneMuonPtGT20GeV.root");
      }
      tree = (TTree*)gDirectory->Get("myTree");

   }
   Init(tree);
}
NtupleBase::~NtupleBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NtupleBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtupleBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NtupleBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   hltResults = 0;
   l1Results = 0;
   gmrPt = 0;
   gmrEta = 0;
   gmrPhi = 0;
   gmrP = 0;
   gmrPx = 0;
   gmrPy = 0;
   gmrPz = 0;
   gmrTheta = 0;
   gmrVx = 0;
   gmrVy = 0;
   gmrVz = 0;
   gmrCharge = 0;
   gmrNDoF = 0;
   gmrChi2 = 0;
   gmrChi2Norm = 0;
   gmrDXY = 0;
   gmrDTheta = 0;
   gmrDPt = 0;
   gmrDEta = 0;
   gmrDPhi = 0;
   gmrDDXY = 0;
   gmrIso03nTracks = 0;
   gmrIso03sumPt = 0;
   gmrDz = 0;
   gmrD0 = 0;
   gmrDsz = 0;
   gmrDDz = 0;
   gmrDD0 = 0;
   gmrDDsz = 0;
   gmrInnerX = 0;
   gmrInnerY = 0;
   gmrInnerZ = 0;
   gmrNumberOfValidHits = 0;
   gmrNumberOfLostHits = 0;
   gmrTrackerHits = 0;
   pytEnergy = 0;
   pytPt = 0;
   pytEta = 0;
   pytPhi = 0;
   pytP = 0;
   pytPx = 0;
   pytPy = 0;
   pytPz = 0;
   pytTheta = 0;
   pytVx = 0;
   pytVy = 0;
   pytVz = 0;
   pytCharge = 0;
   pytMass = 0;
   pytIndex = 0;
   pytNMother = 0;
   pytMotherIndex = 0;
   pytMotherPdgId = 0;
   pytId = 0;
   pytStatus = 0;
   pytNDaughter = 0;
   pytDaughterIndex = 0;
   pytDaughterPdgId = 0;
   jetPt = 0;
   jetEta = 0;
   jetPhi = 0;
   jetP = 0;
   jetPx = 0;
   jetPy = 0;
   jetPz = 0;
   jetTheta = 0;
   jetEnergy = 0;
   jetEt = 0;
   metEta = 0;
   metPhi = 0;
   metPx = 0;
   metPy = 0;
   metPz = 0;
   metSumEt = 0;
   metEt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
   fChain->SetBranchAddress("genEventWeight", &genEventWeight, &b_genEventWeight);
   fChain->SetBranchAddress("genEventScale", &genEventScale, &b_genEventScale);
   fChain->SetBranchAddress("genEventProcID", &genEventProcID, &b_genEventProcID);
   fChain->SetBranchAddress("evtProcessId", &evtProcessId, &b_evtProcessId);
   fChain->SetBranchAddress("evtRun", &evtRun, &b_evtRun);
   fChain->SetBranchAddress("evtNum", &evtNum, &b_evtNum);
   fChain->SetBranchAddress("hltResults", &hltResults, &b_hltResults);
   fChain->SetBranchAddress("l1Results", &l1Results, &b_l1Results);
   fChain->SetBranchAddress("gmrSize", &gmrSize, &b_gmrSize);
   fChain->SetBranchAddress("gmrPt", &gmrPt, &b_gmrPt);
   fChain->SetBranchAddress("gmrEta", &gmrEta, &b_gmrEta);
   fChain->SetBranchAddress("gmrPhi", &gmrPhi, &b_gmrPhi);
   fChain->SetBranchAddress("gmrP", &gmrP, &b_gmrP);
   fChain->SetBranchAddress("gmrPx", &gmrPx, &b_gmrPx);
   fChain->SetBranchAddress("gmrPy", &gmrPy, &b_gmrPy);
   fChain->SetBranchAddress("gmrPz", &gmrPz, &b_gmrPz);
   fChain->SetBranchAddress("gmrTheta", &gmrTheta, &b_gmrTheta);
   fChain->SetBranchAddress("gmrVx", &gmrVx, &b_gmrVx);
   fChain->SetBranchAddress("gmrVy", &gmrVy, &b_gmrVy);
   fChain->SetBranchAddress("gmrVz", &gmrVz, &b_gmrVz);
   fChain->SetBranchAddress("gmrCharge", &gmrCharge, &b_gmrCharge);
   fChain->SetBranchAddress("gmrNDoF", &gmrNDoF, &b_gmrNDoF);
   fChain->SetBranchAddress("gmrChi2", &gmrChi2, &b_gmrChi2);
   fChain->SetBranchAddress("gmrChi2Norm", &gmrChi2Norm, &b_gmrChi2Norm);
   fChain->SetBranchAddress("gmrDXY", &gmrDXY, &b_gmrDXY);
   fChain->SetBranchAddress("gmrDTheta", &gmrDTheta, &b_gmrDTheta);
   fChain->SetBranchAddress("gmrDPt", &gmrDPt, &b_gmrDPt);
   fChain->SetBranchAddress("gmrDEta", &gmrDEta, &b_gmrDEta);
   fChain->SetBranchAddress("gmrDPhi", &gmrDPhi, &b_gmrDPhi);
   fChain->SetBranchAddress("gmrDDXY", &gmrDDXY, &b_gmrDDXY);
   fChain->SetBranchAddress("gmrIso03nTracks", &gmrIso03nTracks, &b_gmrIso03nTracks);
   fChain->SetBranchAddress("gmrIso03sumPt", &gmrIso03sumPt, &b_gmrIso03sumPt);
   fChain->SetBranchAddress("gmrDz", &gmrDz, &b_gmrDz);
   fChain->SetBranchAddress("gmrD0", &gmrD0, &b_gmrD0);
   fChain->SetBranchAddress("gmrDsz", &gmrDsz, &b_gmrDsz);
   fChain->SetBranchAddress("gmrDDz", &gmrDDz, &b_gmrDDz);
   fChain->SetBranchAddress("gmrDD0", &gmrDD0, &b_gmrDD0);
   fChain->SetBranchAddress("gmrDDsz", &gmrDDsz, &b_gmrDDsz);
   fChain->SetBranchAddress("gmrInnerX", &gmrInnerX, &b_gmrInnerX);
   fChain->SetBranchAddress("gmrInnerY", &gmrInnerY, &b_gmrInnerY);
   fChain->SetBranchAddress("gmrInnerZ", &gmrInnerZ, &b_gmrInnerZ);
   fChain->SetBranchAddress("gmrNumberOfValidHits", &gmrNumberOfValidHits, &b_gmrNumberOfValidHits);
   fChain->SetBranchAddress("gmrNumberOfLostHits", &gmrNumberOfLostHits, &b_gmrNumberOfLostHits);
   fChain->SetBranchAddress("gmrTrackerHits", &gmrTrackerHits, &b_gmrTrackerHits);
   fChain->SetBranchAddress("pytSize", &pytSize, &b_pytSize);
   fChain->SetBranchAddress("pytEnergy", &pytEnergy, &b_pytEnergy);
   fChain->SetBranchAddress("pytPt", &pytPt, &b_pytPt);
   fChain->SetBranchAddress("pytEta", &pytEta, &b_pytEta);
   fChain->SetBranchAddress("pytPhi", &pytPhi, &b_pytPhi);
   fChain->SetBranchAddress("pytP", &pytP, &b_pytP);
   fChain->SetBranchAddress("pytPx", &pytPx, &b_pytPx);
   fChain->SetBranchAddress("pytPy", &pytPy, &b_pytPy);
   fChain->SetBranchAddress("pytPz", &pytPz, &b_pytPz);
   fChain->SetBranchAddress("pytTheta", &pytTheta, &b_pytTheta);
   fChain->SetBranchAddress("pytVx", &pytVx, &b_pytVx);
   fChain->SetBranchAddress("pytVy", &pytVy, &b_pytVy);
   fChain->SetBranchAddress("pytVz", &pytVz, &b_pytVz);
   fChain->SetBranchAddress("pytCharge", &pytCharge, &b_pytCharge);
   fChain->SetBranchAddress("pytMass", &pytMass, &b_pytMass);
   fChain->SetBranchAddress("pytIndex", &pytIndex, &b_pytIndex);
   fChain->SetBranchAddress("pytNMother", &pytNMother, &b_pytNMother);
   fChain->SetBranchAddress("pytMotherIndex", &pytMotherIndex, &b_pytMotherIndex);
   fChain->SetBranchAddress("pytMotherPdgId", &pytMotherPdgId, &b_pytMotherPdgId);
   fChain->SetBranchAddress("pytId", &pytId, &b_pytId);
   fChain->SetBranchAddress("pytStatus", &pytStatus, &b_pytStatus);
   fChain->SetBranchAddress("pytNDaughter", &pytNDaughter, &b_pytNDaughter);
   fChain->SetBranchAddress("pytDaughterIndex", &pytDaughterIndex, &b_pytDaughterIndex);
   fChain->SetBranchAddress("pytDaughterPdgId", &pytDaughterPdgId, &b_pytDaughterPdgId);
   fChain->SetBranchAddress("jetSize", &jetSize, &b_jetSize);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetP", &jetP, &b_jetP);
   fChain->SetBranchAddress("jetPx", &jetPx, &b_jetPx);
   fChain->SetBranchAddress("jetPy", &jetPy, &b_jetPy);
   fChain->SetBranchAddress("jetPz", &jetPz, &b_jetPz);
   fChain->SetBranchAddress("jetTheta", &jetTheta, &b_jetTheta);
   fChain->SetBranchAddress("jetEnergy", &jetEnergy, &b_jetEnergy);
   fChain->SetBranchAddress("jetEt", &jetEt, &b_jetEt);
   fChain->SetBranchAddress("metEta", &metEta, &b_metEta);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("metPz", &metPz, &b_metPz);
   fChain->SetBranchAddress("metSumEt", &metSumEt, &b_metSumEt);
   fChain->SetBranchAddress("metEt", &metEt, &b_metEt);
   Notify();
}

Bool_t NtupleBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupleBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupleBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NtupleBase_cxx
