//unit :  dxy(cm)  , qoverp(1/GeV)

/* class NtupleDumper
 *
 * Example of Making Ntuple From CMSSW_RootFiles
 *
 * \author Mingshui Chen
 *EMAIL:  mschen@cern.ch        
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "FWCore/Framework/interface/TriggerNames.h"


class TFile;
class TH1D;

#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH3D.h>
#include "TFile.h"
#include "TText.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include <TLegend.h>//c
#include <Math/VectorUtil.h>
#include <TVector3.h>
#include <utility>

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include <vector>
#include <TTree.h>

class NtupleDumper : public edm::EDAnalyzer {
	public:
		explicit NtupleDumper( const edm::ParameterSet & );
		~NtupleDumper();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void beginJob(const edm::EventSetup& ) ;
		virtual void endJob() ;
		float DeltaPhi(float v1, float v2);
		float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
		TLorentzVector sumMuMu(const reco::Muon &muon1, const reco::Muon &muon2 );

		edm::ESHandle<ParticleDataTable> pdt_;

		bool applyNormalization_;


		bool debug_;


		bool runOnAODSIM_;

		bool runOnCSA07Soup_;
		// try to get the processId

		bool hasHLTTriggerResults_;
		bool hasL1TriggerResults_;

		bool initTriggerNamesService_;

		int minNumMuons_;
	
		bool _DoYouWannaGetJetTag;
		edm::InputTag _JetTagsLabel;


		edm::TriggerNames triggerNames_; 


		double angle ( double x1, double y1, double z1, double x2, double y2, double z2);
		double xyAngle ( double x1, double y1, double x2, double y2) ;

		bool verbose;
		double overallLumi;

		void _fillJet(const edm::Handle<reco::CaloJetCollection> caloJets) ;
		void _fillMet(const edm::Handle<reco::CaloMETCollection> caloMet);
		void _fillGmr(const edm::Handle<reco::MuonCollection> muons);
		void _fillPyt(const edm::Handle<reco::CandidateCollection>  genParticles);
		void _fillL1Trig(const edm::Handle<l1extra::L1ParticleMapCollection> L1PMC);
		void _fillHLTrig(const edm::Handle<edm::TriggerResults> hltTR);
		void _fillEvtInfo(const edm::Event& event);
		void _trigInit();
		void _jetInit();
		void _gmrInit();
		void _pytInit();
		void _metInit();
	private:
		std::string outputFileName_;
		edm::InputTag hepMC_, tracks_, caloTowers_, caloJets_, genJets_, 
			tauTags_, electrons_, muons_, photons_, standAlones_;  
		bool skipGenJets_;

		edm::InputTag csa07EventWeightTag_;
		double _OveralLumiForSoup;
		std::string _CSA07WeightModuleLabel;


		TFile * outputFile_ ;
		bool l1MuonEmulator_;
		bool takeGMTfromL1Extra_;

		TH1F *h_numEvents;

		TTree * myTree;

		// - - - -for Event
		int evtRun;
		int evtNum;
		float evtWeight;
		int evtProcessId;
		float genEventScale;
		float genEventWeight;
		int genEventProcID;
		float genEventFilterEff;

		// - - - -for Trigger
		vector<bool> *hltResults;
		vector<bool> *l1Results;

		// - -- MET
		vector<float> *metEta	    ;
		vector<float> *metPhi	    ;
		vector<float> *metPx	    ;
		vector<float> *metPy	    ;
		vector<float> *metPz	    ;
		vector<float> *metEt	    ;
		vector<float> *metSumEt	    ;

		// - - - -JET
		int jetSize;
		vector<float> *jetPt	    ;
		vector<float> *jetEta	    ;
		vector<float> *jetPhi	    ;
		vector<float> *jetP	    ;
		vector<float> *jetPx	    ;
		vector<float> *jetPy	    ;
		vector<float> *jetPz	    ;
		vector<float> *jetTheta	    ;
		vector<float> *jetEnergy    ;
		vector<float> *jetEt	    ;
		vector<float> *jetTagDiscriminator;

		// - - -  -PYTHIA
		int             pytSize;
		vector <float> *pytEnergy         ;
		vector <float> *pytPt	          ;
		vector <float> *pytEta	          ;
		vector <float> *pytPhi	          ;
		vector <float> *pytP	          ;
		vector <float> *pytPx	          ;
		vector <float> *pytPy	          ;
		vector <float> *pytPz	          ;
		vector <float> *pytTheta          ;
		vector <float> *pytVx	          ;
		vector <float> *pytVy	          ;
		vector <float> *pytVz	          ;
		vector <float> *pytCharge         ;
		vector <float> *pytMass	          ;
		vector <int> *pytIndex          ;
		vector <int> *pytNMother        ;
		vector <int> *pytMotherIndex    ;
		vector <int> *pytMotherPdgId;
		vector <int> *pytNDaughter        ;
		vector <int> *pytDaughterIndex    ;
		vector <int> *pytDaughterPdgId;
		vector <int> *pytId	          ;
		vector <int> *pytStatus         ;

		// - -  - -Deault Global Muon Collection
		int gmrSize ;
		vector<float> *gmrEnergy      ;
		vector<float> *gmrPt	      ;
		vector<float> *gmrEta	      ;
		vector<float> *gmrPhi	      ;
		vector<float> *gmrP	      ;
		vector<float> *gmrPx	      ;
		vector<float> *gmrPy	      ;
		vector<float> *gmrPz	      ;
		vector<float> *gmrTheta	      ;
		vector<float> *gmrVx	      ;
		vector<float> *gmrVy	      ;
		vector<float> *gmrVz	      ;
		vector<float> *gmrCharge      ;
		vector<float> *gmrNDoF	      ;
		vector<float> *gmrChi2	      ;
		vector<float> *gmrChi2Norm    ;
		vector<float> *gmrDXY	      ;
		vector<float> *gmrDTheta      ;
		vector<float> *gmrDPt	      ;
		vector<float> *gmrDEta	      ;
		vector<float> *gmrDPhi	      ;
		vector<float> *gmrDDXY	      ;
		vector<float> *gmrIso03nTracks;
		vector<float> *gmrIso03sumPt  ;

		vector<float> *gmrDz	 ;
		vector<float> *gmrD0	 ;
		vector<float> *gmrDsz	 ;
		vector<float> *gmrDDz	 ;
		vector<float> *gmrDD0	 ;
		vector<float> *gmrDDsz	 ;
		vector<float> *gmrInnerX ;
		vector<float> *gmrInnerY ;
		vector<float> *gmrInnerZ ;

		vector<float> *gmrNumberOfValidHits;
		vector<float> *gmrNumberOfLostHits;
		vector<float> *gmrTrackerHits;


};

