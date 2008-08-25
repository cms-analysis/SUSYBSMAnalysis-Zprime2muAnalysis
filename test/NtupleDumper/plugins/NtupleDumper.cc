#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TLorentzVector.h"
//#include "Math.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include <algorithm>



#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/EDProduct.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertError.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"

#include "NtupleDumper.h"


#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "PhysicsTools/HepMCCandAlgos/interface/CSA07ProcessId.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

using namespace std;
using namespace edm;
using namespace reco;

const double MUMASS = 0.10566 ;

int num2muons = 0;
int num2muonsWithSameSign = 0;
int num2gmts = 0;
int num2gmtsWithSameSign = 0;


NtupleDumper::NtupleDumper( const ParameterSet & cfg ) :
	outputFileName_( cfg.getUntrackedParameter<string>("histoutputFile") ),
	hepMC_      ( cfg.getParameter<InputTag>( "hepMC" ) ),
	tracks_     ( cfg.getParameter<InputTag>( "tracks" ) ),
	caloTowers_ ( cfg.getParameter<InputTag>( "caloTowers" ) ),
	caloJets_   ( cfg.getParameter<InputTag>( "caloJets" ) ),
	genJets_    ( cfg.getParameter<InputTag>( "genJets" ) ),
	tauTags_    ( cfg.getParameter<InputTag>( "tauTags" ) ),
	electrons_  ( cfg.getParameter<InputTag>( "electrons" ) ),
	muons_      ( cfg.getParameter<InputTag>( "muons" ) ),
	standAlones_      ( cfg.getParameter<InputTag>( "standAlones" ) ),
	skipGenJets_( cfg.getUntrackedParameter<bool>( "skipGenJets", false ) 
		    ) {


		photons_ = cfg.getParameter<InputTag>("photons");

		l1MuonEmulator_ = cfg.getUntrackedParameter<bool>("L1MuonEmulator",false);
		takeGMTfromL1Extra_ = cfg.getUntrackedParameter<bool>("takeGMTfromL1Extra",false);

		debug_= cfg.getUntrackedParameter<bool>("debug",false);

		runOnAODSIM_ = cfg.getUntrackedParameter<bool>( "runOnAODSIM", false ) ;

		hasHLTTriggerResults_=cfg.getUntrackedParameter<bool>("hasHLTTriggerResults",false);

		hasL1TriggerResults_=cfg.getUntrackedParameter<bool>("hasL1TriggerResults",false);

		runOnCSA07Soup_ = cfg.getUntrackedParameter<bool>("runOnCSA07Soup",false);

		initTriggerNamesService_=false;

		if(runOnCSA07Soup_ ) 
		{
			csa07EventWeightTag_ = cfg.getParameter<InputTag>("csa07EventWeightTag");
			_OveralLumiForSoup = cfg.getUntrackedParameter<double>("OveralLumiForSoup",100.);
			_CSA07WeightModuleLabel = cfg.getUntrackedParameter<string>("CSA07WeightModuleLabel");
		}

		//	skim events with more then (minNumMuons_ -1 ) muons
		minNumMuons_ = cfg.getUntrackedParameter<int>("minNumMuons",0);

		//- jet tag
		_DoYouWannaGetJetTag = cfg.getUntrackedParameter<bool>("DoYouWannaGetJetTag", false);
		if( _DoYouWannaGetJetTag ) _JetTagsLabel = cfg.getParameter<InputTag>("JetTagsLabel");

	}

void NtupleDumper::analyze(const Event& event, const EventSetup& eventSetup) {


	/// get the particle data table
	ESHandle<ParticleDataTable> pdt;
	eventSetup.getData( pdt );

	// using namespace edm;
	cout << ">>> processing event # " << event.id() <<", time: " << event.time().value()
		<< endl;


	if(debug_==true) std::cout<<"Befour get the genParticles Collection"<<std::endl;


	edm::Handle<reco::CandidateCollection> genParticles;
	event.getByLabel("genParticleCandidates", genParticles);

	/*
	   if(debug_==true && printDecayTree_) 
	   {
	   particleDecayDrawer( genParticles,  eventSetup );
	   particleListDrawer( genParticles,  eventSetup );
	   }
	 */
	if(debug_==true)  std::cout<<"After get the genParticles collection"<<std::endl;
	// Muons
	Handle<MuonCollection> muons;
	event.getByLabel( muons_, muons );
	cout << "====> number of Muons      " << muons->size() << endl;


	if(debug_==true)   std::cout<<"After get the muon collection"<<std::endl;

	// Photons
	Handle<reco::PhotonCollection> photons;
	event.getByLabel(photons_, photons);

	Handle<TrackCollection> tracks;
	event.getByLabel( tracks_, tracks );
	Handle<TrackCollection> standAlones;
	event.getByLabel( standAlones_, standAlones );

	//---calo MET
	edm::Handle<reco::CaloMETCollection> caloMet;
	event.getByLabel("met", caloMet);
	//recoCaloMETs_met_ //only one

	//---------Jet, there are several kinds of jets reconstruction, I'm not sure which one is default
	edm::Handle<reco::CaloJetCollection> caloJets;
	event.getByLabel(caloJets_, caloJets);

	edm::Handle<edm::TriggerResults> hltTR;
	if(hasHLTTriggerResults_) {
		event.getByLabel(edm::InputTag ("TriggerResults::HLT"),hltTR);
		//     try {event.getByLabel(hltTriggerInputTag_,hltTR);} 
		//     catch( cms::Exception& ex ) { LogWarning("NtupleDumper") << "Trigger results: " << triggerInputTag_ << " not found"; }
		//     if (!hltTR.isValid())
		//       throw cms::Exception("ProductNotValid") << "TriggerResults product not valid";
	}

	edm::Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
	edm::Handle<l1extra::L1ParticleMapCollection> L1PMC;
	if(hasL1TriggerResults_) {
		//    event.getByLabel(l1ParticleMapTag_,L1PMC);
		event.getByLabel("l1extraParticleMap",L1PMC);
		//    event.getByLabel(l1GTReadoutRecTag_,L1GTRR);
		event.getByLabel("l1extraParticleMap",L1GTRR);
	}
	//----------Chowder soup weight -----------

	
	_fillEvtInfo(event);




	edm::Handle<reco::JetTagCollection> bTagHandle;
	if(_DoYouWannaGetJetTag) 
	{
		// Get b tag information
		event.getByLabel(_JetTagsLabel, bTagHandle);
	}

	//fill Ntuple---------

	if((int)muons->size()>=minNumMuons_) 
	{

		_metInit();
		_jetInit();
		_pytInit();
		_gmrInit();
		_trigInit();

		_fillMet(caloMet);
		_fillJet(caloJets);
		_fillPyt(genParticles);
		_fillGmr(muons);

		//-------fill trigger results
		if (hasHLTTriggerResults_==true)
		{
			_fillHLTrig(hltTR);
		}

		//if(debug_==true) std::cout<<"test3"<<endl;

		if (hasL1TriggerResults_==true)
		{
			_fillL1Trig(L1PMC);
		}

		// if(debug_==true) std::cout<<"test4"<<endl;
		//----------don't forget fill the tree for every event

		// -- FIXME need a separate function to do this
		if(_DoYouWannaGetJetTag)
		{
			const reco::JetTagCollection & bTags = *(bTagHandle.product());
			// Loop over jets and study b tag info.
			for (int i = 0; i != (int)bTags.size(); ++i) {
				jetTagDiscriminator->push_back(bTags[i].discriminator());
				//jetTagEta->push_back(bTags[i].jet()->eta());
				//jetTagPhi->push_back(bTags[i].jet()->phi());
			}
		}
	}
	myTree->Fill();


}//end Analyze() main function







TLorentzVector NtupleDumper::sumMuMu(const reco::Muon &muon1, const reco::Muon &muon2 ){
	TLorentzVector vmu1, vmu2, vmm; // LorentzVectors for mu+ mu- 
	vmu1.SetPtEtaPhiM(muon1.pt(), muon1.eta(), muon1.phi(), MUMASS);
	vmu2.SetPtEtaPhiM(muon2.pt(), muon2.eta(), muon2.phi(), MUMASS);
	vmm = vmu1 + vmu2;		
	return vmm;
}






float NtupleDumper::DeltaPhi(float v1, float v2)
{ // Computes the correctly normalized phi difference
	// v1, v2 = phi of object 1 aNd 2
	float diff = fabs(v2 - v1);
	float corr = 2*acos(-1.) - diff;
	if (diff < acos(-1.)){ return diff;} else { return corr;} 
}
//------------------------------------------------------------------------------

float NtupleDumper::GetDeltaR(float eta1, float eta2, float phi1, float phi2)
{ // Computes the DeltaR of two objects from their eta and phi values

	return sqrt( (eta1-eta2)*(eta1-eta2) 
			+ DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );

}

double NtupleDumper::angle ( double x1, double y1, double z1, double x2, double y2, double z2) {
	//  * returns three-dimensional Angle between two objects;
	//  * defined via scalar product: 
	//  *   angle = acos((v1 * v2)/(|v1| * |v2|))
	return acos((x1*x2 + y1*y2 + z1*z2)/sqrt((x1*x1 + y1*y1 + z1*z1)*(x2*x2 + y2*y2 + z2*z2)));
}

double NtupleDumper::xyAngle ( double x1, double y1, double x2, double y2) {
	//  * returns two-dimensional Angle between two objects; //in tranvers plane
	//  * defined via scalar product: 
	//  *   angle = acos((v1 * v2)/(|v1| * |v2|))
	return acos((x1*x2 + y1*y2)/sqrt((x1*x1 + y1*y1)*(x2*x2 + y2*y2)));
}


//------------------------------------------------------------------------------
NtupleDumper::~NtupleDumper() {
	delete outputFile_;
}
//------------------------------------------------------------------------------
void NtupleDumper::beginJob( const EventSetup & ) {
	outputFile_  = new TFile( outputFileName_.c_str(), "RECREATE" ) ;

	myTree = new TTree("myTree", "myTree");

	myTree->Branch("evtWeight",&evtWeight,"evtWeight/F");
	myTree->Branch("genEventWeight",&genEventWeight,"genEventWeight/F");
	myTree->Branch("genEventScale",&genEventScale,"genEventScale/F");
	myTree->Branch("genEventFilterEff",&genEventFilterEff,"genEventFilterEff/F");
	myTree->Branch("genEventProcID",&genEventProcID,"genEventProcID/I");
	myTree->Branch("evtProcessId",&evtProcessId,"evtProcessId/I");
	myTree->Branch("evtRun",&evtRun,"evtRun/I");
	myTree->Branch("evtNum",&evtNum,"evtNum/I");

	myTree->Branch("hltResults",&hltResults);
	myTree->Branch("l1Results",&l1Results);

	myTree->Branch("gmrSize",&gmrSize,"gmrSize/I");
	myTree->Branch("gmrPt",&gmrPt);
	myTree->Branch("gmrEta",&gmrEta);
	myTree->Branch("gmrPhi",&gmrPhi);
	myTree->Branch("gmrP",&gmrP);
	myTree->Branch("gmrPx",&gmrPx);
	myTree->Branch("gmrPy",&gmrPy);
	myTree->Branch("gmrPz",&gmrPz);
	myTree->Branch("gmrTheta",&gmrTheta);
	myTree->Branch("gmrVx",&gmrVx);
	myTree->Branch("gmrVy",&gmrVy);
	myTree->Branch("gmrVz",&gmrVz);
	myTree->Branch("gmrCharge",&gmrCharge);
	myTree->Branch("gmrNDoF",&gmrNDoF);
	myTree->Branch("gmrChi2",&gmrChi2);
	myTree->Branch("gmrChi2Norm",&gmrChi2Norm);
	myTree->Branch("gmrDXY",&gmrDXY);
	myTree->Branch("gmrDTheta",&gmrDTheta);
	myTree->Branch("gmrDPt",&gmrDPt);
	myTree->Branch("gmrDEta",&gmrDEta);
	myTree->Branch("gmrDPhi",&gmrDPhi);
	myTree->Branch("gmrDDXY",&gmrDDXY);
	myTree->Branch("gmrIso03nTracks",&gmrIso03nTracks);
	myTree->Branch("gmrIso03sumPt",&gmrIso03sumPt);

	myTree->Branch("gmrDz",&gmrDz);
	myTree->Branch("gmrD0",&gmrD0);
	myTree->Branch("gmrDsz",&gmrDsz);
	myTree->Branch("gmrDDz",&gmrDDz);
	myTree->Branch("gmrDD0",&gmrDD0);
	myTree->Branch("gmrDDsz",&gmrDDsz);
	myTree->Branch("gmrInnerX",&gmrInnerX);
	myTree->Branch("gmrInnerY",&gmrInnerY);
	myTree->Branch("gmrInnerZ",&gmrInnerZ);

	myTree->Branch("gmrNumberOfValidHits",&gmrNumberOfValidHits);
	myTree->Branch("gmrNumberOfLostHits",&gmrNumberOfLostHits);
	myTree->Branch("gmrTrackerHits",&gmrTrackerHits);




	myTree->Branch("pytSize",&pytSize,"pytSize/I");
	myTree->Branch("pytEnergy",&pytEnergy);
	myTree->Branch("pytPt",&pytPt);
	myTree->Branch("pytEta",&pytEta);
	myTree->Branch("pytPhi",&pytPhi);
	myTree->Branch("pytP",&pytP);
	myTree->Branch("pytPx",&pytPx);
	myTree->Branch("pytPy",&pytPy);
	myTree->Branch("pytPz",&pytPz);
	myTree->Branch("pytTheta",&pytTheta);
	myTree->Branch("pytVx",&pytVx);
	myTree->Branch("pytVy",&pytVy);
	myTree->Branch("pytVz",&pytVz);
	myTree->Branch("pytCharge",&pytCharge);
	myTree->Branch("pytMass",&pytMass);
	myTree->Branch("pytIndex",&pytIndex);
	myTree->Branch("pytNMother",&pytNMother);
	myTree->Branch("pytMotherIndex",&pytMotherIndex);
	myTree->Branch("pytMotherPdgId",&pytMotherPdgId);
	myTree->Branch("pytId",&pytId);
	myTree->Branch("pytStatus",&pytStatus);
	myTree->Branch("pytNDaughter",&pytNDaughter);
	myTree->Branch("pytDaughterIndex",&pytDaughterIndex);
	myTree->Branch("pytDaughterPdgId",&pytDaughterPdgId);


	myTree->Branch("jetSize",&jetSize,"jetSize/I");
	myTree->Branch("jetPt",&jetPt);
	myTree->Branch("jetEta",&jetEta);
	myTree->Branch("jetPhi",&jetPhi);
	myTree->Branch("jetP",&jetP);
	myTree->Branch("jetPx",&jetPx);
	myTree->Branch("jetPy",&jetPy);
	myTree->Branch("jetPz",&jetPz);
	myTree->Branch("jetTheta",&jetTheta);
	myTree->Branch("jetEnergy",&jetEnergy);
	myTree->Branch("jetEt",&jetEt);

	myTree->Branch("jetTagDiscriminator",&jetTagDiscriminator);

	myTree->Branch("metEta",&metEta);
	myTree->Branch("metPhi",&metPhi);
	myTree->Branch("metPx",&metPx);
	myTree->Branch("metPy",&metPy);
	myTree->Branch("metPz",&metPz);
	myTree->Branch("metSumEt",&metSumEt);
	myTree->Branch("metEt",&metEt);
	//makeNtuple


}

void NtupleDumper::endJob() {


	outputFile_->Write() ;
	outputFile_->Close() ;

}


void NtupleDumper::_metInit(){
	metEta	  =new vector<float>  ;
	metPhi	  =new vector<float>  ;
	metPx	  =new vector<float>  ;
	metPy	  =new vector<float>  ;
	metPz	  =new vector<float>  ;
	metEt	  =new vector<float>  ;
	metSumEt	  =new vector<float>  ;
}

void NtupleDumper::_jetInit(){
	jetPt	  =new vector<float>  ;
	jetEta	  =new vector<float>  ;
	jetPhi	  =new vector<float>  ;
	jetP	  =new vector<float>  ;
	jetPx	  =new vector<float>  ;
	jetPy	  =new vector<float>  ;
	jetPz	  =new vector<float>  ;
	jetTheta	  =new vector<float>  ;
	jetEnergy  =new vector<float>  ;
	jetEt	  =new vector<float>  ;
	jetTagDiscriminator = new vector<float>;
}
void NtupleDumper::_pytInit(){
	pytEnergy     =new vector<float>;  
	pytPt	      =new vector<float>;  
	pytEta	      =new vector<float>;  
	pytPhi	      =new vector<float>;  
	pytP	      =new vector<float>;  
	pytPx	      =new vector<float>;  
	pytPy	      =new vector<float>;  
	pytPz	      =new vector<float>;  
	pytTheta      =new vector<float>;  
	pytVx	      =new vector<float>;  
	pytVy	      =new vector<float>;  
	pytVz	      =new vector<float>;  
	pytCharge     =new vector<float>;  
	pytMass	      =new vector<float>;  
	pytIndex      =new vector<int>;    
	pytNMother    =new vector<int>;    
	pytMotherIndex=new vector<int>;    
	pytMotherPdgId=new vector<int>;   
	pytNDaughter    =new vector<int>;    
	pytDaughterIndex=new vector<int>;    
	pytDaughterPdgId=new vector<int>;   
	pytId	      =new vector<int>;  
	pytStatus     =new vector<int>;    
}
void NtupleDumper::_gmrInit(){
	gmrEnergy     =new vector<float>;  
	gmrPt	      =new vector<float>;
	gmrEta	      =new vector<float>;
	gmrPhi	      =new vector<float>;
	gmrP	      =new vector<float>;
	gmrPx	      =new vector<float>;
	gmrPy	      =new vector<float>;
	gmrPz	      =new vector<float>;
	gmrTheta	      =new vector<float>;
	gmrVx	      =new vector<float>;
	gmrVy	      =new vector<float>;
	gmrVz	      =new vector<float>;
	gmrCharge      =new vector<float>;
	gmrNDoF	      =new vector<float>;
	gmrChi2	      =new vector<float>;
	gmrChi2Norm    =new vector<float>;
	gmrDXY	      =new vector<float>;
	gmrDTheta      =new vector<float>;
	gmrDPt	      =new vector<float>;
	gmrDEta	      =new vector<float>;
	gmrDPhi	      =new vector<float>;
	gmrDDXY	      =new vector<float>;
	gmrIso03nTracks=new vector<float>;
	gmrIso03sumPt  =new vector<float>;

	gmrDz	  =new vector<float>;
	gmrD0	  =new vector<float>;
	gmrDsz	  =new vector<float>;
	gmrDDz	  =new vector<float>;
	gmrDD0	  =new vector<float>;
	gmrDDsz	  =new vector<float>;
	gmrInnerX =new vector<float>;
	gmrInnerY =new vector<float>;
	gmrInnerZ =new vector<float>;

	gmrNumberOfValidHits = new vector<float>;
	gmrNumberOfLostHits  = new vector<float>;
	gmrTrackerHits  = new vector<float>;
}
void NtupleDumper::_trigInit(){
	hltResults = new vector<bool>;
	l1Results = new vector<bool>;
}

void NtupleDumper::_fillEvtInfo(const edm::Event& event){
	//----fill event info-----//FIXME-//need to fill always the run and num of event
	evtRun = (int)event.id().run();
	evtNum = (int)event.id().event();


	if( !runOnCSA07Soup_ ) return;

		double weight = 1;
		Handle<double> weightHandle;
		event.getByLabel (csa07EventWeightTag_, weightHandle);
		weight = * weightHandle;

		//------ get processid for CSA07
		int processId  = csa07::csa07ProcessId(event, _OveralLumiForSoup, _CSA07WeightModuleLabel );

		float genScale = 0;
		float genWeight =0;
		int genProcID = 0;
		float genFilterEff = -1;

		Handle<double> genWeightHandle;
		event.getByLabel ("genEventWeight", genWeightHandle);
		genWeight = * genWeightHandle;
		Handle<double> genScaleHandle;
		event.getByLabel ("genEventScale", genScaleHandle);
		genScale = * genScaleHandle;
		Handle<int> genProcIDHandle;
		event.getByLabel ("genEventProcID", genProcIDHandle);
		genProcID = * genProcIDHandle;

		bool runOnChowder = false; // check if this is a chowder sample
		if (genProcID == 4) { // it's chowder!
			runOnChowder = true;
			event.getByLabel( _CSA07WeightModuleLabel, "AlpgenProcessID", genProcIDHandle);
			genProcID = *genProcIDHandle;
		}
		if (!runOnChowder) {
			edm::Handle<double> filterEffH;
			event.getByLabel("genEventRunInfo", "FilterEfficiency", filterEffH);
			genFilterEff = *filterEffH;
		}
	//----fill event 
	evtWeight = weight;
	evtProcessId = processId;
	genEventProcID = genProcID;
	genEventScale = genScale;
	genEventWeight = genWeight;
	genEventFilterEff = genFilterEff;
}

void NtupleDumper::_fillHLTrig(const edm::Handle<edm::TriggerResults> hltTR){
	vector<string> m_TrigNames;
	m_TrigNames.clear();
	m_TrigNames.push_back ("HLT1MuonIso"       ) ;
	m_TrigNames.push_back ("HLT1MuonNonIso"    ) ;
	m_TrigNames.push_back ("CandHLT2MuonIso"   ) ;
	m_TrigNames.push_back ("HLT2MuonNonIso"    ) ;

	int TrSize = hltTR->size();
	// if(debug_==true) std::cout<<"hltTR.size "<<TrSize<<endl;

	if(!initTriggerNamesService_){
		initTriggerNamesService_=true;
		triggerNames_.init(*hltTR);
	}

	bool fired ;
	// loop over trigger paths
	for (int tr=0 ; tr < (int)m_TrigNames.size() ; ++tr) {
		int ind = triggerNames_.triggerIndex(m_TrigNames[tr]);
		//if(debug_==true) std::cout<<"hltTR->find()="<<ind<<endl;
		fired = false;                  // 
		if ((ind < TrSize) && (hltTR->accept (ind)))
			fired = true ;
		//if(debug_==true) std::cout<<"fired = "<<fired<<endl;
		hltResults->push_back(fired);
	} // loop over triggers
}

void NtupleDumper::_fillL1Trig(const edm::Handle<l1extra::L1ParticleMapCollection> L1PMC){
	vector<int> m_TrigNames;
	m_TrigNames.clear();
	m_TrigNames.push_back ( l1extra::L1ParticleMap::kSingleMu7 ) ;
	m_TrigNames.push_back ( l1extra::L1ParticleMap::kDoubleMu3 ) ;

	bool fired ;
	// loop over trigger paths
	for (int tr=0 ; tr < (int)m_TrigNames.size() ; ++tr) {
		int ind = m_TrigNames[tr] ; //
		fired=false;
		if(ind <(int) L1PMC->size())
			fired = (bool)(*L1PMC)[ind].triggerDecision();                  // 
		l1Results->push_back(fired);
	} // loop over triggers
}

void NtupleDumper::_fillPyt(const edm::Handle<reco::CandidateCollection>  genParticles){
	//---filled particles number---
	unsigned int npyt = 0;

	int idx  = -1;
	int iMo1 = -1;
	int iDa1 = -1;
	std::vector<const reco::Candidate *> cands_;
	cands_.clear();
	vector<const Candidate *>::const_iterator found = cands_.begin();
	for( CandidateCollection::const_iterator p = genParticles->begin();
			p != genParticles->end(); ++ p ) {
		cands_.push_back( & * p );
	}

	// if(debug_==true) std::cout<<"test5"<<endl;

	for ( CandidateCollection::const_iterator gen = genParticles->begin();
			gen != genParticles->end(); gen++) {

		// if(debug_==true) std::cout<<"test9"<<endl;

		//e(11),mu(13),tau(15),Z(23),W(24),h0,H0(25,25),b(5),t(6)
		if (gen->status() == 3) { 

			//  if(debug_==true) std::cout<<"test8"<<endl;

			pytId->push_back(gen->pdgId());
			pytStatus->push_back(gen->status());
			pytPt->push_back(gen->pt());
			pytPhi->push_back(gen->phi());
			pytEta->push_back(gen->eta());
			pytTheta->push_back(gen->theta());
			pytP->push_back(gen->p());

			//  if(debug_==true) std::cout<<"test15"<<endl;
			pytEnergy->push_back(gen->energy());
			pytPx->push_back(gen->px());
			pytPy->push_back(gen->py());
			pytPz->push_back(gen->pz());

			//  if(debug_==true) std::cout<<"test16"<<endl;
			pytMass->push_back(gen->mass());
			pytVx->push_back(gen->vx());
			pytVy->push_back(gen->vy());
			pytVz->push_back(gen->vz());

			//  if(debug_==true) std::cout<<"test10"<<endl;

			idx = gen - genParticles->begin();
			pytIndex->push_back(idx);


			//  if(debug_==true) std::cout<<"test11"<<endl;
			//two mother? how to deal with that
			pytCharge->push_back(gen->charge());  


			int nMo = gen->numberOfMothers();
			pytNMother->push_back(nMo);
			int nDa = gen->numberOfDaughters();
			pytNDaughter->push_back(nDa);

			//   if(debug_==true) std::cout<<"test12"<<endl;	
			found = find( cands_.begin(), cands_.end(), gen->mother(0) );
			if ( found != cands_.end() ) iMo1 = found - cands_.begin() ;

			//if(debug_==true) std::cout<<"test13_iMo1="<<iMo1<<endl;
			pytMotherIndex->push_back(iMo1);  
			if(iMo1>=(int)genParticles->size()) {
				cout<<"iMo1 index override"<<endl;
				continue;
			}
			int mothId;
			if(iMo1<0) mothId = 0;
			else mothId=(genParticles->begin()+iMo1)->pdgId();
			//  if(debug_==true) std::cout<<"test15"<<endl;
			pytMotherPdgId->push_back(mothId);
			//  if(debug_==true) std::cout<<"test14"<<endl;

			//if (nDa>0 )
			//{
			found = find( cands_.begin(), cands_.end(), gen->daughter(0) );
			if ( found != cands_.end() ) iDa1 = found - cands_.begin() ;

			//if(debug_==true) std::cout<<"test13_iMo1="<<iMo1<<endl;
			pytDaughterIndex->push_back(iDa1);  
			if(iDa1>=(int)genParticles->size()) {
				cout<<"iDa1 index override"<<endl;
				continue;
			}
			int daugId;
			if(iDa1<0) daugId = 0;
			else daugId=(genParticles->begin()+iDa1)->pdgId();
			//  if(debug_==true) std::cout<<"test15"<<endl;
			pytDaughterPdgId->push_back(daugId);
			//}

			npyt ++;

		}//if status ==3

		else {

			// if(debug_==true) std::cout<<"test8"<<endl;
			//e(11),mu(13),tau(15),Z(23),W(24),h0,H0(25,25),b(5),t(6), jpsi(443), pi(+/-211, 111), kaon(+/-321, 311), 
			int id = abs(gen->pdgId());
			if (id==11 || id ==13 || id==15 || id==23 || id==24 || id ==5 || id==6 ||id==32 ||
					id==1 || id==2 || id==3 || id==4 ||id==21 )
				// - - - - - - - - - - - -when i run on chowder(w+j, tt+j, z+j), if i include pi kaon, 
				// - - - - - - - - - - - -the number of pi kaon will be up to 500 particles per event
				// -	-	-	-	it increase the ntuple size to 1MB/100evts, 
				// - - - - - - -for 20M events, it will dump 200GB files
				// - - -	- - 	-  -	 - it's crazy, so don't include them into ntuple
				//		id==211 || id == 111 || id== 321 || id==311 || id == 443)
			{

				if( (id==1 || id==2 || id==3 ||id==21) && gen->pt()<10. ) continue;

				//	  if(debug_==true) std::cout<<"test7"<<endl;
				pytId->push_back(gen->pdgId());
				pytStatus->push_back(gen->status());
				pytPt->push_back(gen->pt());
				pytPhi->push_back(gen->phi());
				pytEta->push_back(gen->eta());
				pytTheta->push_back(gen->theta());
				pytP->push_back(gen->p());
				pytEnergy->push_back(gen->energy());
				pytPx->push_back(gen->px());
				pytPy->push_back(gen->py());
				pytPz->push_back(gen->pz());
				pytMass->push_back(gen->mass());
				pytVx->push_back(gen->vx());
				pytVy->push_back(gen->vy());
				pytVz->push_back(gen->vz());

				idx = gen - genParticles->begin();
				pytIndex->push_back(idx);

				//two mother? how to deal with that
				pytCharge->push_back(gen->charge());   

				int nMo = gen->numberOfMothers();
				pytNMother->push_back(nMo);

				found = find( cands_.begin(), cands_.end(), gen->mother(0) );
				if ( found != cands_.end() ) iMo1 = found - cands_.begin() ;

				pytMotherIndex->push_back(iMo1);   
				int mothId=(genParticles->begin()+iMo1)->pdgId();
				pytMotherPdgId->push_back(mothId);

				int nDa = gen->numberOfDaughters();
				pytNDaughter->push_back(nDa);
				found = find( cands_.begin(), cands_.end(), gen->daughter(0) );
				if ( found != cands_.end() ) iDa1 = found - cands_.begin() ;

				pytDaughterIndex->push_back(iDa1);  
				int daugId;
				if(iDa1<0) daugId = 0;
				else daugId=(genParticles->begin()+iDa1)->pdgId();
				pytDaughterPdgId->push_back(daugId);

				npyt++;
			}//if it's the particle we'd like
		}//else
	}//genparticles collection
	//        if(debug_==true) std::cout<<"test6"<<endl;


	pytSize = npyt;
	//-----------fill pyt particle end
}

void NtupleDumper::_fillGmr(const edm::Handle<reco::MuonCollection> muons){
	//-----------fill gmr start-------
	gmrSize = muons->size(); 

	std::cout<<"gmrSize  and muons->size(): "<<gmrSize<<" "<<muons->size()<<std::endl;

	for( MuonCollection::const_iterator muon= muons->begin();        
			muon != muons->end();  muon++ ) {


		if(muon->isIsolationValid())
		{
			MuonIsolation isoR03 = muon->getIsolationR03();
			gmrIso03nTracks->push_back(isoR03.nTracks);
			gmrIso03sumPt->push_back(isoR03.sumPt);
		}
		else
		{
			gmrIso03nTracks->push_back(-1);
			gmrIso03sumPt->push_back(-1);
		}
		gmrPt->push_back(muon->pt());
		gmrPhi->push_back(muon->phi());
		gmrEta->push_back(muon->eta());
		gmrP->push_back(muon->p());
		gmrEnergy->push_back(muon->energy());
		gmrPx->push_back(muon->px());
		gmrPy->push_back(muon->py());
		gmrPz->push_back(muon->pz());
		gmrTheta->push_back(muon->theta());
		gmrCharge->push_back(muon->charge());
		gmrVx->push_back(muon->vx());
		gmrVy->push_back(muon->vy());
		gmrVz->push_back(muon->vz());

		if(muon->track().isNonnull())
		{
			gmrTrackerHits->push_back(muon->track()->numberOfValidHits());
		}
		else 
		{
			gmrTrackerHits->push_back(0);
		}

		//   if(debug_==true) std::cout<<"test1"<<endl;
		if(muon->combinedMuon().isNonnull())
		{
			gmrDXY->push_back(muon->combinedMuon()->dxy());
			gmrDDXY->push_back(muon->combinedMuon()->dxyError());
			gmrDPhi->push_back(muon->combinedMuon()->phiError());
			gmrDEta->push_back(muon->combinedMuon()->etaError());
			gmrDPt->push_back(muon->combinedMuon()->ptError());
			gmrChi2->push_back(muon->combinedMuon()->chi2());
			gmrNDoF->push_back(muon->combinedMuon()->ndof());
			gmrChi2Norm->push_back(muon->combinedMuon()->normalizedChi2());
			gmrDTheta->push_back(muon->combinedMuon()->thetaError());
			gmrDz->push_back(muon->combinedMuon()->dz());
			gmrD0->push_back(muon->combinedMuon()->d0());
			gmrDsz->push_back(muon->combinedMuon()->dsz());
			gmrDDz->push_back(muon->combinedMuon()->dzError());
			gmrDD0->push_back(muon->combinedMuon()->d0Error());
			gmrDDsz->push_back(muon->combinedMuon()->dszError());
			gmrNumberOfValidHits->push_back(muon->combinedMuon()->numberOfValidHits());
			gmrNumberOfLostHits ->push_back(muon->combinedMuon()->numberOfLostHits ());
		}
		else
		{
			if(muon->track().isNonnull())
			{
				gmrDXY->push_back(muon->track()->dxy());
				gmrDDXY->push_back(muon->track()->dxyError());
				gmrDPhi->push_back(muon->track()->phiError());
				gmrDEta->push_back(muon->track()->etaError());
				gmrDPt->push_back(muon->track()->ptError());
				gmrChi2->push_back(muon->track()->chi2());
				gmrNDoF->push_back(muon->track()->ndof());
				gmrChi2Norm->push_back(muon->track()->normalizedChi2());
				gmrDTheta->push_back(muon->track()->thetaError());
				gmrDz->push_back(muon->track()->dz());
				gmrD0->push_back(muon->track()->d0());
				gmrDsz->push_back(muon->track()->dsz());
				gmrDDz->push_back(muon->track()->dzError());
				gmrDD0->push_back(muon->track()->d0Error());
				gmrDDsz->push_back(muon->track()->dszError());
				gmrNumberOfValidHits->push_back(muon->track()->numberOfValidHits());
				gmrNumberOfLostHits ->push_back(muon->track()->numberOfLostHits ());
			}
			else 
			{
				gmrDXY->push_back(-1);
				gmrDDXY->push_back(-1);
				gmrDPhi->push_back(-1);
				gmrDEta->push_back(-1);
				gmrDPt->push_back(-1);
				gmrChi2->push_back(-1);
				gmrNDoF->push_back(-1);
				gmrChi2Norm->push_back(-1);
				gmrDTheta->push_back(-1);
				gmrDz->push_back(-1);
				gmrD0->push_back(-1);
				gmrDsz->push_back(-1);
				gmrDDz->push_back(-1);
				gmrDD0->push_back(-1);
				gmrDDsz->push_back(-1);
				gmrNumberOfValidHits->push_back(-1);
				gmrNumberOfLostHits ->push_back(-1);
			}
		}

		//*************not exist in AODSIM format-------------
		if(!runOnAODSIM_){
			gmrInnerX->push_back(muon->combinedMuon()->innerPosition().X());
			gmrInnerY->push_back(muon->combinedMuon()->innerPosition().Y());
			gmrInnerZ->push_back(muon->combinedMuon()->innerPosition().Z());
		}
		// if(debug_==true) std::cout<<"test2"<<endl;

		//    gmrQoverPError->push_back(muon->combinedMuon()->qoverpError());
		//    gmrQoverP->push_back(muon->combinedMuon()->qoverp());

	}
	//-------------gmr fill end
}
void NtupleDumper::_fillMet(const edm::Handle<reco::CaloMETCollection> caloMet){
	//-------------met fill start-------
	CaloMETCollection::const_iterator met = caloMet->begin(); 
	metEt->push_back( met->et() );
	metSumEt->push_back( met->sumEt() );
	metPhi->push_back( met->phi() );
	metEta->push_back( met->eta() );
	metPx->push_back( met->px() );
	metPy->push_back( met->py() );
	metPz->push_back( met->pz() );      
	//-----------------met fill end-------
}
void NtupleDumper::_fillJet(const edm::Handle<reco::CaloJetCollection> caloJets) {
	//-----------------jet fill start-----------
	jetSize=caloJets->size();	
	for(CaloJetCollection::const_iterator jet = caloJets->begin();
			jet!=caloJets->end();jet++)
	{
		jetEt->push_back(jet->et());    	
		jetPt->push_back(jet->pt());
		jetPhi->push_back(jet->phi());
		jetEta->push_back(jet->eta());
		jetP->push_back(jet->p());
		jetEnergy->push_back(jet->energy());
		jetPx->push_back(jet->px());
		jetPy->push_back(jet->py());
		jetPz->push_back(jet->pz());
		jetTheta->push_back(jet->theta());
	}
	//------------------jet fill end
}
DEFINE_FWK_MODULE( NtupleDumper );



