// -*- C++ -*-
//
// Package:    QCDBackground
// Class:      QCDBackground
// 
/**\class QCDBackground QCDBackground.cc UserArea/QCDBackground/src/QCDBackground.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Theodore Nicholas Kypreos,8 R-031,+41227675208,
//         Created:  Thu Aug  5 17:04:29 CEST 2010
// $Id: QCDBackground.cc,v 1.1 2010/09/28 10:26:02 kypreos Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
//
// class declaration
//
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"



#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h" 
#include "FWCore/Common/interface/TriggerNames.h" 

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "FWCore/Common/interface/TriggerNames.h"
//#include "UserArea/QCDBackground/interface/AllJetsNtuple.h"

class QCDAnalObject{

  public:
	QCDAnalObject();
	QCDAnalObject(pat::Jet const jet);
	QCDAnalObject(pat::Jet const jet, int const charge);

	bool const IsCloserMuon(pat::Muon const muon);
	bool const HasMuon() const{return _hasMuon;}
	pat::Muon const muon() const {return _muon;}
	pat::Jet const jet() const {return _jet;}
	void KillMuon(){_hasMuon = false;}
	int const charge() const{return _jetCharge;}

	void const Print() const;

//	bool const isValid() const;	
//	QCDAnalObject const CompareObjectMuons(QCDAnalObject& obj, pat::Muon const muon);

  private:

	pat::Muon _muon;
	pat::Jet _jet;	
	bool _hasMuon;
	int _jetCharge;
	double const static _MAX_DELTAR;

	double const GetDeltaR(pat::Muon const muon) const;
	bool const WithinDeltaR(pat::Muon const muon) const;

	


};
double const QCDAnalObject::_MAX_DELTAR = 0.3;



QCDAnalObject::QCDAnalObject(pat::Jet const jet):_hasMuon(false),_jetCharge(0){
	_jet = jet; 
}
QCDAnalObject::QCDAnalObject(pat::Jet const jet, int const charge):_hasMuon(false){
	_jet = jet; 
	_jetCharge = charge;
}
bool const QCDAnalObject::IsCloserMuon(pat::Muon const muon){
	bool isCloser = false;
//	std::cout<<"checking the deltaR... "<<GetDeltaR(muon)<<std::endl;
	if (!_hasMuon && WithinDeltaR(muon)) {
		_muon = muon;
		_hasMuon = true;
		isCloser = true;
	} else if(_hasMuon) {
		double deltaR1 = reco::deltaR(_muon,_jet);
		double deltaR2 = reco::deltaR(muon,_jet);
		if (deltaR2 < deltaR1) {
			_muon = muon;
			isCloser = true;
		}
	}
	return isCloser;
}
//
// method for getting deltaR
//
double const QCDAnalObject::GetDeltaR(pat::Muon const muon) const{
	return reco::deltaR(muon,this->_jet);
}
//
// check if DeltaR is within acceptance
//
bool const QCDAnalObject::WithinDeltaR(pat::Muon const muon) const{
	if (GetDeltaR(muon)<=_MAX_DELTAR) return true;
	else return false;
}
/*
QCDAnalObject QCDAnalObject::CompareObjectMuons(QCDAnalObject& obj, pat::Muon const muon){
	double deltaR1 = this->GetDeltaR(muon);
	double deltaR2 = obj.GetDeltaR(muon);
	if (deltaR1 < deltaR2) {
		obj.KillMuon();
		return *this;
	} else {
		this->KillMuon();
		return obj;
	} 
}
*/
//
//
//
void const QCDAnalObject::Print() const{

	printf("jet: pt = %4.1f\t eta = %+3.2f\t phi = %+3.2f",
		_jet.pt(),_jet.eta(),_jet.phi());

	printf("\t hasMuon?: %d", _hasMuon);
	if (_hasMuon) printf("\tmuon: pt = %4.1f\t eta = %+3.2f\t phi = %+3.2f",
		_muon.pt(),_muon.eta(),_muon.phi());
	
	printf("\n");
}


class QCDBackground : public edm::EDAnalyzer {
   public:
      explicit QCDBackground(const edm::ParameterSet&);
      ~QCDBackground();


  private:
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
      // ----------member data ---------------------------

  	std::string _getFilename;	

	TRandom3* _rand3;

	TFile* _outFile;
	TFile* _inFile;
	TH2F* _h2PtVsPhi;

	TH3F* _h3PVsEtaVsPhi1;
	TH3F* _h3PVsEtaVsPhi2;


	TH1F* hCountTriggersMu9;

	TH1F* hMuonP;
	TH1F* hMuonPt;
	TH1F* hMuonEta;
	TH1F* hMuonPhi;
	TH1F* hMuonD0;
	TH1F* hMuonCharge;

	TH1F* hJetPt;
	TH1F* hJetP;
	TH1F* hJetEta;
	TH1F* hJetPhi;
	TH2F* _h2PVsEta;

	TH1F* hJetPt1;
	TH1F* hJetP1;
	TH1F* hJetEta1;
	TH1F* hJetPhi1;
	TH2F* _h2PVsEta1;

	TH1F* hJetPt1Pos;
	TH1F* hJetP1Pos;
	TH1F* hJetEta1Pos;
	TH1F* hJetPhi1Pos;

	TH1F* hJetPt1Neg;
	TH1F* hJetP1Neg;
	TH1F* hJetEta1Neg;
	TH1F* hJetPhi1Neg;

	TH1F* hJetDeltaEta;
	TH1F* hJetDeltaPhi;
	TH1F* hJetDeltaR;



	TH2F* hJetPtVsMuonPt;
	TH2F* hJetPVsMuonP;
	TH2F* hJetPtVsEta;
	TH2F* hJetPtVsEta1;



	TH1F* hNumSelMuons;
	TH1F* hNumSelJets;



	TH1F* hJetMinDeltaEta;
	TH1F* hJetMinDeltaPhi;
	TH1F* hJetMinDeltaR;


	TH1F* hDiMuonMass;
	TH1F* hDiMuonMassSS;
	TH1F* hDiMuonMassOS;

	TH1F* hDiMuonPOS;
	TH1F* hDiMuonPtOS;
	TH1F* hDiMuonEtaOS;
	TH1F* hDiMuonPhiOS;

	TH1F* hDiMuonPSS;
	TH1F* hDiMuonPtSS;
	TH1F* hDiMuonEtaSS;
	TH1F* hDiMuonPhiSS;




	TH1F* hJetMassAll;
	TH1F* hJetMassAll1;
	TH1F* hJetMassCharged;
	TH1F* hJetMassNoCharge;	


	TH1F* hDiJetDeltaPhi;
	TH1F* hDiJetDeltaPhiNoCharge;

	TTree* _outTree;
	int _run;
	int _event;
	int _lumi;
	float _jetPt;
	float _jetP;
	float _jetEta;
	float _jetPhi;
	float _muonSumPt;
	int _hasMatch;	

	int _hasDiMuon;


  typedef struct {
    int charge;
    float pt;
    float eta;
    float phi;

  } _TrackInfo;
	_TrackInfo _jet1, _jet2;
	_TrackInfo _outJet, _outMuon;
	_TrackInfo _muon1, _muon2;

	TTree* _outTreeMasses;
	TTree* _outTreeMassesMuon;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
QCDBackground::QCDBackground(const edm::ParameterSet& iConfig):_run(0),_event(0),_lumi(0),_jetPt(0.0),_jetP(0.0),_jetEta(0.0),_jetPhi(0.0),_hasMatch(0)

{
  _getFilename  	= iConfig.getUntrackedParameter<std::string>("getFilename", "badger.root");
   //now do what ever initialization is needed
	_rand3 = new TRandom3(0);
}


QCDBackground::~QCDBackground()
{
	_rand3->Delete(); 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
QCDBackground::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//
// trigger checks
//
  	edm::Handle<edm::TriggerResults> triggerResults;
  	iEvent.getByLabel("TriggerResults", triggerResults);
//  	iEvent.getByLabel("patTrigger", triggerResults);

//
// check the di-jet triggers
//
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
	int trigMu9 = triggerNames.triggerIndex("HLT_Mu9");
	int dijet1 = triggerNames.triggerIndex("HLT_DiJetAve15U_8E29");
	int dijet2 = triggerNames.triggerIndex("HLT_DiJetAve30U_8E29");
	if (!triggerResults->accept(dijet1) && !triggerResults->accept(dijet2)) return;
	

	if (!dijet1 && !dijet2) return;

//
// get the run event lumi
//
  	_run 	= (int) iEvent.id().run();
  	_event 	= (int)iEvent.luminosityBlock();
  	_lumi 	= (int)iEvent.id().event();


//
// check if the muon trigger fired
//
	if (triggerResults->accept(trigMu9)) {
		hCountTriggersMu9->Fill(1);
	} else hCountTriggersMu9->Fill(0);



///////////////////////////////////////////////////////////////////////////////
//
// jets
//
///////////////////////////////////////////////////////////////////////////////
	edm::Handle<pat::JetCollection> pfJets;
	iEvent.getByLabel("cleanPatJets", pfJets);
	
	typedef std::vector<QCDAnalObject> AnalObjects;
	AnalObjects anals; anals.clear();

	
	for (pat::JetCollection::const_iterator jet = pfJets->begin(); jet != pfJets->end(); ++jet){
		if (fabs(jet->eta())>2.5) continue;
		int jetCharge = _rand3->Rndm()>0.5?1:-1;
		QCDAnalObject qcdObj(*jet,jetCharge);
		anals.push_back(qcdObj);
	}	
//	std::cout<<"number of jet objects: "<<anals.size()<<std::endl;

///////////////////////////////////////////////////////////////////////////////
//
// muons
//
///////////////////////////////////////////////////////////////////////////////
//	edm::Handle<std::vector<pat::Muon> > muonColl;
	edm::Handle<pat::MuonCollection > muonColl;
//	iEvent.getByLabel("cleanPatMuons", muonColl);
	iEvent.getByLabel("selectedMuons","selectedMuons",muonColl);

//	std::cout<<"size of muon collection: "<<muonColl->size()<<std::endl;
	if (muonColl->size() > 1) _hasDiMuon = 1;
	else _hasDiMuon = 0;


	
///////////////////////////////////////////////////////////////////////////////
//
// find the closest muon to each jet 
//
///////////////////////////////////////////////////////////////////////////////
/*
	for (pat::MuonCollection::const_iterator muon = muonColl->begin(); muon != muonColl->end(); ++muon){

		QCDAnalObject ano();
//		AnalObjects anals; anals.clear();
		for (AnalObjects::iterator anit = anals.begin(); anit!=anals.end(); ++anit){
			bool isCloser = anit->IsCloserMuon(*muon);
			
		}

	}
*/

	pat::MuonCollection selectedMuons= *muonColl;

	AnalObjects selAnals; selAnals.clear();

//	std::cout<<std::endl;
//	std::cout<<"size of anals: "<<anals.size()<<std::endl;	
//	std::cout<<"num matched: "<<selAnals.size()<<std::endl;

	for (AnalObjects::iterator anit = anals.begin(); anit!=anals.end(); ++anit){
		pat::MuonCollection::iterator muon = selectedMuons.begin();
		pat::MuonCollection::iterator closestMuon = selectedMuons.end();
//		anit->Print();
		for ( ;muon != selectedMuons.end(); ++muon){
			bool isCloser = anit->IsCloserMuon(*muon);
			if (isCloser){
//				std::cout<<"this is within range!"<<std::endl;
//				if (anit->HasMuon()) std::cout<<"\t\t i can has muon"<<std::endl;
				closestMuon = muon;
			}
		}
		selAnals.push_back(*anit);	
		if (closestMuon != selectedMuons.end()) {
//			std::cout<<"erasing the matched muon!"<<std::endl;
			selectedMuons.erase(closestMuon);
		}
	}	
//	std::cout<<std::endl;
//	std::cout<<"do some final printing..."<<std::endl;
	for (AnalObjects::const_iterator anit= selAnals.begin(); anit != selAnals.end(); ++anit) {
//		anit->Print();
		_outJet.charge = anit->charge(); 
		_outJet.pt		= anit->jet().pt();
		_outJet.eta		= anit->jet().eta();
		_outJet.phi		= anit->jet().phi();

		
		_hasMatch	= (int)anit->HasMuon();
		
		if (_hasMatch){
			_outMuon.charge = anit->muon().charge(); 
			_outMuon.pt		= anit->muon().pt();
			_outMuon.eta	= anit->muon().eta();
			_outMuon.phi	= anit->muon().phi();
			_muonSumPt		= anit->muon().trackIso();
		}else {
			_outMuon.charge = 0;	
    		_outMuon.pt		= 0;	
    		_outMuon.eta	= 0;	
			_outMuon.phi	= 0;	
			_muonSumPt		= 0;
		}

//	_outTree	->Branch("hasMatch"	, &_hasMatch,   "hasMatch/I"	);	
//	_outTree	->Branch("hasDiMuon", &_hasDiMuon,  "hasDiMuon"	);	
//	_outTree	->Branch("muonSumPt", &_muonSumPt,  "muonSumPt"	);	
		_outTree->Fill();
	}

///////////////////////////////////////////////////////////////////////////////
//
// get jet pairs and make a mass tree  omnomnomnomnom
//
///////////////////////////////////////////////////////////////////////////////


	for (AnalObjects::const_iterator anit1= selAnals.begin(); anit1 != selAnals.end(); ++anit1) {

		_jet1.charge 	= anit1->charge(); 
		_jet1.pt		= anit1->jet().pt();
		_jet1.eta		= anit1->jet().eta();
		_jet1.phi		= anit1->jet().phi();
		for (AnalObjects::const_iterator anit2= anit1; anit2 != selAnals.end(); ++anit2) {
			if (anit1==anit2) continue;			
			_jet2.charge 	= anit2->charge(); 
			_jet2.pt		= anit2->jet().pt();
			_jet2.eta		= anit2->jet().eta();
			_jet2.phi		= anit2->jet().phi();
			
			_outTreeMasses->Fill();	
		}

	}

//
// dump muon ntuple
//
	for (pat::MuonCollection::const_iterator muon1 = selectedMuons.begin(); muon1!=selectedMuons.end();++muon1){
		_muon1.charge 	= muon1->charge(); 
		_muon1.pt		= muon1->pt();
		_muon1.eta		= muon1->eta();
		_muon1.phi		= muon1->phi();
		for (pat::MuonCollection::const_iterator muon2 = muon1; muon2!=selectedMuons.end();++muon2){
			if (muon1 == muon2) continue;
			_muon2.charge 	= muon2->charge(); 
			_muon2.pt		= muon2->pt();
			_muon2.eta		= muon2->eta();
			_muon2.phi		= muon2->phi();
			_outTreeMassesMuon->Fill();	
		}	
	}
//	return;

#ifdef BADGERBADGER 





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
	pat::JetCollection selectedJets;
	typedef std::pair<pat::Jet,int> JetC;
	typedef std::vector<JetC> JetCColl;
	JetCColl selJetsC; selJetsC.clear();

	for (pat::JetCollection::const_iterator jet = pfJets->begin(); jet != pfJets->end(); ++jet){
		if (fabs(jet->eta())>2.5) continue;
		selectedJets.push_back(*jet);
		_h3PVsEtaVsPhi1->Fill(jet->p(),jet->eta(), jet->phi());

		int jetCharge = _rand3->Rndm()>0.5?1:-1;
		selJetsC.push_back(JetC(*jet,jetCharge));

		hJetP->Fill(	jet->p()	);
		hJetPt->Fill(	jet->pt()	);
		hJetEta->Fill(	jet->eta()	);
		hJetPhi->Fill(	jet->phi()	);
		_h2PVsEta	->Fill(jet->p(),jet->eta());	
		hJetPtVsEta ->Fill(jet->p(), jet->eta());
	}

	
//
// get the muons
// 


	pat::MuonCollection badger= *muonColl;
/*
	if (muonColl->size()<2) return;
	std::cout<<std::endl<<std::endl;
	std::cout<<"size of orig: "<<muonColl->size()<<std::endl;
	std::cout<<"size of copy:"<<badger.size()<<std::endl;

	if (badger.empty()) return;

	std::cout<<"1: pt: "<<badger.front().pt()<<std::endl;

	pat::MuonCollection::iterator  badgerit = badger.begin();
	
	badger.erase(badgerit);
	std::cout<<"size of orig: "<<muonColl->size()<<std::endl;
	std::cout<<"size of copy:"<<badger.size()<<std::endl;
	std::cout<<"1: pt: "<<badger.front().pt()<<std::endl;
*/		
/*
	for (pat::MuonCollection::const_iterator badgerit = badger.begin(); badgerit!=badger.end(); ++badgerit){		
			
		for (pat::MuonCollection::const_iterator muon= muonColl->begin(); muon!=muonColl->end(); ++muon){		
			if (badgerit == muon){
				std::cout<<"same it!"<<std::endl;
			} else if (&(*badgerit) == &(*muon)){
				std::cout<<"same thingy!"<<std::endl;
			}else printf("nope\n");
			
		}
	}
*/
	

//	return;
/*
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);

        int trig1 = triggerNames.triggerIndex("HLT_Activity_L1A");
        int trig2 = triggerNames.triggerIndex("HLT_Mu9");

        if (triggerResults->accept(trig1)) std::cout<<triggerNames.triggerName(trig1)<<std::endl;
        else return;
        if (triggerResults->accept(trig2)) std::cout<<triggerNames.triggerName(trig2)<<std::endl;
        else return;
*/



/*
	edm::Handle<reco::PFJetCollection> pfJets;
//	iEvent.getByLabel("ak5PFJets", pfJets);
	iEvent.getByLabel("cleanPatJets", pfJets);

//	std::cout<<"number or ak5PFJets: "<<pfJets->size()<<std::endl;

	edm::Handle<reco::MuonCollection> muonColl;
//	iEvent.getByLabel("muons", muonColl);
//	iEvent.getByLabel("ufMuonSelector", muonColl);
	iEvent.getByLabel("cleanPatMuons", muonColl);
*/

//	return;

//	std::cout<<muonColl->size()<<std::endl;
//	std::vector<reco::Muon> selectedMuons; 
	pat::MuonCollection selectedMuons;

/*
	for (pat::MuonCollection::const_iterator muon = muonColl->begin(); muon != muonColl->end(); ++muon){
		if ( !muon->isGlobalMuon()) continue;
		if ( !muon->isTrackerMuon()) continue;
		reco::TrackRef globalTrack = muon->globalTrack();

		if (globalTrack->normalizedChi2() > 10) continue;
		if (globalTrack->hitPattern().numberOfValidPixelHits()	 < 1	) continue;
		if (globalTrack->hitPattern().numberOfValidTrackerHits() < 10	) continue;
		if (globalTrack->hitPattern().numberOfValidMuonHits() 	 < 1	) continue;
*/
	

	for (pat::MuonCollection::const_iterator muon = muonColl->begin(); muon != muonColl->end(); ++muon){
//		if (muon->trackIso() > 3) continue; 
		reco::TrackRef globalTrack = muon->globalTrack();
		selectedMuons.push_back(*muon);
		hMuonP->Fill(	muon->p()	);
		hMuonPt->Fill(	muon->pt()	);
		hMuonEta->Fill(	muon->eta()	);
		hMuonPhi->Fill(	muon->phi()	);
		hMuonD0->Fill(	globalTrack->dxy()	);
		hMuonCharge->Fill(	muon->charge()	);
	}
//
// get jets
//

	hNumSelJets	->Fill(selectedJets.size());
	hNumSelMuons->Fill(selectedMuons.size());

//
// match jets to muons
//


	for (pat::JetCollection::const_iterator jet = selectedJets.begin(); jet != selectedJets.end(); ++jet){
//		if (selectedMuons.size()>=2) break;
//		bool matchedToMuon = false;
		double const _MinMatch = 0.3;
		double dEta = 99999, dPhi = 99999, delR = 99999;

		for (pat::MuonCollection::const_iterator muon = selectedMuons.begin(); muon != selectedMuons.end(); ++muon){
//			dEta 	=std::min(fabs(dEta), fabs(muon->eta()-jet->eta()));
//			dPhi 	=std::min(fabs(dPhi), fabs(reco::deltaPhi(*muon,*jet)));	
//			delR	=std::min(delR, reco::deltaR(*muon,*jet));	
			if (reco::deltaR(*muon,*jet)<delR) {
				delR	=reco::deltaR(*muon,*jet);	
				dEta 	=muon->eta()-jet->eta();
				dPhi 	=reco::deltaPhi(*muon,*jet);	
			}

			int hasMatch = 0;
			if (reco::deltaR(*muon,*jet)<= _MinMatch) {  
				hasMatch = 1;
				hJetPt1		->Fill(	jet->pt()	);
				hJetP1		->Fill(	jet->p()	);
				hJetEta1	->Fill(	jet->eta()	);
				hJetPhi1	->Fill(	jet->phi()	);
				hJetPtVsEta1 ->Fill(jet->eta(), jet->pt());
				_h2PVsEta1	->Fill(jet->p(),jet->eta());	

				if (muon->charge()>0){
					hJetPt1Pos		->Fill(	jet->pt()	);
					hJetP1Pos		->Fill(	jet->p()	);
					hJetEta1Pos	->Fill(	jet->eta()	);
					hJetPhi1Pos	->Fill(	jet->phi()	);
				} else if (muon->charge()<0){
					hJetPt1Neg		->Fill(	jet->pt()	);
					hJetP1Neg		->Fill(	jet->p()	);
					hJetEta1Neg	->Fill(	jet->eta()	);
					hJetPhi1Neg	->Fill(	jet->phi()	);
				}

				hJetPtVsMuonPt->Fill(muon->pt(), jet->pt());
				hJetPVsMuonP->Fill(muon->p(), jet->p());
//				if (selectedMuons.size()<2)
				_h3PVsEtaVsPhi2->Fill(jet->p(),jet->eta(), jet->phi());
			}
			hJetDeltaEta	->Fill(dEta);	
			hJetDeltaPhi	->Fill(dPhi);	
			hJetDeltaR		->Fill(delR);


			_jetP 		= jet->p();	
			_jetPt 		= jet->pt();	
			_jetEta 	= jet->eta();	
			_jetPhi 	= jet->phi();	
			_hasMatch 	= hasMatch;


			_outTree->Fill();

		}
	}
//
//
//

//	typedef std::vector<JetC> JetCColl;
//	JetCColl selJetsC; selJetsC.clear();

//	for (pat::JetCollection::const_iterator jet = selectedJets.begin(); jet != selectedJets.end(); ++jet){
	for (JetCColl::const_iterator jeti = selJetsC.begin(); jeti != selJetsC.end(); ++jeti){
		pat::Jet const* jet = &(jeti->first);

		TLorentzVector v1, v2, v3;
		v1.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(), 0.1045);
		int	ibin = _h2PtVsPhi->FindBin(jet->p(),jet->eta());
		double weight = _h2PtVsPhi->GetBinContent(ibin);
		if (weight >= 1) continue;

		_jet1.charge = jeti->second;
		_jet1.pt	= jet->pt();
		_jet1.eta	= jet->eta();
		_jet1.phi	= jet->phi();

		for (JetCColl::const_iterator jeti1 = selJetsC.begin(); jeti1 != selJetsC.end(); ++jeti1){
			if (jeti1==jeti) continue;
			pat::Jet const* jet1 = &(jeti1->first);
//			_TrackInfo _jet1, _jet2;
			_jet2.charge = jeti2->second;
			_jet2.pt	= jet2->pt();
			_jet2.eta	= jet2->eta();
			_jet2.phi	= jet2->phi();
			_outTreeMasses->Fill();

			if (jeti->second != jeti1->second) continue;

			double weight1 = weight; 
			v2.SetPtEtaPhiM(jet1->pt(),jet1->eta(),jet1->phi(), 0.1045);
			v3 = v1+v2;
//			ibin = _h2PtVsPhi->FindBin(jet1->p(),jet1->eta(), jet1->phi());
			ibin = _h2PtVsPhi->FindBin(jet1->p(),jet1->eta());
			weight1 *= _h2PtVsPhi->GetBinContent(ibin);
			if (weight1 >= 1) continue;
//			std::cout<<weight1<<std::endl;
			hJetMassAll->Fill(v3.M());
			hJetMassAll1->Fill(v3.M(), weight1);
		}



	}

//
//
//

	for (pat::MuonCollection::const_iterator muon = selectedMuons.begin(); muon != selectedMuons.end(); ++muon){
		TLorentzVector v1, v2, v3;
		v1.SetPtEtaPhiM(muon->pt(),muon->eta(),muon->phi(), 0.1045);
		for (pat::MuonCollection::const_iterator muon1 = muon; muon1 != selectedMuons.end(); ++muon1){
			if (muon1 == muon) continue;
//			if (muon1->charge() == muon->charge()) continue;
			v2.SetPtEtaPhiM(muon1->pt(),muon1->eta(),muon1->phi(), 0.1045);
			v3 = v1+v2;
			hDiMuonMass->Fill(v3.M());
			if (muon1->charge() == muon->charge()) {
				hDiMuonMassSS->Fill(v3.M());
				hDiMuonPSS	->Fill(v3.P());
				hDiMuonPtSS	->Fill(v3.Pt());
				hDiMuonEtaSS->Fill(v3.Eta());
				hDiMuonPhiSS->Fill(v3.Phi());
			} else if (muon1->charge() != muon->charge()) {
				hDiMuonMassOS->Fill(v3.M());
				hDiMuonPOS	->Fill(v3.P());
				hDiMuonPtOS	->Fill(v3.Pt());
				hDiMuonEtaOS->Fill(v3.Eta());
				hDiMuonPhiOS->Fill(v3.Phi());
			}
		}
	}

#endif

/*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/
}


// ------------ method called once each job just before starting event loop  ------------
void 
QCDBackground::beginJob()
{
	_outFile 	= new TFile(_getFilename.c_str(),"recreate");
//	_inFile 	= new TFile("weights1.root", "open");
	_inFile 	= new TFile("weights1xyz.root", "open");
	_h2PtVsPhi = (TH2F*)_inFile->Get("h2PtVsPhi")->Clone();
	_inFile->Close();

	_outFile->cd();
  	_outTree	= new TTree("tree", "myTree");

	_outTree	->Branch("run"		, &_run		,	"run/I"			);	
	_outTree	->Branch("event"	, &_event	,   "event/I"		);	
	_outTree	->Branch("lumi"		, &_lumi	,   "lumi/I"		);	
  	_outTree	->Branch("jet"		, &_outJet, "charge/I:pt/F:eta/F:phi/F");
  	_outTree	->Branch("muon"		, &_outMuon, "charge/I:pt/F:eta/F:phi/F");
	_outTree	->Branch("hasMatch"	, &_hasMatch,   "hasMatch/I"	);	
	_outTree	->Branch("hasDiMuon", &_hasDiMuon,  "hasDiMuon"	);	
	_outTree	->Branch("muonSumPt", &_muonSumPt,  "muonSumPt"	);	


//	_outTree	->Branch("jetPt"	, &_jetPt	,   "jetPt/F"		);	
//	_outTree	->Branch("jetP"		, &_jetP	,   "jetP/F"		);	
//	_outTree	->Branch("jetEta"	, &_jetEta	,   "jetEta/F"		);	
//	_outTree	->Branch("jetPhi"	, &_jetPhi	,   "jetPhi/F"		);	


	_outTreeMasses	= new TTree("treeMasses", "jetMasses");
  	_outTreeMasses->Branch("jet1"	, &_jet1, "charge/I:pt/F:eta/F:phi/F");
  	_outTreeMasses->Branch("jet2"	, &_jet2, "charge/I:pt/F:eta/F:phi/F");
	_outTreeMasses	->Branch("hasDiMuon", &_hasDiMuon,  "hasDiMuon"	);	

	_outTreeMassesMuon	= new TTree("treeMassesMuon", "muonMasses");
  	_outTreeMassesMuon->Branch("muon1"	, &_muon1, "charge/I:pt/F:eta/F:phi/F");
  	_outTreeMassesMuon->Branch("muon2"	, &_muon2, "charge/I:pt/F:eta/F:phi/F");


	hMuonP 		= new TH1F("hMuonP", "Muon p", 600	, 0, 3000);
	hMuonPt 		= new TH1F("hMuonPt", "Muon p_{T}", 600	, 0, 3000);
	hMuonEta		= new TH1F("hMuonEta", "Muon #eta", 32	, -2.5, 2.5);
	hMuonPhi		= new TH1F("hMuonPhi", "Muon #phi", 200	, -M_PI, M_PI);
	hMuonD0			= new TH1F("hMuonD0", "Muon D0", 800, -10, 10);
	hMuonCharge		= new TH1F("hMuonCharge", "Muon charge", 500, -2, 2);




//	Double_t pbins[] = {0,15,25,50,100,500,3000};		


//	_h3PVsEtaVsPhi1= new TH3F("h3PVsEtaVsPhi1", "P Vs #eta Vs. #phi" , 6	, pbins, 32, -2.5, 2.5, 64, -M_PI,M_PI);
//	_h3PVsEtaVsPhi2= new TH3F("h3PVsEtaVsPhi2", "P Vs #eta Vs. #phi" , 6	, pbins	, 32, -2.5, 2.5, 64, -M_PI,M_PI);
	_h3PVsEtaVsPhi1= new TH3F("h3PVsEtaVsPhi1", "P Vs #eta Vs. #phi" , 5	, 5, 1000, 16, -2.5, 2.5, 10, -M_PI,M_PI);
	_h3PVsEtaVsPhi2= new TH3F("h3PVsEtaVsPhi2", "P Vs #eta Vs. #phi" , 5	, 5, 1000, 16, -2.5, 2.5, 10, -M_PI,M_PI);

	hCountTriggersMu9 = new TH1F("hCountTriggersMu9", "count HLT_Mu9 triggers", 2, 0, 2);
	hJetPt 		= new TH1F("hJetPt", "jet p_{T}", 600	, 0, 3000);
	hJetP 		= new TH1F("hJetP", "jet p", 600	, 0, 3000);
	hJetEta		= new TH1F("hJetEta", "jet #eta", 32	, -2.5, 2.5);
	hJetPhi		= new TH1F("hJetPhi", "jet #phi", 200	, -M_PI, M_PI);
	_h2PVsEta	= new TH2F("h2PVsEta", "jet p vs. #eta", 5, 15, 1000, 10, -2.5,2.5);
	_h2PVsEta1	= new TH2F("h2PVsEta1", "jet p vs. #eta", 5, 15, 1000, 10, -2.5,2.5);

//	hJetChi2	= new TH1F("hJetChi2", "jet #chi^{2}", 100, 0, 100);

	hJetP1 			= new TH1F("hJetP1", "jet p", 600	, 0, 3000);
	hJetPt1 		= new TH1F("hJetPt1", "jet p_{T}", 600, 0, 3000);
	hJetEta1		= new TH1F("hJetEta1", "jet #eta", 200, -2.5, 2.5);
	hJetPhi1		= new TH1F("hJetPhi1", "jet #phi", 200, -M_PI, M_PI);

	hJetP1Pos 			= new TH1F("hJetP1Pos", "jet p", 600	, 0, 3000);
	hJetPt1Pos 		= new TH1F("hJetPt1Pos", "jet p_{T}", 600, 0, 3000);
	hJetEta1Pos		= new TH1F("hJetEta1Pos", "jet #eta", 200, -2.5, 2.5);
	hJetPhi1Pos		= new TH1F("hJetPhi1Pos", "jet #phi", 200, -M_PI, M_PI);

	hJetP1Neg 			= new TH1F("hJetP1Neg", "jet p", 600	, 0, 3000);
	hJetPt1Neg 		= new TH1F("hJetPt1Neg", "jet p_{T}", 600, 0, 3000);
	hJetEta1Neg		= new TH1F("hJetEta1Neg", "jet #eta", 200, -2.5, 2.5);
	hJetPhi1Neg		= new TH1F("hJetPhi1Neg", "jet #phi", 200, -M_PI, M_PI);

	hNumSelJets		= new TH1F("hNumSelJets", "num selected jets", 300, 0, 300);
	hNumSelMuons	= new TH1F("hNumSelMuons", "num selected muons", 100, 0, 100);

	hJetDeltaEta	= new TH1F("hJetsDeltaEta"	, "muon-jet #Delta#eta"	, 200, -5, 5);
	hJetDeltaPhi	= new TH1F("hJetsDeltaPhi"	, "muon-jet #Delta#phi"	, 200, -M_PI, M_PI);
	hJetDeltaR		= new TH1F("hJetsDeltaR"	, "muon-jet #Deltar"	, 200, 0, 10);

	hJetMinDeltaEta		= new TH1F("hJetsMinDeltaEta", "muon-jet #Delta#eta"	, 200, -5, 5);
	hJetMinDeltaPhi		= new TH1F("hJetsMinDeltaPhi", "muon-jet #Delta#phi"	, 200, -M_PI, M_PI);
	hJetMinDeltaR		= new TH1F("hJetsMinDeltaR"	, "muon-jet #Deltar"	, 200, 0, 10);


	hDiMuonMass= new TH1F("hDiMuonMass", "di-muon mass", 400, 0, 1000);
	hDiMuonMassSS= new TH1F("hDiMuonMassSS", "same-sign di-muon mass", 400, 0, 1000);
	hDiMuonMassOS= new TH1F("hDiMuonMassOS", "oppo-sign di-muon mass", 400, 0, 1000);

	hDiMuonPSS		= new TH1F("hDiMuonPSS"	, "same-sign di-muon p", 400, 0, 1000);
	hDiMuonPtSS		= new TH1F("hDiMuonPtSS", "same-sign di-muon p_{T}", 400, 0, 1000);
	hDiMuonEtaSS	= new TH1F("hDiMuonEtaSS", "same-sign di-muon #eta", 50, -2.5, 2.5);
	hDiMuonPhiSS	= new TH1F("hDiMuonPhiSS", "same-sign di-muon #phi", 64, -M_PI, M_PI);

	hDiMuonPOS		= new TH1F("hDiMuonPOS"	, "oppo-sign di-muon p", 400, 0, 1000);
	hDiMuonPtOS		= new TH1F("hDiMuonPtOS", "oppo-sign di-muon p_{T}", 400, 0, 1000);
	hDiMuonEtaOS	= new TH1F("hDiMuonEtaOS", "oppo-sign di-muon #eta", 50, -2.5, 2.5);
	hDiMuonPhiOS	= new TH1F("hDiMuonPhiOS", "oppo-sign di-muon #phi", 64, -M_PI, M_PI);

	hJetMassAll		= new TH1F("hJetMassAll", "di-jet mass", 400, 0, 1000);
	hJetMassAll1		= new TH1F("hJetMassAll1", "di-jet mass", 400, 0, 1000);
	hJetMassAll1->Sumw2();
	hJetMassAll->Sumw2();
	hDiMuonMass->Sumw2();
	hDiMuonMassSS->Sumw2();
	hDiMuonMassOS->Sumw2();
	hJetMassNoCharge	= new TH1F("hJetMassNoCharge"	, "di-jet mass", 1000, 0, 1000);
	hJetMassCharged		= new TH1F("hJetMassCharged"	, "di-jet mass", 400, 0, 1000);

	hDiJetDeltaPhi	= new TH1F("hDiJetDeltaPhi", "#Delta#phi between jets", 200, -M_PI, M_PI);
	hDiJetDeltaPhiNoCharge	= new TH1F("hDiJetDeltaPhiNoCharge", "#Delta#phi between jets", 200, -M_PI, M_PI);

	hJetPtVsMuonPt	= new TH2F("hJetPtVsMuonPt"	, "jet p_{T} vs. muon p_{T}", 400, 0, 1000, 400, 0, 1000);
	hJetPVsMuonP	= new TH2F("hJetPVsMuonP"	, "jet p vs. muon p", 400, 0, 1000, 400, 0, 1000);
	hJetPtVsEta 	= new TH2F("hJetPtVsEta" 	, "jet p_{T} vs. #eta", 32, -2.5, 2.5, 400, 0, 1000);
	hJetPtVsEta1	= new TH2F("hJetPtVsEta1"	, "jet p_{T} vs. #eta", 32, -2.5, 2.5, 400, 0, 1000);
/*
	for (int i = 0; i<3; ++i){
		TString s= "q_{jet} = ";
		if (i == 0) s+="0";	
		if (i == 1) s+="+1";	
		if (i == 2) s+="-1";	
		hqJetDeltaEta[i]	= new TH1F(TString::Format("hJetsDeltaEta%d",i)	, TString("muon-jet #Delta#eta")+s	, 200, -5, 5);
		hqJetDeltaPhi[i]	= new TH1F(TString::Format("hJetsDeltaPhi%d",i)	, TString("muon-jet #Delta#phi")+s	, 200, -M_PI, M_PI);
		hqJetDeltaR[i]		= new TH1F(TString::Format("hJetsDeltaR%d"	,i)	, TString("muon-jet #Delta#r")+s	, 200, 0, 10);

//		hqJetMinDeltaEta[i]	= new TH1F(TString::Format("hJetsMinDeltaEta%d",i)	, TString("muon-jet #Delta#eta")+s	, 200, -5, 5);
//		hqJetMinDeltaPhi[i]	= new TH1F(TString::Format("hJetsMinDeltaPhi%d",i)	, TString("muon-jet #Delta#phi")+s	, 200, -M_PI, M_PI);
//		hqJetMinDeltaR[i]	= new TH1F(TString::Format("hJetsMinDeltaR%d",i)	, TString("muon-jet #Delta#r")+s	, 200, 0, 10);
		
	}
*/

//	TH1F* hJetDeltaEta;
//	TH1F* hJetDeltaPhi;
//	TH1F* hJetDeltaR;


}

// ------------ method called once each job just after ending the event loop  ------------
void 
QCDBackground::endJob() {
	_outFile	->cd();

	_outTree->Write();
	_outTreeMasses->Write();
	_outTreeMassesMuon->Write();

	hCountTriggersMu9->Write();
	hMuonP			->Write();
	hMuonPt			->Write();
	hMuonEta		->Write();
	hMuonPhi		->Write();
	hMuonD0			->Write();
	hMuonCharge		->Write();

	hJetP			->Write();
	hJetPt			->Write();
	hJetEta			->Write();
	hJetPhi			->Write();
	_h2PVsEta	->Write();	
	_h2PVsEta1	->Write();

	hJetP1			->Write();
	hJetPt1			->Write();
	hJetEta1		->Write();
	hJetPhi1		->Write();


	hJetP1Pos		->Write();
	hJetPt1Pos		->Write();
	hJetEta1Pos		->Write();
	hJetPhi1Pos		->Write();

	hJetP1Neg		->Write();
	hJetPt1Neg		->Write();
	hJetEta1Neg		->Write();
	hJetPhi1Neg		->Write();

	hJetDeltaR		->Write();
	hJetDeltaEta	->Write();	
	hJetDeltaPhi	->Write();	


	_h3PVsEtaVsPhi1->Write();
	_h3PVsEtaVsPhi2->Write();

/*

	hJetPtVsEta->Write();
	hJetPtVsEta1->Write();
	hNumSelJets		->Write();
	hNumSelMuons	->Write();
//	hJetDeltaEta	->Write();	
//	hJetDeltaPhi	->Write();	
//	hJetDeltaR		->Write();
	hJetMinDeltaEta	->Write();	
	hJetMinDeltaPhi	->Write();	
	hJetMinDeltaR	->Write();
	hDiJetDeltaPhi->Write();
	hDiJetDeltaPhiNoCharge->Write();
*/
/*	for (int i = 0; i < 3;++i){
		hqJetDeltaEta[i]	->Write();	
		hqJetDeltaPhi[i]	->Write();	
		hqJetDeltaR[i]		->Write();

//		hqJetMinDeltaEta[i]	->Write();	
//		hqJetMinDeltaPhi[i]	->Write();	
//		hqJetMinDeltaR[i]	->Write();
	}
*/
	hDiMuonMass->Write();
	hDiMuonMassOS->Write();
	hDiMuonPOS->Write();
	hDiMuonPtOS->Write();
	hDiMuonEtaOS->Write();
	hDiMuonPhiOS->Write();
	hDiMuonMassSS->Write();
	hDiMuonPSS->Write();
	hDiMuonPtSS->Write();
	hDiMuonEtaSS->Write();
	hDiMuonPhiSS->Write();

	hJetMassAll->Write();
	hJetMassAll1->Write();
	hJetMassNoCharge->Write();
	hJetMassCharged->Write();

	hJetPtVsMuonPt->Write();
	hJetPVsMuonP->Write();
	_outFile->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(QCDBackground);
