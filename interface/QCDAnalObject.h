#ifndef QCDAnalObject_h_
#define QCDAnalObject_h_
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

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

#endif
