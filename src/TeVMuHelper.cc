#include <ostream>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TeVMuHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

using namespace reco;
using namespace edm;
using namespace std;

const unsigned TeVMuHelper::tevMuCuts = TeVMuHelper::PT | TeVMuHelper::ISO;
const unsigned TeVMuHelper::heepCuts  = TeVMuHelper::PTOTHER;

const unsigned TeVMuHelper::topSkim      = TeVMuHelper::PT; // | TeVMuHelper::TOPTRIGGER;
const unsigned TeVMuHelper::topElCuts    = TeVMuHelper::ELTIGHT | TeVMuHelper::D0EL | TeVMuHelper::COLLEMU;
const unsigned TeVMuHelper::topMuCuts    = TeVMuHelper::D0 | TeVMuHelper::NSIHITS | TeVMuHelper::CHI2DOF;
const unsigned TeVMuHelper::topLeptonId  = TeVMuHelper::topElCuts | TeVMuHelper::topMuCuts;
const unsigned TeVMuHelper::topIsolation = TeVMuHelper::ISOS;
const unsigned TeVMuHelper::topLepCuts   = TeVMuHelper::topSkim | TeVMuHelper::topLeptonId | TeVMuHelper::topIsolation;

//TeVMuHelper::CutNameMap TeVMuHelper::cutNames;

const bool debug = false;

// JMTBAD don't hardcode inputtag names

TeVMuHelper::TeVMuHelper() :
  _cutMask(0),
  // should take these values from the cfg
  ptCut(20),
  isoCut(10),
  ptOthCut(80),
  chi2dofCut(5),
  d0Cut(0.25),
  d0ElCut(0.04),
  nSiHitsCut(7),
  collinearMuon_dRmax(0.1),
  Sratio_min(0.92)
{
  /*
  if (cutNames.empty()) {
    cutNames[PT]      	 = "PT";
    cutNames[ISO]     	 = "ISO";
    cutNames[CHI2DOF] 	 = "CHI2DOF";
    cutNames[D0]      	 = "D0";
    cutNames[NSIHITS] 	 = "NSIHITS";
    cutNames[PTOTHER] 	 = "PTOTHER";
    cutNames[D0EL]    	 = "D0EL";
    cutNames[COLLEMU] 	 = "COLLEMU";
    cutNames[ISOS]    	 = "ISOS";
    cutNames[ELTIGHT] 	 = "ELTIGHT";
    cutNames[TOPMET]  	 = "TOPMET";
    cutNames[ZVETO]   	 = "ZVETO";
    cutNames[NJETS2]  	 = "NJETS2";
    cutNames[NJETS01] 	 = "NJETS01";
    cutNames[WEAKISO] 	 = "WEAKISO";
    cutNames[TOPTRIGGER] = "TOPTRIGGER";
    }*/
}

void TeVMuHelper::initEvent(const Event& evt) {
  event = &evt;
  METx_uncorrJES = METy_uncorrJES = MET_uncorrJES = -1;
  METx = METy = MET = -1;
  Zvetoed = -1;
  nJets = -1;
  topTriggered = -1;
}

void TeVMuHelper::_eventOK() const {
  if (event == 0)
    throw cms::Exception("TeVMuHelper") << "event object needed but not set!\n";
}
  
void TeVMuHelper::_cacheTrigger() {
  _eventOK();

  edm::Handle<edm::TriggerResults> hltRes;
  event->getByLabel(InputTag("TriggerResults::HLT"), hltRes);
  edm::TriggerNames hltTrigNames;
  hltTrigNames.init(*hltRes);

  vector<string> hltPaths;
  // JMTBAD fix paths
  //hltPaths.push_back("HLT1MuonNonIso");
  //hltPaths.push_back("HLT1ElectronRelaxed");

  bool fired = false;
  for (unsigned i = 0; i < hltPaths.size(); ++i) {
    int ndx = hltTrigNames.triggerIndex(hltPaths[i]);
    fired = fired || hltRes->accept(ndx);
    if (fired) break;
  }

  topTriggered = fired ? 1 : 0;
}

void TeVMuHelper::_cacheMET() {
  _eventOK();

  Handle<View<pat::METType> > METs;
  event->getByLabel("selectedLayer1METs", METs);
  const pat::MET* theMET = toConcretePtr<pat::MET>(METs->at(0));

  // Top group MET cut: do not correct for JES, but correct for
  // muons. They only correct for muons that pass the cuts, while
  // the default MET muon correction uses all muons. They say it is
  // a negligible difference, though, so try using the default.
  METx_uncorrJES = theMET->corEx(pat::MET::uncorrJES);
  METy_uncorrJES = theMET->corEy(pat::MET::uncorrJES);
  MET_uncorrJES  = sqrt(METx_uncorrJES*METx_uncorrJES + METy_uncorrJES*METy_uncorrJES);

  // Go ahead and store the JES-corrected (default) MET, but it isn't
  // used yet.
  METx = theMET->px();
  METy = theMET->py();
  MET  = theMET->pt();
}

void TeVMuHelper::_cacheZvetoed() {
  _eventOK();

  Zvetoed = 0;

  const char* leptonNames[2] = {
    "selectedLayer1Electrons",
    "selectedLayer1Muons"
  };

  for (unsigned z = 0; z < 2; z++) {
    Handle<View<Candidate> > leptons;
    event->getByLabel(leptonNames[z], leptons);
    const unsigned nlep = leptons->size();
    
    for (unsigned i = 0; i < nlep; i++) {
      const Candidate& ilep = (*leptons)[i];

      if (leptonIsCut(ilep, topLepCuts)) continue;
      
      for (unsigned j = i+1; j < nlep; j++) {
	const Candidate& jlep = (*leptons)[j];
	if (leptonIsCut(jlep, topLepCuts)) continue;
	
	// no charge check yet

	double mass = (ilep.p4() + jlep.p4()).mass();
	if (mass > 76 && mass < 106) {
	  Zvetoed = (((1 << i) | (1 << j)) << 1) | z; // an incredibly cryptic code for which pair of leptons caused the veto
	  goto done; // goto-considered-useful
	}
      }
    }
  }
 done:
  ;
}

bool TeVMuHelper::collinearMuon(const GsfElectron* electron) const {
  _eventOK();

  Handle<View<Candidate> > muons;
  event->getByLabel("selectedLayer1Muons", muons);

  for (unsigned i = 0; i < muons->size(); i++) {
    if (reco::deltaR((*muons)[i], *electron) < collinearMuon_dRmax)
      return true;
  }

  return false;
}

double TeVMuHelper::calcIsolationS(const Muon* muon) const {
  if (!muon->isIsolationValid())
    throw cms::Exception("TeVMuHelper") << "muon isolation not valid in passIsolationS()!\n";

  const MuonIsolation& iso = muon->isolationR03();
  return iso.sumPt + iso.hadEt + iso.hoEt + iso.emEt;
}

double TeVMuHelper::calcIsolationS(const GsfElectron* electron) const {
  _eventOK();

  double S = 0;
  const GsfTrackRef& elTk = electron->gsfTrack();

  Handle<TrackCollection> tracks;
  event->getByLabel("ctfWithMaterialTracks", tracks);

  TrackCollection::const_iterator tk = tracks->begin();
  for ( ; tk != tracks->end(); ++tk) {
    double pt = tk->pt();
    if (pt > 1 &&
	fabs(tk->d0()) < 0.1 &&
	tk->numberOfValidHits() >= 7 &&
	fabs(elTk->vz() - tk->vz()) < 0.5) {

      double dR = reco::deltaR(*tk, *elTk);
      if (dR > 0.01 && dR < 0.3) // guess at a veto cone
	S += pt;
    }
  }
  
  return S;
}

bool TeVMuHelper::passIsolationS(double S, double pt) const {
  return pt/(pt+S) > Sratio_min;
}

void TeVMuHelper::_cacheNJets() {
  _eventOK();

  Handle<PatElectronCollection> electrons;
  event->getByLabel("selectedLayer1Electrons", electrons);

  Handle<View<Candidate> > jets;
  event->getByLabel("selectedLayer1Jets", jets);

  nJets = 0;

  for (unsigned j = 0; j < jets->size(); j++) {
    const CandidateBaseRef& jet = jets->refAt(j);
    if (jet->pt() > 30 && fabs(jet->eta()) < 2.4) {
      bool isEl = false;

      PatElectronCollection::const_iterator el = electrons->begin();
      for ( ; el != electrons->end(); ++el) {
	if (reco::deltaR(*el, *jet) < 0.3 && !electronIsCut(&*el, topElCuts | topIsolation)) {
	  isEl = true;
	  break;
	}
      }

      if (!isEl) nJets++;
    }
  }
}

unsigned TeVMuHelper::electronIsCut(const GsfElectron* electron,
				    unsigned cutMask) const {
  unsigned result = 0;
  
  if (cutMask == 0)
    cutMask = _cutMask;

  if ((cutMask & ISOS) && !passIsolationS(calcIsolationS(electron), electron->pt()))
    result |= ISOS;

  if (cutMask & ELTIGHT) {
    const pat::Electron* patEl = toConcretePtr<pat::Electron>(*((Candidate*)electron));
    if (patEl != 0 && patEl->leptonID("tight") == 0)
      result |= ELTIGHT;
  }

  const GsfTrackRef& tk = electron->gsfTrack();
  if (tk.isNonnull()) {
    // Electron d0 cut.
    if ((cutMask & D0EL) && fabs(tk->d0()) > d0ElCut)
      result |= D0EL;
  }
  
  // If any muon is found within a certain dR of the electron, cut.
  if ((cutMask & COLLEMU) && collinearMuon(electron))
    result |= COLLEMU;

  return result;
}

unsigned TeVMuHelper::muonIsCut(const Muon* muon,
				unsigned cutMask) const {
  unsigned result = 0;

  if (cutMask == 0)
    cutMask = _cutMask;

  // Sum pT in cone of dR=0.3 cut (with veto cone of dR=0.01), using
  // all tracks in the isolation definition.
  if ((cutMask & ISO) && muon->isIsolationValid()) {
    double iso = muon->isolationR03().sumPt;
    if (iso > isoCut)
      result |= ISO;
  }

  if ((cutMask & ISOS) && !passIsolationS(calcIsolationS(muon), muon->pt()))
    result |= ISOS;

  const TrackRef& tk = muon->combinedMuon();
  if (tk.isNonnull()) {
    // Cut on chi2/dof.
    if ((cutMask & CHI2DOF) && tk->normalizedChi2() > chi2dofCut)
      result |= CHI2DOF;

    // Cut on d0.
    if ((cutMask & D0) && fabs(tk->d0()) > d0Cut)
      result |= D0;
  }

  // Cut on number of silicon hits for the tracker track.
  const TrackRef& tktk = muon->track();
  if ((cutMask & NSIHITS) && tktk.isNonnull() &&
      tktk->numberOfValidHits() < nSiHitsCut)
    result |= NSIHITS;

  return result & cutMask;
}

unsigned TeVMuHelper::leptonIsCut(const Candidate& lepton,
				  unsigned cutMask) const {
  unsigned result = 0;

  if (cutMask == 0)
    cutMask = _cutMask;

  if ((cutMask & PT) && lepton.pt() <= ptCut)
    result |= PT;

  if ((cutMask & PTOTHER) && lepton.pt() <= ptOthCut)
    result |= PTOTHER;

  int pdgId = abs(lepton.pdgId());

  // Try electron cuts.
  if (pdgId == 11) {
    const GsfElectron* electron = toConcretePtr<GsfElectron>(lepton);
    if (electron != 0)
      result |= electronIsCut(electron, cutMask);
    //else
    //  throw cms::Exception("TeVMuHelper") << "unable to cast lepton of class " << typeid(lepton).name() << " to electron!\n";
  }

  // Try muon cuts.
  if (pdgId == 13) {
    const Muon* muon = toConcretePtr<Muon>(lepton);
    if (muon != 0)
      result |= muonIsCut(muon, cutMask);
    //else
    //  throw cms::Exception("TeVMuHelper") << "unable to cast lepton of class " << typeid(lepton).name() << " to muon!\n";
  }

  return result & cutMask;
}

unsigned TeVMuHelper::dileptonIsCut(const CompositeCandidate& dil,
				    unsigned cutMask) {
  unsigned result = 0;

  if (cutMask == 0)
    cutMask = _cutMask;

  // Top group trigger requirement.
  if (cutMask & TOPTRIGGER) {
    if (topTriggered < 0)
      _cacheTrigger();
    if (topTriggered == 0)
      result |= TOPTRIGGER;
  }

  // Top group's MET cut.
  if (cutMask & TOPMET) {
    if (MET_uncorrJES < 0) _cacheMET(); // i.e. not calculated yet for this event
	  
    int id0 = abs(dileptonDaughter(dil, 0)->pdgId());
    int id1 = abs(dileptonDaughter(dil, 1)->pdgId());

    // same flavor: MET > 30 AND (MET > 0.6 * dil pT OR alpha > 0.25)
    // opp flavor: MET > 20
    if ((id0 != id1 && MET < 20) || (id0 == id1 && MET < 30))
      result |= TOPMET;
    else if (id0 == id1 && MET < 0.6*dil.pt()) {
      // Define the alpha they do in the figure.
      double alpha = (-dil.px()*METx_uncorrJES - dil.py()*METy_uncorrJES) / MET / dil.pt();
      if (alpha < 0.25)
	result |= TOPMET;
    }
  }

  // Top group's Z-veto.
  if (cutMask & ZVETO) {
    if (Zvetoed < 0)
      _cacheZvetoed(); // i.e. not calculated yet for this event
    if (Zvetoed > 0)
      result |= ZVETO;
  }

  // Top group's jet count.
  if (cutMask & NJETS2) {
    if (nJets < 0)
      _cacheNJets();

    if (nJets < 2)
      result |= NJETS2;
  }

  if (cutMask & NJETS01) {
    if (nJets < 0)
      _cacheNJets();

    if (nJets >= 2)
      result |= NJETS01;
  }

  // Any lepton cuts.
  if (cutMask & LEPCUTS) {
    unsigned res0 = leptonIsCut(*dileptonDaughter(dil, 0), cutMask);
    unsigned res1 = leptonIsCut(*dileptonDaughter(dil, 1), cutMask);
    unsigned lepres = res0 | res1;
    
    if (cutMask & WEAKISO) {
      static const unsigned allIso = TeVMuHelper::ISO | TeVMuHelper::ISOS;
      // If weakIso, then allow the dilepton to have one lepton non-isolated.
      
      // If exactly one of cut0 or cut1 has an isolation bit set
      // (assumes that only ISO or only ISOS is in the cutmask to
      // begin with), turn off the isolation cut bit in cut.
      lepres &= ~((res0 & allIso) ^ (res1 & allIso));
    }

    result |= lepres;
  }

  result &= cutMask;

  /*
  if (debug) {
    _eventOK();

    ostringstream out;
    out << "dileptonIsCut(): cutting for:  ";
    for (CutNameMap::const_iterator it = cutNames.begin(); it != cutNames.end(); ++it)
      if (it->first & cutMask) out << it->second << " ";
    out << endl;
    
    out << "METx_uncorrJES: " << METx_uncorrJES << " METy_uncorrJES: " << METy_uncorrJES
	<< " MET_uncorrJES: " << MET_uncorrJES << " Zvetoed: " << Zvetoed << " nJets: " << nJets << endl;

    out << "Leptons in event:\n";

    Handle<PatElectronCollection> electrons;
    event->getByLabel("selectedLayer1Electrons", electrons);
    out << " electrons (" << electrons->size() << "):\n";
    for (PatElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); ++it)
      out << "  E: " << it->energy() << " pt: " << it->pt() << " eta: " << it->eta() << " phi: " << it->phi() << " pid: " << it->pdgId()
	  << " q: " << it->charge() << " isoS: " << calcIsolationS((GsfElectron*)(&*it)) << endl;

    Handle<View<Candidate> muons;
    event->getByLabel("selectedLayer1Muons", muons);
    out << " muons (" << muons->size() << "):\n";
    for (MuonCollection::const_iterator it = muons->begin(); it != muons->end(); ++it)
      out << "  E: " << it->energy() << " pt: " << it->pt() << " eta: " << it->eta() << " phi: " << it->phi() << " pid: " << it->pdgId()
	  << " q: " << it->charge() << " isoS: " << calcIsolationS(&*it) << endl;

    out << "dilepton pt, eta, phi, mass, px, py: "
	<< dil.pt() << " " << dil.eta() << " " << dil.phi() << " " << dil.mass() << " " << dil.px() << " " << dil.py()
	<< "\ndaughter leptons: #, pt, eta, phi, pdgid, charge\n";
    const CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
    out << lep0->pt() << " " << lep0->eta() << " " << lep0->phi() << " " << lep0->pdgId() << " " << lep0->charge() << endl;
    const CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);
    out << lep1->pt() << " " << lep1->eta() << " " << lep1->phi() << " " << lep1->pdgId() << " " << lep1->charge() << endl;
    if (result) {
      out << "was CUT for: ";
      for (CutNameMap::const_iterator it = cutNames.begin(); it != cutNames.end(); ++it)
	if (it->first & result) out << it->second << " ";
      out << endl;
    }
    else
      out << "was NOT cut!\n";
    out << endl;

    LogInfo("TeVMuHelper") << out.str();
  }
  */

  return result;
}

/*
double otherIsolation(const Muon& muon,
		      const double dR,
		      const double dRveto,
		      const double minPt) {
  const TrackRef& trk = muon.combinedMuon();

  // Get the map taking global muon combined tracks to tracker iso
  // deposits.
  Handle<MuIsoDepositAssociationMap> trkIso;
  event->getByLabel("muGlobalIsoDepositCtfTk", trkIso);
  
  // Try to get the MuIsoDeposit.
  MuIsoDeposit depTrk;
  try {
    depTrk  = (*trkIso)[trk];
  } catch (const cms::Exception& e) {
    // If anything goes wrong, return 0 so we don't crash with an
    // exception.
    edm::LogWarning("TeVMuHelper") << "couldn't get MuIsoDeposit! muon::isIsolationValid(): " << muon.isIsolationValid();
    return 0;
  }

  if (minPt == 1.5)
    LogInfo("TMH") << "otherIsolation(): dR=" << dR << " dRveto=" << dRveto << " minPt=" << minPt << " muon pt: " << muon.pt() << " MID::print(): " << depTrk.print();

  // Set up a veto around the muon track with the requested size.
  MuIsoDeposit::Veto veto;
  veto.vetoDir = MuIsoDeposit::Direction(trk->eta(), trk->phi());
  veto.dR = dRveto;
  const MuIsoDeposit::Vetos vetos(1, veto);

  // The pair is (sumPt, nTracks) in the annulus.
  std::pair<double, int> dep =
    depTrk.depositAndCountWithin(dR, vetos, minPt);

  return dep.first;
}
*/
