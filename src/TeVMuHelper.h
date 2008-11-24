#ifndef TEVMUHELPER_H
#define TEVMUHELPER_H

#include <map>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/Event.h"

typedef std::vector<pat::Electron> PatElectronCollection;

class TeVMuHelper {
public:
  // To go in bitfields.
  enum CutResult {
    PASS=0,
    // Least significant word is for per-lepton cuts:
    PT=0x01, ISO=0x02, CHI2DOF=0x04, D0=0x08,
    NSIHITS=0x10, PTOTHER=0x20, D0EL=0x40, COLLEMU=0x80,
    ISOS=0x100, ELTIGHT=0x200,
    LEPCUTS=0xFFFF,
    // Most significant word is for per-event or per-dilepton cuts:
    TOPMET=0x10000, ZVETO=0x20000, NJETS2=0x40000, NJETS01=0x80000,
    WEAKISO=0x100000, TOPTRIGGER=0x200000,
    EVTCUTS=0xFFFF0000,
    ALLCUTS=0xFFFFFFFF
  };

  // Some defined cut choices.

  // AN 2007/038 cuts (pT and isolation)
  static const unsigned tevMuCuts;

  // AN 2008/044 e mu bkg study broken cuts.
  static const unsigned heepCuts;
  
  // Various combinations of AN 2008/015 cuts ("top group").
  static const unsigned topSkim;
  static const unsigned topElCuts;   
  static const unsigned topMuCuts;
  static const unsigned topLeptonId;
  static const unsigned topIsolation;
  static const unsigned topLepCuts;
  static const unsigned topEventCuts;
  static const unsigned topGroupCuts;

  typedef std::map<unsigned, const char*> CutNameMap;
  static CutNameMap cutNames;

  TeVMuHelper();
  void initEvent(const edm::Event& evt);
  void setCutMask(const unsigned mask) { _cutMask = mask; }

  bool collinearMuon(const reco::PixelMatchGsfElectron* electron) const;

  double calcIsolationS(const reco::Muon* muon) const;
  double calcIsolationS(const reco::PixelMatchGsfElectron* electron) const;
  bool passIsolationS(double S, double pt) const;

  unsigned electronIsCut(const reco::PixelMatchGsfElectron* electron,
			 unsigned cutMask=0) const;
  unsigned muonIsCut(const reco::Muon* muon,
		     unsigned cutMask=0) const;
  unsigned leptonIsCut(const reco::Candidate& lepton,
		       unsigned cutMask=0) const;

  unsigned dileptonIsCut(const reco::CompositeCandidate& dil,
			 unsigned cutMask=0);

 private:
  void _eventOK() const;

  void _cacheTrigger();
  void _cacheMET();
  void _cacheZvetoed();
  void _cacheNJets();

  const edm::Event* event;

  unsigned _cutMask;

  // Cut values.
  double ptCut;
  double isoCut;
  double ptOthCut;
  double chi2dofCut;
  double d0Cut;
  double d0ElCut;
  int nSiHitsCut;
  double collinearMuon_dRmax;
  double Sratio_min;

  // Per-event quantity cache.
  double METx_uncorrJES;
  double METy_uncorrJES;
  double MET_uncorrJES;
  double METx;
  double METy;
  double MET;
  int Zvetoed;
  int nJets;
  int topTriggered;
};

#endif
