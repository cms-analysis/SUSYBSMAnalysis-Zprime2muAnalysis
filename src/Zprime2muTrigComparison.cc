#include <string>

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muTrigComparison.h"

using namespace std;

void Zprime2muTrigComparison::analyze(const edm::Event& event, 
				      const edm::EventSetup& eSetup) {
  // Delegate filling our muon vectors to the parent class.
  Zprime2muAnalysis::analyze(event, eSetup);

  compareTrigDecision(event);
}

bool Zprime2muTrigComparison::TriggerTranslator(const string& algo,
						const unsigned int lvl,
						const unsigned int nmu,
						ostream& out) const {
  // See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonHLT
  // https://twiki.cern.ch/twiki/bin/view/CMS/L1TriggerTableHLTExercise
  // L1 Quality codes are still the same, but we are suggested to use the 
  // methods useInSingleMuonTrigger() etc. instead of using the quality
  // code directly:
  // https://twiki.cern.ch/twiki/bin/view/CMS/GMTEmulator

  // Arrays are first-indexed by trigger level then by the number of
  // muons required.
  const unsigned int l = lvl - 1;
  const unsigned int n = nmu - 1;

  // 19/6/2008: CMSSW_2_1_X trigger table
  // these are for L = 2e30
  const double ptMin[3][2] = {
    { 7, 3},
    {12, 3},
    {15, 3}
  };
  const double etaMax[3] = { 2.5, 2.5, 2.5 };
  // 19/6/2008: in CMSSW_2_1_X, requirement is just that the vertex
  // be < 2 cm from the beam axis, i.e. sqrt(vx**2+vy**2) < 2.0
  const double dxy2Max[3] = { 9999, 9999, 2.0*2.0 };
  const double dzMax[3] = { 9999, 9999, 9999 };
  const double nSigmaPt[3] = { 0, 0, 0 };
  // 12/05/2004: a muon must have more than 5 pixel plus silicon hits in
  // the tracker (DAQ TDR, p. 306).
  //const int minHits[3] = { 0, 0, 6 };
  // 28/10/2007: for CMSSW_1_6_0, min hits is 0
  // 19/6/2008: same in 2_1_X
  const int minHits[3] = { 0, 0, 0 };

  // Check if the algorithm string passed in is one we recognized.
  // String value is not currently used, but will be used in the
  // future to differentiate between, e.g., Iso and NonIso triggers.
  if (algo != "L1_SingleMu7" &&
      algo != "L1_DoubleMu3" &&
      algo != "SingleMuNoIsoL2PreFiltered" && 
      algo != "SingleMuNoIsoL3PreFiltered" &&
      algo != "DiMuonNoIsoL2PreFiltered" && 
      algo != "DiMuonNoIsoL3PreFiltered") {
    edm::LogWarning("TriggerTranslator")
      << "+++ unrecognized algorithm " << algo << "! +++";
    return false;
  }

  out << "TriggerTranslator, " << algo << ":\n";

  unsigned int muonsPass = 0;
  //vector<zp2mu::Muon>::const_iterator pmu; //, pmu_prev;

  //vector<double> zvtx;
  //vector<double>::const_iterator pvtx;
  for (unsigned imu = 0; imu < allLeptons[lvl].size(); imu++) {
    const reco::CandidateBaseRef& mu = allLeptons[lvl][imu];
    reco::TrackRef tk, tktk;
    if (lvl < lL3)
      tk = mu->get<reco::TrackRef>();
    else
      tk = mu->get<reco::TrackRef, reco::CombinedMuonTag>();

    // JMTBAD HLTMuonPrefilter cuts not on pt but on what they call ptLx
    // nSigmaPt[lL1] = 0, so ptLx(lL1) = pt(lL1), but safeguard against
    // accidentally setting nSigmaPt[lL1] != 0:
    // 19/6/2008: in 2_1_X the ptLx cut is removed, but keep this code
    // around for good luck.
    double ptLx;
    if (lvl == lL1)
      ptLx = mu->pt();
    else
      ptLx = (1 + nSigmaPt[l]*invPError(&*tk)*mu->p())*mu->pt();

    out << "    mu #" << recLevelHelper.id(mu)
	<< " pt: " << mu->pt() << " ptLx: " << ptLx << " (cut: " << ptMin[l][n] 
	<< ") eta: " << mu->eta() << " (cut: " << etaMax[l] << ")\n";
    
    bool pass = ptLx >= ptMin[l][n] && fabs(mu->eta()) <= etaMax[l];

    // don't bother evaluating the other constraints if they aren't set
    if (dxy2Max[l] != 9999) {
      double vx = mu->vx();
      double vy = mu->vy();
      pass = pass && vx*vx + vy*vy < dxy2Max[l];
    }

    if (dzMax[l] != 9999) {
      pass = pass && fabs(mu->vz()) < dzMax[l];
    }

    if (minHits[l] != 0) {
      pass = pass &&
	tk->hitPattern().numberOfValidTrackerHits() >= minHits[l];
    }

    if (pass) {
      // here check extra stuff depending on the trigger algorithm
      bool passExtra = true;

      if (lvl == lL1) {
	int quality
	  = toConcrete<l1extra::L1MuonParticle>(mu).gmtMuonCand().quality();
	// In single muon trigger, GMT uses only muons with quality > 3.
	// In dimuon trigger, GMT uses muons with quality = 3 and 5-7.
	if ((nmu == 1 && quality < 4) ||
	    (nmu == 2 && (quality < 3 || quality == 4)))
	  passExtra = false;
      }
      /*
      // JMTBAD disabled for now, not currently in HLT
      else if (lvl == lL3 && nmu == 2) {
        for (unsigned jmu = 0; jmu < allLeptons[lvl].size(); jmu++) {
	  const reco::CandidateBaseRef& mu_prev = allLeptons[lvl][jmu];
          // Skip ghost tracks (see p. 308 of DAQ TDR)
	  if (fabs(mu->eta() - mu_prev->eta()) < 0.01 &&
	      fabs(mu->phi() - mu_prev->phi()) < 0.05 &&
	      fabs(mu->pt()  - mu_prev->pt())  < 0.1) passExtra = false;
	}
      }
      */

      if (passExtra)
	muonsPass++;
    }

    if (muonsPass == nmu)
      break;
  }

  bool result = muonsPass >= nmu;

  out << "  " << algo << ": ";
  if (result) out << "pass!";
  else        out << "fail!";
  out << endl;

  return result;
}

void Zprime2muTrigComparison::compareTrigDecision(const edm::Event& event) const {
  const unsigned nPaths = 2;
  static const string trigNames[3][nPaths] = {
    { "L1_SingleMu7", "L1_DoubleMu3" },
    { "SingleMuNoIsoL2PreFiltered", "DiMuonNoIsoL2PreFiltered" },
    { "SingleMuNoIsoL3PreFiltered", "DiMuonNoIsoL3PreFiltered" }
  };

  ostringstream out, xout;
  out << "Homemade decisions:\n";

  for (unsigned int lvl = lL1; lvl <= lL3; lvl++) {
    unsigned homemade_trigbits = 0;

    for (unsigned ipath = 0; ipath < nPaths; ipath++) {
      const string& trigName = trigNames[lvl - lL1][ipath];

      // Try to emulate HLT algorithms.
      bool fired = TriggerTranslator(trigName, lvl, ipath+1, xout);
      
      // If the event passes, set the corresponding bit in trigbits.
      if (fired) homemade_trigbits = homemade_trigbits | (1 << ipath);
      
      out << "  " << setw(35) << trigName << ": decision = " << fired << endl;
    }

    // "Official" muon HLTs are calculated only when corresponding
    // previous levels gave OK, while we calculate the decision for a
    // given level regardless of previous levels' decisions.
    if (lvl >= lL2)
      homemade_trigbits &= trigDecision.getWord(lL1);
    if (lvl >= lL3)
      homemade_trigbits &= trigDecision.getWord(lL2);

    // Compare official and homemade decisions.
    // JMTBAD this warning will also be fired about L2 when running on
    // AOD and on events for which L2 trigbits do not equal L3
    // ones... see storeHLTDecision()
    if (homemade_trigbits != trigDecision.getWord(lvl)) {
      edm::LogWarning("compareTrigDecision")
	<< "+++ Warning: official L" << lvl
	<< " decision disagrees with homemade decision:"
	<< " official: " << trigDecision.getWord(lvl)
	<< " homemade: " << homemade_trigbits << " +++";
      
      if (verbosity == VERBOSITY_NONE)
	dumpEvent(true);
    }
  }

  edm::LogInfo("compareTrigDecision") << xout.str() << endl << out.str();
}

DEFINE_FWK_MODULE(Zprime2muTrigComparison);
