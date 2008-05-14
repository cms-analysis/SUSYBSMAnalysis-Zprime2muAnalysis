//
// Authors: Jason Mumford, Jordan Tucker, Slava Valuev, UCLA
//
// TODO:
//   replace strings with "muon" to "lepton"
//   appropriate eta values for electrons in plots

#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TPad.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TText.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muResolution.h"

using namespace std;

Zprime2muResolution::Zprime2muResolution(const edm::ParameterSet& config) 
  : Zprime2muAnalysis(config), histoFile(0)
{
  outputFile = config.getUntrackedParameter<string>("outputFile");

  // get the parameters specific to the data sample on which we are running
  string dataSet = config.getParameter<string>("dataSet");
  edm::ParameterSet dataSetConfig =
    config.getParameter<edm::ParameterSet>(dataSet);

  peakMass = dataSetConfig.getParameter<double>("peakMass");
  lowerMassWin = dataSetConfig.getParameter<double>("lowerMassWin");
  upperMassWin = dataSetConfig.getParameter<double>("upperMassWin");
  binSize = dataSetConfig.getParameter<int>("binSize");
  maxTrigMass = dataSetConfig.getParameter<double>("maxTrigMass");

  string fn = config.getUntrackedParameter<string>("histoFile");
  useHistosFromFile = config.getUntrackedParameter<bool>("useHistosFromFile");
  if (useHistosFromFile) {
    histoFile = new TFile(fn.c_str(), "READ");
    getHistosFromFile();
  }
  else {
    BookResHistos();
    BookEffHistos();
    BookPtResHistos();
    BookDilResHistos();
    BookChargeResHistos();
  }
}

Zprime2muResolution::~Zprime2muResolution() {
  if (useHistosFromFile) {
    histoFile->Close();
    delete histoFile;
  }
}

void Zprime2muResolution::beginJob(const edm::EventSetup& eSetup) {
  if (useHistosFromFile) {
    DrawResHistos();
    // We're done, so stop running and don't waste time reading events.
    // JMTBAD Is there a better way to stop running than to throw an
    // exception?
    throw cms::Exception("Zprime2muResolution") 
      << "using histos from file; done.\n";
  }
}

void Zprime2muResolution::analyze(const edm::Event& event, 
				  const edm::EventSetup& eSetup) {
  // Delegate filling our muon vectors to the parent class.
  Zprime2muAnalysis::analyze(event, eSetup);

  // Fill resolution histos.
  calcResolution();
}

void Zprime2muResolution::endJob() {
  Zprime2muAnalysis::endJob();

  DrawResHistos();
  //DrawAcceptance();
}

void Zprime2muResolution::BookResHistos() {
  const double pi = TMath::Pi();

  // Acceptance studies.
  GenMassAllEvents = fs->make<TH1F>("GenMassAllEvents",
				    "Gen mass, all events", 24, 200., 5000.);
  GenMassInAccept  = fs->make<TH1F>("GenMassInAccept",
				    "Gen mass, events in acceptance",
				    24, 200., 5000.);
  GenMassAllEvents->Sumw2();  
  GenMassInAccept->Sumw2();

  // Origin of muons
  Origin[0] = fs->make<TH1F>("Origin0","Particle Id of Mother of all mu's", 20, 0, 20);
  Origin[1] = fs->make<TH1F>("Origin1","Particle Id of Mother of opp-sign dilepton mu's",
			     20, 0, 20);

  // Trigger decisions
  const int maxTrigBits = 8;
  TrigResult[0][0] = fs->make<TH1F>("TrigResult00", "Gen, all events",        maxTrigBits, 0., double(maxTrigBits));
  TrigResult[1][0] = fs->make<TH1F>("TrigResult10", "L1 trigger, all events", maxTrigBits, 0., double(maxTrigBits));
  TrigResult[2][0] = fs->make<TH1F>("TrigResult20", "L2 trigger, all events", maxTrigBits, 0., double(maxTrigBits));
  TrigResult[3][0] = fs->make<TH1F>("TrigResult30", "L3 trigger, all events", maxTrigBits, 0., double(maxTrigBits));
  TrigResult[0][1] = fs->make<TH1F>("TrigResult01", "Gen, #eta < 2.4",        maxTrigBits, 0., double(maxTrigBits));
  TrigResult[1][1] = fs->make<TH1F>("TrigResult11","L1 trigger, #eta < 2.4",  maxTrigBits, 0., double(maxTrigBits));
  TrigResult[2][1] = fs->make<TH1F>("TrigResult21","L2 trigger, #eta < 2.4",  maxTrigBits, 0., double(maxTrigBits));
  TrigResult[3][1] = fs->make<TH1F>("TrigResult31","L3 trigger, #eta < 2.4",  maxTrigBits, 0., double(maxTrigBits));
  TrigResult[0][2] = fs->make<TH1F>("TrigResult02", "Gen, #eta < 2.1",        maxTrigBits, 0., double(maxTrigBits));
  TrigResult[1][2] = fs->make<TH1F>("TrigResult12","L1 trigger, #eta < 2.1",  maxTrigBits, 0., double(maxTrigBits));
  TrigResult[2][2] = fs->make<TH1F>("TrigResult22","L2 trigger, #eta < 2.1",  maxTrigBits, 0., double(maxTrigBits));
  TrigResult[3][2] = fs->make<TH1F>("TrigResult32","L3 trigger, #eta < 2.1",  maxTrigBits, 0., double(maxTrigBits));

  TrigMass[0][0] = fs->make<TH1F>("TrigMass00", "Gen mass, Gen, all events", 27, 0., maxTrigMass);
  TrigMass[1][0] = fs->make<TH1F>("TrigMass10", "Gen mass, L1, all events",  27, 0., maxTrigMass);
  TrigMass[2][0] = fs->make<TH1F>("TrigMass20", "Gen mass, L2, all events",  27, 0., maxTrigMass);
  TrigMass[3][0] = fs->make<TH1F>("TrigMass30", "Gen mass, L3, all events",  27, 0., maxTrigMass);
  TrigMass[0][1] = fs->make<TH1F>("TrigMass01", "Gen mass, Gen, #eta < 2.4", 27, 0., maxTrigMass);
  TrigMass[1][1] = fs->make<TH1F>("TrigMass11", "Gen mass, L1, #eta < 2.4",  27, 0., maxTrigMass);
  TrigMass[2][1] = fs->make<TH1F>("TrigMass21", "Gen mass, L2, #eta < 2.4",  27, 0., maxTrigMass);
  TrigMass[3][1] = fs->make<TH1F>("TrigMass31", "Gen mass, L3, #eta < 2.4",  27, 0., maxTrigMass);
  TrigMass[0][2] = fs->make<TH1F>("TrigMass02", "Gen mass, Gen, #eta < 2.1", 27, 0., maxTrigMass);
  TrigMass[1][2] = fs->make<TH1F>("TrigMass12", "Gen mass, L1, #eta < 2.1",  27, 0., maxTrigMass);
  TrigMass[2][2] = fs->make<TH1F>("TrigMass22", "Gen mass, L2, #eta < 2.1",  27, 0., maxTrigMass);
  TrigMass[3][2] = fs->make<TH1F>("TrigMass32", "Gen mass, L3, #eta < 2.1",  27, 0., maxTrigMass);

  // Define errors
  for (int i_rec = 0; i_rec < NUM_REC_LEVELS; i_rec++) {
    for (int j = 0; j < 3; j++) {
      TrigMass[i_rec][j]->Sumw2();
    }
  }

  EventsInAccFailed  = fs->make<TH1F>("EventsInAccFailed", "Events in acc., failed", 5, 0., 5.);
  L1TrigFailSingleMu = fs->make<TH1F>("L1TrigFailSingleMu", "pT, failed 1mu L1", 25, 0., 25.);
  L1TrigFailMu2VsMu1 = fs->make<TH2F>("L1TrigFailMu2VsMu1", "pT mu2 vs pT mu1, failed di-mu L1",
				      25, 0., 25., 25, 0., 25.);
  L1TrigPassSingleMu = fs->make<TH1F>("L1TrigPassSingleMu", "pT, passed 1mu L1", 25, 0., 25.);
  L1TrigPassMu2VsMu1 = fs->make<TH2F>("L1TrigPassMu2VsMu1", "pT mu2 vs pT mu1, passed di-mu L1",
				      25, 0., 25., 25, 0., 25.);

  L2MuonHits    = fs->make<TH1F>("L2MuonHits", "Muon hits, Level-2", 25, 0., 25.);
  L3TrackerHits = fs->make<TH1F>("L3TrackerHits", "Tracker hits (pixel+silicon), Level-3",
				 25, 0., 25.);
  GMRMuonHits[0]= fs->make<TH1F>("GMRMuonHits0", "Muon hits, GMR, barrel",  25, 0., 25.);
  GMRMuonHits[1]= fs->make<TH1F>("GMRMuonHits1", "Muon hits, GMR, overlap", 25, 0., 25.);
  GMRMuonHits[2]= fs->make<TH1F>("GMRMuonHits2", "Muon hits, GMR, endcap",  25, 0., 25.);
  GMRChi2dof[0] = fs->make<TH1F>("GMRChi2dof0", "Chi2/d.o.f., GMR, barrel",  100, 0., 5.);
  GMRChi2dof[1] = fs->make<TH1F>("GMRChi2dof1", "Chi2/d.o.f., GMR, overlap", 100, 0., 5.);
  GMRChi2dof[2] = fs->make<TH1F>("GMRChi2dof2", "Chi2/d.o.f., GMR, endcap",  100, 0., 5.);

  // Muons per event
  NMuons[0][0] = fs->make<TH1F>("NMuons00", "Muons/Event, Gen", 10, 0, 10);
  NMuons[1][0] = fs->make<TH1F>("NMuons10", "Muons/Event, L1",  10, 0, 10);
  NMuons[2][0] = fs->make<TH1F>("NMuons20", "Muons/Event, L2",  10, 0, 10);
  NMuons[3][0] = fs->make<TH1F>("NMuons30", "Muons/Event, L3",  10, 0, 10);
  NMuons[0][1] = fs->make<TH1F>("NMuons01", "Muons/Event (w/ opp-sign dilepton), Gen",
				10, 0, 10);
  NMuons[1][1] = fs->make<TH1F>("NMuons11", "Muons/Event (w/ opp-sign dilepton), L1",
				10, 0, 10);
  NMuons[2][1] = fs->make<TH1F>("NMuons21", "Muons/Event (w/ opp-sign dilepton), L2",
				10, 0, 10);
  NMuons[3][1] = fs->make<TH1F>("NMuons31", "Muons/Event (w/ opp-sign dilepton), L3",
				10, 0, 10);
  NMuons[0][2] = fs->make<TH1F>("NMuons02", "Muons/Event, Gen, triggered", 10, 0., 10.);
  NMuons[1][2] = fs->make<TH1F>("NMuons12", "Muons/Event, L1, triggered",  10, 0., 10.);
  NMuons[2][2] = fs->make<TH1F>("NMuons22", "Muons/Event, L2, triggered",  10, 0., 10.);
  NMuons[3][2] = fs->make<TH1F>("NMuons32", "Muons/Event, L3, triggered",  10, 0., 10.);
  NMuons[0][3] = fs->make<TH1F>("NMuons03", "Muons/Event in acc., Gen, failed", 10, 0., 10.);
  NMuons[1][3] = fs->make<TH1F>("NMuons13", "Muons/Event in acc., L1, failed",  10, 0., 10.);
  NMuons[2][3] = fs->make<TH1F>("NMuons23", "Muons/Event in acc., L2, failed",  10, 0., 10.);
  NMuons[3][3] = fs->make<TH1F>("NMuons33", "Muons/Event in acc., L3, failed",  10, 0., 10.);
  NumDilVsRec  = fs->make<TH1F>("NumDilVsRec", "Events w/ opp-sign dilepton vs Rec Level",
				5, 0,  5);
  SignOfDil[0] = fs->make<TH1F>("SignOfDil0", "Number of opp-sign dileptons, Gen",
				2, 0,  2);
  SignOfDil[1] = fs->make<TH1F>("SignOfDil1", "Number of opp-sign, like-sign dileptons, L1",
				3, 0,  3);
  SignOfDil[2] = fs->make<TH1F>("SignOfDil2", "Number of opp-sign, like-sign dileptons, L2",
				3, 0,  3);
  SignOfDil[3] = fs->make<TH1F>("SignOfDil3", "Number of opp-sign, like-sign dileptons, L3",
				3, 0,  3);

  // Main kinematic variables for all muons
  for (int i = 0; i <= MAX_LEVELS; i++) {
    MuonEta[i] = fs->make<TH1F>(nameHist("MuonEta", i).c_str(), "Eta", 100, -5.,  5. );
    MuonRap[i] = fs->make<TH1F>(nameHist("MuonRap", i).c_str(), "Y",   100, -5.,  5. );
    MuonPhi[i] = fs->make<TH1F>(nameHist("MuonPhi", i).c_str(), "Phi", 100,  -pi, pi);
    if (i==1) {
      MuonPt[i] = fs->make<TH1F>(nameHist("MuonPt", i).c_str(), "Pt", 100,  0.,  150.);
      MuonPz[i] = fs->make<TH1F>(nameHist("MuonPz", i).c_str(), "Pz", 100,  0.,  800.);
      MuonP[i]  = fs->make<TH1F>(nameHist("MuonP", i).c_str(), "P",  100, -5.,  800.);
      MuonPVsEta[i]  = fs->make<TH2F>(nameHist("MuonPVsEta", i).c_str(), "P vs Eta",  100, -6., 6., 100, 0., 1000);
      MuonPtVsEta[i] = fs->make<TH2F>(nameHist("MuonPtVsEta", i).c_str(), "Pt vs Eta", 100, -6., 6., 100, 0., 150);
    }
    else {
      MuonPt[i] = fs->make<TH1F>(nameHist("MuonPt", i).c_str(), "Pt", 100, 0.,    peakMass);
      MuonPz[i] = fs->make<TH1F>(nameHist("MuonPz", i).c_str(), "Pz", 100, 0., 2.*peakMass);
      MuonP[i]  = fs->make<TH1F>(nameHist("MuonP", i).c_str(), "P",  100, 0., 2.*peakMass);
      MuonPVsEta[i]  = fs->make<TH2F>(nameHist("MuonPVsEta", i).c_str(), "P vs Eta",  100, -6., 6., 100, 0., 
				      2.*peakMass);
      MuonPtVsEta[i] = fs->make<TH2F>(nameHist("MuonPtVsEta", i).c_str(), "Pt vs Eta", 100, -6., 6., 100, 0., 
				      peakMass);
    }

    string titl = recLevelHelper.levelName(i) + ", all muons";
    SumPtR03[i][0] = fs->make<TH1F>(nameHist("SumPtR03", i, 0).c_str(),
				    recLevelHelper.levelName(i).c_str(), 30, 0, 30);
    titl = recLevelHelper.levelName(i) + ", w/ opp-sign dimuon";
    SumPtR03[i][1] = fs->make<TH1F>(nameHist("SumPtR03", i, 1).c_str(),
				    titl.c_str(), 30, 0, 30);

    // Main kinematic variables for opp-sign dileptons and
    // muons associated with opp-sign dileptons
    if (i == 1) {
      // AllDilMass[i] = fs->make<TH1F>(nameHist("AllDilMass", i).c_str(), "Opp-sign Dilepton Mass before Q cuts",
      //			       100, 0, 1200);
      DilMass[i]    = fs->make<TH1F>(nameHist("DilMass", i).c_str(), "Opp-sign Dilepton Mass", 100, 0, 1200);
      DilMassVsEta[i] = fs->make<TH2F>(nameHist("DilMassVsEta", i).c_str(), "Dilepton Mass vs Eta", 50, -10., 10.,
				       100, 0, 1200);
      DilMassVsY[i] = fs->make<TH2F>(nameHist("DilMassVsY", i).c_str(), "Dilepton Mass vs Y", 50, -5., 5., 100, 0,
				     1200);
      MuMVsMuP[i][2] = fs->make<TH2F>(nameHist("MuMVsMuP", i, 2).c_str(),"Pt  of mu- vs mu+", 100,0, 150,100, 0,150);
    }
    else {
      // AllDilMass[i] = fs->make<TH1F>(nameHist("AllDilMass", i).c_str(), "Opp-sign Dilepton Mass before Q cuts",
      //		       binSize, lowerMassWin,
      //		       upperMassWin);
      DilMass[i]    = fs->make<TH1F>(nameHist("DilMass", i).c_str(), "Opp-sign Dilepton Mass",
				     binSize, lowerMassWin,
				     upperMassWin);
      DilMassVsEta[i] = fs->make<TH2F>(nameHist("DilMassVsEta", i).c_str(), "Dilepton Mass vs Eta", 50, -10., 10., 
				       50, 0, upperMassWin);
      DilMassVsY[i] = fs->make<TH2F>(nameHist("DilMassVsY", i).c_str(), "Dilepton Mass vs Y", 50, -5., 5.,   50, 0, 
				     upperMassWin);
      MuMVsMuP[i][2] = fs->make<TH2F>(nameHist("MuMVsMuP", i, 2).c_str(),"Pt  of mu- vs mu+",
				      100, 0., peakMass,
				      100, 0., peakMass);
    }
    string tit = "Opp-sign Dilepton Mass, " + recLevelHelper.levelName(i);
    AllDilMass[i] =
      fs->make<TH1F>(nameHist("AllDilMass", i).c_str(), tit.c_str(), 50, 0., upperMassWin);

    MuMVsMuP[i][0] = fs->make<TH2F>(nameHist("MuMVsMuP", i, 0).c_str(),"Eta of mu- vs mu+", 100, -3, 3, 100, -3, 3);
    MuMVsMuP[i][1] = fs->make<TH2F>(nameHist("MuMVsMuP", i, 1).c_str(),"Phi of mu- vs mu+", 100, -pi, pi, 100, -pi, pi);
    for (int j = 0; j < 3; j++) {
      Phi[i][j]      = fs->make<TH1F>(nameHist("Phi", i, j).c_str(), "Phi",       100, -pi, pi);
      PVsEta[i][j]   = fs->make<TProfile>(nameHist("PVsEta", i, j).c_str(), "P vs Eta",  50, -6., 6.);
      PtVsEta[i][j]  = fs->make<TProfile>(nameHist("PtVsEta", i, j).c_str(), "Pt vs Eta", 50, -6., 6.);
      Eta[i][j]      = fs->make<TH1F>(nameHist("Eta", i, j).c_str(), "Eta",          100, -5., 5.);
      Rapidity[i][j] = fs->make<TH1F>(nameHist("Rapidity", i, j).c_str(), "y",            100, -5., 5.);
      if (j != 2) {
	if (i == 1) { // mu+ and mu- at L1
	  Pt[i][j] = fs->make<TH1F>(nameHist("Pt", i, j).c_str(), "Pt", 100, 0., 200.);
	  P[i][j]  = fs->make<TH1F>(nameHist("P", i, j).c_str(), "P",  100, 0., 800.);
	  Pz[i][j] = fs->make<TH1F>(nameHist("Pz", i, j).c_str(), "Pz", 100, 0., 800.);
	}
	else {
	  Pt[i][j] = fs->make<TH1F>(nameHist("Pt", i, j).c_str(), "Pt", 100, 0.,    peakMass);
	  P[i][j]  = fs->make<TH1F>(nameHist("P", i, j).c_str(), "P",  100, 0., 2.*peakMass);
	  Pz[i][j] = fs->make<TH1F>(nameHist("Pz", i, j).c_str(), "Pz", 100, 0., 2.*peakMass);
	}
      }
      else { // Special range for dilepton values.
	if (i==1) {
	  Pt[i][j] = fs->make<TH1F>(nameHist("Pt", i, j).c_str(), "Pt", 100, 0., 200.);
	  P[i][j]  = fs->make<TH1F>(nameHist("P", i, j).c_str(), "P",  100, 0., 800.);
	  Pz[i][j] = fs->make<TH1F>(nameHist("Pz", i, j).c_str(), "Pz", 100, 0., 800.);	
	}
	else {
	  Pt[i][j] = fs->make<TH1F>(nameHist("Pt", i, j).c_str(), "Pt", 100, 0., .5*peakMass);
	  P[i][j]  = fs->make<TH1F>(nameHist("P", i, j).c_str(), "P",
				    100, 0, 1000.+upperMassWin);
	  Pz[i][j] = fs->make<TH1F>(nameHist("Pz", i, j).c_str(), "Pz",
				    100, 0, 1000.+upperMassWin);
	}
      }
    }
  }
  for (int i = 0; i < NUM_REC_LEVELS; i++) {
    ZonDilMass[i] = fs->make<TH1F>(nameHist("ZonDilMass", i).c_str(), "Highest Z mass", 100, 0., 120.);
    ZofDilMass[i] = fs->make<TH1F>(nameHist("ZofDilMass", i).c_str(), "Second  Z mass", 100, 0., 120.);
  }
  Eta4muons = fs->make<TH1F>("Eta4muons", "Eta", 100, -5.,  5.);
  Pt4muons  = fs->make<TH1F>("Pt4muons", "pT",  100,  0.,100.);

  // Differences between different levels of reconstruction
  for (int i = l1; i <= MAX_LEVELS; i++) {
    string tit_eta = recLevelHelper.levelName(i) + " eta - Gen eta";
    string tit_phi = recLevelHelper.levelName(i) + " phi - Gen phi";
    string tit_pt  = recLevelHelper.levelName(i) + " pT - Gen pT";
    string tit_p   = "(" + recLevelHelper.levelName(i) + " P - Gen P)/(Gen P)";
    string tit_ppr = "(" + recLevelHelper.levelName(i) + " P - Gen P)/(Gen P) vs Gen P";
    if (i == l1) {
      EtaRes[i] = fs->make<TH1F>(nameHist("EtaRes", i).c_str(), tit_eta.c_str(), 100, -0.1,  0.1);
      PhiRes[i] = fs->make<TH1F>(nameHist("PhiRes", i).c_str(), tit_phi.c_str(), 100, -0.1,  0.1);
      PtDiff[i] = fs->make<TH1F>(nameHist("PtDiff", i).c_str(), tit_pt.c_str(),  100, 
				 -700.-.2*peakMass,
				 700.+.2*peakMass);
      PRes[i]   = fs->make<TH1F>(nameHist("PRes", i).c_str(), tit_p.c_str(),   100, -1., 1.);
      PResVsP[i]= fs->make<TProfile>(nameHist("PResVsP", i).c_str(), tit_ppr.c_str(), 50,
				     0., upperMassWin, -1., 1.);
    }
    else if (i == l2) {
      EtaRes[i] = fs->make<TH1F>(nameHist("EtaRes", i).c_str(), tit_eta.c_str(), 100, -0.01, 0.01);
      PhiRes[i] = fs->make<TH1F>(nameHist("PhiRes", i).c_str(), tit_phi.c_str(), 100, -0.01, 0.01);
      PtDiff[i] = fs->make<TH1F>(nameHist("PtDiff", i).c_str(), tit_pt.c_str(),  100,
				 -0.4*peakMass, 0.4*peakMass);
      PRes[i]   = fs->make<TH1F>(nameHist("PRes", i).c_str(), tit_p.c_str(),   100, -1., 1.);
      PResVsP[i]= fs->make<TProfile>(nameHist("PResVsP", i).c_str(), tit_ppr.c_str(), 50,
				     0., upperMassWin, -1., 1.);
    }
    else {
      EtaRes[i] = fs->make<TH1F>(nameHist("EtaRes", i).c_str(), tit_eta.c_str(), 100, -0.001,  0.001);
      PhiRes[i] = fs->make<TH1F>(nameHist("PhiRes", i).c_str(), tit_phi.c_str(), 100, -0.0005, 0.0005);
      PtDiff[i] = fs->make<TH1F>(nameHist("PtDiff", i).c_str(), tit_pt.c_str(),  100, 
				 -0.1*peakMass, 0.1*peakMass);
      PRes[i]   = fs->make<TH1F>(nameHist("PRes", i).c_str(), tit_p.c_str(),   100, -0.3, 0.3);
      PResVsP[i]= fs->make<TProfile>(nameHist("PResVsP", i).c_str(), tit_ppr.c_str(), 50,
				     0., upperMassWin, -0.3, 0.3);
    }
  }

  GenPhiResVsPhi[0] = fs->make<TProfile>("GenPhiResVsPhi0", "|L1 phi - Gen phi| vs Gen phi",
					 25, -pi, pi, -0.1,  0.1);
  GenPhiResVsPhi[1] = fs->make<TProfile>("GenPhiResVsPhi1", "|L2 phi - Gen phi| vs Gen phi",
					 25, -pi, pi, -0.1,  0.1);
  GenPhiResVsPhi[2] = fs->make<TProfile>("GenPhiResVsPhi2", "|L3 phi - Gen phi| vs Gen phi",
					 25, -pi, pi, -0.01, 0.01);

  GenInvPtRes[0] = fs->make<TH1F>("GenInvPtRes0", "(L1 1/Pt - Gen 1/Pt)/(Gen 1/Pt)", 100,
				  -.6*peakMass/140.,
				  .6*peakMass/140.);
  GenInvPtRes[1] = fs->make<TH1F>("GenInvPtRes1", "(L2 1/Pt - Gen 1/Pt)/(Gen 1/Pt)", 100, -2, 2);
  GenInvPtRes[2] = fs->make<TH1F>("GenInvPtRes2", "(L3 1/Pt - Gen 1/Pt)/(Gen 1/Pt)",
				  100, -0.3, 0.3);
  GenInvPtResVsPt[0] = fs->make<TProfile>("GenInvPtResVsPt0",
					  "(L1 1/pT - Gen 1/pT)/(Gen 1/pT) vs Gen pT",
					  50, 0., peakMass,
					  -.6*peakMass/140.,
					  .6*peakMass/140.);
  GenInvPtResVsPt[1] = fs->make<TProfile>("GenInvPtResVsPt1",
					  "(L2 1/pT - Gen 1/pT)/(Gen 1/pT) vs Gen pT",
					  50, 0., peakMass, -2., 2.);
  GenInvPtResVsPt[2] = fs->make<TProfile>("GenInvPtResVsPt2",
					  "(L3 1/pT - Gen 1/pT)/(Gen 1/pT) vs Gen pT",
					  50, 0., peakMass, -0.3, 0.3);

  GenInvPRes[0] = fs->make<TH1F>("GenInvPRes0", "(L1 1/P - Gen 1/P)/(Gen 1/P)", 100,
				 -upperMassWin/200.,
				 upperMassWin/200.);
  GenInvPRes[1] = fs->make<TH1F>("GenInvPRes1", "(L2 1/P - Gen 1/P)/(Gen 1/P)", 100, -2.,  2.);
  GenInvPRes[2] = fs->make<TH1F>("GenInvPRes2", "(L3 1/P - Gen 1/P)/(Gen 1/P)", 100, -0.3, 0.3);
  GenPResVsPt[0] = fs->make<TProfile>("GenPResVsPt0", "(L1 P - Gen P)/(Gen P) vs Gen pT",
				      50, 0., .6*peakMass, -1.,  1.);
  GenPResVsPt[1] = fs->make<TProfile>("GenPResVsPt1", "(L2 P - Gen P)/(Gen P) vs Gen pT",
				      50, 0., .6*peakMass, -1.,  1.);
  GenPResVsPt[2] = fs->make<TProfile>("GenPResVsPt2", "(L3 P - Gen P)/(Gen P) vs Gen pT",
				      50, 0., .7*peakMass, -0.3, 0.3);
  GenInvPResVsPt[0] = fs->make<TProfile>("GenInvPResVsPt0","(L1 1/P - Gen 1/P)/(Gen 1/P) vs Gen pT",
					 50, 0., peakMass,
					 -upperMassWin/200.,
					 upperMassWin/200.);
  GenInvPResVsPt[1] = fs->make<TProfile>("GenInvPResVsPt1","(L2 1/P - Gen 1/P)/(Gen 1/P) vs Gen pT",
					 50, 0., peakMass, -2.,  2.);
  GenInvPResVsPt[2] = fs->make<TProfile>("GenInvPResVsPt2","(L3 1/P - Gen 1/P)/(Gen 1/P) vs Gen pT",
					 50, 0., peakMass, -0.3, 0.3);

  GenEtaResScat[0] = fs->make<TH2F>("GenEtaResScat0", "L1 Eta vs Gen Eta", 100, -3, 3, 100, -3, 3);
  GenEtaResScat[1] = fs->make<TH2F>("GenEtaResScat1", "L2 Eta vs Gen Eta", 100, -3, 3, 100, -3, 3);
  GenEtaResScat[2] = fs->make<TH2F>("GenEtaResScat2", "L3 Eta vs Gen Eta", 100, -3, 3, 100, -3, 3);
  GenPhiResScat[0]=fs->make<TH2F>("GenPhiResScat0", "L1 Phi vs Gen Phi", 100, -pi, pi, 100, -pi, pi);
  GenPhiResScat[1]=fs->make<TH2F>("GenPhiResScat1", "L2 Phi vs Gen Phi", 100, -pi, pi, 100, -pi, pi);
  GenPhiResScat[2]=fs->make<TH2F>("GenPhiResScat2", "L3 Phi vs Gen Phi", 100, -pi, pi, 100, -pi, pi);
  GenPtResScat[0]  = fs->make<TH2F>("GenPtResScat0", "L1 Pt vs Gen Pt",
				    100, 0., peakMass, 100, 0., 200.);
  GenPtResScat[1]  = fs->make<TH2F>("GenPtResScat1", "L2 Pt vs Gen Pt",
				    100, 0., peakMass,
				    100, 0., peakMass);
  GenPtResScat[2]  = fs->make<TH2F>("GenPtResScat2", "L3 Pt vs Gen Pt",
				    100, 0., peakMass,
				    100, 0., peakMass);
  AllDilMassRes = fs->make<TH1F>("AllDilMassRes", "L3 Mass - Gen Mass, before Q cuts", 100, 
				 -.2*peakMass,
				 .2*peakMass);
  GenDilMassRes[0] = fs->make<TH1F>("GenDilMassRes0", "L1 Mass - Gen Mass", 100,
				    -peakMass, peakMass);
  GenDilMassRes[1] = fs->make<TH1F>("GenDilMassRes1", "L2 Mass - Gen Mass", 100,
				    -peakMass, peakMass);
  GenDilMassRes[2] = fs->make<TH1F>("GenDilMassRes2", "L3 Mass - Gen Mass", 100,
				    -.2*peakMass,
				    .2*peakMass);
  GenDilMassFrRes[0] = fs->make<TH1F>("GenDilMassFrRes0","(L1 Mass-Gen Mass)/(Gen Mass)",100,-1.,1.);
  GenDilMassFrRes[1] = fs->make<TH1F>("GenDilMassFrRes1","(L2 Mass-Gen Mass)/(Gen Mass)",100,-1.,1.);
  GenDilMassFrRes[2] = fs->make<TH1F>("GenDilMassFrRes2","(L3 Mass-Gen Mass)/(Gen Mass)",
				      100, -0.2, 0.2);
  GenDilMassResScat[0] = fs->make<TProfile>("GenDilMassResScat0", "L1-Gen Mass vs Gen Mass", 25, 0.,
					    upperMassWin, 0., 
					    peakMass*peakMass, " ");
  GenDilMassResScat[1] = fs->make<TProfile>("GenDilMassResScat1", "L2-Gen Mass vs Gen Mass", 25, 0.,
					    upperMassWin, 0., 
					    peakMass*peakMass, " ");
  GenDilMassResScat[2] = fs->make<TProfile>("GenDilMassResScat2", "L3-Gen Mass vs Gen Mass", 25, 0.,
					    upperMassWin, 0., 
					    0.04*peakMass*peakMass, " ");
  GenDilMassFrResScat[0] = fs->make<TProfile>("GenDilMassFrResScat0",
					      "(L1-Gen Mass)/(Gen Mass) vs Gen Mass",
					      25, 0., upperMassWin, 
					      -1., 1., " ");
  GenDilMassFrResScat[1] = fs->make<TProfile>("GenDilMassFrResScat1",
					      "(L2-Gen Mass)/(Gen Mass) vs Gen Mass",
					      25, 0., upperMassWin, 
					      -1., 1., " ");
  GenDilMassFrResScat[2] = fs->make<TProfile>("GenDilMassFrResScat2",
					      "(L3-Gen Mass)/(Gen Mass) vs Gen Mass",
					      25, 0., upperMassWin, 
					      -0.2, 0.2, " ");

  L1EtaRes[0] = fs->make<TH1F>("L1EtaRes0", "L2 Eta - L1 Eta", 100, -0.1, 0.1);
  L1EtaRes[1] = fs->make<TH1F>("L1EtaRes1", "L3 Eta - L1 Eta", 100, -0.1, 0.1);
  L1PhiRes[0] = fs->make<TH1F>("L1PhiRes0", "L2 Phi - L1 Phi", 100, -0.1, 0.1);
  L1PhiRes[1] = fs->make<TH1F>("L1PhiRes1", "L3 Phi - L1 Phi", 100, -0.1, 0.1);
  L1PtDiff[0] = fs->make<TH1F>("L1PtDiff0", "L2 Pt - L1 Pt",   100, -200.,
			       -200.+peakMass);
  L1PtDiff[1] = fs->make<TH1F>("L1PtDiff1", "L3 Pt - L1 Pt",   100, -150., 
			       -150.+peakMass);

  L1EtaResScat[0] = fs->make<TH2F>("L1EtaResScat0", "L2 Eta vs L1 Eta", 100, -3, 3, 100, -3, 3);
  L1EtaResScat[1] = fs->make<TH2F>("L1EtaResScat1", "L3 Eta vs L1 Eta", 100, -3, 3, 100, -3, 3);
  L1PhiResScat[0] = fs->make<TH2F>("L1PhiResScat0", "L2 Phi vs L1 Phi", 100, -pi, pi, 100, -pi, pi);
  L1PhiResScat[1] = fs->make<TH2F>("L1PhiResScat1", "L3 Phi vs L1 Phi", 100, -pi, pi, 100, -pi, pi);
  L1PtResScat[0]  = fs->make<TH2F>("L1PtResScat0", "L2 Pt vs L1 Pt", 100, 0., 200.,
				   100, 0., peakMass);
  L1PtResScat[1]  = fs->make<TH2F>("L1PtResScat1", "L3 Pt vs L1 Pt", 100, 0., 200.,
				   100, 0., peakMass);

  L2EtaRes = fs->make<TH1F>("L2EtaRes", "L3 Eta - L2 Eta", 100, -0.025, 0.025);
  L2PhiRes = fs->make<TH1F>("L2PhiRes", "L3 Phi - L2 Phi", 100, -0.05,  0.05);
  L2PtDiff = fs->make<TH1F>("L2PtDiff", "L3 Pt - L2 Pt",   100,
			    -.6*peakMass, .6*peakMass);
  L2EtaResScat = fs->make<TH2F>("L2EtaResScat", "L3 Eta vs L2 Eta", 100, -3,  3,  100, -3, 3);
  L2PhiResScat = fs->make<TH2F>("L2PhiResScat", "L3 Phi vs L2 Phi", 100,  -pi, pi,  100, -pi, pi);
  L2PtResScat  = fs->make<TH2F>("L2PtResScat", "L3 Pt vs L2 Pt", 100, 0., peakMass,
				100, 0., peakMass);
}

void Zprime2muResolution::BookEffHistos() {
  const double pi = TMath::Pi();
  string tit;
  for (int i = lgen; i <= MAX_LEVELS; i++) {
    tit = "Gen eta, " + recLevelHelper.levelName(i) + " muons";
    EffVsEta[i] = fs->make<TH1F>(nameHist("EffVsEta", i).c_str(), tit.c_str(), 50, -2.5, 2.5);
    tit = "Gen phi, " + recLevelHelper.levelName(i) + " muons";
    EffVsPhi[i] = fs->make<TH1F>(nameHist("EffVsPhi", i).c_str(), tit.c_str(), 63,  -pi, pi);
    tit = "Gen pT, "  + recLevelHelper.levelName(i) + " muons";
    EffVsPt[i]  = fs->make<TH1F>(nameHist("EffVsPt", i).c_str(), tit.c_str(), 35,  0., 3500.);
  }

  RecMass[0][0] = fs->make<TH1F>("RecMass00", "Gen mass, gen, all",  27, 0., maxTrigMass);
  RecMass[0][1] = fs->make<TH1F>("RecMass01", "Gen mass, gen, in acceptance", 27, 0., maxTrigMass);
  RecMass[0][2] = fs->make<TH1F>("RecMass02", "placeholder",  27, 0., maxTrigMass);
  RecMass[1][0] = fs->make<TH1F>("RecMass10", "Gen mass, L3, 2mu",   27, 0., maxTrigMass);
  RecMass[1][1] = fs->make<TH1F>("RecMass11", "Gen mass, L3, dimu",  27, 0., maxTrigMass);
  RecMass[1][2] = fs->make<TH1F>("RecMass12", "Gen mass, L3, 2mu, w/ cuts",  27, 0., maxTrigMass);
  RecMass[2][0] = fs->make<TH1F>("RecMass20", "Gen mass, GMR, 2mu",  27, 0., maxTrigMass);
  RecMass[2][1] = fs->make<TH1F>("RecMass21", "Gen mass, GMR, dimu", 27, 0., maxTrigMass);
  RecMass[2][2] = fs->make<TH1F>("RecMass22", "Gen mass, GMR, 2mu, w/ cuts",  27, 0., maxTrigMass);
  RecMass[3][0] = fs->make<TH1F>("RecMass30", "Gen mass, OPT, 2mu",  27, 0., maxTrigMass);
  RecMass[3][1] = fs->make<TH1F>("RecMass31", "Gen mass, OPT, dimu", 27, 0., maxTrigMass);
  RecMass[3][2] = fs->make<TH1F>("RecMass32", "Gen mass, OPT, 2mu, w/ cuts",  27, 0., maxTrigMass);

  const string locsubs[4] = {"Barrel","Overlap","Endcap","Outside"};
  vector<string> loctitles;
  for (int i = 0; i < 4; i++)
    for (int j = i; j < 4; j++) {
      loctitles.push_back(locsubs[i] + "-" + locsubs[j]);
    }
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 10; j++) {
      string title;
      if (i == 0) title = loctitles[j] + ", 2mu";
      else title = loctitles[j] + ", dimu";
      RecMassByLoc[i][j] =
	fs->make<TH1F>(nameHist("RecMassByLoc", i, j).c_str(), title.c_str(),
		       27, 0., maxTrigMass);
    }
  }

  // Define errors
  for (int i = lgen; i <= MAX_LEVELS; i++) {
    EffVsEta[i]->Sumw2();
    EffVsPhi[i]->Sumw2();
    EffVsPt[i]->Sumw2();
  }
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      RecMass[i][j]->Sumw2();
    }
  }

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 10; j++)
      RecMassByLoc[i][j]->Sumw2();
}

void Zprime2muResolution::BookPtResHistos() {
  ResidualPt[0] = fs->make<TProfile>("ResidualPt0", "(delta 1/pT)/(1/pT) vs Eta, pT < 10",
				     12, 0., 2.4, -10., 10.); // very long tails
  ResidualPt[1] = fs->make<TProfile>("ResidualPt1", "(delta 1/pT)/(1/pT) vs Eta, pT < 25",
				     12, 0., 2.4, -10., 10.);
  ResidualPt[2] = fs->make<TProfile>("ResidualPt2", "(delta 1/pT)/(1/pT) vs Eta, pT < 50",
				     12, 0., 2.4, -10., 10.);
  ResidualPt[3] = fs->make<TProfile>("ResidualPt3", "(delta 1/pT)/(1/pT) vs Eta, pT < 75",
				     12, 0., 2.4, -10., 10.);
  ResidualPt[4] = fs->make<TProfile>("ResidualPt4", "(delta 1/pT)/(1/pT) vs Eta, pT < 100",
				     12, 0., 2.4, -10., 10.);
  ResidualPt[5] = fs->make<TProfile>("ResidualPt5", "(delta 1/pT)/(1/pT) vs Eta, pT < 1000",
				     12, 0., 2.4, -10., 10.);

  // 1/p resolution for GMR and OPT, separately for barrel, overlap and endcap.
  MuInvPRes[0][0] = fs->make<TH1F>("MuInvPRes00", "1/P res, GMR, barrel",  100, -0.3, 0.3);
  MuInvPRes[0][1] = fs->make<TH1F>("MuInvPRes01", "1/P res, GMR, overlap", 100, -0.3, 0.3);
  MuInvPRes[0][2] = fs->make<TH1F>("MuInvPRes02", "1/P res, GMR, endcap",  100, -0.3, 0.3);

  MuInvPRes[1][0] = fs->make<TH1F>("MuInvPRes10", "1/P res, OPT, barrel",  100, -0.3, 0.3);
  MuInvPRes[1][1] = fs->make<TH1F>("MuInvPRes11", "1/P res, OPT, overlap", 100, -0.3, 0.3);
  MuInvPRes[1][2] = fs->make<TH1F>("MuInvPRes12", "1/P res, OPT, endcap",  100, -0.3, 0.3);

  // 1/pT resolution for GMR and OPT, separately for barrel, overlap and endcap
  MuInvPtRes[0][0] = fs->make<TH1F>("MuInvPtRes00", "1/pT res, GMR, barrel",  100, -0.3, 0.3);
  MuInvPtRes[0][1] = fs->make<TH1F>("MuInvPtRes01", "1/pT res, GMR, overlap", 100, -0.3, 0.3);
  MuInvPtRes[0][2] = fs->make<TH1F>("MuInvPtRes02", "1/pT res, GMR, endcap",  100, -0.3, 0.3);

  MuInvPtRes[1][0] = fs->make<TH1F>("MuInvPtRes10", "1/pT res, OPT, barrel",  100, -0.3, 0.3);
  MuInvPtRes[1][1] = fs->make<TH1F>("MuInvPtRes11", "1/pT res, OPT, overlap", 100, -0.3, 0.3);
  MuInvPtRes[1][2] = fs->make<TH1F>("MuInvPtRes12", "1/pT res, OPT, endcap",  100, -0.3, 0.3);

  // pT res for all fits
  MuPtDiff[0] = fs->make<TH1F>("MuPtDiff0",  "L3 pT - Gen pT", 100, -1500., 1500.);
  MuPtDiff[1] = fs->make<TH1F>("MuPtDiff1", "OPT pT - Gen pT", 100, -1500., 1500.);
  TotInvPtRes[0] = fs->make<TH1F>("TotInvPtRes0", 
				  "(L3 1/pT - Gen 1/pT)/(Gen 1/pT)", 100, -.5, .5);
  TotInvPtRes[1] = fs->make<TH1F>("TotInvPtRes1", "(Tracker Only 1/pT - Gen 1/pT)/(Gen 1/pT)", 
				  100, -.5, .5);
  TotInvPtRes[2] = 
    fs->make<TH1F>("TotInvPtRes2", "(Tracker + 1 mu-station 1/pT - Gen 1/pT)/(Gen 1/pT)", 
		   100, -.5, .5);
  //TotInvPtRes[3] = fs->make<TH1F>("TotInvPtRes3", "(ABCM 1/pT - Gen 1/pT)/(Gen 1/pT)", 
  //		    100, -.5, .5);
  TotInvPtRes[3] = fs->make<TH1F>("TotInvPtRes3", "(Outer 1/pT - Gen 1/pT)/(Gen 1/pT)", 
				  100, -.5, .5);
  // 1/pT pull for all fits
  TotInvPtPull[0] = fs->make<TH1F>("TotInvPtPull0", "L3 1/pT Pull", 100, -10., 10.);
  TotInvPtPull[1] = fs->make<TH1F>("TotInvPtPull1", "Tracker Only 1/pT Pull", 100, -10., 10.);
  TotInvPtPull[2] = fs->make<TH1F>("TotInvPtPull2", "Tracker + 1 mu-station 1/pT Pull", 
				   100, -10., 10.);
  //TotInvPtPull[3] = fs->make<TH1F>("TotInvPtPull3", "ABCM 1/pT Pull", 100, -10., 10.);
  TotInvPtPull[3] = fs->make<TH1F>("TotInvPtPull3", "Outer 1/pT Pull", 100, -10., 10.);

  // pT res for all fits (split up by barrel and endcap)
  InvPtRes[0][0] = fs->make<TH1F>("InvPtRes00", "eta < 1.04, L3", 100, -.3, .3); 
  InvPtRes[0][1] = fs->make<TH1F>("InvPtRes01", "eta > 1.04, L3", 100, -.3, .3); 
  InvPtRes[1][0] = fs->make<TH1F>("InvPtRes10", "eta < 1.04, Tracker Only", 100, -.3, .3); 
  InvPtRes[1][1] = fs->make<TH1F>("InvPtRes11", "eta > 1.04, Tracker Only", 100, -.3, .3); 
  InvPtRes[2][0] = fs->make<TH1F>("InvPtRes20", "eta < 1.04, Tracker + 1 mu-station",
				  100, -.3, .3); 
  InvPtRes[2][1] = fs->make<TH1F>("InvPtRes21", "eta > 1.04, Tracker + 1 mu-station", 
				  100, -.3, .3); 
  InvPtRes[3][0] = fs->make<TH1F>("InvPtRes30", "eta < 1.04, GMR", 100, -.3, .3); 
  InvPtRes[3][1] = fs->make<TH1F>("InvPtRes31", "eta > 1.04, GMR", 100, -.3, .3); 

  // 1/pT pull for all fits (split up by barrel and endcap)
  InvPtPull[0][0] = fs->make<TH1F>("InvPtPull00", "eta < 1.04, L3", 100, -10., 10.); 
  InvPtPull[0][1] = fs->make<TH1F>("InvPtPull01", "eta > 1.04, L3", 100, -10., 10.); 
  InvPtPull[1][0] = fs->make<TH1F>("InvPtPull10", "eta < 1.04, Tracker Only", 100, -10., 10.); 
  InvPtPull[1][1] = fs->make<TH1F>("InvPtPull11", "eta > 1.04, Tracker Only", 100, -10., 10.);
  InvPtPull[2][0] = fs->make<TH1F>("InvPtPull20", "eta < 1.04, Tracker + 1 mu-station", 
				   100, -10., 10.); 
  InvPtPull[2][1] = fs->make<TH1F>("InvPtPull21", "eta > 1.04, Tracker + 1 mu-station", 
				   100, -10., 10.); 
  InvPtPull[3][0] = fs->make<TH1F>("InvPtPull30", "eta < 1.04, GMR", 100, -10., 10.); 
  InvPtPull[3][1] = fs->make<TH1F>("InvPtPull31", "eta > 1.04, GMR", 100, -10., 10.); 
}

void Zprime2muResolution::BookDilResHistos(){
  // Dilepton mass spectra:
  //  - dilepton only (i=0);
  //  - dilepton plus nearby photon candidates (i=1).
  double mass_min = lowerMassWin;
  double mass_max = upperMassWin;
  for (int i = lgen; i <= MAX_LEVELS; i++) {
    if (i == l1 || i == l2) mass_min = 0.;
    else                    mass_min = lowerMassWin;
    for (int j = 0; j < 2; j++) {
      DilMassComp[i][j] =
	fs->make<TH1F>(nameHist("DilMassComp", i, j).c_str(), recLevelHelper.levelName(i).c_str(),
		       binSize, mass_min, mass_max);
    }
  }

  // Invariant mass resolutions:
  //  - dilepton mass vs mass reconstructed from GEANT muons (0, 3);
  //  - dilepton mass w.r.t. the true (PYTHIA) resonance mass (1, 4);
  //  - dilepton+photon(s) mass w.r.t. the true (PYTHIA) resonance mass (2, 5);
  //     - all events (0-2);
  //     - events at the Z' mass peak (3-5).
  for (int j = 0; j < 6; j++)
    DilMassRes[0][j] = 0;
  for (int i = l1; i <= MAX_LEVELS; i++) {
    for (int j = 0; j < 6; j++) {
      if (i == l1 || i == l2) {
	DilMassRes[i][j] = fs->make<TH1F>(nameHist("DilMassRes", i, j).c_str(), recLevelHelper.levelName(i).c_str(), 100, -1.,  1.);
      }
      else {
	DilMassRes[i][j] = fs->make<TH1F>(nameHist("DilMassRes", i, j).c_str(), recLevelHelper.levelName(i).c_str(), 100, -0.3, 0.3);
      }
    }
  }
  for (int j = 0; j < 6; j++) {
    DilMassRes[MAX_LEVELS+1][j] = fs->make<TH1F>(nameHist("DilMassRes", MAX_LEVELS+1, j).c_str(), "OPT, D0 method",
						 100, -0.3, 0.3);
  }

  // Dilepton pT resolution
  DilPtRes[0] = 0;
  for (int i = l1; i <= MAX_LEVELS; i++) {
    DilPtRes[i] = fs->make<TH1F>(nameHist("DilPtRes", i).c_str(), recLevelHelper.levelName(i).c_str(), 100, -1., 4.);
  }

  MuPVsMuM[0] = fs->make<TH2F>("MuPVsMuM0", "pT Rel Error for dM > 0.1", 100, 0., 1.,
			       100, 0., 1.);
  MuPVsMuM[1] = fs->make<TH2F>("MuPVsMuM1", "pT Rel Error for dM < 0.1", 100, 0., 1., 
			       100, 0., 1.);
}

void Zprime2muResolution::BookChargeResHistos() {
  for (int i = l1; i <= MAX_LEVELS; i++) {
    string tit = recLevelHelper.levelName(i) + " charge - Gen charge";
    QRes[i] = fs->make<TH1F>(nameHist("QRes", i).c_str(), tit.c_str(), 7, -3.5, 3.5);
  }

  QResVsPt[0][0] = fs->make<TH1F>("QResVsPt00", "L1 pT, right charge", 50, 0., 150.);
  QResVsPt[1][0] = fs->make<TH1F>("QResVsPt10", "L2 pT, right charge", 50, 0., 500.);
  QResVsPt[2][0] = fs->make<TH1F>("QResVsPt20", "L3 pT, right charge", 50, 0., 500.);
  QResVsPt[0][1] = fs->make<TH1F>("QResVsPt01", "L1 pT, wrong charge", 50, 0., 150.);
  QResVsPt[1][1] = fs->make<TH1F>("QResVsPt11", "L2 pT, wrong charge", 50, 0., 500.);
  QResVsPt[2][1] = fs->make<TH1F>("QResVsPt21", "L3 pT, wrong charge", 50, 0., 500.);
  QResVsInvPt[0][0] = fs->make<TH1F>("QResVsInvPt00", "1/L1pT, right charge", 50, 0., 0.04);
  QResVsInvPt[1][0] = fs->make<TH1F>("QResVsInvPt10", "1/L2pT, right charge", 50, 0., 0.04);
  QResVsInvPt[2][0] = fs->make<TH1F>("QResVsInvPt20", "1/L3pT, right charge", 50, 0., 0.04);
  QResVsInvPt[0][1] = fs->make<TH1F>("QResVsInvPt01", "1/L1pT, wrong charge", 50, 0., 0.04);
  QResVsInvPt[1][1] = fs->make<TH1F>("QResVsInvPt11", "1/L2pT, wrong charge", 50, 0., 0.04);
  QResVsInvPt[2][1] = fs->make<TH1F>("QResVsInvPt21", "1/L3pT, wrong charge", 50, 0., 0.04);
  QResVsP[0][0]  = fs->make<TH1F>("QResVsP00", "L1 P, right charge", 50, 0., 750.);
  QResVsP[1][0]  = fs->make<TH1F>("QResVsP10", "L2 P, right charge", 50, 0., 
				  upperMassWin);
  QResVsP[2][0]  = fs->make<TH1F>("QResVsP20", "L3 P, right charge", 50, 0., 
				  upperMassWin);
  QResVsP[0][1]  = fs->make<TH1F>("QResVsP01", "L1 P, wrong charge", 50, 0., 750.);
  QResVsP[1][1]  = fs->make<TH1F>("QResVsP11", "L2 P, wrong charge", 50, 0., 
				  upperMassWin);
  QResVsP[2][1]  = fs->make<TH1F>("QResVsP21", "L3 P, wrong charge", 50, 0., 
				  upperMassWin);
  QResVsInvP[0][0] = fs->make<TH1F>("QResVsInvP00", "1/L1P, right charge", 50, 0., 0.02);
  QResVsInvP[1][0] = fs->make<TH1F>("QResVsInvP10", "1/L2P, right charge", 50, 0., 0.02);
  QResVsInvP[2][0] = fs->make<TH1F>("QResVsInvP20", "1/L3P, right charge", 50, 0., 0.02);
  QResVsInvP[0][1] = fs->make<TH1F>("QResVsInvP01", "1/L1P, wrong charge", 50, 0., 0.02);
  QResVsInvP[1][1] = fs->make<TH1F>("QResVsInvP11", "1/L2P, wrong charge", 50, 0., 0.02);
  QResVsInvP[2][1] = fs->make<TH1F>("QResVsInvP21", "1/L3P, wrong charge", 50, 0., 0.02);

  // Define errors
  for (int i_rec = l1; i_rec <= MAX_LEVELS; i_rec++) {
    QRes[i_rec]->Sumw2();
  }
  for (int i_rec = 0; i_rec < NUM_REC_LEVELS-1; i_rec++) {
    for (int j = 0; j < 2; j++) {
      QResVsPt[i_rec][j]->Sumw2();
      QResVsInvPt[i_rec][j]->Sumw2();
      QResVsP[i_rec][j]->Sumw2();
      QResVsInvP[i_rec][j]->Sumw2();
    }
  }
}

//-----------------------------------------------------------------------------
//         Fill resolution histograms for L1, L2, and L3 values
//-----------------------------------------------------------------------------
void Zprime2muResolution::calcResolution(const bool debug) {
  bool accept2_4 = false, accept2_1 = false;
  double gen_mass = -999.;
  unsigned idi;
  LeptonRefVector::const_iterator plep;
  reco::CompositeCandidateCollection::const_iterator pdi;

  if (allDileptons[lgen].size() > 0) {
    // JMTBAD fill the histograms for every dil in the event?
    // Fill histograms for acceptance studies.
    gen_mass = allDileptons[lgen][0].mass(); // (in GeV)
    GenMassAllEvents->Fill(gen_mass);
    // Check whether both muons are inside the full eta coverage.
    if (numDaughtersInAcc(allDileptons[lgen][0], ETA_CUT) >= 2)
      GenMassInAccept->Fill(gen_mass);
  }

  // Get origin of generated leptons.
  for (plep = allLeptons[lgen].begin(); plep != allLeptons[lgen].end(); plep++)
    Origin[0]->Fill(getOrigin(motherId(*plep)));

  // some additional plots for improved diagnostics
  fillEffHistos();
  fillPtResHistos(debug);
  fillChargeResHistos(debug);

  for (int i_rec = 0; i_rec < NUM_REC_LEVELS; i_rec++) {
    // Total number of muons found in each level of reconstruction.
    NMuons[i_rec][0]->Fill((double)allLeptons[i_rec].size());
  }

  for (int i_rec = 0; i_rec <= MAX_LEVELS; i_rec++) {
    // Fill y, eta, phi, p, pt, pz, pt vs eta, p vs eta plots for all muons.
    fillMuonHistos(i_rec, debug);
  }

  // Main loop over generated, L1, L2, L3 and off-line information.
  for (int i_rec = 0; i_rec <= MAX_LEVELS; i_rec++) {
    const reco::CompositeCandidateCollection& dileptons = getDileptons(i_rec);

    if (i_rec == 0) {
      if (dileptons.size() == 1) { // one di-muon at generation
	gen_mass = dileptons[0].mass()/1000.; // (in TeV)
	// Check whether both muons are inside the full eta coverage.
	if (numDaughtersInAcc(dileptons[0], ETA_CUT) >= 2) {
	  accept2_4 = true;
	  // At least one muon is in the limited muon trigger acceptance.
	  if (numDaughtersInAcc(dileptons[0], TRIGGER_ETA_CUT[l1]) >= 1)
	    accept2_1 = true;
	}
      }
      else {
	if (!doingHiggs && useGen)
	  edm::LogWarning("Zprime2muResolution")
	    << "+++ Warning in calcResolution: " << dileptons.size()
	    << " generated dimuons found! +++\n";
      }
    }

    // Trigger decision
    bool decision = false;
    if (i_rec <= l3) {
      decision = passTrigger(i_rec);
      if (debug)
	LogTrace("Zprime2muResolution")
	  << "irec " << i_rec << " passTrig " << decision
	  << " word " << trigWord[i_rec];
      //if (accept2_4 && gen_mass > 2.6) 
      TrigResult[i_rec][0]->Fill((double)trigWord[i_rec]);
      if (decision) TrigMass[i_rec][0]->Fill(gen_mass);
      if (accept2_4) {
	TrigResult[i_rec][1]->Fill((double)trigWord[i_rec]);
	//if (gen_mass > 2.6)
	if (decision) TrigMass[i_rec][1]->Fill(gen_mass);
      }
      if (accept2_1) {
	TrigResult[i_rec][2]->Fill((double)trigWord[i_rec]);
	if (decision) TrigMass[i_rec][2]->Fill(gen_mass);
      }

      // Check what happens to useful events (both mu's in acceptance)
      if (accept2_4 && !decision) {
	EventsInAccFailed->Fill((double)i_rec);
	NMuons[i_rec][3]->Fill((double)allLeptons[i_rec].size());
      }

      // More info for L1
      if (i_rec == l1) {
	if (decision) {
	  if (allLeptons[i_rec].size() == 1)
	    L1TrigPassSingleMu->Fill(allLeptons[i_rec][0]->pt());
	  else if (allLeptons[i_rec].size() == 2)
	    L1TrigPassMu2VsMu1->Fill(allLeptons[i_rec][0]->pt(),
				     allLeptons[i_rec][1]->pt());
	}
	else {
	  if (allLeptons[i_rec].size() == 1)
	    L1TrigFailSingleMu->Fill(allLeptons[i_rec][0]->pt());
	  else if (allLeptons[i_rec].size() == 2)
	    L1TrigFailMu2VsMu1->Fill(allLeptons[i_rec][0]->pt(),
				     allLeptons[i_rec][1]->pt());
	}
      }
    }

    if (i_rec == l2) {
      // Number of muon hits at Level-2
      for (plep = allLeptons[l2].begin(); plep != allLeptons[l2].end(); plep++)
	L2MuonHits->Fill(nHits(*plep, HITS_MU));
    }
    else if (i_rec == l3) {
      // Number of tracker (silicon + pixel) hits at Level-3
      for (plep = allLeptons[l3].begin(); plep != allLeptons[l3].end(); plep++)
	L3TrackerHits->Fill(nHits(*plep, HITS_TRK));
    }

    // Number of muon hits and chi2/d.o.f. for off-line (GMR) muons.
    // Unfortunately, the number of muon hits includes RPC hits, and there
    // seems to be no way to subtract them.
    if (i_rec == lgmr) {
      for (plep = allLeptons[lgmr].begin(); plep != allLeptons[lgmr].end(); plep++) {
	int muHits = nHits(*plep, HITS_MU);
	const reco::Track* tk = getMainTrack(*plep);
	if (tk == 0) continue;
	double chi2dof = tk->chi2()/tk->ndof();
	double eta = fabs((*plep)->eta());
	if (eta < 0.9) {
	  GMRMuonHits[0]->Fill(muHits); // barrel
	  GMRChi2dof[0]->Fill(chi2dof);
	}
	else if (eta < 1.2) {
	  GMRMuonHits[1]->Fill(muHits); // overlap
	  GMRChi2dof[1]->Fill(chi2dof);
	}
	else {
	  GMRMuonHits[2]->Fill(muHits); // endcap
	  GMRChi2dof[2]->Fill(chi2dof);
	}
      }
    }

    // Cut on the Trigger bits
    if (i_rec <= l3 && cutTrig[i_rec] && !decision) break;

    // Total number of muons for events passing the Trigger.
    if (i_rec <= l3) NMuons[i_rec][2]->Fill((double)allLeptons[i_rec].size());

    // Fill histos with opposite sign and same sign dileptons
    if (i_rec > lgen && i_rec <= l3) fillSignOfDilepton(i_rec, debug);

    // If an opposite-sign dilepton is found, fill histos.  A dilepton
    // should always be found at the generated level.
    for (idi = 0, pdi = dileptons.begin(); pdi != dileptons.end(); pdi++, idi++) {

      //if (i_rec == 3) dumpDilEvents(2, i_rec, *pdi);

      // Fill mass histograms for all dileptons found, before we apply
      // the "off-line" track quality cuts.
      if (i_rec <= l3) AllDilMass[i_rec]->Fill(pdi->mass());
      if (i_rec == l3) {
	if (dileptons.size() == allDileptons[lgen].size())
	  // JMTBAD rely on ordering
	  AllDilMassRes->
	    Fill(pdi->mass() - allDileptons[lgen][idi].mass());
	if (useGen && allDileptons[lgen].size() <= 0)
	  edm::LogWarning("Zprime2muResolution")
	    << "+++ Warning: no dilepton in the MC! +++\n";
      }

      // Check the "off-line" track quality and apply the cuts
      int qcut_mum = 0, qcut_mup = 0; // quality cut number
      if (DO_QCUTS ?
	  dilQCheck(*pdi, QSEL, qcut_mum, qcut_mup) : true) {
  
	// Keep track of the number of opposite-sign dileptons found
	// at each level of reconstruction.
	if (i_rec == 0) {
	  NumDilVsRec->Fill(0.);
	  SignOfDil[0]->Fill(0.);
	  // Num Gen dileptons with both muons inside the eta coverage.
	  if (accept2_4) {
	    NumDilVsRec->Fill(1.);
	    SignOfDil[0]->Fill(1.);
	  }
	  // Get origin of generated muons.
	  for (unsigned ilep = 0; ilep < pdi->numberOfDaughters(); ilep++)
	    Origin[1]->Fill(getOrigin(motherId(dileptonDaughter(*pdi, ilep))));
	}
	else if (i_rec <= l3) NumDilVsRec->Fill(i_rec+1);
	
	// Plot sum Pt in cone of dR=0.3 for dilepton leptons
	if (!doingElectrons && i_rec >= lgmr && passTrigger())
	  for (unsigned ilep = 0; ilep < pdi->numberOfDaughters(); ilep++) {
	    const reco::Muon& mu
	      = toConcrete<reco::Muon>(dileptonDaughter(*pdi, ilep));
	    if (mu.isIsolationValid())
	      SumPtR03[i_rec][1]->Fill(mu.getIsolationR03().sumPt);
	  }

	if (debug) {
	  LogTrace("Zprime2muResolution") << recLevelHelper.levelName(i_rec) << " values:";
	  LogTrace("Zprime2muResolution")
	    << "#|Charge |   Eta   |   Phi   |    P    |"
	    << "    Pt   |    Pz   |   Rap   |  Mass  ";
	  LogTrace("Zprime2muResolution")
	    << "-------------------------------------------------"
	    << "------------------------------";
	}

	// Get number of muons in events with dilepton.
	if (i_rec <= l3)
	  NMuons[i_rec][1]->Fill((double)allLeptons[i_rec].size());

	// Dilepton invariant mass.
	DilMass[i_rec]->Fill(pdi->mass());
	DilMassVsEta[i_rec]->Fill(pdi->eta(), pdi->mass());
	DilMassVsY[i_rec]->Fill(pdi->rapidity(), pdi->mass());

	// L1 and L2 electrons have no charge (they are just
	// superclusters); skip the below which depends on being able
	// to find the daughter by charge.
	if (doingElectrons && (i_rec == l1 || i_rec == l2))
	  continue;

	const reco::CandidateBaseRef& mum = dileptonDaughterByCharge(*pdi, -1);
	const reco::CandidateBaseRef& mup = dileptonDaughterByCharge(*pdi, +1);

	// A few mu- vs mu+ scatter plots
	MuMVsMuP[i_rec][0]->Fill(mup->eta(), mum->eta());
	MuMVsMuP[i_rec][1]->Fill(mup->phi(), mum->phi()); // JMTBAD phi
	MuMVsMuP[i_rec][2]->Fill(mup->pt(),  mum->pt());

	// Loop over "particle types": 0 stands for mu-, 1 - for mu+,
	// 2 - for opposite sign dilepton.
	for (int i_part = 0; i_part < 3; i_part++) {
	  const reco::Candidate* cand = 0;
	  // Split histos into mu-, and mu+, and dilepton.
	  if (i_part == 0)      cand = &*mum;
	  else if (i_part == 1) cand = &*mup;
	  else if (i_part == 2) cand = &*pdi;

	  if (debug) {
	    std::ostringstream strstrm;
	    strstrm << setprecision(5) << i_part << "| ";
	    if (i_part < 2)
	      strstrm << setw(4) << cand->charge() << "  | ";
	    else
	      strstrm << "      | ";
	    strstrm << endl << cand->p4();
	    LogTrace("Zprime2muAnalysis") << strstrm.str();
	  }

	  Eta[i_rec][i_part]->Fill(cand->eta());
	  Phi[i_rec][i_part]->Fill(cand->phi());
	  P[i_rec][i_part]->Fill(cand->p());
	  Pt[i_rec][i_part]->Fill(cand->pt());
	  Pz[i_rec][i_part]->Fill(fabs(cand->pz()));
	  Rapidity[i_rec][i_part]->Fill(cand->rapidity());
	  PVsEta[i_rec][i_part]->Fill(cand->eta(), cand->p());
	  PtVsEta[i_rec][i_part]->Fill(cand->eta(), cand->pt());

	  //int id = pdi->id();
	  // Calculate the resolution of each rec muon.
	  if (i_rec > 0 && i_part < 2) {
	    // Compare with generated values.
	    if (dileptons.size() == allDileptons[lgen].size()) {
	      const reco::CandidateBaseRef& murec = i_part == 0 ? mum : mup;
	      int charge = i_part == 0 ? -1 : 1;
	      const reco::CandidateBaseRef& mugen =
		dileptonDaughterByCharge(allDileptons[lgen][idi], charge);

	       // Make sure the tracks are reasonably close together.
	      if (matchTracks(murec, mugen)) {
		// Fill resolution histos.
		double pres = murec->p()/mugen->p() - 1.;
		EtaRes[i_rec]->Fill(murec->eta()-mugen->eta());
		PhiRes[i_rec]->Fill(murec->phi()-mugen->phi());
		PtDiff[i_rec]->Fill(murec->pt()-mugen->pt());
		PRes[i_rec]->Fill(pres);
		PResVsP[i_rec]->Fill(mugen->p(), pres*pres);
		if (i_rec <= l3) {
		  GenPhiResVsPhi[i_rec-1]->Fill(mugen->phi(),
						fabs(murec->phi()-mugen->phi()));
		  GenInvPtRes[i_rec-1]->Fill(((1./murec->pt())-(1./mugen->pt()))
					     /(1./mugen->pt()));
		  //GenInvPtResVsPt[i_rec-1]->Fill(mugen->pt(),
		  //      (fabs(1./murec->pt()-1./mugen->pt()))/(1./mugen->pt()));
		  GenInvPtResVsPt[i_rec-1]->Fill(mugen->pt(),
						 (fabs(murec->pt()-mugen->pt()))/(mugen->pt()));
		  GenInvPRes[i_rec-1]->Fill(mugen->p()/murec->p() - 1.);
		  GenPResVsPt[i_rec-1]->Fill(mugen->pt(), fabs(pres));
		  GenInvPResVsPt[i_rec-1]->Fill(mugen->pt(),
						fabs(mugen->p()/murec->p()-1.));
		  GenEtaResScat[i_rec-1]->Fill(mugen->eta(), murec->eta());
		  GenPhiResScat[i_rec-1]->Fill(mugen->phi(), murec->phi());
		  GenPtResScat[i_rec-1]->Fill(mugen->pt(),   murec->pt());
		}
	      }
	    }

	    if ((i_rec == l2 || i_rec == l3) &&
		dileptons.size() == allDileptons[l1].size() &&
		!doingElectrons) {
	      // Compare rec 2 and rec 3 values with rec 1 values.
	      // When doingElectrons the dileptonDaughterByCharge will
	      // not work, since L1 and L2 electrons as implemented do
	      // not have charges (they are just superclusters).
	      const reco::CandidateBaseRef& murec = i_part == 0 ? mum : mup;
	      int charge = i_part == 0 ? -1 : 1;
	      const reco::CandidateBaseRef& mul1 =
		dileptonDaughterByCharge(allDileptons[l1][idi], charge);

	      if (matchTracks(murec, mul1)) {
		// Fill resolution histos.
		L1EtaRes[i_rec-2]->Fill(murec->eta()-mul1->eta());
		L1PhiRes[i_rec-2]->Fill(murec->phi()-mul1->phi());
		L1PtDiff[i_rec-2]->Fill(murec->pt() -mul1->pt());
		L1EtaResScat[i_rec-2]->Fill(mul1->eta(), murec->eta());
		L1PhiResScat[i_rec-2]->Fill(mul1->phi(), murec->phi());
		L1PtResScat[i_rec-2]->Fill(mul1->pt(),   murec->pt());
	      }
	      if (i_rec == l3 && dileptons.size() == allDileptons[l2].size()) {
		// Compare rec 2 with rec 3 values.
		const reco::CandidateBaseRef& murec = i_part == 0 ? mum : mup;
		int charge = i_part == 0 ? -1 : 1;
		const reco::CandidateBaseRef& mul2 =
		  dileptonDaughterByCharge(allDileptons[l2][idi], charge);

		if (matchTracks(murec, mul2)) {
		  L2EtaRes->Fill(murec->eta()-mul2->eta());
		  L2PhiRes->Fill(murec->phi()-mul2->phi());
		  L2PtDiff->Fill(murec->pt() -mul2->pt());
		  L2EtaResScat->Fill(mul2->eta(), murec->eta());
		  L2PhiResScat->Fill(mul2->phi(), murec->phi());
		  L2PtResScat->Fill(mul2->pt(),   murec->pt());
		}
	      }
	    }
	  }
	  // Opposite-sign dileptons
	  else if ((i_rec > 0 && i_rec <= l3) && i_part == 2) {
	    if (dileptons.size() == allDileptons[lgen].size()) {
	      double genmass = allDileptons[lgen][idi].mass();
	      double dM   = cand->mass() - genmass;
	      double dMoM = dM/genmass;
	      GenDilMassRes[i_rec-1]->Fill(dM);
	      GenDilMassFrRes[i_rec-1]->Fill(dMoM);
	      GenDilMassResScat[i_rec-1]->Fill(genmass, dM*dM);
	      GenDilMassFrResScat[i_rec-1]->Fill(genmass, dMoM*dMoM);
	    }
	  }
	}
      }
    }
  }

  // Dilepton tests and resolution
  fillDilResHistos(debug);

  if (doingHiggs) { //SAMPLE_INDEX == kH4mu130) { // H -> ZZ* -> 4mu.
    if (allDileptons[lgen].size() != 2) // generated
      edm::LogWarning("Zprime2muResolution") 
	<< "+++ Warning: " << allDileptons[lgen].size()
	<< " generated dilepton(s) were found +++\n"; 
    else {
      ZonDilMass[0]->Fill(allDileptons[lgen][0].mass());
      ZofDilMass[0]->Fill(allDileptons[lgen][1].mass());
    }
    for (pdi = allDileptons[lgen].begin(); pdi != allDileptons[lgen].end();
	 pdi++) {
      const reco::CandidateBaseRef& lep0 = dileptonDaughter(*pdi, 0);
      const reco::CandidateBaseRef& lep1 = dileptonDaughter(*pdi, 1);
      Eta4muons->Fill(lep0->eta());
      Eta4muons->Fill(lep1->eta());
      Pt4muons->Fill(lep0->pt());
      Pt4muons->Fill(lep1->pt());
    }

    // For events passing the trigger
    if (passTrigger()) {
      if (allDileptons[l3].size() >= 2) { // HLT
	ZonDilMass[3]->Fill(allDileptons[l3][0].mass());
	ZofDilMass[3]->Fill(allDileptons[l3][1].mass());
      }
    }
  }
}

void Zprime2muResolution::fillEffHistos() {
  // Efficiency to reconstruct muons from Z' decays at various trigger 
  // levels and by various off-line reconstructors.  L1/HLT efficiencies
  // are not included.
  unsigned ilep, jlep;
  double gen_eta, gen_phi, gen_pt;

  // Efficiency for single muons
  for (ilep = 0; ilep < allLeptons[lgen].size(); ilep++) {
    const reco::CandidateBaseRef& lep = allLeptons[lgen][ilep];
    if (isResonance(motherId(lep))) {
      gen_eta = lep->eta();
      gen_phi = lep->phi();
      gen_pt  = lep->pt();

      EffVsEta[0]->Fill(gen_eta);
      EffVsPhi[0]->Fill(gen_phi);
      EffVsPt[0]->Fill(gen_pt);

      for (int i = l1; i < MAX_LEVELS; i++) {
	const int closestId = id(closestLepton(lep, i));
	// if (passTrigger(i) && closestId >= 0) {
	if (closestId >= 0) {
	  EffVsEta[i]->Fill(gen_eta);
	  EffVsPhi[i]->Fill(gen_phi);
	  EffVsPt[i]->Fill(gen_pt);
	}
      }

      for (jlep = 0; jlep < bestLeptons.size(); jlep++) {
	const reco::CandidateBaseRef& best = bestLeptons[jlep];
	const int rec_level = recLevel(best);
	const reco::CandidateBaseRef& match = closestLepton(best, rec_level);
	const int closestId = id(closestLepton(lep, rec_level)); 
	if (closestId >= 0 && closestId == id(match)) {
	  EffVsEta[MAX_LEVELS]->Fill(gen_eta);
	  EffVsPhi[MAX_LEVELS]->Fill(gen_phi);
	  EffVsPt[MAX_LEVELS]->Fill(gen_pt);
	  //break; // some gen muons can be matched to > 1 rec muon
	}
      }
    }
  }

  // Dimuon reconstruction efficiencies
  double gen_mass = -999.;
  int wheredi = 0;
  if (allDileptons[lgen].size() == 1) { // one generated dimuon
    gen_mass = allDileptons[lgen][0].mass()/1000.; // (in TeV)
    wheredi = int(whereIsDilepton(allDileptons[lgen][0]));
    RecMass[0][0]->Fill(gen_mass);

    // Events with both muons inside the full eta coverage, and passing
    // the trigger
    if (passTrigger() &&
	numDaughtersInAcc(allDileptons[lgen][0], ETA_CUT) >= 2)
      RecMass[0][1]->Fill(gen_mass);
  }

  if (passTrigger()) {
    if (allLeptons[l3].size() > 1)     // At least two muons found by L3
      RecMass[1][0]->Fill(gen_mass);
    if (allDileptons[l3].size() > 0)   // Opposite-sign dimuon found by L3
      RecMass[1][1]->Fill(gen_mass);
    if (allLeptons[lgmr].size() > 1)   // At least two muons found by GMR
      RecMass[2][0]->Fill(gen_mass);
    if (allDileptons[lgmr].size() > 0) // Opposite-sign dimuon found by GMR
      RecMass[2][1]->Fill(gen_mass);
    if (bestLeptons.size() > 1)        // At least two muons found by OPT
      RecMass[3][0]->Fill(gen_mass);
    if (bestDileptons.size() > 0) {    // Opposite-sign dimuon found by OPT
      RecMass[3][1]->Fill(gen_mass);
      RecMassByLoc[1][wheredi]->Fill(gen_mass);
    }

    // Count the number of leptons that pass cuts at OPT level
    int passCut = 0;
    for (int z = 1; z <= 3; z++) {
      passCut = 0;
      const LeptonRefVector& leps
	= z < 3 ? allLeptons[z == 1 ? l3 : lgmr] : bestLeptons;
      LeptonRefVector::const_iterator plep;
      for (plep = leps.begin(); plep != leps.end(); plep++) {
	if (!leptonIsCut(*plep)) passCut++;
	if (passCut > 1) break;
      }

      if (passCut > 1)
	RecMass[z][2]->Fill(gen_mass);
    }
      
    // Fill this also separately for barrel, endcap, overlap region
    // combinations, at least two leptons passing cuts
    if (passCut > 1)
      RecMassByLoc[0][wheredi]->Fill(gen_mass);
  }
}

void Zprime2muResolution::fillPtResHistos(const bool debug) {
  // Function to fill histos of pt resolutions.
  double gen_pt, gen_p, gen_eta, residual, pt_pull, tmperr;
  unsigned ilep;

  // Inverse momentum and 1/pT resolution for GMR and OPT, separately for
  // barrel, endcap, and their overlap.
  if (passTrigger()) { // event passed the trigger
    for (int i = 0; i < 2; i++) {
      const LeptonRefVector& leptons = i == 0 ? allLeptons[lgmr] : bestLeptons;
      for (ilep = 0; ilep < leptons.size(); ilep++) {
	const reco::CandidateBaseRef& lep = leptons[ilep];
	// Find the generated muon closest to the reconstructed muon.
	const reco::CandidateBaseRef& gen_lep = matchedLepton(lep, lgen);
	if (!gen_lep.isNull()) {
	  gen_eta = gen_lep->eta();
	  gen_p   = gen_lep->p();
	  gen_pt  = gen_lep->pt();
	  double invpres = gen_p/lep->p() - 1.;
	  double invptres = gen_pt/lep->pt() - 1.;
	  // if (gen_p > 1700. && gen_p < 2300.) {
	  if (fabs(gen_eta) < 0.9) { // barrel
	    MuInvPRes[i][0]->Fill(invpres);
	    MuInvPtRes[i][0]->Fill(invptres); 
	  }
	  else if (fabs(gen_eta) < 1.2) { // overlap
	    MuInvPRes[i][1]->Fill(invpres);
	    MuInvPtRes[i][1]->Fill(invptres);
	  }
	  else { // endcap
	    MuInvPRes[i][2]->Fill(invpres);
	    MuInvPtRes[i][2]->Fill(invptres);
	  }
	  // }
	}
      }
    }
  }

  // First plot histos of 1/pT resolution vs eta for a given muon.  Here
  // residual = ((1/PtL3)-(1/pT[0]))/(1/pT[0]).  We want to compare this
  // with the histogram given in pg. 27 of the Muon TDR.
  for (ilep = 0; ilep < bestLeptons.size(); ilep++) {
    const reco::CandidateBaseRef& lep = bestLeptons[ilep];
    if (lep->pt() < 100.) continue;
    // Find the generated muon closest to the reconstructed muon.
    const reco::CandidateBaseRef& gen_lep = matchedLepton(lep, lgen);
    if (!gen_lep.isNull()) {
      gen_pt  = gen_lep->pt();
      gen_eta = gen_lep->eta();
      //if ((gen_pt > 800) && (gen_pt < 1200)) {
	//if ((fabs(gen_eta) > 1.2) && (fabs(gen_eta) < 2.1)) {
	//if ((fabs(gen_eta) > 0.9) && (fabs(gen_eta) < 1.2)) {
	//if (fabs(gen_eta) < 0.9) {
	MuPtDiff[1]->Fill(lep->pt() - gen_pt);
	//}
      //}
    }
  }

  for (ilep = 0; ilep < allLeptons[l3].size(); ilep++) {
    const reco::CandidateBaseRef& lep = allLeptons[l3][ilep];
    // matchStudy(lep);

    // Only look at events where reconstructed pT > 100 GeV.
    if (lep->pt() < 100.) continue;

    // Find the generated muon closest to the L3 muon.
    const reco::CandidateBaseRef& gen_lep = matchedLepton(lep, lgen);
    if (!gen_lep.isNull()) {
      // Calculation of 1/pT residual.
      gen_pt  = gen_lep->pt();
      gen_eta = gen_lep->eta();
      residual = (1./lep->pt() - 1./gen_pt)/(1./gen_pt);
      //if ((gen_pt > 800) && (gen_pt < 1200)) {
	//if ((fabs(gen_eta) > 1.2) && (fabs(gen_eta) < 2.1)) {
	//if ((fabs(gen_eta) > 0.9) && (fabs(gen_eta) < 1.2)) {
        //if (fabs(gen_eta) < 0.9) {
      MuPtDiff[0]->Fill(lep->pt() - gen_pt);
        //}
      //}

      // Split histograms into gen_pt < 10, 25, 50, 75, 100, 1000.
      int idh = -1;
      if (gen_pt < 10.)        idh = 0;
      else if (gen_pt < 25.)   idh = 1;
      else if (gen_pt < 50.)   idh = 2;
      else if (gen_pt < 75.)   idh = 3;
      else if (gen_pt < 100.)  idh = 4;
      else if (gen_pt < 1000.) idh = 5;
      if (idh >= 0) ResidualPt[idh]->Fill(abs(gen_eta), abs(residual));

      if (debug) {
	LogTrace("Zprime2muResolution")
	  << " l3 muon = "         << id(lep)
	  << " closest genmuon = " << id(gen_lep);
	LogTrace("Zprime2muResolution")
	  << "   Ptgen: "  << gen_pt  << " Pt3: " << lep->pt()
	  << " Etagen: "   << gen_eta << " Residual: " << residual;
      }
    }
  }

  // Find the 1/pT residual and pulls for all alternative fits to L3.
  // Use same number of muons for all fits so that we can directly compare
  // the change in RMS, sigma and pulls between fits.  Loop over all of L3 
  // muons first.
  if (!passTrigger()) return;

  for (ilep = 0; ilep < allLeptons[l3].size(); ilep++) {
    const reco::CandidateBaseRef& lep = allLeptons[l3][ilep];
    // Only look at events where reconstructed pT > 100 GeV.
    if (lep->pt() < 100.) continue;

    // Get generated muon closest to the L3 muon.  This gen muon will be used
    // for comparison to all fits.
    const reco::CandidateBaseRef& gen_lep = matchedLepton(lep, lgen);

    // Find corresponding tracker-only, FMS and GMR tracks
    const reco::CandidateBaseRef& tk_lep = matchedLepton(lep, ltk);
    const reco::CandidateBaseRef& fms_lep = matchedLepton(lep, lfms);
    const reco::CandidateBaseRef& gmr_lep = matchedLepton(lep, lgmr);

    // Fill only for muons that have tracker-only and FMS tracks, and
    // are matched to the generated tracks
    if (!gen_lep.isNull() && !tk_lep.isNull() && !fms_lep.isNull()) {
      gen_pt = gen_lep->pt();
      int bar_end;

      // Residual and pulls for on-line and off-line reconstructions.
      // separate each fit by barrel or endcap
      if (abs(lep->eta()) < ENDCAP_BARREL_CUT) bar_end = 0;
      else bar_end = 1;
      if (lep->pt() > 0.) {
	residual  = (1./lep->pt() - 1./gen_pt)/(1./gen_pt);
	TotInvPtRes[0]->Fill(residual);
	InvPtRes[0][bar_end]->Fill(residual);
	if ((tmperr = invPtError(lep)) > 0.) {
	  pt_pull = (1./lep->pt() - 1./gen_pt)/tmperr;
	  TotInvPtPull[0]->Fill(pt_pull);
	  InvPtPull[0][bar_end]->Fill(pt_pull);
	}
      }

      // Residual and pulls for tracker-only fit.
      if (abs(tk_lep->eta()) < ENDCAP_BARREL_CUT) bar_end = 0;
      else bar_end = 1;
      if (tk_lep->pt() > 0.) {
	residual  = (1./tk_lep->pt() - 1./gen_pt)/(1./gen_pt);
	TotInvPtRes[1]->Fill(residual);
	InvPtRes[1][bar_end]->Fill(residual);
	if ((tmperr = invPtError(tk_lep)) > 0.) {
	  pt_pull = (1./tk_lep->pt() - 1./gen_pt)/tmperr;
	  TotInvPtPull[1]->Fill(pt_pull);
	  InvPtPull[1][bar_end]->Fill(pt_pull);
	}
      }

      if (!gmr_lep.isNull()) {
      // Residual and pulls for internal-seeded GMR.
	if (abs(gmr_lep->eta()) < ENDCAP_BARREL_CUT) bar_end = 0;
	else bar_end = 1;
	if (gmr_lep->pt() > 0.) {
	  residual  = (1./gmr_lep->pt() - 1./gen_pt)/(1./gen_pt);
	  InvPtRes[3][bar_end]->Fill(residual);
	  if ((tmperr = invPtError(gmr_lep)) > 0.) {
	    pt_pull = (1./gmr_lep->pt() - 1./gen_pt)/tmperr;
	    InvPtPull[3][bar_end]->Fill(pt_pull);
	  }
	}
      }

#if 0 // JMTBAD no forward measurement 
      // pT measurement at outer surface of tracker.
      // (Should not be the default)
      if (lep->forwardPt() > 0.) {
	residual  = (1./lep->forwardPt()-1./gen_pt)/(1./gen_pt);
	TotInvPtRes[3]->Fill(residual);
	if (lep->errForwardInvPt() > 0.) {
	  pt_pull = (1./lep->forwardPt()-1./gen_pt)/lep->errForwardInvPt();
	  TotInvPtPull[3]->Fill(pt_pull);
	}
      }
#endif

      // Residual and pulls using the innermost muon hit for
      // pT measurement.  This is alleged to be the best measurement
      // of pT in the detector.
      if (abs(fms_lep->eta()) < ENDCAP_BARREL_CUT) bar_end = 0;
      else bar_end = 1;
      if (fms_lep->pt() > 0.) {
	residual  = (1./fms_lep->pt() - 1./gen_pt)/(1./gen_pt);
	TotInvPtRes[2]->Fill(residual);
	InvPtRes[2][bar_end]->Fill(residual);
	if ((tmperr = invPtError(fms_lep)) > 0.) {
	  pt_pull = (1./fms_lep->pt() - 1./gen_pt)/tmperr;
	  TotInvPtPull[2]->Fill(pt_pull);
	  InvPtPull[2][bar_end]->Fill(pt_pull);
	}
      }
    }
  }
}

// This function makes histos of the difference between charge assignment
// for various levels of reconstruction.
void Zprime2muResolution::fillChargeResHistos(const bool debug) {
  // Loop over all reconstructed levels. L1 and L2 "electrons" are
  // superclusters, and do not have a charge, so skip them.
  for (int rec = doingElectrons ? l3 : l1; rec <= MAX_LEVELS; rec++) {
    const LeptonRefVector& leptons = getLeptons(rec);
    for (unsigned ilep = 0; ilep < leptons.size(); ilep++) {
      // Find charge of closest gen muon and store difference in histogram
      const reco::CandidateBaseRef& lep = leptons[ilep];
      const reco::CandidateBaseRef& gen_lep = matchedLepton(lep, lgen);
      if (!gen_lep.isNull()) {
	int deltaQ = lep->charge() - gen_lep->charge();
	QRes[rec]->Fill(deltaQ);
	if (rec <= l3 && rec < MAX_LEVELS) {
	  double pt = lep->pt();
	  double p  = lep->p();
	  int idh;
	  if      (deltaQ == 0)      idh =  0; // correct charge assignment
	  else if (abs(deltaQ) == 2) idh =  1; // wrong assignment
	  else throw cms::Exception("Zprime2muResolution")
	    << "+++ Impossible deltaQ in fillChargeResHistos() +++\n";
	  QResVsPt[rec-1][idh]->Fill(pt);
	  QResVsInvPt[rec-1][idh]->Fill(1./pt);
	  QResVsP[rec-1][idh]->Fill(p);
	  QResVsInvP[rec-1][idh]->Fill(1./p);
	}
	if (debug) {
	  LogTrace("Zprime2muResolution")
	    << "Charge at rec level " << rec << ": " << lep->charge()
	    << ", Closest Gen: " << gen_lep->charge();
	}
      }
    }
  }
}

void Zprime2muResolution::fillMuonHistos(const int rec, const bool debug) {
  // Function for filling eta, y, phi, p, pz, pt, pt vs eta, p vs eta
  // for all muons for all levels of reconstruction.
  // Inputs: rec   = level of reconstruction
  //         debug = print statements

  if (rec > MAX_LEVELS) {
    throw cms::Exception("Zprime2muResolution")
      << "+++ Unknown rec. level = " << rec << " in fillMuonHistos() +++\n";
  }

  if (debug) {
    ostringstream out;
    out << "****************************************\n";
    out << "In fillMuonHistos: rec = " << rec << endl;
    out << "#|Charge |   Eta   |   Phi   |    P    |"
	<< "    Pt   |    Pz   |   Rap   |  Mass  " << endl;
    out << "-------------------------------------------------"
	<< "------------------------------";
    LogTrace("Zprime2muResolution") << out.str();
  }

  const LeptonRefVector& leptons = getLeptons(rec);
  
  for (unsigned ilep = 0; ilep < leptons.size(); ilep++) {
    const reco::CandidateBaseRef& lep = leptons[ilep];

    MuonEta[rec]->Fill(lep->eta());
    MuonRap[rec]->Fill(lep->rapidity());
    MuonPhi[rec]->Fill(lep->phi());
    MuonPt[rec]->Fill(lep->pt());
    MuonPz[rec]->Fill(abs(lep->pz()));
    MuonP[rec]->Fill(lep->p());
    MuonPtVsEta[rec]->Fill(lep->eta(), lep->pt());
    MuonPVsEta[rec]->Fill( lep->eta(), lep->p() );

    if (!doingElectrons && rec >= lgmr) {
      const reco::Muon& mu = toConcrete<reco::Muon>(lep);
      if (mu.isIsolationValid())
	SumPtR03[rec][0]->Fill(mu.getIsolationR03().sumPt);
    }

    if (debug)
      LogTrace("Zprime2muResolution") << setprecision(5) << id(lep) << "| "
				      << setw(4) << lep->charge() << "  | "
				      << endl << lep->p4();
  }
}

void Zprime2muResolution::fillSignOfDilepton(const int rec, const bool debug) {
  // Small function for filling histogram which contains number of opposite
  // sign dileptons, + like-sign dileptons, and - like-sign dileptons.  The
  // routine is called for each level of reconstruction.
  // Inputs: rec          = level of reconstruction
  //         debug        = true for print statements

  int total_charge = 0;
  int nLeptons = int(allLeptons[rec].size());
  LeptonRefVector::const_iterator plep;
  ostringstream out;

  // Only enter routine if there are 2 or more muons.
  // Add up total charge of all muons.
  if (nLeptons < 2) return;

  for (plep = allLeptons[rec].begin(); plep != allLeptons[rec].end(); plep++)
    total_charge += (*plep)->charge();
  
  if (debug)
    out << "fillSignOfDilepton: rec level = " << rec
	<< " numLeptons = " << nLeptons
	<< " total_charge = " << total_charge << "; ";
  
  // If total charge is equal to number of muons, then positive dilepton
  // could be found
  if (total_charge == nLeptons) {
    SignOfDil[rec]->Fill(2.);
    if (debug) out << "positive same-sign dilepton found";
  }
  // If total charge is equal to -(number of muons), then negative dilepton
  // could be found
  else if (total_charge == -nLeptons) {
    SignOfDil[rec]->Fill(1.);
    if (debug) out << "negative same-sign dilepton found";
  }
  // Otherwise opposite sign dilepton was found
  else {
    SignOfDil[rec]->Fill(0.);
    if (debug) out << "opposite-sign dilepton found";
  }

  if (debug) LogTrace("fillSignOfDilepton") << out.str();
}

void Zprime2muResolution::fillDilResHistos(const bool debug) {
  // Histograms to compare dilepton mass resolution between various fits.
  unsigned int n_gen = allDileptons[lgen].size();

  LorentzVector genDimuV, recDimuV, genResV, recResV;

  const double mass_min = peakMass - 0.1*peakMass;
  const double mass_max = peakMass + 0.1*peakMass;

  // Search for dileptons at all trigger levels and for all off-line
  // reconstruction methods.
  for (int rec = 0; rec <= MAX_LEVELS; rec++) {
    // Check whether the event passed the trigger.
    if (rec >= l1 && rec <= l3) {
      if (cutTrig[rec] && !passTrigger(rec)) return;
    }

    const reco::CompositeCandidateCollection& dileptons = getDileptons(rec);
    unsigned int n_dil = dileptons.size();

    // Highest mass reconstructed at various trigger levels and by various
    // off-line fitting algorithms
    if (n_dil > 0) {
      DilMassComp[rec][0]->Fill(dileptons[0].mass());
      DilMassComp[rec][1]->Fill(resV(rec, 0).mass());
    }

    // Only continue on to calculate resolutions if considering
    // reconstructed dileptons.
    if (rec == lgen) continue;

    for (unsigned int i_dil = 0; i_dil < n_dil; i_dil++) {
      // Very poor match of dileptons; needs to be improved!
      if (n_dil == n_gen) {
	// Invariant mass resolution calculated relative to the mass
	// reconstructed from GEANT muons.
	genDimuV = allDileptons[lgen][i_dil].p4();
	recDimuV = dileptons[i_dil].p4();
	genResV  = resV(lgen, i_dil);
	recResV  = resV(rec,  i_dil);
	double geant_mass  = genDimuV.mass();
	double pythia_mass = genResV.mass();
	double dil_mass    = recDimuV.mass();
	double res_mass    = recResV.mass();

	// Extra smearing; only for tests!
	// Modify random generator seed if required.  The seed is set to
	// the current machine clock.
	//gRandom->SetSeed(0);
	//double epsil = 0.042;
	//double gsmear = gRandom->Gaus(0., epsil);
	//dil_mass *= (1.+gsmear);
	//res_mass *= (1.+gsmear);

	double dil_geant_resol  = (dil_mass - geant_mass)/geant_mass;
	double dil_pythia_resol = (dil_mass - pythia_mass)/pythia_mass;
	double res_pythia_resol = (res_mass - pythia_mass)/pythia_mass;

	// All events
	DilMassRes[rec][0]->Fill(dil_geant_resol);
	DilMassRes[rec][1]->Fill(dil_pythia_resol);
	DilMassRes[rec][2]->Fill(res_pythia_resol);

	// Events within +/-10% of the Z' mass peak
	if (geant_mass > mass_min  && geant_mass < mass_max) {
	  DilMassRes[rec][3]->Fill(dil_geant_resol);
	}
	if (pythia_mass > mass_min && pythia_mass < mass_max) {
	  DilMassRes[rec][4]->Fill(dil_pythia_resol);
	  DilMassRes[rec][5]->Fill(res_pythia_resol);
	}

	// An attempt to use method used in search for RS Gravitons by D0
	// (hep-ex/0505018, p.4).  They assign both muons the same value
	// of transverse momentum based on the weighted average (in 1/pT)
	// of their individual pT's.  They say this resulted in 30% decrease
	// in the RMS of the inv. mass distribution.
	if (rec == MAX_LEVELS) {
	  // Calculate reweighted pT.
	  const reco::CandidateBaseRef& lep0
	    = dileptonDaughter(dileptons[i_dil], 0);
	  const reco::CandidateBaseRef& lep1
	    = dileptonDaughter(dileptons[i_dil], 1);
	  double pTp = lep0->pt();
	  double pTm = lep1->pt();
	  double wp  = 1./pTp/(1./pTp + 1./pTm);
	  double wm  = 1./pTm/(1./pTp + 1./pTm);
	  double pTw = pTp*wp + pTm*wm;
	  // LogTrace("Zprime2muResolution")
          //   << "pTp = " << pTp << " pTm = " << pTm
	  //   << " wp = " << wp << " wm = " << wm
	  //   << " pTw = " << pTw << " " << 2./(1./pTp + 1./pTm);

	  // Construct new dimuon and calculate its mass.
	  LorentzVector vlep0, vlep1, vdil;
	  SetP4M(vlep0, pTw, lep0->phi(), lep0->p(), lep0->theta(), leptonMass);
	  SetP4M(vlep1, pTw, lep1->phi(), lep1->p(), lep1->theta(), leptonMass);
	  vdil = vlep0 + vlep1;
	  double D0mass = vdil.mass();
	  // LogTrace("Zprime2muResolution") << " Gen mass = " << geant_mass
	  //      << " rec mass = " << recDimuV.M()
	  //      << " D0 mass = " << D0mass;
	  double D0dil_geant_resol = (D0mass - geant_mass)/geant_mass;
	  DilMassRes[rec+1][0]->Fill(D0dil_geant_resol);
	}

	// The 1/pT resolution for the dilepton has a very strange shape.
	// I see 2 peaks, with the lower peak cutoff on the left at -1.
	// This is because the pT of the dilepton is a falling distribution
	// starting at zero.
	// double residual =
	//    (1./recDimuV.Pt() - 1./genDimuV.Pt())/(1./genDimuV.Pt());

	// The pT resolution is nicer; there is only one peak, but it is
	// still antisymmetric.
	DilPtRes[rec]->Fill((recDimuV.pt() - genDimuV.pt())/genDimuV.pt());
      }
    }
  }

  // Studies of pT errors for muons in poorly-reconstructed and
  // well-reconstructed dimuons.  It was used to determine quality cuts by
  // looking at events in the tail of Drell-Yan dimuon invariant mass
  // distribution.
  int qcut_mum = 0, qcut_mup = 0;
  if (bestDileptons.size() > 0) {
    if (DO_QCUTS ? 
	dilQCheck(bestDileptons[0], QSEL, qcut_mum, qcut_mup) : true) {
      // Require the reconstructed dilepton invariant mass be above a
      // certain mass value ("tail" of Drell-Yan distributions).
      double recm = resV(MAX_LEVELS, 0).mass(); // highest mass dilepton
      if (recm >= mass_min) {
	const reco::CandidateBaseRef& mum
	  = dileptonDaughterByCharge(bestDileptons[0], -1);
	const reco::CandidateBaseRef& mup
	  = dileptonDaughterByCharge(bestDileptons[0], +1);
	double err_pt_neg = ptError(mum)/mum->pt();
	double err_pt_pos = ptError(mup)/mup->pt();

	// Fill histos of relative pT error for mu+ vs mu- for two types of
	// dileptons.
	if (allDileptons[lgen].size() > 0) {
	  double genm = resV(lgen, 0).mass();

	  if ((recm-genm)/genm > 0.1) {
	    MuPVsMuM[0]->Fill(err_pt_neg, err_pt_pos);
	    if (debug) dumpDilQuality();
	  }
	  else {
	    MuPVsMuM[1]->Fill(err_pt_neg, err_pt_pos);
	  }
	}
      }
    }
  }
}

int Zprime2muResolution::getOrigin(const int motherId) {
  int origin = 0;
  // Group mother id into broad "origin" types defined similarly to
  // MuonSimtrackAnalyser.cc in ORCA.
  enum muonOrigin {undef=0, pi=1, K=2, KL=3, eta=4, rho=5,
		   c=6, b=7, tau=8, Z=10, W=11, H=12, Zprime=13, G=14};
  //edm::ESHandle<edm::ParticleDataTable> pdt;
  //evSetup.getData(pdt);
  switch (abs(motherId)) {
  case 211: origin = pi;  break;          // pi+/-
  case 321: origin = K;   break;          // K+/-
  case 311: origin = KL;  break;          // K0
  case 221: origin = eta; break;          // eta
  case 113: origin = rho; break;          // rho0
  case 411: case 421: case 431: case 443: // D+/-, D0, D_s+/-, j/psi
    origin = c; break;
  case 4122: case 100443:                 // Lambda_c+, psi'
    origin = c; break;
  case 511: case 521: case 531: case 553: // B0, B+/-, B_s0, Y
    origin = b; break;
  case 5122: case 100553:                 // Lambda_b0, Y'
    origin = b; break;
  case 15: origin = tau; break;           // tau
    // origin = 9 is blank
  case 23: origin = Z; break;             // Z0
  case 24: origin = W; break;             // W+/-
  case 25: origin = H; break;             // SM H
  case 32: origin = Zprime; break;        // Z'
  case 5000039: origin = G; break;        // G*
    // rest blank
  }
    
  return origin;
}

void Zprime2muResolution::getHistosFromFile() {
  TDirectory* histoDir = 0;
  histoFile->GetObject("Zprime2muResolution", histoDir);
  if (histoDir == 0)
    throw cms::Exception("Zprime2muResolution")
      << "getHistosFromFile failed; no TDirectory named Zprime2muResolution"
      << " in the file!\n";

  // JMTBAD there is undoubtedly a better way to do this...
  histoDir->GetObject("GenMassAllEvents", GenMassAllEvents);
  histoDir->GetObject("GenMassInAccept", GenMassInAccept);
  histoDir->GetObject("EventsInAccFailed", EventsInAccFailed);
  histoDir->GetObject("L1TrigFailSingleMu", L1TrigFailSingleMu);
  histoDir->GetObject("L1TrigFailMu2VsMu1", L1TrigFailMu2VsMu1);
  histoDir->GetObject("L1TrigPassSingleMu", L1TrigPassSingleMu);
  histoDir->GetObject("L1TrigPassMu2VsMu1", L1TrigPassMu2VsMu1);
  histoDir->GetObject("L2MuonHits", L2MuonHits);
  for (int i = 0; i < 3; i++) {
    histoDir->GetObject(nameHist("GMRMuonHits", i).c_str(), GMRMuonHits[i]);
    histoDir->GetObject(nameHist("GMRChi2dof", i).c_str(), GMRChi2dof[i]);
  }
  histoDir->GetObject("L3TrackerHits", L3TrackerHits);
  histoDir->GetObject("NumDilVsRec", NumDilVsRec);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      histoDir->GetObject(nameHist("RecMass", i, j).c_str(), RecMass[i][j]);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 10; j++)
      histoDir->GetObject(nameHist("RecMassByLoc", i,j).c_str(), RecMassByLoc[i][j]);
  for (int i = 0; i <= MAX_LEVELS; i++) {
    histoDir->GetObject(nameHist("EffVsEta", i).c_str(), EffVsEta[i]);
    histoDir->GetObject(nameHist("EffVsPhi", i).c_str(), EffVsPhi[i]);
    histoDir->GetObject(nameHist("EffVsPt", i).c_str(), EffVsPt[i]);
  }
  for (int i = 0; i < NUM_REC_LEVELS; i++) {
    for (int j = 0; j < 3; j++) {
      histoDir->GetObject(nameHist("TrigResult", i, j).c_str(), TrigResult[i][j]);
      histoDir->GetObject(nameHist("TrigMass", i, j).c_str(), TrigMass[i][j]);
    }
    for (int j = 0; j < 4; j++)
      histoDir->GetObject(nameHist("NMuons", i, j).c_str(), NMuons[i][j]);
    histoDir->GetObject(nameHist("SignOfDil", i).c_str(), SignOfDil[i]);
    histoDir->GetObject(nameHist("ZonDilMass", i).c_str(), ZonDilMass[i]);
    histoDir->GetObject(nameHist("ZofDilMass", i).c_str(), ZofDilMass[i]);
  }
  for (int i = 0; i <= MAX_LEVELS; i++) {
    histoDir->GetObject(nameHist("AllDilMass", i).c_str(), AllDilMass[i]);
    histoDir->GetObject(nameHist("DilMass", i).c_str(), DilMass[i]);
    histoDir->GetObject(nameHist("DilMassVsEta", i).c_str(), DilMassVsEta[i]);
    histoDir->GetObject(nameHist("DilMassVsY", i).c_str(), DilMassVsY[i]);
    histoDir->GetObject(nameHist("MuonEta", i).c_str(), MuonEta[i]);
    histoDir->GetObject(nameHist("MuonPhi", i).c_str(), MuonPhi[i]);
    histoDir->GetObject(nameHist("MuonRap", i).c_str(), MuonRap[i]);
    histoDir->GetObject(nameHist("MuonP", i).c_str(), MuonP[i]);
    histoDir->GetObject(nameHist("MuonPt", i).c_str(), MuonPt[i]);
    histoDir->GetObject(nameHist("MuonPz", i).c_str(), MuonPz[i]);
    histoDir->GetObject(nameHist("MuonPVsEta", i).c_str(), MuonPVsEta[i]);
    histoDir->GetObject(nameHist("MuonPtVsEta", i).c_str(), MuonPtVsEta[i]);
    for (int j = 0; j < 3; j++) {
      histoDir->GetObject(nameHist("Eta", i, j).c_str(), Eta[i][j]);
      histoDir->GetObject(nameHist("Phi", i, j).c_str(), Phi[i][j]);
      histoDir->GetObject(nameHist("Rapidity", i, j).c_str(), Rapidity[i][j]);
      histoDir->GetObject(nameHist("P", i, j).c_str(), P[i][j]);
      histoDir->GetObject(nameHist("Pt", i, j).c_str(), Pt[i][j]);
      histoDir->GetObject(nameHist("Pz", i, j).c_str(), Pz[i][j]);
      histoDir->GetObject(nameHist("PVsEta", i, j).c_str(), PVsEta[i][j]);
      histoDir->GetObject(nameHist("PtVsEta", i, j).c_str(), PtVsEta[i][j]);
      histoDir->GetObject(nameHist("MuMVsMuP", i, j).c_str(), MuMVsMuP[i][j]);
    }
    for (int j = 0; j < 2; j++)
      histoDir->GetObject(nameHist("SumPtR03", i, j).c_str(), SumPtR03[i][j]);
  }
  histoDir->GetObject("Eta4muons", Eta4muons);
  histoDir->GetObject("Pt4muons", Pt4muons);
  for (int k = l1; k <= MAX_LEVELS; k++) {
    histoDir->GetObject(nameHist("EtaRes", k).c_str(), EtaRes[k]);
    histoDir->GetObject(nameHist("PhiRes", k).c_str(), PhiRes[k]);
    histoDir->GetObject(nameHist("PtDiff", k).c_str(), PtDiff[k]);
    histoDir->GetObject(nameHist("PRes", k).c_str(), PRes[k]);
    histoDir->GetObject(nameHist("PResVsP", k).c_str(), PResVsP[k]);
  }
  for (int k=0; k<3; k++) {
    histoDir->GetObject(nameHist("GenPhiResVsPhi", k).c_str(), GenPhiResVsPhi[k]);
    histoDir->GetObject(nameHist("GenInvPtRes", k).c_str(), GenInvPtRes[k]);
    histoDir->GetObject(nameHist("GenInvPtResVsPt", k).c_str(), GenInvPtResVsPt[k]);
    histoDir->GetObject(nameHist("GenInvPRes", k).c_str(), GenInvPRes[k]);
    histoDir->GetObject(nameHist("GenPResVsPt", k).c_str(), GenPResVsPt[k]);
    histoDir->GetObject(nameHist("GenInvPResVsPt", k).c_str(), GenInvPResVsPt[k]);
    histoDir->GetObject(nameHist("GenEtaResScat", k).c_str(), GenEtaResScat[k]);
    histoDir->GetObject(nameHist("GenPhiResScat", k).c_str(), GenPhiResScat[k]);
    histoDir->GetObject(nameHist("GenPtResScat", k).c_str(), GenPtResScat[k]);
    histoDir->GetObject(nameHist("GenDilMassRes", k).c_str(), GenDilMassRes[k]);
    histoDir->GetObject(nameHist("GenDilMassFrRes", k).c_str(), GenDilMassFrRes[k]);
    histoDir->GetObject(nameHist("GenDilMassResScat", k).c_str(), GenDilMassResScat[k]);
    histoDir->GetObject(nameHist("GenDilMassFrResScat", k).c_str(), GenDilMassFrResScat[k]);
  }
  histoDir->GetObject("AllDilMassRes", AllDilMassRes);
  for (int l=0; l<2; l++) {
    histoDir->GetObject(nameHist("Origin", l).c_str(), Origin[l]);
    histoDir->GetObject(nameHist("L1EtaRes", l).c_str(), L1EtaRes[l]);
    histoDir->GetObject(nameHist("L1PhiRes", l).c_str(), L1PhiRes[l]);
    histoDir->GetObject(nameHist("L1PtDiff", l).c_str(), L1PtDiff[l]);
    histoDir->GetObject(nameHist("L1EtaResScat", l).c_str(), L1EtaResScat[l]);
    histoDir->GetObject(nameHist("L1PhiResScat", l).c_str(), L1PhiResScat[l]);
    histoDir->GetObject(nameHist("L1PtResScat", l).c_str(), L1PtResScat[l]);
  }
  histoDir->GetObject("L2EtaRes", L2EtaRes);
  histoDir->GetObject("L2PhiRes", L2PhiRes);
  histoDir->GetObject("L2PtDiff", L2PtDiff);
  histoDir->GetObject("L2EtaResScat", L2EtaResScat);
  histoDir->GetObject("L2PhiResScat", L2PhiResScat);
  histoDir->GetObject("L2PtResScat", L2PtResScat);
  for (int i = 0; i < 2; i++) {
    histoDir->GetObject(nameHist("MuPtDiff", i).c_str(), MuPtDiff[i]);
    for (int j = 0; j < 3; j++) {
      histoDir->GetObject(nameHist("MuInvPRes", i, j).c_str(), MuInvPRes[i][j]);
      histoDir->GetObject(nameHist("MuInvPtRes", i, j).c_str(), MuInvPtRes[i][j]);
    }
  }
  for (int i_pt = 0; i_pt < 6; i_pt++)
    histoDir->GetObject(nameHist("ResidualPt", i_pt).c_str(), ResidualPt[i_pt]);
  for (int i = 0; i < 4; i++) {
    histoDir->GetObject(nameHist("TotInvPtRes", i).c_str(), TotInvPtRes[i]);
    histoDir->GetObject(nameHist("TotInvPtPull", i).c_str(), TotInvPtPull[i]);
    for (int j = 0; j < 2; j++) {
      histoDir->GetObject(nameHist("InvPtRes", i, j).c_str(), InvPtRes[i][j]);
      histoDir->GetObject(nameHist("InvPtPull", i, j).c_str(), InvPtPull[i][j]);
    }
  }
  for (int i = 0; i <= MAX_LEVELS; i++) {
    for (int j = 0; j < 2; j++)
      histoDir->GetObject(nameHist("DilMassComp", i, j).c_str(), DilMassComp[i][j]);
    if (i > 0) {
      for (int j = 0; j < 6; j++)
	histoDir->GetObject(nameHist("DilMassRes", i, j).c_str(), DilMassRes[i][j]);
      histoDir->GetObject(nameHist("DilPtRes", i).c_str(), DilPtRes[i]);
    }
  }
  for (unsigned int j = 0; j < 2; j++)
    histoDir->GetObject(nameHist("MuPVsMuM", j).c_str(), MuPVsMuM[j]);
  for (int i = l1; i <= MAX_LEVELS; i++)
    histoDir->GetObject(nameHist("QRes", i).c_str(), QRes[i]);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      histoDir->GetObject(nameHist("QResVsPt", i, j).c_str(), QResVsPt[i][j]);
      histoDir->GetObject(nameHist("QResVsInvPt", i, j).c_str(), QResVsInvPt[i][j]);
      histoDir->GetObject(nameHist("QResVsP", i, j).c_str(), QResVsP[i][j]);
      histoDir->GetObject(nameHist("QResVsInvP", i, j).c_str(), QResVsInvP[i][j]);
    }
  }
}

void Zprime2muResolution::DrawResHistos() {
  const string str_part[3] = {
    "(Opp-sign Dil) Mu-", "(Opp-sign Dil) Mu+", "Opp-sign Dilepton"
  };

  TText t;
  t.SetTextFont(12);
  t.SetTextSize(.03);

  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 700);
  TPostScript *ps = new TPostScript(outputFile.c_str(), 111);

  const int NUM_PAGES = 90;
  TPad *pad[NUM_PAGES];
  for (int i_page = 0; i_page <= NUM_PAGES; i_page++)
    pad[i_page] = new TPad("","", .05, .05, .95, .93);

  int page = 0;
  ostringstream strpage;
  string tit;
  TPaveLabel *title;

  // Origin of Muons.
  const int nx = 20; // mother codes from MuonReconstructionNtuple.cc
  char *mother[nx] = {"  ","pi","K ","K0","eta","rho","c ","b ","tau","  ",
		      "Z ","W ","H ","Z'","G*","  ","  ","  ","  ","  "};
  for (int ibin = 0; ibin < nx; ibin++) {
    Origin[0]->GetXaxis()->SetBinLabel(ibin+1,mother[ibin]);
    Origin[1]->GetXaxis()->SetBinLabel(ibin+1,mother[ibin]);
  }

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  title = new TPaveLabel(0.1,0.94,0.9,0.98,"Mother Particle Id of Muons");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  gStyle->SetOptLogy(1);
  pad[page]->Draw();
  pad[page]->Divide(1,2);
  pad[page]->cd(1);  Origin[0]->Draw();
  pad[page]->cd(2);  Origin[1]->Draw();
  gStyle->SetOptLogy(0);
  c1->Update();

  // Number of muons in each level.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Number of Muons");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2, 4);
  for (int i=0; i<4; i++) {
    pad[page]->cd(2*i+1);  gPad->SetLogy(1);  NMuons[i][0]->Draw();
    pad[page]->cd(2*i+2);  gPad->SetLogy(1);  NMuons[i][1]->Draw();
  }
  // pad[page]->cd(9);  NumDilVsRec->Draw();
  c1->Update();

  // Trigger bits
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Trigger Bits");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  for (int irec = 1; irec < NUM_REC_LEVELS; irec++) {
    pad[page]->cd(2*irec-1);  gPad->SetLogy(1);  TrigResult[irec][0]->Draw();
    pad[page]->cd(2*irec);                       NMuons[irec][2]->Draw();
  }
  c1->Update();

  // Muons failing and passing L1 Global Trigger
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon pT at L1");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,2);
  pad[page]->cd(1);  L1TrigFailSingleMu->Draw();
  pad[page]->cd(2);  L1TrigPassSingleMu->Draw();
  pad[page]->cd(3);  L1TrigFailMu2VsMu1->Draw("box");
  pad[page]->cd(4);  L1TrigPassMu2VsMu1->Draw("box");
  c1->Update();

  // Total number of muon and tracker hits (pixel plus silicon) at Level-3
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, "Muon and Tracker hits");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(1,2);
  pad[page]->cd(1);  L2MuonHits->Draw();
  pad[page]->cd(2);  L3TrackerHits->Draw();
  c1->Update();

  // Total number of muon hits at GMR
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, "Muon hits");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,2);
  for (int i = 0; i < 3; i++) {
    pad[page]->cd(i+1);  GMRMuonHits[i]->Draw();
  }
  TH1F* GMRMuonHitsTot = (TH1F*)GMRMuonHits[0]->Clone();
  GMRMuonHitsTot->SetNameTitle("GMRMuonHitsTot", "Muon hits, GMR, total");
  GMRMuonHitsTot->Add(GMRMuonHits[0], GMRMuonHits[1], 1., 1.);
  GMRMuonHitsTot->Add(GMRMuonHitsTot, GMRMuonHits[2], 1., 1.);
  pad[page]->cd(4);  GMRMuonHitsTot->Draw();
  c1->Update();
  delete GMRMuonHitsTot;

  // Chi^2/d.o.f. at GMR
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98, "Chi^2/d.o.f.");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,2);
  for (int i = 0; i < 3; i++) {
    pad[page]->cd(i+1);  GMRChi2dof[i]->Draw();
  }
  TH1F* GMRChi2dofTot = (TH1F*)GMRChi2dof[0]->Clone();
  GMRChi2dofTot->SetNameTitle("GMRChi2dofTot", "Chi2/d.o.f., GMR, total");
  GMRChi2dofTot->Add(GMRChi2dof[0], GMRChi2dof[1], 1., 1.);
  GMRChi2dofTot->Add(GMRChi2dofTot, GMRChi2dof[2], 1., 1.);
  pad[page]->cd(4);  GMRChi2dofTot->Draw();
  c1->Update();
  delete GMRChi2dofTot;

  string histtitle;
  int nbins;
  Stat_t f_bin, ent_bin, err_bin;
  // some temp histograms for plotting ratios
  TH1F* ratio[2][6];

#ifdef LATER
  // Off-line reconstruction efficiencies for muons
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Off-line Efficiency for Muons from Z' Decays");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  gStyle->SetOptStat(110010);
  pad[page]->Draw();
  pad[page]->Divide(1,3);
  TH1F *effeta, *effphi, *effpt;
  // Off-line muon reconstruction efficiency as a function of true eta.
  // Only muons produced in Z' decays are used; L1/HLT efficiency is not
  // included.
  pad[page]->cd(1);  gPad->SetGrid(1);
  effeta = (TH1F*)EffVsEta[3]->Clone();
  effeta->SetNameTitle("effeta", "L3 muon efficiency vs eta");
  effeta->Divide(EffVsEta[3], EffVsEta[0], 1., 1., "B");
  effeta->SetMinimum(0.50);  effeta->SetMaximum(1.02);
  effeta->Draw();  effeta->Draw("samehisto");
  // Off-line muon reconstruction efficiency as a function of true phi.
  pad[page]->cd(2);  gPad->SetGrid(1);
  effphi = (TH1F*)EffVsPhi[3]->Clone();
  effphi->SetNameTitle("effphi", "L3 muon efficiency vs phi");
  effphi->Divide(EffVsPhi[3], EffVsPhi[0], 1., 1., "B");
  effphi->SetMinimum(0.50);  effphi->SetMaximum(1.02);
  effphi->Draw();  effphi->Draw("samehisto");
  // Off-line muon reconstruction efficiency as a function of true pT.
  pad[page]->cd(3);  gPad->SetGrid(1);
  effpt = (TH1F*)EffVsPt[3]->Clone();
  effpt->SetNameTitle("effpt", "L3 muon efficiency vs pT");
  effpt->Divide(EffVsPt[3], EffVsPt[0], 1., 1., "B");
  effpt->SetMinimum(0.50);  effpt->SetMaximum(1.02);
  effpt->Draw();  effpt->Draw("samehisto");
  c1->Update();
  delete effeta;
  delete effphi;
  delete effpt;
#endif

  // Reconstruction efficiencies for muons from Z' decays as a function of eta
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Efficiency for Muons from Z' Decays");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  gStyle->SetOptStat(10);
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  TH1F *effvseta[MAX_LEVELS+1];
  // Muon reconstruction efficiency (different reconstructors) as a
  // function of true eta. Only muons produced in Z' decays are used;
  // L1/HLT efficiency is not included.
  for (int irec = 1; irec <= MAX_LEVELS; irec++) {
    pad[page]->cd(irec);  gPad->SetGrid(1);
    effvseta[irec] = (TH1F*)EffVsEta[irec]->Clone();
    histtitle = recLevelHelper.levelName(irec) + " muon efficiency vs eta";
    effvseta[irec]->SetTitle(histtitle.c_str());
    effvseta[irec]->Divide(EffVsEta[irec], EffVsEta[0], 1., 1., "B");
    effvseta[irec]->SetMinimum(0.60);  effvseta[irec]->SetMaximum(1.02);
    effvseta[irec]->Draw();  effvseta[irec]->Draw("samehisto");
  }
  c1->Update();
  for (int irec = 1; irec <= MAX_LEVELS; irec++)
    delete effvseta[irec];

  // Trigger efficiency as a function of dimuon mass
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Trigger Efficiency");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(3,4);
  TH1F *rat[NUM_REC_LEVELS][3];

  gStyle->SetOptStat(111111);
  for (int j = 0; j < 3; j++) {
    pad[page]->cd(j+1);  gPad->SetGrid(1);
    TrigMass[0][j]->Draw("histo");
  }
  c1->Update();

  gStyle->SetOptStat(110010);
  for (int irec = l1; irec < NUM_REC_LEVELS; irec++) {
    for (int j = 0; j < 3; j++) {
      pad[page]->cd(3*irec+j+1);  gPad->SetGrid(1);
      rat[irec][j] = (TH1F*)TrigMass[irec][j]->Clone();
      if (j == 0)
	histtitle = "Trigger eff. (all events) vs mass, ";
      else if (j == 1)
	histtitle = "Trigger eff. (#eta < 2.4) vs mass, ";
      else if (j == 2)
	histtitle = "Trigger eff. (#eta < 2.1) vs mass, ";
      histtitle += recLevelHelper.levelName(irec);
      rat[irec][j]->SetTitle(histtitle.c_str());
      rat[irec][j]->Divide(TrigMass[irec][j], TrigMass[0][j], 1., 1., "B");
      rat[irec][j]->SetMinimum(0.89); rat[irec][j]->SetMaximum(1.01);
      rat[irec][j]->Draw();  rat[irec][j]->Draw("samehisto");
    }
  }
  c1->Update();
  for (int irec = l1; irec < NUM_REC_LEVELS; irec++)
    for (int j = 0; j < 3; j++)
      delete rat[irec][j];

  // Reconstruction efficiencies as a function of dimuon mass
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Dimuon Reconstruction Efficiency");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  gStyle->SetOptStat(110010);
  pad[page]->Draw();
  pad[page]->Divide(3,4);
  TH1F *efftot[4], *eff2mu[4], *effdimu[4], *eff2mucuts[4];
  string level;
  for (int i = 1; i < 4; i++) {
    if (i == 1)      level = "L3";
    else if (i == 2) level = "GMR";
    else if (i == 3) level = "OPT";

    // Total - acceptance, trigger and "off-line" reconstruction - efficiency.
    pad[page]->cd(i);  gPad->SetGrid(1);
    efftot[i] = (TH1F*)RecMass[i][1]->Clone();
    histtitle = "Total efficiency vs mass, " + level;
    efftot[i]->SetTitle(histtitle.c_str());
    efftot[i]->Divide(RecMass[i][1], RecMass[0][0], 1., 1., "B");
    efftot[i]->SetMinimum(0.50);  efftot[i]->SetMaximum(1.02);
    efftot[i]->Draw();  efftot[i]->Draw("samehisto");

    // "Off-line" efficiency to find at least two muons, relative to events
    // in acceptance and accepted by the L1/HLT triggers.
    pad[page]->cd(i+3);  gPad->SetGrid(1);
    eff2mu[i] = (TH1F*)RecMass[i][0]->Clone();
    histtitle = "Two-mu efficiency vs mass, " + level;
    eff2mu[i]->SetTitle(histtitle.c_str());
    eff2mu[i]->Divide(RecMass[i][0], RecMass[0][1], 1., 1., "B");
    eff2mu[i]->SetMinimum(0.70);  eff2mu[i]->SetMaximum(1.02);
    eff2mu[i]->Draw();  eff2mu[i]->Draw("samehisto");

    // "Off-line" efficiency to find opposite-sign dimuon, relative to events
    // in acceptance and accepted by the L1/HLT triggers, and with two muons
    // passing cuts.
    pad[page]->cd(i+6);  gPad->SetGrid(1);
    eff2mucuts[i] = (TH1F*)RecMass[i][1]->Clone();
    histtitle = "Dimuon w/ cuts efficiency vs mass, " + level;
    eff2mucuts[i]->SetTitle(histtitle.c_str());
    eff2mucuts[i]->Divide(RecMass[i][1], RecMass[i][2], 1., 1., "B");
    eff2mucuts[i]->SetMinimum(0.70);  eff2mucuts[i]->SetMaximum(1.02);
    eff2mucuts[i]->Draw();  eff2mucuts[i]->Draw("samehisto");

    // "Off-line" efficiency to find opposite-sign dimuon, relative to events
    // in acceptance and accepted by the L1/HLT triggers.
    pad[page]->cd(i+9);  gPad->SetGrid(1);
    effdimu[i] = (TH1F*)RecMass[i][1]->Clone();
    histtitle = "Dimuon efficiency vs mass, " + level;
    effdimu[i]->SetTitle(histtitle.c_str());
    effdimu[i]->Divide(RecMass[i][1], RecMass[0][1], 1., 1., "B");
    effdimu[i]->SetMinimum(0.70);  effdimu[i]->SetMaximum(1.02);
    effdimu[i]->Draw();  effdimu[i]->Draw("samehisto");
  }
  c1->Update();
  for (int i = 1; i < 4; i++) {
    delete efftot[i];
    delete eff2mu[i];
    delete eff2mucuts[i];
    delete effdimu[i];
  }

  // Opposite-sign efficiencies as a function of dimuon mass,
  // separated by where the muons are
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Opposite-sign Reconstruction Efficiency (OPT)");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  gStyle->SetOptStat(110010);
  pad[page]->Draw();
  pad[page]->Divide(2,5);
  TH1F *effdimusep[10];
  ostringstream out;
  out << "Opp-sign efficiency:\n";
  for (int i = 0; i < 10; i++) {
    pad[page]->cd(i+1);  gPad->SetGrid(1);
    effdimusep[i] = (TH1F*)RecMassByLoc[1][i]->Clone();
    out << RecMassByLoc[1][i]->GetTitle()
	<< ": " << RecMassByLoc[1][i]->Integral()
	<< "/"  << RecMassByLoc[0][i]->Integral() << " = " 
	<< RecMassByLoc[1][i]->Integral()/RecMassByLoc[0][i]->Integral()
	<< endl;
    effdimusep[i]->Divide(RecMassByLoc[1][i], RecMassByLoc[0][i], 1., 1., "B");
    effdimusep[i]->SetMinimum(0.80);  effdimusep[i]->SetMaximum(1.02);
    effdimusep[i]->Draw();  effdimusep[i]->Draw("samehisto");
  }
  edm::LogInfo("Zprime2muResolution") << out.str();
  c1->Update();
  for (int i = 0; i < 10; i++)
    delete effdimusep[i];

  // Number of Dileptons with opposite sign and like sign
  char *genlabel[3] = {"total #","# found w/ eta cut"};
  for (int i_bin = 0; i_bin < 2; i_bin++) {
    SignOfDil[0]->GetXaxis()->SetBinLabel(i_bin+1,genlabel[i_bin]);
  }
  char *sign[3] = {"opposite","same neg","same pos"};
  for (int i_hist = 1; i_hist < NUM_REC_LEVELS; i_hist++) {
    for (int j_bin = 0; j_bin < 3; j_bin++) {
      SignOfDil[i_hist]->GetXaxis()->SetBinLabel(j_bin+1,sign[j_bin]);
    }
  }
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Sign of Dileptons");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  gStyle->SetOptStat(111111);
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  pad[page]->cd(1);  SignOfDil[0]->Draw();
  pad[page]->cd(2);  EventsInAccFailed->Draw();
  for (int i_rec = 1; i_rec < NUM_REC_LEVELS; i_rec++) {
    pad[page]->cd(2*i_rec+1);  SignOfDil[i_rec]->Draw();
    pad[page]->cd(2*i_rec+2);  NMuons[i_rec][3]->Draw();
  }
  c1->Update();

  // Draw Eta, Phi, Rap, momentum for all mu's for all rec levels
  for (int i=0; i<=MAX_LEVELS; i++) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = recLevelHelper.levelName(i) + " Values, all Muons";
    delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,4);
    pad[page]->cd(1);  MuonEta[i]->Draw();
    pad[page]->cd(2);  MuonPhi[i]->Draw();
    pad[page]->cd(3);  MuonRap[i]->Draw();
    pad[page]->cd(4);  if(i==0)         gPad->SetLogy(1); MuonP[i]->Draw();
    pad[page]->cd(5);  if(i==0 || i==1) gPad->SetLogy(1); MuonPt[i]->Draw();
    pad[page]->cd(6);  if(i==0)         gPad->SetLogy(1); MuonPz[i]->Draw();
    pad[page]->cd(7);  MuonPVsEta[i]->Draw();
    pad[page]->cd(8);  MuonPtVsEta[i]->Draw();
    c1->Update();
  }

  // Draw Eta, Phi, Rap, momentum for all levels Mu+ and Mu- associated
  // with an opp-sign dilepton.
  for (int i=0; i<=MAX_LEVELS; i++) {
    for (int j=0; j<2; j++) {
      ps->NewPage();
      c1->Clear();
      c1->cd(0);
      tit = str_part[j] + " " + recLevelHelper.levelName(i) + " Values";
      delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
      title->SetFillColor(10);
      title->Draw();
      strpage << "- " << (++page) << " -";
      t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
      pad[page]->Draw();
      pad[page]->Divide(2,4);
      pad[page]->cd(1);  Eta[i][j]->Draw();
      pad[page]->cd(2);  Phi[i][j]->Draw();
      pad[page]->cd(3);  Rapidity[i][j]->Draw();
      pad[page]->cd(4);  P[i][j]->Draw();
      pad[page]->cd(5);  if(i == 1) gPad->SetLogy(1);  Pt[i][j]->Draw();
      pad[page]->cd(6);  Pz[i][j]->Draw();
      pad[page]->cd(7);  PVsEta[i][j]->Draw();
      pad[page]->cd(8);  PtVsEta[i][j]->Draw();
      c1->Update();
    }
  }

  // Opposite sign dileptons, at Gen, L1, L2 and L3
  for (int i_rec = 0; i_rec <= MAX_LEVELS; i_rec++) {
    // Eta, Phi, Rap and momentum
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = str_part[2] + " " + recLevelHelper.levelName(i_rec) + " Values";
    delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,4);
    pad[page]->cd(1);  Eta[i_rec][2]->Draw();
    pad[page]->cd(2);  Phi[i_rec][2]->Draw();
    pad[page]->cd(3);  Rapidity[i_rec][2]->Draw();
    pad[page]->cd(4);  P[i_rec][2]->Draw();
    pad[page]->cd(5);  Pt[i_rec][2]->Draw();
    pad[page]->cd(6);  Pz[i_rec][2]->Draw();
    pad[page]->cd(7);  PVsEta[i_rec][2]->Draw();
    pad[page]->cd(8);  PtVsEta[i_rec][2]->Draw();
    c1->Update();

    // Mass, Mass vs Eta, Mass vs Phi, and mu+ vs mu- scatter plots
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    if (i_rec != 3)
      pad[page]->Divide(2,3);
    else
      pad[page]->Divide(2,4);

    int zone = 1;
    if (i_rec == 3) {
      // If off-line, first plot the mass distribution before the quality cuts
      pad[page]->cd(zone++);  AllDilMass[i_rec]->Draw();  zone++;
    }
    pad[page]->cd(zone++);  DilMass[i_rec]->Draw();

    // Various fits to mass distributions
    //binnedMassFits(i_rec);

    pad[page]->cd(zone++);  DilMassVsEta[i_rec]->Draw();
    pad[page]->cd(zone++);  DilMassVsY[i_rec]->Draw();
    pad[page]->cd(zone++);  MuMVsMuP[i_rec][0]->Draw();
    pad[page]->cd(zone++);  MuMVsMuP[i_rec][1]->Draw();
    pad[page]->cd(zone++);  MuMVsMuP[i_rec][2]->Draw();
    c1->Update();
  }

  // Eta Resolution in each Reconstruction level.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon Eta Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);  GenEtaResScat[0]->Draw();
  pad[page]->cd(2);  GenEtaResScat[1]->Draw();
  pad[page]->cd(3);  GenEtaResScat[2]->Draw();
  pad[page]->cd(4);  L1EtaResScat[0]->Draw();
  pad[page]->cd(5);  L1EtaResScat[1]->Draw();
  pad[page]->cd(6);  L2EtaResScat->Draw();
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon Eta Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i_rec = l1; i_rec <= MAX_LEVELS; i_rec++) {
    pad[page]->cd(i_rec);
    EtaRes[i_rec]->Draw();  EtaRes[i_rec]->Fit("gaus","Q");
  }
  c1->Update();

  // Phi Resolution in each level.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon Phi Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);  GenPhiResScat[0]->Draw();
  pad[page]->cd(2);  GenPhiResScat[1]->Draw();
  pad[page]->cd(3);  GenPhiResScat[2]->Draw();
  pad[page]->cd(4);  L1PhiResScat[0]->Draw();
  pad[page]->cd(5);  L1PhiResScat[1]->Draw();
  pad[page]->cd(6);  L2PhiResScat->Draw();
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon Phi Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i_rec = l1; i_rec <= MAX_LEVELS; i_rec++) {
    pad[page]->cd(i_rec);
    PhiRes[i_rec]->Draw();  PhiRes[i_rec]->Fit("gaus","Q");
  }
  c1->Update();

  // Pt Resolution in each level.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon pT Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);  GenPtResScat[0]->Draw();
  pad[page]->cd(2);  GenPtResScat[1]->Draw();
  pad[page]->cd(3);  GenPtResScat[2]->Draw();
  pad[page]->cd(4);  L1PtResScat[0]->Draw();
  pad[page]->cd(5);  L1PtResScat[1]->Draw();
  pad[page]->cd(6);  L2PtResScat->Draw();
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon pT Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i_rec = l1; i_rec <= MAX_LEVELS; i_rec++) {
    pad[page]->cd(i_rec);
    PtDiff[i_rec]->Draw();  PtDiff[i_rec]->Fit("gaus","Q");
  }
  c1->Update();

  // 1/Pt resolution
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon 1/pT Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);  GenInvPtRes[0]->Draw();  GenInvPtRes[0]->Fit("gaus","Q");
  pad[page]->cd(2);  GenInvPtResVsPt[0]->Draw();
  pad[page]->cd(3);  GenInvPtRes[1]->Draw();  GenInvPtRes[1]->Fit("gaus","Q");
  pad[page]->cd(4);  GenInvPtResVsPt[1]->Draw();
  pad[page]->cd(5);  GenInvPtRes[2]->Draw();  GenInvPtRes[2]->Fit("gaus","Q");
  pad[page]->cd(6);  GenInvPtResVsPt[2]->Draw();
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon 1/pT Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 2; i++) {
      pad[page]->cd(2*j+i+1);  MuInvPtRes[i][j]->Draw();
      MuInvPtRes[i][j]->Fit("gaus","Q");
    }
  }
  TH1F* MuInvPtResGMRTot = (TH1F*)MuInvPtRes[0][0]->Clone();
  MuInvPtResGMRTot->SetNameTitle("MuInvPtResGMRTot", "1/pT res, GMR, total");
  MuInvPtResGMRTot->Add(MuInvPtRes[0][0], MuInvPtRes[0][1], 1., 1.);
  MuInvPtResGMRTot->Add(MuInvPtResGMRTot, MuInvPtRes[0][2], 1., 1.);
  pad[page]->cd(7);
  MuInvPtResGMRTot->Draw();  MuInvPtResGMRTot->Fit("gaus","Q");
  TH1F* MuInvPtResOPTTot = (TH1F*)MuInvPtRes[1][0]->Clone();
  MuInvPtResOPTTot->SetNameTitle("MuInvPtResOPTTot", "1/pT res, OPT, total");
  MuInvPtResOPTTot->Add(MuInvPtRes[1][0], MuInvPtRes[1][1], 1., 1.);
  MuInvPtResOPTTot->Add(MuInvPtResOPTTot, MuInvPtRes[1][2], 1., 1.);
  pad[page]->cd(8);
  MuInvPtResOPTTot->Draw();  MuInvPtResOPTTot->Fit("gaus","Q");
  c1->Update();
  delete MuInvPtResGMRTot;
  delete MuInvPtResOPTTot;

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "1/pT Resolution For Different Fits");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      pad[page]->cd(2*i+j+1);
      InvPtRes[i][j]->Draw();  InvPtRes[i][j]->Fit("gaus","Q");
    }
  }
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"1/pT Pull For Different Fits");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      pad[page]->cd(2*i+j+1);
      InvPtPull[i][j]->Draw();  InvPtPull[i][j]->Fit("gaus","Q");
    }
  }
  c1->Update();

  // P resolution
  TH1F *PResVsPMod[MAX_LEVELS+1];
  for (int i_rec = l1; i_rec <= MAX_LEVELS; i_rec++) {
    histtitle = "Sqrt(Var((" + recLevelHelper.levelName(i_rec) + "-Gen P)/Gen P)) vs Gen P";
    nbins = PResVsP[i_rec]->GetNbinsX();
    PResVsPMod[i_rec] =
      fs->make<TH1F>(nameHist("PResVsPMod", i_rec).c_str(),
		     histtitle.c_str(), nbins,
		     PResVsP[i_rec]->GetXaxis()->GetXmin(),
		     PResVsP[i_rec]->GetXaxis()->GetXmax());
    for (int ibin = 1; ibin <= nbins; ibin++) {
      f_bin   = PResVsP[i_rec]->GetBinContent(ibin);
      ent_bin = PResVsP[i_rec]->GetBinEntries(ibin);
      if (f_bin > 0.) {
	f_bin = sqrt(f_bin);
	if (ent_bin > 0.) {err_bin = f_bin/sqrt(2.*ent_bin);}
	else              {err_bin = 0.;}
      }
      else    {f_bin = 0.; err_bin = 0.;}
      PResVsPMod[i_rec]->SetBinContent(ibin, f_bin);
      PResVsPMod[i_rec]->SetBinError(ibin, err_bin);
    }
  }

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,"Muon Momentum Resolution (I)");
  title->SetFillColor(10);  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  for (int i_rec = l1; i_rec <= l3; i_rec++) {
    pad[page]->cd(2*i_rec-1);  PRes[i_rec]->Draw();
    PRes[i_rec]->Fit("gaus","Q");
    pad[page]->cd(2*i_rec);    PResVsPMod[i_rec]->Draw();
  }
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,"Muon Momentum Resolution (II)");
  title->SetFillColor(10);  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,5);
  for (int i_rec = lgmr; i_rec <= MAX_LEVELS; i_rec++) {
    pad[page]->cd(2*i_rec-7);  PRes[i_rec]->Draw();
    PRes[i_rec]->Fit("gaus","Q");
    pad[page]->cd(2*i_rec-6);  PResVsPMod[i_rec]->Draw();
  }
  c1->Update();
  // The following has to be after c1->Update();
  for (int i_rec = l1; i_rec <= MAX_LEVELS; i_rec++)
    delete PResVsPMod[i_rec];

  // 1/P resolution
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon 1/P Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);  GenInvPRes[0]->Draw();  GenInvPRes[0]->Fit("gaus","Q");
  pad[page]->cd(2);  GenInvPResVsPt[0]->Draw();
  pad[page]->cd(3);  GenInvPRes[1]->Draw();  GenInvPRes[1]->Fit("gaus","Q");
  pad[page]->cd(4);  GenInvPResVsPt[1]->Draw();
  pad[page]->cd(5);  GenInvPRes[2]->Draw();  GenInvPRes[2]->Fit("gaus","Q");
  pad[page]->cd(6);  GenInvPResVsPt[2]->Draw();
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon 1/P Resolution");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 2; i++) {
      pad[page]->cd(2*j+i+1);  MuInvPRes[i][j]->Draw();
      MuInvPRes[i][j]->Fit("gaus","Q");
    }
  }
  TH1F* MuInvPResGMRTot = (TH1F*)MuInvPRes[0][0]->Clone();
  MuInvPResGMRTot->SetNameTitle("MuInvPResGMRTot", "1/P res, GMR, total");
  MuInvPResGMRTot->Add(MuInvPRes[0][0], MuInvPRes[0][1], 1., 1.);
  MuInvPResGMRTot->Add(MuInvPResGMRTot, MuInvPRes[0][2], 1., 1.);
  pad[page]->cd(7);
  MuInvPResGMRTot->Draw();  MuInvPResGMRTot->Fit("gaus","Q");
  TH1F* MuInvPResOPTTot = (TH1F*)MuInvPRes[1][0]->Clone();
  MuInvPResOPTTot->SetNameTitle("MuInvPResOPTTot", "1/P res, OPT, total");
  MuInvPResOPTTot->Add(MuInvPRes[1][0], MuInvPRes[1][1], 1., 1.);
  MuInvPResOPTTot->Add(MuInvPResOPTTot, MuInvPRes[1][2], 1., 1.);
  pad[page]->cd(8);
  MuInvPResOPTTot->Draw();  MuInvPResOPTTot->Fit("gaus","Q");
  c1->Update();
  delete MuInvPResGMRTot;
  delete MuInvPResOPTTot;

  // Invariant mass
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			 "Invariant Mass Distributions, Dilepton Only");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,5);
  for (int i = 0; i < MAX_LEVELS+1; i++) {
    pad[page]->cd(i+1);
    DilMassComp[i][0]->Draw();
    DilMassComp[i][0]->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass");
    DilMassComp[i][0]->GetYaxis()->SetTitle("Entries");
    DilMassComp[i][0]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			 "Invariant Mass Distributions, Dilepton Only (log scale)");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,5);
  for (int i = 0; i < MAX_LEVELS+1; i++) {
    pad[page]->cd(i+1);
    gPad->SetLogy(1);
    DilMassComp[i][0]->Draw();
    DilMassComp[i][0]->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass");
    DilMassComp[i][0]->GetYaxis()->SetTitle("Entries");
    DilMassComp[i][0]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			 "Invariant Mass Distributions, Including Photons");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,5);
  for (int i = 0; i < MAX_LEVELS+1; i++) {
    pad[page]->cd(i+1);
    DilMassComp[i][1]->Draw();
    DilMassComp[i][1]->GetXaxis()->SetTitle("Resonance mass");
    DilMassComp[i][1]->GetYaxis()->SetTitle("Entries");
    DilMassComp[i][1]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
			 "Invariant Mass Distributions, Including Photons (log scale)");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,5);
  for (int i = 0; i < MAX_LEVELS+1; i++) {
    pad[page]->cd(i+1);
    gPad->SetLogy(1);
    DilMassComp[i][1]->Draw();
    DilMassComp[i][1]->GetXaxis()->SetTitle("Resonance mass");
    DilMassComp[i][1]->GetYaxis()->SetTitle("Entries");
    DilMassComp[i][1]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  // Mass Resolution in each level.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Invariant Mass Resolution, Dilepton Only");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  pad[page]->cd(1);  GenDilMassRes[0]->Draw();
  GenDilMassRes[0]->Fit("gaus","Q");
  pad[page]->cd(2);  GenDilMassFrRes[0]->Draw();
  GenDilMassFrRes[0]->Fit("gaus","Q");
  pad[page]->cd(3);  GenDilMassRes[1]->Draw();
  GenDilMassRes[1]->Fit("gaus","Q");
  pad[page]->cd(4);  GenDilMassFrRes[1]->Draw();
  GenDilMassFrRes[1]->Fit("gaus","Q");
  pad[page]->cd(5);  AllDilMassRes->Draw();
  AllDilMassRes->Fit("gaus","Q");
  pad[page]->cd(7);  GenDilMassRes[2]->Draw();
  GenDilMassRes[2]->Fit("gaus","Q");
  pad[page]->cd(8);  GenDilMassFrRes[2]->Draw();
  GenDilMassFrRes[2]->Fit("gaus","Q");
  c1->Update();   

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Invariant Mass Resolution, Dilepton Only");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  // All this mess below is just to take square root of the content of every
  // bin (leave errors for later).  Clone() does not work since we have to
  // make TH1F from TProfile.
  TH1F *MassResScat[NUM_REC_LEVELS-1], *MassFrResScat[NUM_REC_LEVELS-1];
  for (int i_rec = 0; i_rec < NUM_REC_LEVELS-1; i_rec++) {
    ostringstream histname;
    histname << "Sqrt(Var(L" << i_rec + 1 << "-Gen Mass)) vs Gen Mass";
    nbins = GenDilMassResScat[i_rec]->GetNbinsX();
    MassResScat[i_rec] =
      fs->make<TH1F>(nameHist("MassResScat", i_rec).c_str(),
		     histname.str().c_str(), nbins,
		     GenDilMassResScat[i_rec]->GetXaxis()->GetXmin(),
		     GenDilMassResScat[i_rec]->GetXaxis()->GetXmax());
    for (int ibin = 1; ibin <= nbins; ibin++) {
      f_bin   = GenDilMassResScat[i_rec]->GetBinContent(ibin);
      ent_bin = GenDilMassResScat[i_rec]->GetBinEntries(ibin);
      if (f_bin > 0.) {
	f_bin = sqrt(f_bin);
	if (ent_bin > 0.) {err_bin = f_bin/sqrt(2.*ent_bin);}
	else              {err_bin = 0.;}
      }
      else    {f_bin = 0.; err_bin = 0.;}
      MassResScat[i_rec]->SetBinContent(ibin, f_bin);
      MassResScat[i_rec]->SetBinError(ibin, err_bin);
    }
    pad[page]->cd(2*i_rec+1);  MassResScat[i_rec]->Draw();

    histtitle = "Sqrt(Var((" + recLevelHelper.levelName(i_rec+1) +
      "-Gen Mass)/Gen Mass)) vs Gen Mass";
    nbins = GenDilMassFrResScat[i_rec]->GetNbinsX();
    MassFrResScat[i_rec] =
      fs->make<TH1F>(nameHist("MassFrResScat", i_rec).c_str(),
		     histtitle.c_str(), nbins,
		     GenDilMassFrResScat[i_rec]->GetXaxis()->GetXmin(),
		     GenDilMassFrResScat[i_rec]->GetXaxis()->GetXmax());
    for (int ibin = 1; ibin <= nbins; ibin++) {
      f_bin   = GenDilMassFrResScat[i_rec]->GetBinContent(ibin);
      ent_bin = GenDilMassFrResScat[i_rec]->GetBinEntries(ibin);
      if (f_bin > 0.) {
	f_bin = sqrt(f_bin);
	if (ent_bin > 0.) {err_bin = f_bin/sqrt(2.*ent_bin);}
	else              {err_bin = 0.;}
      }
      else    {f_bin = 0.; err_bin = 0.;}
      MassFrResScat[i_rec]->SetBinContent(ibin, f_bin);
      MassFrResScat[i_rec]->SetBinError(ibin, err_bin);
    }
    pad[page]->cd(2*i_rec+2);  MassFrResScat[i_rec]->Draw();
  }
  c1->Update();
  for (int i_rec = 0; i_rec < NUM_REC_LEVELS-1; i_rec++) {
    delete MassResScat[i_rec];
    delete MassFrResScat[i_rec];
  }

  // Dilepton mass resolution relative to the mass reconstructed from GEANT
  // muons.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Dilepton Mass Resolution, Muons Only");
  title->SetFillColor(10);  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i = l2; i <= MAX_LEVELS; i++) {
    pad[page]->cd(i-1);
    DilMassRes[i][0]->Draw();  DilMassRes[i][0]->Fit("gaus","Q");
    DilMassRes[i][0]->GetXaxis()->SetTitle("(Rec mass - Gen mass)/Gen mass");
    DilMassRes[i][0]->GetYaxis()->SetTitle("Entries");
    DilMassRes[i][0]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  // Dilepton mass resolution with respect to the true resonance mass from
  // PYTHIA.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Dilepton Mass Resolution w.r.t. True Resonance Mass");
  title->SetFillColor(10);  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i = l2; i <= MAX_LEVELS; i++) {
    pad[page]->cd(i-1);
    DilMassRes[i][1]->Draw();  DilMassRes[i][1]->Fit("gaus","Q");
    DilMassRes[i][1]->GetXaxis()->SetTitle("(Rec mass - Gen mass)/Gen mass");
    DilMassRes[i][1]->GetYaxis()->SetTitle("Entries");
    DilMassRes[i][1]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  // Resonance (dilepton+photons) mass resolution with respect to the
  // true resonance mass from PYTHIA.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Resonance Mass Resolution w.r.t. True Resonance Mass");
  title->SetFillColor(10);  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i = l2; i <= MAX_LEVELS; i++) {
    pad[page]->cd(i-1);
    DilMassRes[i][2]->Draw();  DilMassRes[i][2]->Fit("gaus","Q");
    DilMassRes[i][2]->GetXaxis()->SetTitle("(Rec mass - Gen mass)/Gen mass");
    DilMassRes[i][2]->GetYaxis()->SetTitle("Entries");
    DilMassRes[i][2]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  // Dilepton mass resolution at the Z' mass, muons only.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Dilepton Mass Resolution at Z' Mass, Muons Only");
  title->SetFillColor(10);  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i = l2; i <= MAX_LEVELS; i++) {
    pad[page]->cd(i-1);
    DilMassRes[i][3]->Draw();  DilMassRes[i][3]->Fit("gaus","Q");
    DilMassRes[i][3]->GetXaxis()->SetTitle("(Rec mass - Gen mass)/Gen mass");
    DilMassRes[i][3]->GetYaxis()->SetTitle("Entries");
    DilMassRes[i][3]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  // Dilepton mass resolution at the Z' mass, w.r.t. the true resonance mass.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Dilepton Mass Resolution at Z' Mass, w.r.t. True Mass");
  title->SetFillColor(10);  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i = l2; i <= MAX_LEVELS; i++) {
    pad[page]->cd(i-1);
    DilMassRes[i][4]->Draw();  DilMassRes[i][4]->Fit("gaus","Q");
    DilMassRes[i][4]->GetXaxis()->SetTitle("(Rec mass - Gen mass)/Gen mass");
    DilMassRes[i][4]->GetYaxis()->SetTitle("Entries");
    DilMassRes[i][4]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  // Resonance (dilepton+photons) mass resolution at Z' mass with respect to 
  // the true resonance mass from PYTHIA.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1, 0.94, 0.9, 0.98,
				       "Resonance Mass Resolution at Z' Mass, w.r.t. True Mass");
  title->SetFillColor(10);  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i = l2; i <= MAX_LEVELS; i++) {
    pad[page]->cd(i-1);
    DilMassRes[i][5]->Draw();  DilMassRes[i][5]->Fit("gaus","Q");
    DilMassRes[i][5]->GetXaxis()->SetTitle("(Rec mass - Gen mass)/Gen mass");
    DilMassRes[i][5]->GetYaxis()->SetTitle("Entries");
    DilMassRes[i][5]->SetTitleOffset(1.2, "Y");
  }
  c1->Update();

  // Charge misassignment in each level.
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,"Muon Charge Assignment");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  gStyle->SetOptStat(1111);
  gStyle->SetOptLogy(1);
  pad[page]->Draw();
  pad[page]->Divide(2,4);
  for (int i_rec = l1; i_rec <= MAX_LEVELS; i_rec++) {
    QRes[i_rec]->SetMinimum(1.);  QRes[i_rec]->SetMaximum(3000.);
    pad[page]->cd(i_rec);  QRes[i_rec]->Draw("hist");
  }
  /* for (int ibin = 1; ibin <= QRes[3]->GetNbinsX(); ibin++) {
     LogTrace("Zprime2muResolution") << "ibin = "     << ibin
     << " entries = " << QRes[3]->GetBinContent(ibin);
     } */
  gStyle->SetOptStat(111111);
  gStyle->SetOptLogy(0);
  c1->Update();

  // Charge misassignment vs pT and 1/pT
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,
				       "Charge Misassignment vs pT and 1/pT");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);
  ratio[0][0] = (TH1F*)QResVsPt[0][0]->Clone();
  ratio[0][0]->SetNameTitle("QResVsPtRatio00", "(wrong Q)/(right Q) vs L1 pT");
  ratio[0][0]->Divide(QResVsPt[0][1], QResVsPt[0][0], 1., 1.);
  ratio[0][0]->Draw();
  pad[page]->cd(2);
  ratio[0][1] = (TH1F*)QResVsInvPt[0][0]->Clone();
  ratio[0][1]->SetNameTitle("QResVsPtRatio01", "(wrong Q)/(right Q) vs L1 1/pT");
  ratio[0][1]->Divide(QResVsInvPt[0][1], QResVsInvPt[0][0], 1., 1.);
  ratio[0][1]->Draw();
  pad[page]->cd(3);
  ratio[0][2] = (TH1F*)QResVsPt[1][0]->Clone();
  ratio[0][2]->SetNameTitle("QResVsPtRatio02", "(wrong Q)/(right Q) vs L2 pT");
  ratio[0][2]->Divide(QResVsPt[1][1], QResVsPt[1][0], 1., 1.); 
  ratio[0][2]->Draw();
  pad[page]->cd(4);
  ratio[0][3] = (TH1F*)QResVsInvPt[1][0]->Clone();
  ratio[0][3]->SetNameTitle("QResVsPtRatio03", "(wrong Q)/(right Q) vs L2 1/pT");
  ratio[0][3]->Divide(QResVsInvPt[1][1], QResVsInvPt[1][0], 1., 1.);
  ratio[0][3]->Draw();
  pad[page]->cd(5);
  ratio[0][4] = (TH1F*)QResVsPt[2][0]->Clone();
  ratio[0][4]->SetNameTitle("QResVsPtRatio04", "(wrong Q)/(right Q) vs L3 pT");
  ratio[0][4]->Divide(QResVsPt[2][1], QResVsPt[2][0], 1., 1.); 
  ratio[0][4]->Draw();
  pad[page]->cd(6);
  ratio[0][5] = (TH1F*)QResVsInvPt[2][0]->Clone();
  ratio[0][5]->SetNameTitle("QResVsPtRatio05", "(wrong Q)/(right Q) vs L3 1/pT");
  ratio[0][5]->Divide(QResVsInvPt[2][1], QResVsInvPt[2][0], 1., 1.);
  ratio[0][5]->Draw();
  c1->Update();

  // Charge misassignment vs P and 1/P
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,
				       "Charge Misassignment vs P and 1/P");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  pad[page]->Draw();
  pad[page]->Divide(2,3);
  pad[page]->cd(1);
  ratio[1][0] = (TH1F*)QResVsP[0][0]->Clone();
  ratio[1][0]->SetNameTitle("QResVsPtRatio10", "(wrong Q)/(right Q) vs L1 P");
  ratio[1][0]->Divide(QResVsP[0][1], QResVsP[0][0], 1., 1.); 
  ratio[1][0]->Draw();
  pad[page]->cd(2);
  ratio[1][1] = (TH1F*)QResVsInvP[0][0]->Clone();
  ratio[1][1]->SetNameTitle("QResVsPtRatio11", "(wrong Q)/(right Q) vs L1 1/P");
  ratio[1][1]->Divide(QResVsInvP[0][1], QResVsInvP[0][0], 1., 1.);
  ratio[1][1]->Draw();
  pad[page]->cd(3);
  ratio[1][2] = (TH1F*)QResVsP[1][0]->Clone();
  ratio[1][2]->SetNameTitle("QResVsPtRatio12", "(wrong Q)/(right Q) vs L2 P");
  ratio[1][2]->Divide(QResVsP[1][1], QResVsP[1][0], 1., 1.); 
  ratio[1][2]->Draw();
  pad[page]->cd(4);
  ratio[1][3] = (TH1F*)QResVsInvP[1][0]->Clone();
  ratio[1][3]->SetNameTitle("QResVsPtRatio13", "(wrong Q)/(right Q) vs L2 1/P");
  ratio[1][3]->Divide(QResVsInvP[1][1], QResVsInvP[1][0], 1., 1.);
  ratio[1][3]->Draw();
  pad[page]->cd(5);
  ratio[1][4] = (TH1F*)QResVsP[2][0]->Clone();
  ratio[1][4]->SetNameTitle("QResVsPtRatio14", "(wrong Q)/(right Q) vs L3 P");
  ratio[1][4]->Divide(QResVsP[2][1], QResVsP[2][0], 1., 1.);
  ratio[1][4]->SetMinimum(-0.1); 
  ratio[1][4]->SetMaximum(0.5); 
  ratio[1][4]->Draw();
  pad[page]->cd(6);
  ratio[1][5] = (TH1F*)QResVsInvP[2][0]->Clone();
  ratio[1][5]->SetNameTitle("QResVsPtRatio15", "(wrong Q)/(right Q) vs L3 1/P");
  ratio[1][5]->Divide(QResVsInvP[2][1], QResVsInvP[2][0], 1., 1.);
  ratio[1][5]->SetMinimum(-0.1); 
  ratio[1][5]->SetMaximum(0.5); 
  ratio[1][5]->Draw();
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,
				       "#Sigma Pt (#DeltaR < 0.3)");
  title->SetFillColor(10);
  title->Draw();
  strpage << "- " << (++page) << " -";
  t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
  gStyle->SetOptLogy(1);
  pad[page]->Draw();
  pad[page]->Divide(2,5);
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 2; j++) {
      pad[page]->cd(i*2+j+1);
      SumPtR03[i+lgmr][j]->Draw();
    }
  gStyle->SetOptLogy(0);
  c1->Update();

  if (doingHiggs) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = "H(130) --> ZZ* -> 4 mu";
    delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(2,2);
    pad[page]->cd(1);  ZonDilMass[0]->Draw();
    pad[page]->cd(2);  ZofDilMass[0]->Draw();
    pad[page]->cd(3);  ZonDilMass[3]->Draw();
    pad[page]->cd(4);  ZofDilMass[3]->Draw();
    c1->Update();

    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = "H(130) --> ZZ* -> 4 mu";
    delete title; title = new TPaveLabel(0.1,0.94,0.9,0.98,tit.c_str());
    title->SetFillColor(10);
    title->Draw();
    strpage << "- " << (++page) << " -";
    t.DrawText(.9, .02, strpage.str().c_str());  strpage.str("");
    pad[page]->Draw();
    pad[page]->Divide(1,2);
    pad[page]->cd(1);  Eta4muons->Draw();
    pad[page]->cd(2);  Pt4muons->Draw();
    c1->Update();
  }

  // Close postscript file.
  ps->Close();

  delete title;
  delete ps;
  delete c1;

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      delete ratio[i][j];
}

void Zprime2muResolution::DrawAcceptance() {
  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 320);
  TPostScript *eps = new TPostScript("accept_vs_mass.eps", 113);
  eps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  TH1F* accept;
  accept = (TH1F*)GenMassAllEvents->Clone();
  accept->SetTitle("");
  accept->Divide(GenMassInAccept, GenMassAllEvents, 1., 1.);
  for (int ibin = 1; ibin <= accept->GetNbinsX(); ibin++) {
    edm::LogInfo("Zprime2muResolution")
      << "ibin = " << ibin
      << " mass = " << accept->GetXaxis()->GetBinCenter(ibin)
      << " acceptance = " << accept->GetBinContent(ibin);
  }
  accept->SetMinimum(0.0);
  accept->SetMaximum(1.05);
  accept->GetXaxis()->SetTitle("M#mu^{+}#mu^{-} (GeV)");
  accept->GetYaxis()->SetTitle("Acceptance");
  accept->Draw("histo");
  c1->Update();
  eps->Close();
  delete accept;
  delete eps;
  delete c1;
}
