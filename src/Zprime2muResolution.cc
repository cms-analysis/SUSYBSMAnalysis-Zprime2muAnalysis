//
// Authors: Jason Mumford, Jordan Tucker, Slava Valuev, UCLA
//

#include <string>
#include <vector>

#include "TCanvas.h"
#include "TPad.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TText.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muResolution.h"

using namespace std;

Zprime2muResolution::Zprime2muResolution(const edm::ParameterSet& config) 
  : Zprime2muAnalysis(config)
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

  BookResHistos();
  BookEffHistos();
  BookPtResHistos();
  BookDilResHistos();
  BookChargeResHistos();
}

Zprime2muResolution::~Zprime2muResolution() {
  DeleteHistos();
}

void Zprime2muResolution::BookResHistos() {
  // Origin of muons
  Origin[0] = new TH1F("","Particle Id of Mother of all mu's", 20, 0, 20);
  Origin[1] = new TH1F("","Particle Id of Mother of opp-sign dilepton mu's",
		       20, 0, 20);

  // Trigger decisions
  TrigResult[0][0] = new TH1F("", "Gen, all events",        5, 0., 5.);
  TrigResult[1][0] = new TH1F("", "L1 trigger, all events", 5, 0., 5.);
  TrigResult[2][0] = new TH1F("", "L2 trigger, all events", 5, 0., 5.);
  TrigResult[3][0] = new TH1F("", "L3 trigger, all events", 5, 0., 5.);
  TrigResult[0][1] = new TH1F("", "Gen, #eta < 2.4",        5, 0., 5.);
  TrigResult[1][1] = new TH1F("","L1 trigger, #eta < 2.4",  5, 0., 5.);
  TrigResult[2][1] = new TH1F("","L2 trigger, #eta < 2.4",  5, 0., 5.);
  TrigResult[3][1] = new TH1F("","L3 trigger, #eta < 2.4",  5, 0., 5.);
  TrigResult[0][2] = new TH1F("", "Gen, #eta < 2.1",        5, 0., 5.);
  TrigResult[1][2] = new TH1F("","L1 trigger, #eta < 2.1",  5, 0., 5.);
  TrigResult[2][2] = new TH1F("","L2 trigger, #eta < 2.1",  5, 0., 5.);
  TrigResult[3][2] = new TH1F("","L3 trigger, #eta < 2.1",  5, 0., 5.);

  TrigMass[0][0] = new TH1F("", "Gen mass, Gen, all events", 27, 0., 5.4);
  TrigMass[1][0] = new TH1F("", "Gen mass, L1, all events",  27, 0., 5.4);
  TrigMass[2][0] = new TH1F("", "Gen mass, L2, all events",  27, 0., 5.4);
  TrigMass[3][0] = new TH1F("", "Gen mass, L3, all events",  27, 0., 5.4);
  TrigMass[0][1] = new TH1F("", "Gen mass, Gen, #eta < 2.4", 27, 0., 5.4);
  TrigMass[1][1] = new TH1F("", "Gen mass, L1, #eta < 2.4",  27, 0., 5.4);
  TrigMass[2][1] = new TH1F("", "Gen mass, L2, #eta < 2.4",  27, 0., 5.4);
  TrigMass[3][1] = new TH1F("", "Gen mass, L3, #eta < 2.4",  27, 0., 5.4);
  TrigMass[0][2] = new TH1F("", "Gen mass, Gen, #eta < 2.1", 27, 0., 5.4);
  TrigMass[1][2] = new TH1F("", "Gen mass, L1, #eta < 2.1",  27, 0., 5.4);
  TrigMass[2][2] = new TH1F("", "Gen mass, L2, #eta < 2.1",  27, 0., 5.4);
  TrigMass[3][2] = new TH1F("", "Gen mass, L3, #eta < 2.1",  27, 0., 5.4);

  // Define errors
  for (int i_rec = 0; i_rec < NUM_REC_LEVELS; i_rec++) {
    for (int j = 0; j < 3; j++) {
      TrigMass[i_rec][j]->Sumw2();
    }
  }

  EventsInAccFailed  = new TH1F("", "Events in acc., failed", 5, 0., 5.);
  L1TrigFailSingleMu = new TH1F("", "pT, failed 1mu L1", 25, 0., 25.);
  L1TrigFailMu2VsMu1 = new TH2F("", "pT mu2 vs pT mu1, failed di-mu L1",
				25, 0., 25., 25, 0., 25.);
  L1TrigPassSingleMu = new TH1F("", "pT, passed 1mu L1", 25, 0., 25.);
  L1TrigPassMu2VsMu1 = new TH2F("", "pT mu2 vs pT mu1, passed di-mu L1",
				25, 0., 25., 25, 0., 25.);

  L2MuonHits    = new TH1F("", "Muon hits, Level-2", 25, 0., 25.);
  L3TrackerHits = new TH1F("", "Tracker hits (pixel+silicon), Level-3",
			   25, 0., 25.);
  GMRMuonHits[0]= new TH1F("", "Muon hits, GMR, barrel",  25, 0., 25.);
  GMRMuonHits[1]= new TH1F("", "Muon hits, GMR, overlap", 25, 0., 25.);
  GMRMuonHits[2]= new TH1F("", "Muon hits, GMR, endcap",  25, 0., 25.);
  GMRChi2dof[0] = new TH1F("", "Chi2/d.o.f., GMR, barrel",  100, 0., 5.);
  GMRChi2dof[1] = new TH1F("", "Chi2/d.o.f., GMR, overlap", 100, 0., 5.);
  GMRChi2dof[2] = new TH1F("", "Chi2/d.o.f., GMR, endcap",  100, 0., 5.);

  // Muons per event
  NMuons[0][0] = new TH1F("", "Muons/Event, Gen", 10, 0, 10);
  NMuons[1][0] = new TH1F("", "Muons/Event, L1",  10, 0, 10);
  NMuons[2][0] = new TH1F("", "Muons/Event, L2",  10, 0, 10);
  NMuons[3][0] = new TH1F("", "Muons/Event, L3",  10, 0, 10);
  NMuons[0][1] = new TH1F("", "Muons/Event (w/ opp-sign dilepton), Gen",
			  10, 0, 10);
  NMuons[1][1] = new TH1F("", "Muons/Event (w/ opp-sign dilepton), L1",
			  10, 0, 10);
  NMuons[2][1] = new TH1F("", "Muons/Event (w/ opp-sign dilepton), L2",
			  10, 0, 10);
  NMuons[3][1] = new TH1F("", "Muons/Event (w/ opp-sign dilepton), L3",
			  10, 0, 10);
  NMuons[0][2] = new TH1F("", "Muons/Event, Gen, triggered", 10, 0., 10.);
  NMuons[1][2] = new TH1F("", "Muons/Event, L1, triggered",  10, 0., 10.);
  NMuons[2][2] = new TH1F("", "Muons/Event, L2, triggered",  10, 0., 10.);
  NMuons[3][2] = new TH1F("", "Muons/Event, L3, triggered",  10, 0., 10.);
  NMuons[0][3] = new TH1F("", "Muons/Event in acc., Gen, failed", 10, 0., 10.);
  NMuons[1][3] = new TH1F("", "Muons/Event in acc., L1, failed",  10, 0., 10.);
  NMuons[2][3] = new TH1F("", "Muons/Event in acc., L2, failed",  10, 0., 10.);
  NMuons[3][3] = new TH1F("", "Muons/Event in acc., L3, failed",  10, 0., 10.);
  NumDilVsRec  = new TH1F("", "Events w/ opp-sign dilepton vs Rec Level",
			   5, 0,  5);
  SignOfDil[0] = new TH1F("", "Number of opp-sign dileptons, Gen",
			  2, 0,  2);
  SignOfDil[1] = new TH1F("", "Number of opp-sign, like-sign dileptons, L1",
			  3, 0,  3);
  SignOfDil[2] = new TH1F("", "Number of opp-sign, like-sign dileptons, L2",
			  3, 0,  3);
  SignOfDil[3] = new TH1F("", "Number of opp-sign, like-sign dileptons, L3",
			  3, 0,  3);

  // Main kinematic variables for all muons
  for (int i = 0; i < NUM_REC_LEVELS; i++) {
    MuonEta[i] = new TH1F("", "Eta", 100, -5.,  5. );
    MuonRap[i] = new TH1F("", "Y",   100, -5.,  5. );
    MuonPhi[i] = new TH1F("", "Phi", 100,  0.,  6.3);
    if (i==1) {
      MuonPt[i] = new TH1F("", "Pt", 100,  0.,  150.);
      MuonPz[i] = new TH1F("", "Pz", 100,  0.,  800.);
      MuonP[i]  = new TH1F("", "P",  100, -5.,  800.);
      MuonPVsEta[i]  = new TH2F("", "P vs Eta",  100, -6., 6., 100, 0., 1000);
      MuonPtVsEta[i] = new TH2F("", "Pt vs Eta", 100, -6., 6., 100, 0., 150);
    }
    else {
      MuonPt[i] = new TH1F("", "Pt", 100, 0.,    peakMass);
      MuonPz[i] = new TH1F("", "Pz", 100, 0., 2.*peakMass);
      MuonP[i]  = new TH1F("", "P",  100, 0., 2.*peakMass);
      MuonPVsEta[i]  = new TH2F("", "P vs Eta",  100, -6., 6., 100, 0., 
				2.*peakMass);
      MuonPtVsEta[i] = new TH2F("", "Pt vs Eta", 100, -6., 6., 100, 0., 
				   peakMass);
    }

    // Main kinematic variables for opp-sign dileptons and
    // muons associated with opp-sign dileptons
    if (i == 1) {
      // AllDilMass[i] = new TH1F("", "Opp-sign Dilepton Mass before Q cuts",
      //			       100, 0, 1200);
      DilMass[i]    = new TH1F("", "Opp-sign Dilepton Mass", 100, 0, 1200);
      DilMassVsEta[i] = new TH2F("", "Dilepton Mass vs Eta", 50, -10., 10.,
				 100, 0, 1200);
      DilMassVsY[i] = new TH2F("", "Dilepton Mass vs Y", 50, -5., 5., 100, 0,
			       1200);
      MuMVsMuP[i][2] = new TH2F("","Pt  of mu- vs mu+", 100,0, 150,100, 0,150);
    }
    else {
      // AllDilMass[i] = new TH1F("", "Opp-sign Dilepton Mass before Q cuts",
      //		       binSize, lowerMassWin,
      //		       upperMassWin);
      DilMass[i]    = new TH1F("", "Opp-sign Dilepton Mass",
			       binSize, lowerMassWin,
			       upperMassWin);
      DilMassVsEta[i] = new TH2F("", "Dilepton Mass vs Eta", 50, -10., 10., 
				 50, 0, upperMassWin);
      DilMassVsY[i] = new TH2F("", "Dilepton Mass vs Y", 50, -5., 5.,   50, 0, 
			       upperMassWin);
      MuMVsMuP[i][2] = new TH2F("","Pt  of mu- vs mu+",
				100, 0., peakMass,
				100, 0., peakMass);
    }
    string tit = "Opp-sign Dilepton Mass, " + str_level[i];
    AllDilMass[i] =
      new TH1F("", tit.c_str(), 50, 0., upperMassWin);

    MuMVsMuP[i][0] = new TH2F("","Eta of mu- vs mu+", 100, -3, 3, 100, -3, 3);
    MuMVsMuP[i][1] = new TH2F("","Phi of mu- vs mu+", 100, 0, 6.28,100,0,6.28);
    for (int j = 0; j < 3; j++) {
      Phi[i][j]      = new TH1F("", "Phi",       100, 0., 6.3);
      PVsEta[i][j]   = new TProfile("", "P vs Eta",  50, -6., 6.);
      PtVsEta[i][j]  = new TProfile("", "Pt vs Eta", 50, -6., 6.);
      Eta[i][j]      = new TH1F("", "Eta",          100, -5., 5.);
      Rapidity[i][j] = new TH1F("", "y",            100, -5., 5.);
      if (j != 2) {
	if (i == 1) { // mu+ and mu- at L1
	  Pt[i][j] = new TH1F("", "Pt", 100, 0., 200.);
	  P[i][j]  = new TH1F("", "P",  100, 0., 800.);
	  Pz[i][j] = new TH1F("", "Pz", 100, 0., 800.);
	}
	else {
	  Pt[i][j] = new TH1F("", "Pt", 100, 0.,    peakMass);
	  P[i][j]  = new TH1F("", "P",  100, 0., 2.*peakMass);
	  Pz[i][j] = new TH1F("", "Pz", 100, 0., 2.*peakMass);
	}
      }
      else { // Special range for dilepton values.
	if (i==1) {
	  Pt[i][j] = new TH1F("", "Pt", 100, 0., 200.);
	  P[i][j]  = new TH1F("", "P",  100, 0., 800.);
	  Pz[i][j] = new TH1F("", "Pz", 100, 0., 800.);	
	}
	else {
	  Pt[i][j] = new TH1F("", "Pt", 100, 0., .5*peakMass);
	  P[i][j]  = new TH1F("", "P",
			      100, 0, 1000.+upperMassWin);
	  Pz[i][j] = new TH1F("", "Pz",
			      100, 0, 1000.+upperMassWin);
	}
      }
    }
    ZonDilMass[i] = new TH1F("", "Highest Z mass", 100, 0., 120.);
    ZofDilMass[i] = new TH1F("", "Second  Z mass", 100, 0., 120.);
  }
  Eta4muons = new TH1F("", "Eta", 100, -5.,  5.);
  Pt4muons  = new TH1F("", "pT",  100,  0.,100.);

  // Differences between different levels of reconstruction
  for (int i = l1; i <= MAX_LEVELS; i++) {
    string tit_eta = str_level[i] + " eta - Gen eta";
    string tit_phi = str_level[i] + " phi - Gen phi";
    string tit_pt  = str_level[i] + " pT - Gen pT";
    string tit_p   = "(" + str_level[i] + " P - Gen P)/(Gen P)";
    string tit_ppr = "(" + str_level[i] + " P - Gen P)/(Gen P) vs Gen P";
    if (i == l1) {
      EtaRes[i] = new TH1F("", tit_eta.c_str(), 100, -0.1,  0.1);
      PhiRes[i] = new TH1F("", tit_phi.c_str(), 100, -0.1,  0.1);
      PtDiff[i] = new TH1F("", tit_pt.c_str(),  100, 
			   -700.-.2*peakMass,
			   700.+.2*peakMass);
      PRes[i]   = new TH1F("", tit_p.c_str(),   100, -1., 1.);
      PResVsP[i]= new TProfile("", tit_ppr.c_str(), 50,
			       0., upperMassWin, -1., 1.);
    }
    else if (i == l2) {
      EtaRes[i] = new TH1F("", tit_eta.c_str(), 100, -0.01, 0.01);
      PhiRes[i] = new TH1F("", tit_phi.c_str(), 100, -0.01, 0.01);
      PtDiff[i] = new TH1F("", tit_pt.c_str(),  100,
			-0.4*peakMass, 0.4*peakMass);
      PRes[i]   = new TH1F("", tit_p.c_str(),   100, -1., 1.);
      PResVsP[i]= new TProfile("", tit_ppr.c_str(), 50,
			       0., upperMassWin, -1., 1.);
    }
    else {
      EtaRes[i] = new TH1F("", tit_eta.c_str(), 100, -0.001,  0.001);
      PhiRes[i] = new TH1F("", tit_phi.c_str(), 100, -0.0005, 0.0005);
      PtDiff[i] = new TH1F("", tit_pt.c_str(),  100, 
		        -0.1*peakMass, 0.1*peakMass);
      PRes[i]   = new TH1F("", tit_p.c_str(),   100, -0.3, 0.3);
      PResVsP[i]= new TProfile("", tit_ppr.c_str(), 50,
			       0., upperMassWin, -0.3, 0.3);
    }
  }

  GenPhiResVsPhi[0] = new TProfile("", "|L1 phi - Gen phi| vs Gen phi",
				   25, 0., 6.3, -0.1,  0.1);
  GenPhiResVsPhi[1] = new TProfile("", "|L2 phi - Gen phi| vs Gen phi",
				   25, 0., 6.3, -0.1,  0.1);
  GenPhiResVsPhi[2] = new TProfile("", "|L3 phi - Gen phi| vs Gen phi",
				   25, 0., 6.3, -0.01, 0.01);

  GenInvPtRes[0] = new TH1F("", "(L1 1/Pt - Gen 1/Pt)/(Gen 1/Pt)", 100,
			    -.6*peakMass/140.,
			     .6*peakMass/140.);
  GenInvPtRes[1] = new TH1F("", "(L2 1/Pt - Gen 1/Pt)/(Gen 1/Pt)", 100, -2, 2);
  GenInvPtRes[2] = new TH1F("", "(L3 1/Pt - Gen 1/Pt)/(Gen 1/Pt)",
			    100, -0.3, 0.3);
  GenInvPtResVsPt[0] = new TProfile("",
                        "(L1 1/pT - Gen 1/pT)/(Gen 1/pT) vs Gen pT",
				    50, 0., peakMass,
				    -.6*peakMass/140.,
				     .6*peakMass/140.);
  GenInvPtResVsPt[1] = new TProfile("",
                        "(L2 1/pT - Gen 1/pT)/(Gen 1/pT) vs Gen pT",
				    50, 0., peakMass, -2., 2.);
  GenInvPtResVsPt[2] = new TProfile("",
                        "(L3 1/pT - Gen 1/pT)/(Gen 1/pT) vs Gen pT",
				    50, 0., peakMass, -0.3, 0.3);

  GenInvPRes[0] = new TH1F("", "(L1 1/P - Gen 1/P)/(Gen 1/P)", 100,
			   -upperMassWin/200.,
			    upperMassWin/200.);
  GenInvPRes[1] = new TH1F("", "(L2 1/P - Gen 1/P)/(Gen 1/P)", 100, -2.,  2.);
  GenInvPRes[2] = new TH1F("", "(L3 1/P - Gen 1/P)/(Gen 1/P)", 100, -0.3, 0.3);
  GenPResVsPt[0] = new TProfile("", "(L1 P - Gen P)/(Gen P) vs Gen pT",
				50, 0., .6*peakMass, -1.,  1.);
  GenPResVsPt[1] = new TProfile("", "(L2 P - Gen P)/(Gen P) vs Gen pT",
				50, 0., .6*peakMass, -1.,  1.);
  GenPResVsPt[2] = new TProfile("", "(L3 P - Gen P)/(Gen P) vs Gen pT",
				50, 0., .7*peakMass, -0.3, 0.3);
  GenInvPResVsPt[0] = new TProfile("","(L1 1/P - Gen 1/P)/(Gen 1/P) vs Gen pT",
				   50, 0., peakMass,
				-upperMassWin/200.,
				 upperMassWin/200.);
  GenInvPResVsPt[1] = new TProfile("","(L2 1/P - Gen 1/P)/(Gen 1/P) vs Gen pT",
				   50, 0., peakMass, -2.,  2.);
  GenInvPResVsPt[2] = new TProfile("","(L3 1/P - Gen 1/P)/(Gen 1/P) vs Gen pT",
				   50, 0., peakMass, -0.3, 0.3);

  GenEtaResScat[0] = new TH2F("", "L1 Eta vs Gen Eta", 100, -3, 3, 100, -3, 3);
  GenEtaResScat[1] = new TH2F("", "L2 Eta vs Gen Eta", 100, -3, 3, 100, -3, 3);
  GenEtaResScat[2] = new TH2F("", "L3 Eta vs Gen Eta", 100, -3, 3, 100, -3, 3);
  GenPhiResScat[0]=new TH2F("", "L1 Phi vs Gen Phi", 100, 0, 6.3, 100, 0, 6.3);
  GenPhiResScat[1]=new TH2F("", "L2 Phi vs Gen Phi", 100, 0, 6.3, 100, 0, 6.3);
  GenPhiResScat[2]=new TH2F("", "L3 Phi vs Gen Phi", 100, 0, 6.3, 100, 0, 6.3);
  GenPtResScat[0]  = new TH2F("", "L1 Pt vs Gen Pt",
			      100, 0., peakMass, 100, 0., 200.);
  GenPtResScat[1]  = new TH2F("", "L2 Pt vs Gen Pt",
			      100, 0., peakMass,
			      100, 0., peakMass);
  GenPtResScat[2]  = new TH2F("", "L3 Pt vs Gen Pt",
			      100, 0., peakMass,
			      100, 0., peakMass);
  AllDilMassRes[0] = AllDilMassRes[1] = 0;
  AllDilMassRes[2] = new TH1F("", "L3 Mass - Gen Mass, before Q cuts", 100, 
			      -.2*peakMass,
			       .2*peakMass);
  GenDilMassRes[0] = new TH1F("", "L1 Mass - Gen Mass", 100,
			      -peakMass, peakMass);
  GenDilMassRes[1] = new TH1F("", "L2 Mass - Gen Mass", 100,
			      -peakMass, peakMass);
  GenDilMassRes[2] = new TH1F("", "L3 Mass - Gen Mass", 100,
			      -.2*peakMass,
			       .2*peakMass);
  GenDilMassFrRes[0] = new TH1F("","(L1 Mass-Gen Mass)/(Gen Mass)",100,-1.,1.);
  GenDilMassFrRes[1] = new TH1F("","(L2 Mass-Gen Mass)/(Gen Mass)",100,-1.,1.);
  GenDilMassFrRes[2] = new TH1F("","(L3 Mass-Gen Mass)/(Gen Mass)",
				100, -0.2, 0.2);
  GenDilMassResScat[0] = new TProfile("", "L1-Gen Mass vs Gen Mass", 25, 0.,
				      upperMassWin, 0., 
			     peakMass*peakMass, " ");
  GenDilMassResScat[1] = new TProfile("", "L2-Gen Mass vs Gen Mass", 25, 0.,
			              upperMassWin, 0., 
			     peakMass*peakMass, " ");
  GenDilMassResScat[2] = new TProfile("", "L3-Gen Mass vs Gen Mass", 25, 0.,
				      upperMassWin, 0., 
			0.04*peakMass*peakMass, " ");
  GenDilMassFrResScat[0] = new TProfile("",
                            "(L1-Gen Mass)/(Gen Mass) vs Gen Mass",
				   25, 0., upperMassWin, 
					-1., 1., " ");
  GenDilMassFrResScat[1] = new TProfile("",
                            "(L2-Gen Mass)/(Gen Mass) vs Gen Mass",
				   25, 0., upperMassWin, 
					-1., 1., " ");
  GenDilMassFrResScat[2] = new TProfile("",
                            "(L3-Gen Mass)/(Gen Mass) vs Gen Mass",
				   25, 0., upperMassWin, 
					-0.2, 0.2, " ");

  L1EtaRes[0] = new TH1F("", "L2 Eta - L1 Eta", 100, -0.1, 0.1);
  L1EtaRes[1] = new TH1F("", "L3 Eta - L1 Eta", 100, -0.1, 0.1);
  L1PhiRes[0] = new TH1F("", "L2 Phi - L1 Phi", 100, -0.1, 0.1);
  L1PhiRes[1] = new TH1F("", "L3 Phi - L1 Phi", 100, -0.1, 0.1);
  L1PtDiff[0] = new TH1F("", "L2 Pt - L1 Pt",   100, -200.,
			 -200.+peakMass);
  L1PtDiff[1] = new TH1F("", "L3 Pt - L1 Pt",   100, -150., 
			 -150.+peakMass);

  L1EtaResScat[0] = new TH2F("", "L2 Eta vs L1 Eta", 100, -3, 3, 100, -3, 3);
  L1EtaResScat[1] = new TH2F("", "L3 Eta vs L1 Eta", 100, -3, 3, 100, -3, 3);
  L1PhiResScat[0] = new TH2F("", "L2 Phi vs L1 Phi", 100, 0, 6.3, 100, 0, 6.3);
  L1PhiResScat[1] = new TH2F("", "L3 Phi vs L1 Phi", 100, 0, 6.3, 100, 0, 6.3);
  L1PtResScat[0]  = new TH2F("", "L2 Pt vs L1 Pt", 100, 0., 200.,
			     100, 0., peakMass);
  L1PtResScat[1]  = new TH2F("", "L3 Pt vs L1 Pt", 100, 0., 200.,
			     100, 0., peakMass);

  L2EtaRes = new TH1F("", "L3 Eta - L2 Eta", 100, -0.025, 0.025);
  L2PhiRes = new TH1F("", "L3 Phi - L2 Phi", 100, -0.05,  0.05);
  L2PtDiff = new TH1F("", "L3 Pt - L2 Pt",   100,
		      -.6*peakMass, .6*peakMass);
  L2EtaResScat = new TH2F("", "L3 Eta vs L2 Eta", 100, -3,  3,  100, -3, 3);
  L2PhiResScat = new TH2F("", "L3 Phi vs L2 Phi", 100,  0, 6.3,  100, 0, 6.3);
  L2PtResScat  = new TH2F("", "L3 Pt vs L2 Pt", 100, 0., peakMass,
			  100, 0., peakMass);
}

void Zprime2muResolution::BookEffHistos() {
  string tit;
  for (int i = lgen; i <= MAX_LEVELS; i++) {
    tit = "Gen eta, " + str_level[i] + " muons";
    EffVsEta[i] = new TH1F("", tit.c_str(), 50, -2.5, 2.5);
    tit = "Gen phi, " + str_level[i] + " muons";
    EffVsPhi[i] = new TH1F("", tit.c_str(), 63,  0.,  6.3);
    tit = "Gen pT, "  + str_level[i] + " muons";
    EffVsPt[i]  = new TH1F("", tit.c_str(), 35,  0., 3500.);
  }

  RecMass[0][0] = new TH1F("", "Gen mass, gen, all",  27, 0., 5.4);
  RecMass[0][1] = new TH1F("", "Gen mass, gen, in acceptance", 27, 0., 5.4);
  RecMass[1][0] = new TH1F("", "Gen mass, L3, 2mu",   27, 0., 5.4);
  RecMass[1][1] = new TH1F("", "Gen mass, L3, dimu",  27, 0., 5.4);
  RecMass[2][0] = new TH1F("", "Gen mass, GMR, 2mu",  27, 0., 5.4);
  RecMass[2][1] = new TH1F("", "Gen mass, GMR, dimu", 27, 0., 5.4);
  RecMass[3][0] = new TH1F("", "Gen mass, TMR, 2mu",  27, 0., 5.4);
  RecMass[3][1] = new TH1F("", "Gen mass, TMR, dimu", 27, 0., 5.4);

  // Define errors
  for (int i = lgen; i <= MAX_LEVELS; i++) {
    EffVsEta[i]->Sumw2();
    EffVsPhi[i]->Sumw2();
    EffVsPt[i]->Sumw2();
  }
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      RecMass[i][j]->Sumw2();
    }
  }
}

void Zprime2muResolution::BookPtResHistos() {
  ResidualPt[0] = new TProfile("", "(delta 1/pT)/(1/pT) vs Eta, pT < 10",
			       12, 0., 2.4, -10., 10.); // very long tails
  ResidualPt[1] = new TProfile("", "(delta 1/pT)/(1/pT) vs Eta, pT < 25",
			       12, 0., 2.4, -10., 10.);
  ResidualPt[2] = new TProfile("", "(delta 1/pT)/(1/pT) vs Eta, pT < 50",
			       12, 0., 2.4, -10., 10.);
  ResidualPt[3] = new TProfile("", "(delta 1/pT)/(1/pT) vs Eta, pT < 75",
			       12, 0., 2.4, -10., 10.);
  ResidualPt[4] = new TProfile("", "(delta 1/pT)/(1/pT) vs Eta, pT < 100",
			       12, 0., 2.4, -10., 10.);
  ResidualPt[5] = new TProfile("", "(delta 1/pT)/(1/pT) vs Eta, pT < 1000",
			       12, 0., 2.4, -10., 10.);

  // 1/p resolution for GMR and TMR, separately for barrel, overlap and endcap.
  MuInvPRes[0][0] = new TH1F("", "1/P res, GMR, barrel",  100, -0.3, 0.3);
  MuInvPRes[0][1] = new TH1F("", "1/P res, GMR, overlap", 100, -0.3, 0.3);
  MuInvPRes[0][2] = new TH1F("", "1/P res, GMR, endcap",  100, -0.3, 0.3);

  MuInvPRes[1][0] = new TH1F("", "1/P res, TMR, barrel",  100, -0.3, 0.3);
  MuInvPRes[1][1] = new TH1F("", "1/P res, TMR, overlap", 100, -0.3, 0.3);
  MuInvPRes[1][2] = new TH1F("", "1/P res, TMR, endcap",  100, -0.3, 0.3);

  // 1/pT resolution for GMR and TMR, separately for barrel, overlap and endcap
  MuInvPtRes[0][0] = new TH1F("", "1/pT res, GMR, barrel",  100, -0.3, 0.3);
  MuInvPtRes[0][1] = new TH1F("", "1/pT res, GMR, overlap", 100, -0.3, 0.3);
  MuInvPtRes[0][2] = new TH1F("", "1/pT res, GMR, endcap",  100, -0.3, 0.3);

  MuInvPtRes[1][0] = new TH1F("", "1/pT res, TMR, barrel",  100, -0.3, 0.3);
  MuInvPtRes[1][1] = new TH1F("", "1/pT res, TMR, overlap", 100, -0.3, 0.3);
  MuInvPtRes[1][2] = new TH1F("", "1/pT res, TMR, endcap",  100, -0.3, 0.3);

  // pT res for all fits
  MuPtDiff[0] = new TH1F("",  "L3 pT - Gen pT", 100, -1500., 1500.);
  MuPtDiff[1] = new TH1F("", "TMR pT - Gen pT", 100, -1500., 1500.);
  TotInvPtRes[0] = new TH1F("", 
			     "(L3 1/pT - Gen 1/pT)/(Gen 1/pT)", 100, -.5, .5);
  TotInvPtRes[1] = new TH1F("", "(Tracker Only 1/pT - Gen 1/pT)/(Gen 1/pT)", 
			     100, -.5, .5);
  TotInvPtRes[2] = 
    new TH1F("", "(Tracker + 1 mu-station 1/pT - Gen 1/pT)/(Gen 1/pT)", 
	     100, -.5, .5);
  //TotInvPtRes[3] = new TH1F("", "(ABCM 1/pT - Gen 1/pT)/(Gen 1/pT)", 
  //		    100, -.5, .5);
  TotInvPtRes[3] = new TH1F("", "(Outer 1/pT - Gen 1/pT)/(Gen 1/pT)", 
			    100, -.5, .5);
  // 1/pT pull for all fits
  TotInvPtPull[0] = new TH1F("", "L3 1/pT Pull", 100, -10., 10.);
  TotInvPtPull[1] = new TH1F("", "Tracker Only 1/pT Pull", 100, -10., 10.);
  TotInvPtPull[2] = new TH1F("", "Tracker + 1 mu-station 1/pT Pull", 
			     100, -10., 10.);
  //TotInvPtPull[3] = new TH1F("", "ABCM 1/pT Pull", 100, -10., 10.);
  TotInvPtPull[3] = new TH1F("", "Outer 1/pT Pull", 100, -10., 10.);

  // pT res for all fits (split up by barrel and endcap)
  InvPtRes[0][0] = new TH1F("eta < 1.04", "L3", 100, -.3, .3); 
  InvPtRes[0][1] = new TH1F("eta > 1.04", "L3", 100, -.3, .3); 
  InvPtRes[1][0] = new TH1F("eta < 1.04", "Tracker Only", 100, -.3, .3); 
  InvPtRes[1][1] = new TH1F("eta > 1.04", "Tracker Only", 100, -.3, .3); 
  InvPtRes[2][0] = new TH1F("eta < 1.04", "Tracker + 1 mu-station",
			    100, -.3, .3); 
  InvPtRes[2][1] = new TH1F("eta > 1.04", "Tracker + 1 mu-station", 
			    100, -.3, .3); 
  InvPtRes[3][0] = new TH1F("eta < 1.04", "GMR", 100, -.3, .3); 
  InvPtRes[3][1] = new TH1F("eta > 1.04", "GMR", 100, -.3, .3); 

  // 1/pT pull for all fits (split up by barrel and endcap)
  InvPtPull[0][0] = new TH1F("eta < 1.04", "L3", 100, -10., 10.); 
  InvPtPull[0][1] = new TH1F("eta > 1.04", "L3", 100, -10., 10.); 
  InvPtPull[1][0] = new TH1F("eta < 1.04", "Tracker Only", 100, -10., 10.); 
  InvPtPull[1][1] = new TH1F("eta > 1.04", "Tracker Only", 100, -10., 10.);
  InvPtPull[2][0] = new TH1F("eta < 1.04", "Tracker + 1 mu-station", 
				100, -10., 10.); 
  InvPtPull[2][1] = new TH1F("eta > 1.04", "Tracker + 1 mu-station", 
				100, -10., 10.); 
  InvPtPull[3][0] = new TH1F("eta < 1.04", "GMR", 100, -10., 10.); 
  InvPtPull[3][1] = new TH1F("eta > 1.04", "GMR", 100, -10., 10.); 
}

void Zprime2muResolution::BookDilResHistos(){
  // Dilepton mass spectra:
  //  - dilepton only (i=0);
  //  - dileptom plus nearby photon candidates (i=1).
  int nbins = binSize;
  double mass_min = lowerMassWin;
  double mass_max = upperMassWin;
  for (int i = lgen; i <= MAX_LEVELS; i++) {
    for (int j = 0; j < 2; j++) {
      if (i == l1 || i == l2) {
	DilMassComp[i][j] = new TH1F("", str_level[i].c_str(), nbins,
				     0., mass_max);
      }
      else {
	DilMassComp[i][j] = new TH1F("", str_level[i].c_str(), nbins,
				     mass_min, mass_max);
      }
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
	DilMassRes[i][j] = new TH1F("", str_level[i].c_str(), 100, -1.,  1.);
      }
      else {
	DilMassRes[i][j] = new TH1F("", str_level[i].c_str(), 100, -0.3, 0.3);
      }
    }
  }
  for (int j = 0; j < 6; j++) {
    DilMassRes[MAX_LEVELS+1][j] = new TH1F("", "TMR, D0 method",
					   100, -0.3, 0.3);
  }

  // Dilepton pT resolution
  DilPtRes[0] = 0;
  for (int i = l1; i <= MAX_LEVELS; i++) {
    DilPtRes[i] = new TH1F("", str_level[i].c_str(), 100, -1., 4.);
  }

  MuPVsMuM[0] = new TH2F("", "pT Rel Error for dM > 0.1", 100, 0., 1.,
			 100, 0., 1.);
  MuPVsMuM[1] = new TH2F("", "pT Rel Error for dM < 0.1", 100, 0., 1., 
			 100, 0., 1.);
}

void Zprime2muResolution::BookChargeResHistos() {
  for (int i = l1; i <= MAX_LEVELS; i++) {
    string tit = str_level[i] + " charge - Gen charge";
    QRes[i] = new TH1F("", tit.c_str(), 7, -3.5, 3.5);
  }

  QResVsPt[0][0] = new TH1F("", "L1 pT, right charge", 50, 0., 150.);
  QResVsPt[1][0] = new TH1F("", "L2 pT, right charge", 50, 0., 500.);
  QResVsPt[2][0] = new TH1F("", "L3 pT, right charge", 50, 0., 500.);
  QResVsPt[0][1] = new TH1F("", "L1 pT, wrong charge", 50, 0., 150.);
  QResVsPt[1][1] = new TH1F("", "L2 pT, wrong charge", 50, 0., 500.);
  QResVsPt[2][1] = new TH1F("", "L3 pT, wrong charge", 50, 0., 500.);
  QResVsInvPt[0][0] = new TH1F("", "1/L1pT, right charge", 50, 0., 0.04);
  QResVsInvPt[1][0] = new TH1F("", "1/L2pT, right charge", 50, 0., 0.04);
  QResVsInvPt[2][0] = new TH1F("", "1/L3pT, right charge", 50, 0., 0.04);
  QResVsInvPt[0][1] = new TH1F("", "1/L1pT, wrong charge", 50, 0., 0.04);
  QResVsInvPt[1][1] = new TH1F("", "1/L2pT, wrong charge", 50, 0., 0.04);
  QResVsInvPt[2][1] = new TH1F("", "1/L3pT, wrong charge", 50, 0., 0.04);
  QResVsP[0][0]  = new TH1F("", "L1 P, right charge", 50, 0., 750.);
  QResVsP[1][0]  = new TH1F("", "L2 P, right charge", 50, 0., 
			    upperMassWin);
  QResVsP[2][0]  = new TH1F("", "L3 P, right charge", 50, 0., 
			    upperMassWin);
  QResVsP[0][1]  = new TH1F("", "L1 P, wrong charge", 50, 0., 750.);
  QResVsP[1][1]  = new TH1F("", "L2 P, wrong charge", 50, 0., 
			    upperMassWin);
  QResVsP[2][1]  = new TH1F("", "L3 P, wrong charge", 50, 0., 
			    upperMassWin);
  QResVsInvP[0][0] = new TH1F("", "1/L1P, right charge", 50, 0., 0.02);
  QResVsInvP[1][0] = new TH1F("", "1/L2P, right charge", 50, 0., 0.02);
  QResVsInvP[2][0] = new TH1F("", "1/L3P, right charge", 50, 0., 0.02);
  QResVsInvP[0][1] = new TH1F("", "1/L1P, wrong charge", 50, 0., 0.02);
  QResVsInvP[1][1] = new TH1F("", "1/L2P, wrong charge", 50, 0., 0.02);
  QResVsInvP[2][1] = new TH1F("", "1/L3P, wrong charge", 50, 0., 0.02);

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
  bool   accept2_4 = false, accept2_1 = false;
  double phival, gen_mass = -999.;
  TLorentzVector tempV;
  vector<zp2mu::Muon>::const_iterator   pmu;
  vector<zp2mu::DiMuon>::const_iterator pdi;
  vector<zp2mu::DiMuon> diMuon;

  // Get origin of generated muons.
  for (pmu = allMuons[0].begin(); pmu != allMuons[0].end(); pmu++) {
    Origin[0]->Fill(getOrigin(pmu->genMotherId()));
  }

  // some additional plots for improved diagnostics
  fillEffHistos();
  fillPtResHistos(debug);
  fillChargeResHistos(debug);

  for (int i_rec = 0; i_rec < NUM_REC_LEVELS; i_rec++) {
    // Total number of muons found in each level of reconstruction.
    NMuons[i_rec][0]->Fill((double)allMuons[i_rec].size());

    // Fill y, eta, phi, p, pt, pz, pt vs eta, p vs eta plots for all muons.
    fillMuonHistos(i_rec, debug);
  }

  // Main loop over generated, L1, L2, L3 and off-line information.
  for (int i_rec = 0; i_rec <= MAX_LEVELS; i_rec++) {

    if (i_rec < MAX_LEVELS) diMuon = allDiMuons[i_rec];
    else if (i_rec == MAX_LEVELS) diMuon = bestDiMuons;

    if (i_rec == 0) {
      if (diMuon.size() == 1) { // one di-muon at generation
	gen_mass = diMuon[0].dimuV().M()/1000.; // (in TeV)
	// Check whether both muons are inside the full eta coverage.
	if (fabs(diMuon[0].muPlus().eta())  < ETA_CUT &&
	    fabs(diMuon[0].muMinus().eta()) < ETA_CUT) {
	  accept2_4 = true;
	  // At least one muon is in the limited muon trigger acceptance.
	  if (fabs(diMuon[0].muPlus().eta())  < TRIGGER_ETA_CUT[l1] ||
	      fabs(diMuon[0].muMinus().eta()) < TRIGGER_ETA_CUT[l1]) {
	    accept2_1 = true;
	  }
	}
      }
      else {
	if (!doingHiggs)
	  edm::LogWarning("Zprime2muResolution")
	    << "+++ Warning in calcResolution: " << diMuon.size()
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
	NMuons[i_rec][3]->Fill((double)allMuons[i_rec].size());
      }

      // More info for L1
      if (i_rec == l1) {
	if (decision) {
	  if (allMuons[i_rec].size() == 1)
	    L1TrigPassSingleMu->Fill(allMuons[i_rec][0].pt());
	  else if (allMuons[i_rec].size() == 2)
	    L1TrigPassMu2VsMu1->Fill(allMuons[i_rec][0].pt(),
				     allMuons[i_rec][1].pt());
	}
	else {
	  if (allMuons[i_rec].size() == 1)
	    L1TrigFailSingleMu->Fill(allMuons[i_rec][0].pt());
	  else if (allMuons[i_rec].size() == 2)
	    L1TrigFailMu2VsMu1->Fill(allMuons[i_rec][0].pt(),
				     allMuons[i_rec][1].pt());
	}
      }
    }

    // Number of muon hits at Level-2
    if (i_rec == l2) {
      for (pmu = allMuons[l2].begin(); pmu != allMuons[l2].end(); pmu++) {
	L2MuonHits->Fill(pmu->nmuHits());
      }
    }
    // Number of tracker (silicon + pixel) hits at Level-3
    if (i_rec == l3) {
      for (pmu = allMuons[l3].begin(); pmu != allMuons[l3].end(); pmu++) {
	L3TrackerHits->Fill(pmu->npixHits() + pmu->nsilHits());
      }
    }
    // Number of muon hits and chi2/d.o.f. for off-line (GMR) muons.
    // Unfortunately, the number of muon hits includes RPC hits, and there
    // seems to be no way to subtract them.
    if (i_rec == lgmr) {
      for (pmu = allMuons[lgmr].begin(); pmu != allMuons[lgmr].end(); pmu++) {
	int muHits = pmu->nmuHits();
	double chi2dof = pmu->chi2()/pmu->dof();
	if (fabs(pmu->eta()) < 0.9) {
	  GMRMuonHits[0]->Fill(muHits); // barrel
	  GMRChi2dof[0]->Fill(chi2dof);
	}
	else if (fabs(pmu->eta()) < 1.2) {
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
    if (i_rec <= l3) NMuons[i_rec][2]->Fill((double)allMuons[i_rec].size());

    // Fill histos with opposite sign and same sign dileptons
    if (i_rec > lgen && i_rec <= l3) fillSignOfDilepton(i_rec, debug);

    // If an opposite-sign dilepton is found, fill histos.  A dilepton
    // should always be found at the generated level.
    for (pdi = diMuon.begin(); pdi != diMuon.end(); pdi++) {

      //if (i_rec == 3) dumpDilEvents(2, i_rec, *pdi);

      // Fill mass histograms for all dileptons found, before we apply
      // the "off-line" track quality cuts.
      if (i_rec <= l3) AllDilMass[i_rec]->Fill(pdi->dimuV().M());
      if (i_rec == l3) {
	if (diMuon.size() == allDiMuons[0].size())
	  AllDilMassRes[i_rec-1]->
	    Fill(pdi->dimuV().M()-allDiMuons[0][pdi->id()].dimuV().M());
	if (allDiMuons[0].size() <= 0)
	  edm::LogWarning("Zprime2muResolution")
	    << "+++ Warning: no dilepton in the MC! +++\n";
      }

      // Check the "off-line" track quality and apply the cuts
      int qcut_mum = 0, qcut_mup = 0; // quality cut number
      if (DO_QCUTS ?
	  dilQCheck(*pdi, QSEL, qcut_mum, qcut_mup, debug) : true) {
  
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
	  Origin[1]->Fill(getOrigin(pdi->muPlus().genMotherId()));
	  Origin[1]->Fill(getOrigin(pdi->muMinus().genMotherId()));
	}
	else if (i_rec <= l3) NumDilVsRec->Fill(i_rec+1);

	if (debug) {
	  LogTrace("Zprime2muResolution") << str_level[i_rec] << " values:";
	  LogTrace("Zprime2muResolution")
	    << "#|Charge |   Eta   |   Phi   |    P    |"
	    << "    Pt   |    Pz   |   Rap   |  Mass  ";
	  LogTrace("Zprime2muResolution")
	    << "-------------------------------------------------"
	    << "------------------------------";
	}

	// Get number of muons in events with dilepton.
	if (i_rec <= l3)
	  NMuons[i_rec][1]->Fill((double)allMuons[i_rec].size());

	// Dilepton invariant mass.
	if (i_rec <= l3) {
	  DilMass[i_rec]->Fill(pdi->dimuV().M());
	  DilMassVsEta[i_rec]->Fill(pdi->dimuV().Eta(),    pdi->dimuV().M());
	  DilMassVsY[i_rec]->Fill(pdi->dimuV().Rapidity(), pdi->dimuV().M());
	}

	TLorentzVector vmum, vmup;
	zp2mu::Muon mum = pdi->muMinus();
	zp2mu::Muon mup = pdi->muPlus();
	vmum.SetPtEtaPhiM(mum.pt(), mum.eta(), mum.phi(), MUMASS);
	vmup.SetPtEtaPhiM(mup.pt(), mup.eta(), mup.phi(), MUMASS);

	// A few mu- vs mu+ scatter plots
	if (i_rec <= l3) {
	  MuMVsMuP[i_rec][0]->Fill(vmup.Eta(), vmum.Eta());
	  MuMVsMuP[i_rec][1]->Fill(mup.phi(),  mum.phi());
	  MuMVsMuP[i_rec][2]->Fill(vmup.Pt(),  vmum.Pt());
	}

	// Loop over "particle types": 0 stands for mu-, 1 - for mu+,
	// 2 - for opposite sign dilepton.
	for (int i_part = 0; i_part < 3; i_part++) {

	  // Split histos into mu-, and mu+, and dilepton.
	  if (i_part == 0)      tempV = vmum;
	  else if (i_part == 1) tempV = vmup;
	  else if (i_part == 2) tempV = pdi->dimuV();

	  if (debug) {
	    std::ostringstream strstrm;
	    strstrm << setprecision(5) << i_part << "| ";
	    if (i_part == 0)
	      strstrm << setw(4) << mum.charge() << "  | ";
	    else if (i_part == 1)
	      strstrm << setw(4) << mup.charge() << "  | ";
	    else
	      strstrm                        << "      | ";
	    strstrm << endl << tempV;
	    LogTrace("Zprime2muAnalysis") << strstrm.str();
	  }

	  // Shift phi range from (-pi; pi) to (0; 2*pi).
	  if (tempV.Phi() < 0.)
	    phival = tempV.Phi() + 2.*TMath::Pi();
	  else
	    phival = tempV.Phi();

	  if (i_rec <= l3) {
	    Eta[i_rec][i_part]->Fill(tempV.Eta());
	    Phi[i_rec][i_part]->Fill(phival);
	    P[i_rec][i_part]->Fill(tempV.P());
	    Pt[i_rec][i_part]->Fill(tempV.Pt());
	    Pz[i_rec][i_part]->Fill(abs(tempV.Pz()));
	    Rapidity[i_rec][i_part]->Fill(tempV.Rapidity());
	    PVsEta[i_rec][i_part]->Fill(tempV.Eta(),  tempV.P());
	    PtVsEta[i_rec][i_part]->Fill(tempV.Eta(), tempV.Pt());
	  }

	  int id = pdi->id();
	  // Calculate the resolution of each rec muon.
	  if (i_rec > 0 && i_part < 2) {
	    zp2mu::Muon murec, mugen, mul1, mul2;

	    // Compare with generated values.
	    if (diMuon.size() == allDiMuons[0].size()) {
	      if (i_part == 0) {
		murec = mum;
		mugen = allDiMuons[0][id].muMinus();
	      }
	      else {
		murec = mup;
		mugen = allDiMuons[0][id].muPlus();
	      }
	      // Make sure the tracks are reasonably close together.
	      if (matchTracks(murec, mugen)) {
		// Fill resolution histos.
		double pres = murec.p()/mugen.p() - 1.;
		EtaRes[i_rec]->Fill(murec.eta()-mugen.eta());
		PhiRes[i_rec]->Fill(murec.phi()-mugen.phi());
		PtDiff[i_rec]->Fill(murec.pt()-mugen.pt());
		PRes[i_rec]->Fill(pres);
		PResVsP[i_rec]->Fill(mugen.p(), pres*pres);
		if (i_rec <= l3) {
		  GenPhiResVsPhi[i_rec-1]->Fill(mugen.phi(),
						fabs(murec.phi()-mugen.phi()));
		  GenInvPtRes[i_rec-1]->Fill(((1./murec.pt())-(1./mugen.pt()))
					     /(1./mugen.pt()));
		  //GenInvPtResVsPt[i_rec-1]->Fill(mugen.pt(),
		  //      (fabs(1./murec.pt()-1./mugen.pt()))/(1./mugen.pt()));
		  GenInvPtResVsPt[i_rec-1]->Fill(mugen.pt(),
						 (fabs(murec.pt()-mugen.pt()))/(mugen.pt()));
		  GenInvPRes[i_rec-1]->Fill(mugen.p()/murec.p() - 1.);
		  GenPResVsPt[i_rec-1]->Fill(mugen.pt(), fabs(pres));
		  GenInvPResVsPt[i_rec-1]->Fill(mugen.pt(),
						fabs(mugen.p()/murec.p()-1.));
		  GenEtaResScat[i_rec-1]->Fill(mugen.eta(), murec.eta());
		  GenPhiResScat[i_rec-1]->Fill(mugen.phi(), murec.phi());
		  GenPtResScat[i_rec-1]->Fill(mugen.pt(),   murec.pt());
		}
	      }
	    }

	    if ((i_rec == l2 || i_rec == l3) &&
		diMuon.size() == allDiMuons[1].size()) {
	      // Compare rec 2 and rec 3 values with rec 1 values.
	      if (i_part == 0) {
		murec = mum;
		mul1  = allDiMuons[1][id].muMinus();
	      }
	      else {
		murec = mup;
		mul1  = allDiMuons[1][id].muPlus();
	      }
	      if (matchTracks(murec, mul1)) {
		// Fill resolution histos.
		L1EtaRes[i_rec-2]->Fill(murec.eta()-mul1.eta());
		L1PhiRes[i_rec-2]->Fill(murec.phi()-mul1.phi());
		L1PtDiff[i_rec-2]->Fill(murec.pt() -mul1.pt());
		L1EtaResScat[i_rec-2]->Fill(mul1.eta(), murec.eta());
		L1PhiResScat[i_rec-2]->Fill(mul1.phi(), murec.phi());
		L1PtResScat[i_rec-2]->Fill(mul1.pt(),   murec.pt());
	      }
	      if (i_rec == l3 && diMuon.size() == allDiMuons[2].size()) {
		// Compare rec 2 with rec 3 values.
		if (i_part == 0) {
		  murec = mum;
		  mul2  = allDiMuons[2][id].muMinus();
		}
		else {
		  murec = mup;
		  mul2  = allDiMuons[2][id].muPlus();
		}
		if (matchTracks(murec, mul2)) {
		  L2EtaRes->Fill(murec.eta()-mul2.eta());
		  L2PhiRes->Fill(murec.phi()-mul2.phi());
		  L2PtDiff->Fill(murec.pt() -mul2.pt());
		  L2EtaResScat->Fill(mul2.eta(), murec.eta());
		  L2PhiResScat->Fill(mul2.phi(), murec.phi());
		  L2PtResScat->Fill(mul2.pt(),   murec.pt());
		}
	      }
	    }
	  }
	  // Opposite-sign dileptons
	  else if ((i_rec > 0 && i_rec <= l3) && i_part == 2) {
	    if (diMuon.size() == allDiMuons[0].size()) {
	      double genmass = allDiMuons[0][id].dimuV().M();
	      double dM   = tempV.M() - genmass;
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
    if (allDiMuons[0].size() != 2) // generated
      edm::LogWarning("Zprime2muResolution") 
	<< "+++ Warning: " << allDiMuons[0].size()
	<< " generated dilepton(s) were found +++\n"; 
    else {
      ZonDilMass[0]->Fill(allDiMuons[0][0].dimuV().M());
      ZofDilMass[0]->Fill(allDiMuons[0][1].dimuV().M());
    }
    for (pdi = allDiMuons[0].begin(); pdi != allDiMuons[0].end(); pdi++) {
      Eta4muons->Fill(pdi->muPlus().eta());
      Eta4muons->Fill(pdi->muMinus().eta());
      Pt4muons->Fill(pdi->muPlus().pt());
      Pt4muons->Fill(pdi->muMinus().pt());
    }

    // For events passing the trigger
    if (passTrigger()) {
      if (allDiMuons[l3].size() >= 2) { // HLT
	ZonDilMass[3]->Fill(allDiMuons[l3][0].dimuV().M());
	ZofDilMass[3]->Fill(allDiMuons[l3][1].dimuV().M());
      }
    }
  }
}

void Zprime2muResolution::fillEffHistos() {
  // Efficiency to reconstruct muons from Z' decays at various trigger 
  // levels and by various off-line reconstructors.  L1/HLT efficiencies
  // are not included.
  vector<zp2mu::Muon>::const_iterator pmu;
  double gen_eta, gen_phi, gen_pt;

  // Efficiency for single muons
  for (pmu = allMuons[lgen].begin(); pmu != allMuons[lgen].end(); pmu++) {
    if (pmu->genMotherId() == 32) {
      gen_eta = pmu->eta();
      gen_phi = pmu->phi();
      gen_pt  = pmu->pt();

      EffVsEta[0]->Fill(gen_eta);
      EffVsPhi[0]->Fill(gen_phi);
      EffVsPt[0]->Fill(gen_pt);

      for (int i = l1; i < MAX_LEVELS; i++) {
	// if (passTrigger(i) && pmu->closestId(i) >= 0) {
	if (pmu->closestId(i) >= 0) {
	  EffVsEta[i]->Fill(gen_eta);
	  EffVsPhi[i]->Fill(gen_phi);
	  EffVsPt[i]->Fill(gen_pt);
	}
      }
      for (vector<zp2mu::Muon>::const_iterator pbest = bestMuons.begin();
	   pbest != bestMuons.end(); pbest++) {
	int rec_level = pbest->recLevel();
	zp2mu::Muon recmu = muonRef(rec_level, pbest->closestId(rec_level));
	if (pmu->closestId(rec_level) >= 0 &&
	    pmu->closestId(rec_level) == recmu.id()) {
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
  if (allDiMuons[lgen].size() == 1) { // one generated dimuon
    gen_mass = allDiMuons[lgen][0].dimuV().M()/1000.; // (in TeV)
    RecMass[0][0]->Fill(gen_mass);

    // Events with both muons inside the full eta coverage, and passing
    // the trigger
    if (fabs(allDiMuons[lgen][0].muPlus().eta())  < ETA_CUT &&
	fabs(allDiMuons[lgen][0].muMinus().eta()) < ETA_CUT) {
      if (passTrigger()) RecMass[0][1]->Fill(gen_mass);
    }
  }

  if (passTrigger()) {
    if (allMuons[l3].size() > 1)     // At least two muons found by L3
      RecMass[1][0]->Fill(gen_mass);
    if (allDiMuons[l3].size() > 0)   // Opposite-sign dimuon found by L3
      RecMass[1][1]->Fill(gen_mass);
    if (allMuons[lgmr].size() > 1)   // At least two muons found by GMR
      RecMass[2][0]->Fill(gen_mass);
    if (allDiMuons[lgmr].size() > 0) // Opposite-sign dimuon found by GMR
      RecMass[2][1]->Fill(gen_mass);
    if (bestMuons.size() > 1)        // At least two muons found by TMR
      RecMass[3][0]->Fill(gen_mass);
    if (bestDiMuons.size() > 0)      // Opposite-sign dimuon found by TMR
      RecMass[3][1]->Fill(gen_mass);
  }
}

void Zprime2muResolution::fillPtResHistos(const bool debug) {
  // Function to fill histos of pt resolutions.
  int idh, gen_id, tracker_id, gmr_id, fms_id, bar_end;
  double gen_pt, gen_p, gen_eta, residual, pt_pull;
  zp2mu::Muon genmu, tkmu, gmrmu, fmsmu;
  vector<zp2mu::Muon> muons;
  vector<zp2mu::Muon>::const_iterator pmu;

  // Inverse momentum and 1/pT resolution for GMR and TMR, separately for
  // barrel, endcap, and their overlap.
  if (passTrigger()) { // event passed the trigger
    for (int i = 0; i < 2; i++) {
      if (i == 0)      muons = allMuons[lgmr];
      else if (i == 1) muons = bestMuons;
      for (pmu = muons.begin(); pmu != muons.end(); pmu++) {
	gen_id = pmu->matchedId(lgen); // Generated muon closest to rec muon
	if (gen_id != -999) {
	  genmu = muonRef(lgen, gen_id);

	  gen_eta = genmu.eta();
	  gen_p   = genmu.p();
	  gen_pt  = genmu.pt();
	  double invpres = gen_p/pmu->p() - 1.;
	  double invptres = gen_pt/pmu->pt() - 1.;
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
  for (pmu = bestMuons.begin(); pmu != bestMuons.end(); pmu++) {
    if (pmu->pt() < 100.) continue;
    gen_id = pmu->matchedId(lgen);
    if (gen_id != -999) {
      genmu   = muonRef(lgen, gen_id);
      gen_pt  = genmu.pt();
      gen_eta = genmu.eta();
      //if ((gen_pt > 800) && (gen_pt < 1200)) {
	//if ((fabs(gen_eta) > 1.2) && (fabs(gen_eta) < 2.1)) {
	//if ((fabs(gen_eta) > 0.9) && (fabs(gen_eta) < 1.2)) {
	//if (fabs(gen_eta) < 0.9) {
	MuPtDiff[1]->Fill(pmu->pt() - gen_pt);
	//}
      //}
    }
  }

  for (pmu = allMuons[l3].begin(); pmu != allMuons[l3].end(); pmu++) {

    // matchStudy(*pmu);

    // Only look at events where reconstructed pT > 100 GeV.
    if (pmu->pt() < 100.) continue;

    // Get index of the generated muon closest to the L3 muon.  The function
    // will return -999 if track was not found.
    gen_id = pmu->matchedId(lgen);
    if (gen_id != -999) {
      genmu = muonRef(lgen, gen_id);

      // Calculation of 1/pT residual.
      gen_eta  = genmu.eta();
      gen_pt   = genmu.pt();
      residual = (1./pmu->pt() - 1./gen_pt)/(1./gen_pt);
      //if ((gen_pt > 800) && (gen_pt < 1200)) {
	//if ((fabs(gen_eta) > 1.2) && (fabs(gen_eta) < 2.1)) {
	//if ((fabs(gen_eta) > 0.9) && (fabs(gen_eta) < 1.2)) {
        //if (fabs(gen_eta) < 0.9) {
      MuPtDiff[0]->Fill(pmu->pt() - gen_pt);
        //}
      //}

      // Split histograms into gen_pt < 10, 25, 50, 75, 100, 1000.
      idh = -1;
      if (gen_pt < 10.)        idh = 0;
      else if (gen_pt < 25.)   idh = 1;
      else if (gen_pt < 50.)   idh = 2;
      else if (gen_pt < 75.)   idh = 3;
      else if (gen_pt < 100.)  idh = 4;
      else if (gen_pt < 1000.) idh = 5;
      if (idh >= 0) ResidualPt[idh]->Fill(abs(gen_eta), abs(residual));

      if (debug) {
	LogTrace("Zprime2muResolution")
	  << " l3 muon = "         << pmu->id()
	  << " closest genmuon = " << genmu.id();
	LogTrace("Zprime2muResolution")
	  << "   Ptgen: "  << gen_pt  << " Pt3: " << pmu->pt()
	  << " Etagen: "   << gen_eta << " Residual: " << residual;
      }
    }
  }

  // Find the 1/pT residual and pulls for all alternative fits to L3.
  // Use same number of muons for all fits so that we can directly compare
  // the change in RMS, sigma and pulls between fits.  Loop over all of L3 
  // muons first.
  if (!passTrigger()) return;

  for (pmu = allMuons[l3].begin(); pmu != allMuons[l3].end(); pmu++) {
    // Only look at events where reconstructed pT > 100 GeV.
    if (pmu->pt() < 100.) continue;

    // Get generated muon closest to the L3 muon.  This gen muon will be used
    // for comparison to all fits.
    gen_id = pmu->matchedId(lgen);

    // Find corresponding tracker-only, FMS and GMR tracks
    tracker_id = pmu->matchedId(ltk);
    fms_id     = pmu->matchedId(lfms);
    gmr_id     = pmu->matchedId(lgmr);

    // Fill only for muons that have tracker-only and FMS tracks, and
    // are matched to the generated tracks
    if (gen_id != -999 && tracker_id != -999 && fms_id != -999) {

      genmu  = muonRef(lgen, gen_id);
      gen_pt = genmu.pt();

      // Residual and pulls for on-line and off-line reconstructions.
      // separate each fit by barrel or endcap
      if (abs(pmu->eta()) < ENDCAP_BARREL_CUT) bar_end = 0;
      else bar_end = 1;
      if (pmu->pt() > 0.) {
	residual  = (1./pmu->pt() - 1./gen_pt)/(1./gen_pt);
	TotInvPtRes[0]->Fill(residual);
	InvPtRes[0][bar_end]->Fill(residual);
	if (pmu->errInvPt() > 0.) {
	  pt_pull = (1./pmu->pt() - 1./gen_pt)/pmu->errInvPt();
	  TotInvPtPull[0]->Fill(pt_pull);
	  InvPtPull[0][bar_end]->Fill(pt_pull);
	}
      }

      // Residual and pulls for tracker-only fit.
      tkmu = muonRef(ltk, tracker_id);
      if (abs(tkmu.eta()) < ENDCAP_BARREL_CUT) bar_end = 0;
      else bar_end = 1;
      if (tkmu.pt() > 0.) {
	residual  = (1./tkmu.pt() - 1./gen_pt)/(1./gen_pt);
	TotInvPtRes[1]->Fill(residual);
	InvPtRes[1][bar_end]->Fill(residual);
	if (tkmu.errInvPt() > 0.) {
	  pt_pull = (1./tkmu.pt() - 1./gen_pt)/tkmu.errInvPt();
	  TotInvPtPull[1]->Fill(pt_pull);
	  InvPtPull[1][bar_end]->Fill(pt_pull);
	}
      }

      // Residual and pulls for internal-seeded GMR.
      if (gmr_id != -999) {
	gmrmu = muonRef(lgmr, gmr_id);
	if (abs(gmrmu.eta()) < ENDCAP_BARREL_CUT) bar_end = 0;
	else bar_end = 1;
	if (gmrmu.pt() > 0.) {
	  residual  = (1./gmrmu.pt() - 1./gen_pt)/(1./gen_pt);
	  InvPtRes[3][bar_end]->Fill(residual);
	  if (gmrmu.errInvPt() > 0.) {
	    pt_pull = (1./gmrmu.pt() - 1./gen_pt)/gmrmu.errInvPt();
	    InvPtPull[3][bar_end]->Fill(pt_pull);
	  }
	}
      }

      // pT measurement at outer surface of tracker.
      // (Should not be the default)
      if (pmu->forwardPt() > 0.) {
	residual  = (1./pmu->forwardPt()-1./gen_pt)/(1./gen_pt);
	TotInvPtRes[3]->Fill(residual);
	if (pmu->errForwardInvPt() > 0.) {
	  pt_pull = (1./pmu->forwardPt()-1./gen_pt)/pmu->errForwardInvPt();
	  TotInvPtPull[3]->Fill(pt_pull);
	}
      }

      // Residual and pulls using the innermost muon hit for
      // pT measurement.  This is alleged to be the best measurement
      // of pT in the detector.
      fmsmu = muonRef(lfms, fms_id);
      if (abs(fmsmu.eta()) < ENDCAP_BARREL_CUT) bar_end = 0;
      else bar_end = 1;
      if (fmsmu.pt() > 0.) {
	residual  = (1./fmsmu.pt() - 1./gen_pt)/(1./gen_pt);
	TotInvPtRes[2]->Fill(residual);
	InvPtRes[2][bar_end]->Fill(residual);
	if (fmsmu.errInvPt() > 0.) {
	  pt_pull = (1./fmsmu.pt() - 1./gen_pt)/fmsmu.errInvPt();
	  TotInvPtPull[2]->Fill(pt_pull);
	  InvPtPull[2][bar_end]->Fill(pt_pull);
	}
      }
    }
  }
}

void Zprime2muResolution::fillChargeResHistos(const bool debug) {
  // This function makes histos of the difference between charge assignment
  // for various levels of reconstruction.
  int gen_id;
  int deltaQ, idh;
  zp2mu::Muon genmu;
  vector<zp2mu::Muon>::const_iterator pmu;

  // Loop over all reconstructed levels.
  for (int rec = l1; rec < MAX_LEVELS; rec++) {

    // Loop over all muons at a given level
    for (pmu = allMuons[rec].begin(); pmu != allMuons[rec].end(); pmu++) {

      // Find charge of closest gen muon and store difference in histogram
      gen_id = pmu->matchedId(lgen);
      if (gen_id != -999) {
	genmu = muonRef(lgen, gen_id);
	deltaQ = pmu->charge() - genmu.charge();
	QRes[rec]->Fill(deltaQ);
	if (rec <= l3) {
	  double pt = pmu->pt();
	  double p  = pmu->p();
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
	    << "Charge at rec level " << rec << ": " << pmu->charge()
	    << ", Closest Gen: " << genmu.charge();
	}
      }
    }
  }

  // Same for TMR muons
  for (pmu = bestMuons.begin(); pmu != bestMuons.end(); pmu++) {
    gen_id = pmu->matchedId(lgen);
    if (gen_id != -999) {
      genmu = muonRef(lgen, gen_id);
      deltaQ = pmu->charge() - genmu.charge();
      QRes[MAX_LEVELS]->Fill(deltaQ);
#ifdef IF_NEEDED
      double pt = pmu->pt();
      double p  = pmu->p();
      if      (deltaQ == 0)      idh =  0; // correct charge assignment
      else if (abs(deltaQ) == 2) idh =  1; // wrong assignment
      else throw cms::Exception("Zprime2muResolution")
	<< "+++ Impossible deltaQ in fillChargeResHistos() +++\n";
      QResVsPt[MAX_LEVELS][idh]->Fill(pt);
      QResVsInvPt[MAX_LEVELS][idh]->Fill(1./pt);
      QResVsP[MAX_LEVELS][idh]->Fill(p);
      QResVsInvP[MAX_LEVELS][idh]->Fill(1./p);
#endif
      if (debug) {
	LogTrace("Zprime2muResolution")
	  << "Charge at for best muons: " << pmu->charge()
	  << ", Closest Gen: " << genmu.charge();
      }
    }
  }
}

void Zprime2muResolution::fillMuonHistos(const int rec, const bool debug) {
  // Function for filling eta, y, phi, p, pz, pt, pt vs eta, p vs eta
  // for all muons for all levels of reconstruction.
  // Inputs: rec   = level of reconstruction
  //         debug = print statements

  TLorentzVector lVect;
  double phival;

  if (rec >= NUM_REC_LEVELS) {
    throw cms::Exception("Zprime2muResolution")
      << "+++ Unknown rec. level = " << rec << " in fillMuonHistos() +++";
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

  // Loop over all muons
  vector<zp2mu::Muon>::const_iterator pmu;
  for (pmu = allMuons[rec].begin(); pmu != allMuons[rec].end(); pmu++) {

    // Use TLorentzVector for easier access to variables such as Pz.
    lVect.SetPtEtaPhiM(pmu->pt(), pmu->eta(), pmu->phi(), MUMASS);

    // Special circumstances for phi
    if (lVect.Phi() < 0.)
      phival = lVect.Phi() + 2.*TMath::Pi();
    else
      phival = lVect.Phi();

    MuonEta[rec]->Fill(lVect.Eta());
    MuonRap[rec]->Fill(lVect.Rapidity());
    MuonPhi[rec]->Fill(phival);
    MuonPt[rec]->Fill(lVect.Pt());
    MuonPz[rec]->Fill(abs(lVect.Pz()));
    MuonP[rec]->Fill(lVect.P());
    MuonPtVsEta[rec]->Fill(lVect.Eta(), lVect.Pt());
    MuonPVsEta[rec]->Fill( lVect.Eta(), lVect.P() );

    if (debug)
      LogTrace("Zprime2muResolution") << setprecision(5) << pmu->id() << "| "
				      << setw(4) << pmu->charge() << "  | "
				      << endl << lVect;
  }
}

void Zprime2muResolution::fillSignOfDilepton(const int rec, const bool debug) {
  // Small function for filling histogram which contains number of opposite
  // sign dileptons, + like-sign dileptons, and - like-sign dileptons.  The
  // routine is called for each level of reconstruction.
  // Inputs: rec          = level of reconstruction
  //         debug        = true for print statements

  int total_charge = 0;
  int nMuons = 0;
  vector<zp2mu::Muon>::const_iterator pmu;

  // Only enter routine if there are 2 or more muons.
  // Add up total charge of all muons.
  if (allMuons[rec].size() > 1) {
    for (pmu = allMuons[rec].begin(); pmu != allMuons[rec].end(); pmu++) {
      total_charge += pmu->charge();
      nMuons++;
    }

    if (debug) {
      	LogTrace("Zprime2muResolution")
	  << "In fillSignOfDilepton routine: rec level = " << rec 
	  << " numMuons = " << nMuons << " total_charge = " << total_charge;
    }

    if (nMuons < 2) return;

    // If total charge is equal to number of muons, then positive dilepton
    // could be found
    if (total_charge == nMuons) {
      SignOfDil[rec]->Fill(2.);
      if (debug)
 	LogTrace("Zprime2muResolution") << "positive same-sign dilepton found";
    }

    // If total charge is equal to -(number of muons), then negative dilepton
    // could be found
    else if (total_charge == -nMuons) {
      SignOfDil[rec]->Fill(1.);
      if (debug)
 	LogTrace("Zprime2muResolution") << "negative same-sign dilepton found";
    }

    // Otherwise opposite sign dilepton was found
    else {
      SignOfDil[rec]->Fill(0.);
      if (debug)
	LogTrace("Zprime2muResolution") << "opposite-sign dilepton found";
    }
  }
}

void Zprime2muResolution::fillDilResHistos(const bool debug) {
  // Histograms to compare dilepton mass resolution between various fits.
  unsigned int n_gen = 0;
  vector<zp2mu::DiMuon> diMuons;
  TLorentzVector genDimuV, recDimuV, genResV, recResV;

  double mass_min = peakMass - 0.1*peakMass;
  double mass_max = peakMass + 0.1*peakMass;

  // Search for dileptons at all trigger levels and for all off-line
  // reconstruction methods.
  for (int rec = 0; rec <= MAX_LEVELS; rec++) {

    // Check whether the event passed the trigger.
    if (rec >= l1 && rec <= l3) {
      if (cutTrig[rec] && !passTrigger(rec)) return;
    }

    if      (rec <  MAX_LEVELS) diMuons = allDiMuons[rec];
    else if (rec == MAX_LEVELS) diMuons = bestDiMuons;

    unsigned int n_dil = diMuons.size();

    // Highest mass reconstructed at various trigger levels and by various
    // off-line fitting algorithms
    if (n_dil > 0) {
      DilMassComp[rec][0]->Fill(diMuons[0].dimuV().M());
      DilMassComp[rec][1]->Fill(diMuons[0].resV().M());
    }

    if (rec == lgen) {
      n_gen = n_dil;
      continue;
    }

    for (unsigned int i_dil = 0; i_dil < n_dil; i_dil++) {
      // Very poor match of dileptons; needs to be improved!
      if (n_dil == n_gen) {
	// Invariant mass resolution calculated relative to the mass
	// reconstructed from GEANT muons.
	genDimuV = allDiMuons[0][i_dil].dimuV();
	recDimuV = diMuons[i_dil].dimuV();
	genResV  = allDiMuons[0][i_dil].resV();
	recResV  = diMuons[i_dil].resV();
	double geant_mass  = genDimuV.M();
	double pythia_mass = genResV.M();
	double dil_mass    = recDimuV.M();
	double res_mass    = recResV.M();

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
	  // Calculated reweighted pT.
	  zp2mu::Muon mup = diMuons[i_dil].muPlus();
	  zp2mu::Muon mum = diMuons[i_dil].muMinus();
	  double pTp = mup.pt();
	  double pTm = mum.pt();
	  double wp  = 1./pTp/(1./pTp + 1./pTm);
	  double wm  = 1./pTm/(1./pTp + 1./pTm);
	  double pTw = pTp*wp + pTm*wm;
	  mup.setPt(pTw);
	  mum.setPt(pTw);
	  // LogTrace("Zprime2muResolution")
          //   << "pTp = " << pTp << " pTm = " << pTm
	  //   << " wp = " << wp << " wm = " << wm
	  //   << " pTw = " << pTw << " " << 2./(1./pTp + 1./pTm);

	  // Construct new dimuon and calculate its mass.
	  TLorentzVector vmup, vmum, vdil;
	  vmup.SetPtEtaPhiM(mup.pt(), mup.eta(), mup.phi(), MUMASS);
	  vmum.SetPtEtaPhiM(mum.pt(), mum.eta(), mum.phi(), MUMASS);
	  vdil = vmup + vmum;
	  zp2mu::DiMuon thisDiMu;
	  thisDiMu.fill(true, i_dil, rec);
	  thisDiMu.setMuPlus(mup);
	  thisDiMu.setMuMinus(mum);
	  thisDiMu.setDimuV(vdil);
	  double D0mass = thisDiMu.dimuV().M();
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
	DilPtRes[rec]->Fill((recDimuV.Pt() - genDimuV.Pt())/genDimuV.Pt());
      }
    }
  }

  // Studies of pT errors for muons in poorly-reconstructed and
  // well-reconstructed dimuons.  It was used to determine quality cuts by
  // looking at events in the tail of Drell-Yan dimuon invariant mass
  // distribution.
  int qcut_mum = 0, qcut_mup = 0;
  if (bestDiMuons.size() > 0) {
    if (DO_QCUTS ? 
	dilQCheck(bestDiMuons[0], QSEL, qcut_mum, qcut_mup, debug) : true) {

      // Require the reconstructed dilepton invariant mass be above a
      // certain mass value ("tail" of Drell-Yan distributions).
      double recm = bestDiMuons[0].resV().M(); // highest mass dilepton
      if (recm >= mass_min) {

	zp2mu::Muon mum = bestDiMuons[0].muMinus();
	zp2mu::Muon mup = bestDiMuons[0].muPlus();
	double err_pt_neg = mum.errPt()/mum.pt();
	double err_pt_pos = mup.errPt()/mup.pt();

	// Fill histos of relative pT error for mu+ vs mu- for two types of
	// dileptons.
	if (allDiMuons[lgen].size() > 0) {
	  double genm = allDiMuons[lgen][0].resV().M();

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

void Zprime2muResolution::DeleteHistos(){
  // ResHistos
  delete EventsInAccFailed;
  delete L1TrigFailSingleMu;
  delete L1TrigFailMu2VsMu1;
  delete L1TrigPassSingleMu;
  delete L1TrigPassMu2VsMu1;
  delete L2MuonHits;
  for (int i = 0; i < 3; i++) {
    delete GMRMuonHits[i];
    delete GMRChi2dof[i];
  }
  delete L3TrackerHits;
  delete NumDilVsRec;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 2; j++)
      delete RecMass[i][j];
  for (int i = 0; i <= MAX_LEVELS; i++) {
    delete EffVsEta[i];
    delete EffVsPhi[i];
    delete EffVsPt[i];
  }
  for (int i = 0; i < NUM_REC_LEVELS; i++) {
    for (int j = 0; j < 3; j++) {
      delete TrigResult[i][j];
      delete TrigMass[i][j];
    }
    delete AllDilMass[i];
    delete DilMass[i];
    delete DilMassVsEta[i];
    delete DilMassVsY[i];
    delete SignOfDil[i];
    delete MuonEta[i];
    delete MuonPhi[i];
    delete MuonRap[i];
    delete MuonP[i];
    delete MuonPt[i];
    delete MuonPz[i];
    delete MuonPVsEta[i];
    delete MuonPtVsEta[i];
    for (int j = 0; j < 4; j++)
      delete NMuons[i][j];
    for (int j = 0; j < 3; j++) {
      delete Eta[i][j];
      delete Phi[i][j];
      delete Rapidity[i][j];
      delete P[i][j];
      delete Pt[i][j];
      delete Pz[i][j];
      delete PVsEta[i][j];
      delete PtVsEta[i][j];
      delete MuMVsMuP[i][j];
    }
    delete ZonDilMass[i];
    delete ZofDilMass[i];
  }
  delete Eta4muons;
  delete Pt4muons;
  for (int k = l1; k <= MAX_LEVELS; k++) {
    delete EtaRes[k];
    delete PhiRes[k];
    delete PtDiff[k];
    delete PRes[k];
    delete PResVsP[k];
  }
  for (int k=0; k<3; k++) {
    delete GenPhiResVsPhi[k];
    delete GenInvPtRes[k];
    delete GenInvPtResVsPt[k];
    delete GenInvPRes[k];
    delete GenPResVsPt[k];
    delete GenInvPResVsPt[k];
    delete GenEtaResScat[k];
    delete GenPhiResScat[k];
    delete GenPtResScat[k];
    delete AllDilMassRes[k];
    delete GenDilMassRes[k];
    delete GenDilMassFrRes[k];
    delete GenDilMassResScat[k];
    delete GenDilMassFrResScat[k];
  }
  for (int l=0; l<2; l++) {
    delete Origin[l];
    delete L1EtaRes[l];
    delete L1PhiRes[l];
    delete L1PtDiff[l];
    delete L1EtaResScat[l];
    delete L1PhiResScat[l];
    delete L1PtResScat[l];
  }
  delete L2EtaRes;
  delete L2PhiRes;
  delete L2PtDiff;
  delete L2EtaResScat;
  delete L2PhiResScat;
  delete L2PtResScat;
  // end of ResHistos

  // PtRes histos
  for (int i = 0; i < 2; i++) {
    delete MuPtDiff[i];
    for (int j = 0; j < 3; j++) {
      delete MuInvPRes[i][j];
      delete MuInvPtRes[i][j];
    }
  }
  for (int i_pt = 0; i_pt < 6; i_pt++)
    delete ResidualPt[i_pt];
  for (int i = 0; i < 4; i++) {
    delete TotInvPtRes[i];
    delete TotInvPtPull[i];
    for (int j = 0; j < 2; j++) {
      delete InvPtRes[i][j];
      delete InvPtPull[i][j];
    }
  }

  // DilRes histos
  for (int i = 0; i <= MAX_LEVELS; i++) {
    for (int j = 0; j < 2; j++)
      delete DilMassComp[i][j];
    for (int j = 0; j < 6; j++)
      delete DilMassRes[i][j];
    delete DilPtRes[i];
  }

  for (unsigned int j = 0; j < 2; j++)
    delete MuPVsMuM[j];
}

void Zprime2muResolution::DrawResHistos(){
  const string str_part[3] = {
    "(Opp-sign Dil) Mu-", "(Opp-sign Dil) Mu+", "Opp-sign Dilepton"
  };

  TText t;
  t.SetTextFont(12);
  t.SetTextSize(.03);

  TCanvas *c1 = new TCanvas("c1", "", 0, 0, 500, 700);
  TPostScript *ps = new TPostScript(outputFile.c_str(), 111);

  const int NUM_PAGES = 70;
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
  GMRMuonHitsTot->SetNameTitle("", "Muon hits, GMR, total");
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
  GMRChi2dofTot->SetNameTitle("", "Chi2/d.o.f., GMR, total");
  GMRChi2dofTot->Add(GMRChi2dof[0], GMRChi2dof[1], 1., 1.);
  GMRChi2dofTot->Add(GMRChi2dofTot, GMRChi2dof[2], 1., 1.);
  pad[page]->cd(4);  GMRChi2dofTot->Draw();
  c1->Update();
  delete GMRChi2dofTot;

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
  effeta->SetNameTitle("", "L3 muon efficiency vs eta");
  effeta->Divide(EffVsEta[3], EffVsEta[0], 1., 1., "B");
  effeta->SetMinimum(0.50);  effeta->SetMaximum(1.02);
  effeta->Draw();  effeta->Draw("samehisto");
  // Off-line muon reconstruction efficiency as a function of true phi.
  pad[page]->cd(2);  gPad->SetGrid(1);
  effphi = (TH1F*)EffVsPhi[3]->Clone();
  effphi->SetNameTitle("", "L3 muon efficiency vs phi");
  effphi->Divide(EffVsPhi[3], EffVsPhi[0], 1., 1., "B");
  effphi->SetMinimum(0.50);  effphi->SetMaximum(1.02);
  effphi->Draw();  effphi->Draw("samehisto");
  // Off-line muon reconstruction efficiency as a function of true pT.
  pad[page]->cd(3);  gPad->SetGrid(1);
  effpt = (TH1F*)EffVsPt[3]->Clone();
  effpt->SetNameTitle("", "L3 muon efficiency vs pT");
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
  string histtitle;
  // Muon reconstruction efficiency (different reconstructors) as a
  // function of true eta. Only muons produced in Z' decays are used;
  // L1/HLT efficiency is not included.
  for (int irec = 1; irec <= MAX_LEVELS; irec++) {
    pad[page]->cd(irec);  gPad->SetGrid(1);
    effvseta[irec] = (TH1F*)EffVsEta[irec]->Clone();
    histtitle = str_level[irec] + " muon efficiency vs eta";
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
      histtitle += str_level[irec];
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
  pad[page]->Divide(3,3);
  TH1F *efftot[4], *eff2mu[4], *effdimu[4];
  string level;
  for (int i = 1; i < 4; i++) {
    if (i == 1)      level = "L3";
    else if (i == 2) level = "GMR";
    else if (i == 3) level = "TMR";

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
    // in acceptance and accepted by the L1/HLT triggers.
    pad[page]->cd(i+6);  gPad->SetGrid(1);
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
    delete effdimu[i];
  }

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
  for (int i=0; i<NUM_REC_LEVELS; i++) {
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = str_level[i] + " Values, all Muons";
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
  for (int i=0; i<NUM_REC_LEVELS; i++) {
    for (int j=0; j<2; j++) {
      ps->NewPage();
      c1->Clear();
      c1->cd(0);
      tit = str_part[j] + " " + str_level[i] + " Values";
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
  for (int i_rec = 0; i_rec < NUM_REC_LEVELS; i_rec++) {
    // Eta, Phi, Rap and momentum
    ps->NewPage();
    c1->Clear();
    c1->cd(0);
    tit = str_part[2] + " " + str_level[i_rec] + " Values";
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
  MuInvPtResGMRTot->SetNameTitle("", "1/pT res, GMR, total");
  MuInvPtResGMRTot->Add(MuInvPtRes[0][0], MuInvPtRes[0][1], 1., 1.);
  MuInvPtResGMRTot->Add(MuInvPtResGMRTot, MuInvPtRes[0][2], 1., 1.);
  pad[page]->cd(7);
  MuInvPtResGMRTot->Draw();  MuInvPtResGMRTot->Fit("gaus","Q");
  TH1F* MuInvPtResTMRTot = (TH1F*)MuInvPtRes[1][0]->Clone();
  MuInvPtResTMRTot->SetNameTitle("", "1/pT res, TMR, total");
  MuInvPtResTMRTot->Add(MuInvPtRes[1][0], MuInvPtRes[1][1], 1., 1.);
  MuInvPtResTMRTot->Add(MuInvPtResTMRTot, MuInvPtRes[1][2], 1., 1.);
  pad[page]->cd(8);
  MuInvPtResTMRTot->Draw();  MuInvPtResTMRTot->Fit("gaus","Q");
  c1->Update();
  delete MuInvPtResGMRTot;
  delete MuInvPtResTMRTot;

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
  int nbins;
  Stat_t f_bin, ent_bin, err_bin;
  TH1F *PResVsPMod[MAX_LEVELS+1];
  for (int i_rec = l1; i_rec <= MAX_LEVELS; i_rec++) {
    histtitle = "Sqrt(Var((" + str_level[i_rec] + "-Gen P)/Gen P)) vs Gen P";
    nbins = PResVsP[i_rec]->GetNbinsX();
    PResVsPMod[i_rec] = new TH1F("dp", histtitle.c_str(), nbins,
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
  MuInvPResGMRTot->SetNameTitle("", "1/P res, GMR, total");
  MuInvPResGMRTot->Add(MuInvPRes[0][0], MuInvPRes[0][1], 1., 1.);
  MuInvPResGMRTot->Add(MuInvPResGMRTot, MuInvPRes[0][2], 1., 1.);
  pad[page]->cd(7);
  MuInvPResGMRTot->Draw();  MuInvPResGMRTot->Fit("gaus","Q");
  TH1F* MuInvPResTMRTot = (TH1F*)MuInvPRes[1][0]->Clone();
  MuInvPResTMRTot->SetNameTitle("", "1/P res, TMR, total");
  MuInvPResTMRTot->Add(MuInvPRes[1][0], MuInvPRes[1][1], 1., 1.);
  MuInvPResTMRTot->Add(MuInvPResTMRTot, MuInvPRes[1][2], 1., 1.);
  pad[page]->cd(8);
  MuInvPResTMRTot->Draw();  MuInvPResTMRTot->Fit("gaus","Q");
  c1->Update();
  delete MuInvPResGMRTot;
  delete MuInvPResTMRTot;

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
  pad[page]->cd(5);  AllDilMassRes[2]->Draw();
  AllDilMassRes[2]->Fit("gaus","Q");
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
    MassResScat[i_rec] = new TH1F("", histname.str().c_str(), nbins,
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

    histtitle = "Sqrt(Var((" + str_level[i_rec+1] +
      "-Gen Mass)/Gen Mass)) vs Gen Mass";
    nbins = GenDilMassFrResScat[i_rec]->GetNbinsX();
    MassFrResScat[i_rec] = new TH1F("dm", histtitle.c_str(), nbins,
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
  pad[page]->cd(1);  TH1F *rat1_1 = (TH1F*)QResVsPt[0][0]->Clone();
  rat1_1->SetNameTitle("", "(wrong Q)/(right Q) vs L1 pT");
  rat1_1->Divide(QResVsPt[0][1], QResVsPt[0][0], 1., 1.);  rat1_1->Draw();
  pad[page]->cd(2);  TH1F *rat1_2 = (TH1F*)QResVsInvPt[0][0]->Clone();
  rat1_2->SetNameTitle("", "(wrong Q)/(right Q) vs L1 1/pT");
  rat1_2->Divide(QResVsInvPt[0][1], QResVsInvPt[0][0], 1., 1.); rat1_2->Draw();
  pad[page]->cd(3);  TH1F *rat1_3 = (TH1F*)QResVsPt[1][0]->Clone();
  rat1_3->SetNameTitle("", "(wrong Q)/(right Q) vs L2 pT");
  rat1_3->Divide(QResVsPt[1][1], QResVsPt[1][0], 1., 1.);  rat1_3->Draw();
  pad[page]->cd(4);  TH1F *rat1_4 = (TH1F*)QResVsInvPt[1][0]->Clone();
  rat1_4->SetNameTitle("", "(wrong Q)/(right Q) vs L2 1/pT");
  rat1_4->Divide(QResVsInvPt[1][1], QResVsInvPt[1][0], 1., 1.); rat1_4->Draw();
  pad[page]->cd(5);  TH1F *rat1_5 = (TH1F*)QResVsPt[2][0]->Clone();
  rat1_5->SetNameTitle("", "(wrong Q)/(right Q) vs L3 pT");
  rat1_5->Divide(QResVsPt[2][1], QResVsPt[2][0], 1., 1.);  rat1_5->Draw();
  pad[page]->cd(6);  TH1F *rat1_6 = (TH1F*)QResVsInvPt[2][0]->Clone();
  rat1_6->SetNameTitle("", "(wrong Q)/(right Q) vs L3 1/pT");
  rat1_6->Divide(QResVsInvPt[2][1], QResVsInvPt[2][0], 1., 1.); rat1_6->Draw();
  c1->Update();
  delete rat1_1;
  delete rat1_2;
  delete rat1_3;
  delete rat1_4;
  delete rat1_5;
  delete rat1_6;

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
  pad[page]->cd(1);  TH1F *rat2_1 = (TH1F*)QResVsP[0][0]->Clone();
  rat2_1->SetNameTitle("", "(wrong Q)/(right Q) vs L1 P");
  rat2_1->Divide(QResVsP[0][1], QResVsP[0][0], 1., 1.);  rat2_1->Draw();
  pad[page]->cd(2);  TH1F *rat2_2 = (TH1F*)QResVsInvP[0][0]->Clone();
  rat2_2->SetNameTitle("", "(wrong Q)/(right Q) vs L1 1/P");
  rat2_2->Divide(QResVsInvP[0][1], QResVsInvP[0][0], 1., 1.); rat2_2->Draw();
  pad[page]->cd(3);  TH1F *rat2_3 = (TH1F*)QResVsP[1][0]->Clone();
  rat2_3->SetNameTitle("", "(wrong Q)/(right Q) vs L2 P");
  rat2_3->Divide(QResVsP[1][1], QResVsP[1][0], 1., 1.);  rat2_3->Draw();
  pad[page]->cd(4);  TH1F *rat2_4 = (TH1F*)QResVsInvP[1][0]->Clone();
  rat2_4->SetNameTitle("", "(wrong Q)/(right Q) vs L2 1/P");
  rat2_4->Divide(QResVsInvP[1][1], QResVsInvP[1][0], 1., 1.); rat2_4->Draw();
  pad[page]->cd(5);  TH1F *rat2_5 = (TH1F*)QResVsP[2][0]->Clone();
  rat2_5->SetNameTitle("", "(wrong Q)/(right Q) vs L3 P");
  rat2_5->Divide(QResVsP[2][1], QResVsP[2][0], 1., 1.);
  rat2_5->SetMinimum(-0.1);  rat2_5->SetMaximum(0.5);  rat2_5->Draw();
  pad[page]->cd(6);  TH1F *rat2_6 = (TH1F*)QResVsInvP[2][0]->Clone();
  rat2_6->SetNameTitle("", "(wrong Q)/(right Q) vs L3 1/P");
  rat2_6->Divide(QResVsInvP[2][1], QResVsInvP[2][0], 1., 1.);
  rat2_6->SetMinimum(-0.1);  rat2_6->SetMaximum(0.5);  rat2_6->Draw();
  c1->Update();
  /*
  delete rat2_1;
  delete rat2_2;
  delete rat2_3;
  delete rat2_4;
  delete rat2_5;
  delete rat2_6;
  */

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
  delete ps;
  delete c1;
}

void Zprime2muResolution::analyze(const edm::Event& event, 
				  const edm::EventSetup& eSetup) {
  // delegate filling our muon vectors to the parent class
  Zprime2muAnalysis::analyze(event, eSetup);
  // fill resolution histos
  calcResolution();
}

void Zprime2muResolution::beginJob(const edm::EventSetup& eSetup) {
}

void Zprime2muResolution::endJob() {
  DrawResHistos();
}
