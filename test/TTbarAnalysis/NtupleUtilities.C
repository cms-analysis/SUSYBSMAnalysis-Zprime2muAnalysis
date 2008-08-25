#ifndef NtupleUtilities_C
#define NtupleUtilities_C
/*
 * =====================================================================================
 *
 *       Filename:  NtupleUtilities.C
 *
 *    Description:  some utilities for RootMaker
 *
 *        Version:  1.0
 *        Created:  05/19/2008 10:42:32 AM CEST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  IHEP, Beijing
 *
 * =====================================================================================
 */
#include "TMath.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TROOT.h"
#include <vector>
float DeltaPhi(float v1, float v2)
{ // Computes the correctly normalized phi difference
	// v1, v2 = phi of object 1 and 2
	float diff = fabs(v2 - v1);
	float corr = 2*acos(-1.) - diff;
	if (diff < acos(-1.)){ return diff;} else { return corr;}

	//  cout<<"enter DeltaPhi"<<endl;
}
//------------------------------------------------------------------------------

float GetDeltaR(float eta1, float eta2, float phi1, float phi2)
{ // Computes the DeltaR of two objects from their eta and phi values

	return sqrt( (eta1-eta2)*(eta1-eta2)
			+ DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
	//cout<<"enter getDeltaR"<<endl;
}

// the numbers' meaning is described in https://twiki.cern.ch/twiki/bin/view/CMS/CSA07ProcessId
int GetProcessId( int genEventProcID, Float_t genEventScale, Float_t genEventFilterEff, Float_t evtWeight, Float_t overallLumi = 100.)
{
	// ---------------------------------NOTE :::::::::::::::::------------------------------------
	// the numbers' meaning is described in https://twiki.cern.ch/twiki/bin/view/CMS/CSA07ProcessId
	// ---------------------------------::::::::::::::::------------------------------------
	// get process Id
	int procId = genEventProcID;
	// get generator event scale
	double ptHat = genEventScale;
	// get generated filter efficiency
	// chowder is -1,  not available for alpgen samples
	double filterEff = genEventFilterEff;
	// get csa07 weight
	double weight = evtWeight * 1000/overallLumi;
	// get the csa07 process id
	int csa07ProcId;
	// chowder processes
	if (procId == 1000)
		csa07ProcId = 0;  // W+0jet		 (weight ~ 5.15)
	else if (procId == 1001 && weight < 1.03)
		csa07ProcId = 1;  // W+1jet 0<pTW<100   (weight ~ 1.02)
	else if (procId == 1001 && weight > 1.03)
		csa07ProcId = 2;  // W+1jet 100<pTW<300 (weight ~ 1.04)
	else if (procId == 1002 && weight > 1.)
		csa07ProcId = 3;  // W+2jet 0<pTW<100   (weight ~ 1.07)
	else if (procId == 1002 && weight < 1.)
		csa07ProcId = 4;  // W+2jet 100<pTW<300 (weight ~ 0.78)
	else if (procId == 1003 && weight > 1.)
		csa07ProcId = 5;  // W+3jet 0<pTW<100   (weight ~ 1.67)
	else if (procId == 1003 && weight < 1.)
		csa07ProcId = 6;  // W+3jet 100<pTW<300 (weight ~ 0.91)
	else if (procId == 1004 && weight > 0.96)
		csa07ProcId = 7;  // W+4jet 0<pTW<100   (weight ~ 0.98)
	else if (procId == 1004 && weight < 0.96)
		csa07ProcId = 8;  // W+4jet 100<pTW<300 (weight ~ 0.95)
	else if (procId == 1005 && weight > 1.)
		csa07ProcId = 9;  // W+5jet 0<pTW<100   (weight ~ 1.35)
	else if (procId == 1005 && weight < 1.)
		csa07ProcId = 10; // W+5jet 100<pTW<300 (weight ~ 0.90)
	else if (procId == 2000)
		csa07ProcId = 11; // Z+0jet		  (weight ~ 1.38)
	else if (procId == 2001 && weight > 0.9)
		csa07ProcId = 12; // Z+1jet 0<pTZ<100   (weight ~ 0.98)
	else if (procId == 2001 && weight < 0.9)
		csa07ProcId = 13; // Z+1jet 100<pTZ<300 (weight ~ 0.83)
	else if (procId == 2002 && weight > 0.9)
		csa07ProcId = 14; // Z+2jet 0<pTZ<100   (weight ~ 0.93)
	else if (procId == 2002 && weight < 0.9)
		csa07ProcId = 15; // Z+2jet 100<pTZ<300 (weight ~ 0.80)
	else if (procId == 2003 && weight > 0.9)
		csa07ProcId = 16; // Z+3jet 0<pTZ<100   (weight ~ 0.94)
	else if (procId == 2003 && weight < 0.9)
		csa07ProcId = 17; // Z+3jet 100<pTZ<300 (weight ~ 0.53)
	else if (procId == 2004 && weight < 0.5)
		csa07ProcId = 18; // Z+4jet 0<pTZ<100   (weight ~ 0.42)
	else if (procId == 2004 && weight > 0.5)
		csa07ProcId = 19; // Z+4jet 100<pTZ<300 (weight ~ 0.63)
	else if (procId == 2005 && weight < 0.8)
		csa07ProcId = 20; // Z+5jet 0<pTZ<100   (weight ~ 0.72)
	else if (procId == 2005 && weight > 0.8)
		csa07ProcId = 21; // Z+5jet 100<pTZ<300 (weight ~ 0.85)
	else if (procId == 3000)
		csa07ProcId = 22; // tt+0jet
	else if (procId == 3001)
		csa07ProcId = 23; // tt+1jet
	else if (procId == 3002)
		csa07ProcId = 24; // tt+2jet
	else if (procId == 3003)
		csa07ProcId = 25; // tt+3jet
	else if (procId == 3004)
		csa07ProcId = 26; // tt+4jet
	// gumbo processes
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(15 < ptHat && ptHat < 20) && filterEff == 1.)
		csa07ProcId = 28; // QCD_Pt_15_20
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(20 < ptHat && ptHat < 30) && filterEff == 1.)
		csa07ProcId = 29; // QCD_Pt_20_30
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(30 < ptHat && ptHat < 50) && filterEff == 1.)
		csa07ProcId = 30; // QCD_Pt_30_50
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(50 < ptHat && ptHat < 80) && filterEff == 1.)
		csa07ProcId = 31; // QCD_Pt_50_80
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(80 < ptHat && ptHat < 120) && filterEff == 1.)
		csa07ProcId = 32; // QCD_Pt_80_120
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(120 < ptHat && ptHat < 170) && filterEff == 1.)
		csa07ProcId = 33; // QCD_Pt_120_170
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(170 < ptHat && ptHat < 230) && filterEff == 1.)
		csa07ProcId = 34; // QCD_Pt_170_230
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(230 < ptHat && ptHat < 300) && filterEff == 1.)
		csa07ProcId = 35; // QCD_Pt_230_300
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(300 < ptHat && ptHat < 380) && filterEff == 1.)
		csa07ProcId = 36; // QCD_Pt_300_380
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(380 < ptHat && ptHat < 470) && filterEff == 1.)
		csa07ProcId = 37; // QCD_Pt_380_470
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(470 < ptHat && ptHat < 600) && filterEff == 1.)
		csa07ProcId = 38; // QCD_Pt_470_600
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(600 < ptHat && ptHat < 800) && filterEff == 1.)
		csa07ProcId = 39; // QCD_Pt_600_800
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(800 < ptHat && ptHat < 1000) && filterEff == 1.)
		csa07ProcId = 40; // QCD_Pt_800_1000
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(1000 < ptHat && ptHat < 1400) && filterEff == 1.)
		csa07ProcId = 41; // QCD_Pt_1000_1400
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(1400 < ptHat && ptHat < 1800) && filterEff == 1.)
		csa07ProcId = 42; // QCD_Pt_1400_1800
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(1800 < ptHat && ptHat < 2200) && filterEff == 1.)
		csa07ProcId = 43; // QCD_Pt_1800_2200
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(2200 < ptHat && ptHat < 2600) && filterEff == 1.)
		csa07ProcId = 44; // QCD_Pt_2200_2600
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(2600 < ptHat && ptHat < 3000) && filterEff == 1.)
		csa07ProcId = 45; // QCD_Pt_2600_3000
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(3000 < ptHat && ptHat < 3500) && filterEff == 1.)
		csa07ProcId = 46; // QCD_Pt_3000_3500
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(3500 < ptHat) && filterEff == 1.)
		csa07ProcId = 47; // QCD_Pt_3500_inf
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68 || procId == 92 || procId == 93 || procId == 94 || procId == 95) &&
			(filterEff == 1.))
		csa07ProcId = 27; // min bias - MUST COME AFTER ALL QCD
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(ptHat < 15))
		csa07ProcId = 48; // PhotonJets_Pt_0_15
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(15 < ptHat && ptHat < 20))
		csa07ProcId = 49; // PhotonJets_Pt_15_20
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(20 < ptHat && ptHat < 30))
		csa07ProcId = 50; // PhotonJets_Pt_20_30
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(30 < ptHat && ptHat < 50))
		csa07ProcId = 51; // PhotonJets_Pt_30_50
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(50 < ptHat && ptHat < 80))
		csa07ProcId = 52; // PhotonJets_Pt_50_80
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(80 < ptHat && ptHat < 120))
		csa07ProcId = 53; // PhotonJets_Pt_80_120
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(120 < ptHat && ptHat < 170))
		csa07ProcId = 54; // PhotonJets_Pt_120_170
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(170 < ptHat && ptHat < 300))
		csa07ProcId = 55; // PhotonJets_Pt_170_300
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(300 < ptHat && ptHat < 500))
		csa07ProcId = 56; // PhotonJets_Pt_300_500
	else if ((procId == 14 || procId == 18 || procId == 29) &&
			(500 < ptHat && ptHat < 7000))
		csa07ProcId = 57; // PhotonJets_Pt_500_7000
	else if (procId == 102 || procId == 123 || procId == 124)
		csa07ProcId = 58; // higgs signal
	else if (procId == 141)
		csa07ProcId = 59; // Zprime signal
	// stew processes
	else if (filterEff == 0.00013)
		csa07ProcId = 60; // B(bar)->JPsi
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68) &&
			(filterEff == 1. || filterEff == 0.964))
		csa07ProcId = 61; // QCD_Pt_0_15
	else if ((460 < procId && procId < 480) &&
			(0 < ptHat && ptHat < 20))
		csa07ProcId = 62; // Bottomonium Pt_0_20
	else if ((460 < procId && procId < 480) &&
			(20 < ptHat))
		csa07ProcId = 63; // Bottomonium Pt_20_inf
	else if ((420 < procId && procId < 440) &&
			(0 < ptHat && ptHat < 20))
		csa07ProcId = 64; // Charmonium Pt_0_20
	else if ((420 < procId && procId < 440) &&
			(20 < ptHat))
		csa07ProcId = 65; // Charmonium Pt_20_inf
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68 || procId == 95) &&
			(filterEff == 0.00019))
		csa07ProcId = 66; // bbe Pt_5_50
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68 || procId == 95) &&
			(filterEff == 0.0068))
		csa07ProcId = 67; // bbe Pt_50_170
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68 || procId == 95) &&
			(filterEff == 0.0195))
		csa07ProcId = 68; // bbe Pt_170_up
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68 || procId == 95) &&
			(filterEff == 0.0097))
		csa07ProcId = 69; // ppEleX
	else if ((procId == 11 || procId == 12 || procId == 13 || procId == 28 || procId == 53 || procId == 68 || procId == 95) &&
			(filterEff == 0.0008))
		csa07ProcId = 70; // ppMuX
	else {
		csa07ProcId = -1; // unknown process
		cout<<"---------------Unknown CSA07 process:-----------"
			<< "OUCH! Unknown CSA07 process with: \n"
			<< "  ID	    : " << procId    << "\n"
			<< "  scale     : " << ptHat     << "\n"
			<< "  filter eff: " << filterEff << "\n"
			<< "  weight    : " << weight    << "\n";
	}
	return csa07ProcId;
}

// - --  deprecated method
int GetProcess(int pytSize, vector<int> *pytStatus, vector<int> *pytId, vector<int> *pytNMother) 
{
	// NOTE: this function is only supported to be used for Choder soup 
	// FIXME there is an issue for minbias 
	// Gumbo
	// 1: dijets, 2:minbias, 5: higgs, 6:Zprime, 7:photon+jets

	// Chowder
	// 20: w+0jet, 21:w+1jet, 22:w+2jets, 23:w+3jets, 24:w+4jets ....
	// 30: z+0jet, 31:z+1jet, 32:z+2jets, 33:z+3jets, 34:z+4jets ....
	// 40: tt+0jet, 41:tt+1jet, 42:tt+2jets, 43:tt+3jets, 44:tt+4jets ....

	// Stew
	// ...................................................................
	// kk runing on stew,
	// because a bug which was in NtupleDumper.cc not exist in the version with which kk running
	// kk, you can use the "evtProcessId" instead of this function
	// the numbers' meaning is described in https://twiki.cern.ch/twiki/bin/view/CMS/CSA07ProcessId

	//	// di-bosons, ww, wz, zz, +jets , tw-----------supposed to be as signal, not soup, please keep in mind 

	// 99 : unknown
	int PrimaryId = 0;// 2:w, 3:z, 4:tt, 5:higgs 6:zprime, 7:photon

	bool isDijetsEvent = false;
	int numPartonJets = 0;
	int numPartonWithNMother2 = 0;
	for(int i = 0; i<pytSize; i++)
	{
		if((*pytStatus)[i]!=3) continue;
		int id = abs ((*pytId)[i]);
		int nmother = (*pytNMother)[i];

		if( nmother==2 && id==23 ) PrimaryId = 3;
		if( nmother==2 && id==24 ) PrimaryId = 2;
		if( nmother==2 && id==6  ) PrimaryId = 4;
		if( nmother==2 && id==25 ) PrimaryId = 5;
		if( nmother==2 && id==32 ) PrimaryId = 6;
		if( nmother==2 && id==22 ) PrimaryId = 7;

		if(i==6) 
		{
			if ( id==1 || id==2 || id==3 || id==4 || id==5 || id==21 )
			{
				if((*pytNMother)[i]==2 ) numPartonJets++;
			}
		}
		if(i==7) 
		{
			if ( id==1 || id==2 || id==3 || id==4 || id==5 || id==21 )
			{
				if((*pytNMother)[i]==2) numPartonJets++;
			}
		}
		if( ( (*pytStatus)[i]==3 ) && ( (*pytNMother)[i]==2 ) )
			numPartonWithNMother2 ++;
	}
	if( numPartonJets == 2 && numPartonWithNMother2 ==2 ) isDijetsEvent = true;
	if(isDijetsEvent == true) return 1 ;
	if(PrimaryId == 2) return PrimaryId*10 + (numPartonWithNMother2-1);
	else if(PrimaryId == 3) return PrimaryId*10 + (numPartonWithNMother2-1);
	else if(PrimaryId == 4) return PrimaryId*10 + (numPartonWithNMother2-2);
	else if(PrimaryId == 5) return 5;
	else if(PrimaryId == 6) return 6;
	else if(PrimaryId == 7) return 7;
	return 99;
}
void printPYTInfo(
//printPYTInfo(pytSize, pytPt, pytEta, pytPhi, pytP, pytPx, pytPy, pytPz, pytMass, pytIndex, pytNMother, pytMotherIndex, pytId, pytStatus);
		Int_t           pytSize,
		//		vector<float>   *pytEnergy,
		vector<float>   *pytPt,
		vector<float>   *pytEta,
		vector<float>   *pytPhi,
		vector<float>   *pytP,
		vector<float>   *pytPx,
		vector<float>   *pytPy,
		vector<float>   *pytPz,
		//		vector<float>   *pytTheta,
		//		vector<float>   *pytVx,
		//		vector<float>   *pytVy,
		//		vector<float>   *pytVz,
		//		vector<float>   *pytCharge,
		vector<float>   *pytMass,
		vector<int>     *pytIndex,
		vector<int>     *pytNMother,
		vector<int>     *pytMotherIndex,
		vector<int>     *pytId,
		vector<int>     *pytStatus
		)
{
	//printf("======> EVENT:%lld <======\n", entry);

	printf(" pytSize = %3d\n",pytSize);

	printf(" idx  | pdgID |Stat| Moth |nMo|    pt       eta     phi   |     px         py         pz        m        p       |\n");

	for(int i=0; i<pytSize; i++){
		printf(" %4d | %5d | %2d | %4d | %2d| %7.3f %10.3f %6.3f | %10.3f %10.3f %10.3f %8.3f %10.3f |\n",
				(*pytIndex)[i],
				(*pytId)[i],
				(*pytStatus)[i],
				(*pytMotherIndex)[i],
				(*pytNMother)[i],
				(*pytPt)[i],
				(*pytEta)[i],
				(*pytPhi)[i],
				(*pytPx)[i],
				(*pytPy)[i],
				(*pytPz)[i],
				(*pytMass)[i],
				(*pytP)[i]
			  );
	}
}

void printGMRInfo(
		Int_t           gmrSize,
		vector<float>   *gmrPt,
		vector<float>   *gmrEta,
		vector<float>   *gmrPhi,
		vector<float>   *gmrP,
		vector<float>   *gmrPx,
		vector<float>   *gmrPy,
		vector<float>   *gmrPz,
		//		vector<float>   *gmrTheta,
		vector<float>   *gmrVx,
		vector<float>   *gmrVy,
		vector<float>   *gmrVz,
		vector<float>   *gmrCharge,
		vector<float>   *gmrNDoF,
		//		vector<float>   *gmrChi2,
		vector<float>   *gmrChi2Norm,
		vector<float>   *gmrDXY,
		//		vector<float>   *gmrDTheta,
		//		vector<float>   *gmrDPt,
		//		vector<float>   *gmrDEta,
		//		vector<float>   *gmrDPhi,
		//		vector<float>   *gmrDDXY,
		vector<float>   *gmrIso03nTracks,
		vector<float>   *gmrIso03sumPt
		)
{
	printf(" gmrSize = %3d\n",gmrSize);
	printf("    pt          eta     phi   |     px         py         pz         p      |  Ndof NormChi2 IsoNtracks IsoSumPt Charge |\n");

	for(Int_t i=0; i<gmrSize; i++){
		printf(" %10.3f %10.3f %6.3f | %10.3f %10.3f %10.3f %10.3f | %5.1f %8.2f %10.1f %8.3f %6.1f  |\n",
				(*gmrPt)[i],
				(*gmrEta)[i],
				(*gmrPhi)[i],
				(*gmrPx)[i],
				(*gmrPy)[i],
				(*gmrPz)[i],
				(*gmrP)[i],
				(*gmrNDoF)[i],
				(*gmrChi2Norm)[i],
				(*gmrIso03nTracks)[i],
				(*gmrIso03sumPt)[i],
				(*gmrCharge)[i]
			  );
		//      cout<<"sumPt "<<(*gmrIso03sumPt)[i]<<endl;
	}
	printf("     Vx         Vy         Vz     |     Dxy        Dz         D0     | DisOf1stHitToVtx |\n");

	for(Int_t i=0; i<gmrSize; i++){

		Float_t gmrDist = sqrt( ((*gmrVx))[i]*((*gmrVx)[i]) + ((*gmrVy)[i])*((*gmrVy)[i]) + ((*gmrVz)[i])*((*gmrVz)[i]) );

		printf(" %10.3f %10.3f %10.3f | %10.3f %10.3f %10.3f | %15.3f |\n",
				(*gmrVy)[i],
				(*gmrVx)[i],
				(*gmrVz)[i],
				(*gmrDXY)[i],
				(*gmrDXY)[i],
				(*gmrDXY)[i],
				gmrDist
			  );
		//      cout<<"sumPt "<<(*gmrIso03sumPt)[i]<<endl;
		//,
	}
}
void printMETInfo(
		vector<float> *metEt
		)
{
	printf(" MET et = %8.3f\n", (*metEt)[0]);
}
void printJETInfo(
		Int_t           jetSize,
		vector<float>   *jetPt,
		vector<float>   *jetEta,
		vector<float>   *jetPhi,
		vector<float>   *jetP,
		vector<float>   *jetPx,
		vector<float>   *jetPy,
		vector<float>   *jetPz
		//		vector<float>   *jetTheta,
		//		vector<float>   *jetEnergy,
		//		vector<float>   *jetEt
		)
{
	printf(" jetSize = %3d\n",jetSize);
	printf(" idx  |    pt          eta     phi   |     px         py         pz         p      |\n");

	for(Int_t i=0; i<jetSize; i++){
		printf(" %4d | %10.3f %10.3f %6.3f | %10.3f %10.3f %10.3f %10.3f |\n",
				i,
				(*jetPt)[i],
				(*jetEta)[i],
				(*jetPhi)[i],
				(*jetPx)[i],
				(*jetPy)[i],
				(*jetPz)[i],
				(*jetP)[i]
			  );
	}

}
//-- METHOD -- gluon decay mode
int GetGluonDecayMode(int gluonIndex, int pytSize, vector<int>* pytIndex, vector<int>* pytStatus, vector<int>* pytId, vector<int>* pytMotherIndex )
{
	//- - - - 0-unknown, 1-gluon->c, 2-gluon->b, 3-gluon->b and c, 4-gluon-> !b and !c
	int gluonDecayMode = 0; 
	bool boolG2B=false;
	bool boolG2C=false;
	for(int iPYT=gluonIndex+1; iPYT<pytSize; iPYT++)
	{
		if((*pytStatus)[iPYT]!=2) continue;
		if(abs((*pytId)[iPYT])==4 && (*pytMotherIndex)[iPYT]==(*pytIndex)[gluonIndex]) boolG2C=true;
		if(abs((*pytId)[iPYT])==5 && (*pytMotherIndex)[iPYT]==(*pytIndex)[gluonIndex]) boolG2B=true;
		/*
		   if((*pytId)[iPYT]==21) gluonDecayMode = GetGluonDecayMode(iPYT, pytSize, pytIndex, pytStatus, pytId, pytMotherIndex);  
		//cout<<"sub gluonDecayMode = "<< gluonDecayMode<<endl;
		if(gluonDecayMode==1 || gluonDecayMode==3) boolG2C = true;
		if(gluonDecayMode==2 || gluonDecayMode==3) boolG2B = true;
		 */
	}
	if (boolG2C && !boolG2B) {gluonDecayMode=1;}
	else if (!boolG2C && boolG2B){gluonDecayMode=2;}
	else if (boolG2C && boolG2B) {gluonDecayMode=3;}
	else{gluonDecayMode=4;}
	return gluonDecayMode;
}
// --METHOD-- Determine if the true muon is prompt 
// NOTE: Z->tau tau, W->tau nu
// -- FIXME  b->c with w, w->mu nu , this mu is not prompt mu
bool IsItPromptMuon(int muonIndex, vector<int>* pytIndex, vector<int>*pytStatus, vector<int>* pytMotherIndex, vector<int>* pytMotherPdgId)
{
	//--- here muonIndex is the ordered number in stored pytParticles
	if(muonIndex > (int)pytIndex->size() || muonIndex < 0 ) cout<<"Muon Index out of range"<<endl;
	bool isprompt = false;
	for(int i = 2; i< muonIndex; i++)
	{
		if((*pytStatus)[muonIndex] >= 2 ) return true;
		else 
		{
			if( ( abs((*pytMotherPdgId)[muonIndex]) == 24 || (*pytMotherPdgId)[muonIndex]==23 ) ) return true;  
			else 
			{
				int muonMotherOriginalIndex = (*pytMotherIndex)[muonIndex];
				for(int j=2; j<muonIndex ; j++)	
				{
					if(muonMotherOriginalIndex == (*pytIndex)[j]) 
					{
						if (abs((*pytMotherPdgId)[j]) ==24 || (*pytMotherPdgId)[j]==23) return true;
					}
				}
			}
		}
	}
	return isprompt;	
}
bool IsItPromptMuon(int muonIndex, vector<Float_t>* pytVx, vector<Float_t>* pytVy)
{
	if(muonIndex > (int)pytVx->size() || muonIndex < 0 ) cout<<"Muon Index out of range"<<endl;
	bool isprompt = false;
	float x = 	(*pytVx)[muonIndex];
	float y =  (*pytVy)[muonIndex];
	float dxy = sqrt(x*x+y*y);
	//- - hard coded value -- FIXME, need evaluate from the impact parameter distribution between prompt and non-prompt muons  
	if( dxy < 0.005 ) isprompt = true; // in cm
	return isprompt;
}
// --METHOD-- look for a true muon to match the gmr muon, return the index of stored pytParticles 
int FindTrueMuonMatchedToThisGMRMuon(int muonIndex, vector<Float_t>* gmrEta, vector<Float_t>* gmrPhi, 
		int pytSize, vector<int>*pytId, vector<Float_t>*pytEta, vector<Float_t>*pytPhi) 
{
	int trueMuonIndex = -1; // -1 : unmatched
	Float_t dr = 999;
	for(int i=0; i<pytSize; i++)
	{
		if(abs((*pytId)[i])!=13) continue;
		Float_t tmpdr = GetDeltaR((*gmrEta)[muonIndex], (*pytEta)[i], (*gmrPhi)[muonIndex], (*pytPhi)[i]);
		if (tmpdr < dr ) { dr = tmpdr; trueMuonIndex = i ; }
	}
	// FIXME here is a hard-coded cut value for matching
	if( dr < 0.3 ) return trueMuonIndex;
	return -1;
}
// --METHOD-- look for a reco jet to match the gmr muon, return the index of stored jets
int FindJetMatchedToThisGMRMuon(int muonIndex, vector<Float_t>* gmrEta, vector<Float_t>* gmrPhi, 
		int jetSize, vector<Float_t>*jetEta, vector<Float_t>*jetPhi, vector<Float_t>*jetPt, double minJetPt) 
{
	int jetIndex = -1; // -1 : unmatched
	Float_t dr = 999;
	for(int i=0; i<jetSize; i++)
	{
		float jetpt = (*jetPt)[i];
		if( jetpt < minJetPt ) continue;
		Float_t tmpdr = GetDeltaR((*gmrEta)[muonIndex], (*jetEta)[i], (*gmrPhi)[muonIndex], (*jetPhi)[i]);
		if (tmpdr < dr ) { dr = tmpdr; jetIndex = i ; }
	}
	// FIXME here is a hard-coded cut value for matching
	// if( dr < 0.5 ) return jetIndex;
	return jetIndex;
}

int FindNearestPartonToThisJet(int jetIndex, vector<Float_t>* jetEta, vector<Float_t>* jetPhi, 
		int pytSize, vector<int>*pytId, vector<Float_t>*pytEta, vector<Float_t>*pytPhi) 
{
	int partonIndex = -1; // -1 : unmatched
	Float_t dr = 999;
	for(int i=0; i<pytSize; i++)
	{
		int id = abs((*pytId)[i]);
		if( id!= 1 && id!=2 && id!=3 && id!=4 && id!=5 && id!=21 ) continue;
		Float_t tmpdr = GetDeltaR((*jetEta)[jetIndex], (*pytEta)[i], (*jetPhi)[jetIndex], (*pytPhi)[i]);
		if (tmpdr < dr ) { dr = tmpdr; partonIndex = i ; }
	}
	// FIXME here is a hard-coded cut value for matching
	if( dr < 0.3 ) return partonIndex;
	return -1;
}
// --METHOD-- look for a parton to match the jet, return the index of stored pytParticles 
//  http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc?revision=1.9&view=markup
//
// Algorithmic Definition: 
// Output: define one associatedParton
// Loop on all particle.
// A particle is a parton if its daughter is a string(code=92) or a cluster(code=93) 
// If (parton is within the cone defined by coneSizeToAssociate) then:
//           if (parton is a b)                                   then associatedParton is the b
//      else if (associatedParton =! b and parton is a c)         then associatedParton is the c
//      else if (associatedParton =! b and associatedParton =! c) then associatedParton is the one with the highest pT
// associatedParton can be -1 --> no partons in the cone
// True Flavour of the jet --> flavour of the associatedParton
//
// ToDo: if more than one b(c) in the cone --> the selected parton is not always the one with highest pT
//
int FindPartonMatchedToThisJet(int jetIndex, vector<Float_t>* jetEta, vector<Float_t>* jetPhi, 
		int pytSize, vector<int>*pytId, vector<Float_t>*pytEta, vector<Float_t>*pytPhi, vector<Float_t>* pytPt, 
		vector<int> *pytNDaughter, vector<int>*pytDaughterPdgId) 
{

	int tempParticle = -1;
	int tempPartonHighestPt = -1;
	int tempNearest = -1;
	float maxPt = 0;
	float minDr = 1000;
	for( int m = 0; m != pytSize; ++ m ) {
		if( (*pytNDaughter)[m] > 0  && ( (*pytDaughterPdgId)[m] == 91 || (*pytDaughterPdgId)[m] == 92 ) ) {
			double dist =  GetDeltaR((*jetEta)[jetIndex], (*pytEta)[m], (*jetPhi)[jetIndex], (*pytPhi)[m]);
			// -- hard coded value
			if( dist <= 0.3 ) {
				if( dist < minDr ) {
					minDr = dist;
					tempNearest = m;
				}
				if( tempParticle == -1 && ( abs( (*pytId)[m] ) == 4 )  ) tempParticle = m;
				if(                         abs( (*pytId)[m] ) == 5    ) tempParticle = m;
				if( (*pytPt)[m] > maxPt ) {
					maxPt = (*pytPt)[m];
					tempPartonHighestPt = m;
				}
			}
		}
	}
	//theHeaviest = tempParticle;
	//theHardest  = tempPartonHighestPt;
	//theNearest = tempNearest;
	if ( tempParticle == -1 ) tempParticle = tempPartonHighestPt;
	// -- if still no parton matched, then use my own simple matcher
	if ( tempParticle == -1 ) tempParticle = FindNearestPartonToThisJet(jetIndex, jetEta, jetPhi, pytSize, pytId, pytEta, pytPhi); 
	return tempParticle;
}
// --METHOD-- look for a parton to match the gmr muon, return the index of stored pytParticles 
int FindPartonMatchedToThisGMRMuon(int muonIndex, vector<Float_t>* gmrEta, vector<Float_t>* gmrPhi, 
		int pytSize, vector<int>*pytId, vector<Float_t>*pytEta, vector<Float_t>*pytPhi) 
{
	int partonIndex = -1; // -1 : unmatched
	Float_t dr = 999;
	for(int i=0; i<pytSize; i++)
	{
		int id = abs((*pytId)[i]);
		if( id!= 1 && id!=2 && id!=3 && id!=4 && id!=5 && id!=21 ) continue;
		Float_t tmpdr = GetDeltaR((*gmrEta)[muonIndex], (*pytEta)[i], (*gmrPhi)[muonIndex], (*pytPhi)[i]);
		if (tmpdr < dr ) { dr = tmpdr; partonIndex = i ; }
	}
	// FIXME here is a hard-coded cut value for matching
	if( dr < 0.3 ) return partonIndex;
	return -1;
}
TString GetPtHat(Float_t pthat) // pthat = genEventScale
{
	if (pthat > 0 && pthat < 15)  { return "0~15";} 
	if (pthat > 15 && pthat < 20) { return "15~20";}
	if (pthat > 20 && pthat < 30) { return "20~30";}
	if (pthat > 30 && pthat < 50) { return "30~50";}
	if (pthat > 50 && pthat < 80) { return "50~80";}
	if (pthat > 80 && pthat < 120){ return "80~120";}
	if (pthat > 120 && pthat < 170) { return "120~170";}
	if (pthat > 170 && pthat < 230) { return "170~230";}
	if (pthat > 230 && pthat < 300) { return "230~300";}
	if (pthat > 300 && pthat < 380) { return "300~380";}
	if (pthat > 380 && pthat < 470) { return "380~470";}
	if (pthat > 470 && pthat < 600) { return "470~600";}
	if (pthat > 600 && pthat < 800) { return "600~800";}
	if (pthat > 800 && pthat < 1000){ return "800~1000";}
	if (pthat > 1000 && pthat < 1400) { return "1000~1400";}
	if (pthat > 1400 && pthat < 1800) { return "1400~1800";}
	if (pthat > 1800 && pthat < 2200) { return "1800~2200";}
	if (pthat > 2200 && pthat < 2600) { return "2200~2600";}
	if (pthat > 2600 && pthat < 3000) { return "2600~3000";}
	if (pthat > 3000 && pthat < 3500) { return "3000~3500";}
	if (pthat > 3500) { return "3500~inf";}
	return "unknown";
}
#endif
