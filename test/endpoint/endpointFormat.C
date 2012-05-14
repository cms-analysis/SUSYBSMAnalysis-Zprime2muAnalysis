#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

typedef struct {

	float k;
	float eta;
	float phi;
	float ptErr;

} _OutType;

void endpointFormat() {

	TFile infile1("tmp1.root","open");
	TFile infile2("tmp2.root","open");

	TTree* intree1 = (TTree*)infile1.Get("tree");
	TTree* intree2 = (TTree*)infile2.Get("tree");
		

	_OutType outTrk;	


	TFile outfile("endpoint.root","recreate");
	TTree* outTree1 = new TTree("muonsTrk","muonsTrk");
	TTree* outTree2 = new TTree("muonsOpt","muonsOpt");

	
	float k,eta,phi;
	intree1->SetBranchAddress("k",&k);
	intree1->SetBranchAddress("eta",&eta);
	intree1->SetBranchAddress("phi",&phi);

	outTree1->Branch("trk1",&outTrk,"k/F:eta:phi:ptErr");

	int numEntries = intree1->GetEntries();
	std::cout<<"num entries: "<<numEntries<<std::endl;

	outTrk.ptErr=1.;

	for (int jEntry =0; jEntry < numEntries; ++jEntry) {
		intree1->GetEntry(jEntry);
		outTrk.k=k;
		outTrk.eta=eta;
		outTrk.phi=phi;

		outTree1->Fill();
	}

	intree2->SetBranchAddress("k",&k);
	intree2->SetBranchAddress("eta",&eta);
	intree2->SetBranchAddress("phi",&phi);
	outTree2->Branch("trk1",&outTrk,"k/F:eta:phi:ptErr");

	numEntries = intree1->GetEntries();
	std::cout<<"num entries: "<<numEntries<<std::endl;

	outTrk.ptErr=1.;

	for (int jEntry =0; jEntry < numEntries; ++jEntry) {
		intree2->GetEntry(jEntry);
		outTrk.k=k;
		outTrk.eta=eta;
		outTrk.phi=phi;

		outTree2->Fill();
	}

	outTree1->Write();
	outTree2->Write();
	


/*
	_outTree->Branch("trk1",&_outTrk1, "k/F:eta:phi:ptErr");
	_outTree->Branch("trk2",&_outTrk2, "k/F:eta:phi:ptErr");
*/

}


