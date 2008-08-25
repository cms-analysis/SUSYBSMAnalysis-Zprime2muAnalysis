{
gInterpreter->LoadMacro("NtupleBase.C+");
gInterpreter->LoadMacro("TTbarAnalysis.C+");

TString s_input           = "/data/1b/mschen/MuonChowder1.root";
//TString s_input           = "/data/1b/mschen/CSA07Zmumu100pb07.root";
		//s_input = "";
//TString s_output          = "rootfiles/Zmumu.root"; // default is test.root
TString s_output          = "rootfiles/Chowder.root"; // default is test.root
Long64_t entryToStart     =      0; // default is         0
Long64_t entriesToAnalyze =      -1;// default is        -1

if (true)
{
  TTbarAnalysis *analyzer = new TTbarAnalysis();
  analyzer->Initialize(s_input);
  analyzer->Loop(s_output,entryToStart,entriesToAnalyze);
  delete analyzer;
}

//gInterpreter->ExecuteMacro("plotJetPerformance.C\(s_output\)");
}
