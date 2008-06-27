{
	gInterpreter->LoadMacro("NtupleBase.C+");
	gInterpreter->LoadMacro("FakeMuonGenerator.cc+");
	gInterpreter->LoadMacro("FakeMuonAnalyzer.C+");
	gInterpreter->LoadMacro("FakeMuonParameterizer.C+");

	if(false)
	{
		FakeMuonParameterizer *parameterizer = new FakeMuonParameterizer();
		parameterizer->Initialize("");
		parameterizer->Loop();
			delete parameterizer;
	}
	if(false)
	{
		FakeMuonAnalyzer *analyzer = new FakeMuonAnalyzer();
		analyzer->Initialize("");
		analyzer->Loop() ;
		delete analyzer;
	}
}
