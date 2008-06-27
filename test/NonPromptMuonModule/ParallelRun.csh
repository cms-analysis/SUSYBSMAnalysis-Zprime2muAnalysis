#! /bin/csh

set THREAD=0
set TOTAL_ENTRIES=50000
set STARTING_ENTRY=0
while ($THREAD < 2)
cat >! tmp_$THREAD.C << EOF
{
	gInterpreter->LoadMacro("NtupleBase.C+");
	gInterpreter->LoadMacro("FakeMuonParameterizer.C+");
	if(true)
	{
		FakeMuonParameterizer *parameterizer = new FakeMuonParameterizer();
		parameterizer->Initialize("");
		parameterizer->Loop("DijetsParameterizing_$THREAD.root", $STARTING_ENTRY, $TOTAL_ENTRIES);
		delete parameterizer;
	}
}
EOF
nohup root -l -q tmp_$THREAD.C >! log_$THREAD &
@ THREAD = ($THREAD + 1)
@ STARTING_ENTRY = ($THREAD * $TOTAL_ENTRIES)
sleep 2
end


set size=`ll DijetsParameterizing_1.root | awk '{print $5}'`
while ( $size <= 1000 )
echo "Not done"
sleep 1
set size=`ll DijetsParameterizing_1.root | awk '{print $5}'`
end

#hadd DijetsParameterizing.root DijetsParameterizing_*.root
#
#cat > tmp_$THREAD.C << EOF
#{
#	gInterpreter->LoadMacro("NtupleBase.C+");
#	gInterpreter->LoadMacro("FakeMuonGenerator.cc+");
#	gInterpreter->LoadMacro("FakeMuonAnalyzer.C+");
#	if(true)
#	{
#		FakeMuonAnalyzer *analyzer = new FakeMuonAnalyzer();
#		analyzer->Initialize("");
#		analyzer->Loop("DijetsPredObsv_$THREAD.root", $STARTING_ENTRY, $TOTAL_ENTRIES) ;
#		delete analyzer;
#	}
#}
#EOF
