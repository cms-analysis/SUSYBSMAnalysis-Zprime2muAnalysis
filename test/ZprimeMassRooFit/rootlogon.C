{
// Loading FloridaStyle
gInterpreter->ExecuteMacro("utils/FloridaStyle.C");
gInterpreter->LoadMacro("utils/Utils.C+");

// Loading RooFit libraries
gSystem->Load("libRooFit");
using namespace RooFit;

gInterpreter->LoadMacro("utils/RooUtils.C+");

// Loading fitting functions
gInterpreter->LoadMacro("functions/RooDrellYanC0.cxx+");
gInterpreter->LoadMacro("functions/RooDrellYanC1.cxx+");
gInterpreter->LoadMacro("functions/RooDrellYanC2.cxx+");
gInterpreter->LoadMacro("functions/RooQCD.cxx+");
gInterpreter->LoadMacro("functions/RooTTbar.cxx+");
}
