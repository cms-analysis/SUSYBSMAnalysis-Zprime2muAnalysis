import FWCore.ParameterSet.Config as cms
from Configuration.Generator.PythiaUESettings_cfi import pythiaUESettingsBlock

source = cms.Source("EmptySource")

generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    comEnergy = cms.double(10000.0),
    crossSection = cms.untracked.double(5.46e-2),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring(
            "MSEL        = 0   ! user defined processes",
            "MSUB(1)     = 1   ! ff -> gamma*/Z0/Z'",
            "MSTP(43)    = 3   ! complete Z0/gamma* interference",
            "CKIN(1)     = 500 ! min sqrt(s hat) (GeV)",
            "CKIN(2)     = -1  ! (no) max sqrt(s hat) (GeV)",
            "MDME(174,1) = 0   !Z decay into d dbar",        
            "MDME(175,1) = 0   !Z decay into u ubar",
            "MDME(176,1) = 0   !Z decay into s sbar",
            "MDME(177,1) = 0   !Z decay into c cbar",
            "MDME(178,1) = 0   !Z decay into b bbar",
            "MDME(179,1) = 0   !Z decay into t tbar",
            "MDME(182,1) = 0   !Z decay into e- e+",
            "MDME(183,1) = 0   !Z decay into nu_e nu_ebar",
            "MDME(184,1) = 1   !Z decay into mu- mu+",
            "MDME(185,1) = 0   !Z decay into nu_mu nu_mubar",
            "MDME(186,1) = 0   !Z decay into tau- tau+",
            "MDME(187,1) = 0   !Z decay into nu_tau nu_taubar",
        ),
        parameterSets = cms.vstring(
            "pythiaUESettings", 
            "processParameters"
        )
    )
)

ProductionFilterSequence = cms.Sequence(generator)
