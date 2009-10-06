import FWCore.ParameterSet.Config as cms
from Configuration.Generator.PythiaUESettings_cfi import pythiaUESettingsBlock

source = cms.Source("EmptySource")

generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    comEnergy = cms.double(10000.0),
    crossSection = cms.untracked.double(1.228e3),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring(
            'MSEL        = 11 ', 
            'MDME(174,1) = 0    !Z0 decay into d dbar', 
            'MDME(175,1) = 0    !Z0 decay into u ubar', 
            'MDME(176,1) = 0    !Z0 decay into s sbar', 
            'MDME(177,1) = 0    !Z0 decay into c cbar', 
            'MDME(178,1) = 0    !Z0 decay into b bbar', 
            'MDME(179,1) = 0    !Z0 decay into t tbar', 
            'MDME(182,1) = 0    !Z0 decay into e- e+', 
            'MDME(183,1) = 0    !Z0 decay into nu_e nu_ebar', 
            'MDME(184,1) = 1    !Z0 decay into mu- mu+', 
            'MDME(185,1) = 0    !Z0 decay into nu_mu nu_mubar', 
            'MDME(186,1) = 0    !Z0 decay into tau- tau+', 
            'MDME(187,1) = 0    !Z0 decay into nu_tau nu_taubar', 
            'CKIN(1)     = 40.',
            'CKIN(2)     = -1.'
        ),
        parameterSets = cms.vstring(
            'pythiaUESettings', 
            'processParameters'
        )
    )
)

ProductionFilterSequence = cms.Sequence(generator)
