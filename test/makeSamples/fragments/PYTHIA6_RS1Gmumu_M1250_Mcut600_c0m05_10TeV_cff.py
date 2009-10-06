import FWCore.ParameterSet.Config as cms
from Configuration.Generator.PythiaUESettings_cfi import pythiaUESettingsBlock

source = cms.Source("EmptySource")

generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    comEnergy = cms.double(10000.0),
    crossSection = cms.untracked.double(1e99),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring(
            'MSEL = 0                    ! (D=1) to select between full user control (0, then use MSUB) and some preprogrammed alternative',
            'PMAS(347,1) = 1250.         ! mass of RS Graviton',
            'PARP(50) = 0.27             ! 0.054 == c = k/M_Pl = 0.01  =>  0.27 == c = 0.05',
            'MSUB(391) = 1               ! q qbar -> G* ',
            'MSUB(392) = 1               ! g g -> G*',
            'MDME(4084,1)=0            ! d dbar',
            'MDME(4085,1)=0            ! u ubar',
            'MDME(4086,1)=0            ! s sbar',
            'MDME(4087,1)=0            ! c cbar',
            'MDME(4088,1)=0            ! b bbar',
            'MDME(4089,1)=0            ! t tbar',
            'MDME(4090,1)=-1            ! bprime bprimebar',
            'MDME(4091,1)=-1            ! tprime tprimebar',
            'MDME(4092,1)=0            ! e+ e-',
            'MDME(4093,1)=0            ! nu_e nu_ebar',
            'MDME(4094,1)=1            ! mu- mu+',
            'MDME(4095,1)=0            ! nu_mu nu_mubar',
            'MDME(4096,1)=0            ! tau- tau+',
            'MDME(4097,1)=0            ! nu_tau  nu_taubar',
            'MDME(4098,1)=-1            ! tauprime- tauprime+ ',
            'MDME(4099,1)=-1            ! nuprime_tau nuprime_taubar ',
            'MDME(4100,1)=0            ! g  g  ',
            'MDME(4101,1)=0            ! gamma gamma ',
            'MDME(4102,1)=0            ! Z Z',
            'MDME(4103,1)=0            ! W W',
            "CKIN(1)     = 600  ! lower invariant mass cutoff (GeV)",
            "CKIN(2)     = -1   ! no upper invariant mass cutoff",
            'CKIN(3)=-1.          ! minimum pt hat for hard interactions',
            'CKIN(4)=-1.          ! maximum pt hat for hard interactions',
        ),
        parameterSets = cms.vstring(
            "pythiaUESettings", 
            "processParameters"
        )
    )
)

ProductionFilterSequence = cms.Sequence(generator)
