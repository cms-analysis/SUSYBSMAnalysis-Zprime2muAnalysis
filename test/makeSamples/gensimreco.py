# cmsDriver.py Configuration/GenProduction/python/PYTHIA6_EXOTICA_ZPrimeSSMMuMu_M500_7TeV_cff.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT:GRun,RAW2DIGI,L1Reco,RECO --conditions START311_V1G1::All --datatier GEN-SIM-RECO --eventcontent RECOSIM --customise Configuration/GenProduction/geometry38X_customize.py,Configuration/GlobalRuns/reco_TLR_311X.py --cust_function customiseGeometry,customisePPMC --processName REDIGI311X --no_exec
# and then hand edited
import sys, os, FWCore.ParameterSet.Config as cms

process = cms.Process('REDIGI311X')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.MixingModule.mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('reco.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('generation_step'))
)

process.GlobalTag.globaltag = 'START311_V2::All'

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    crossSection = cms.untracked.double(2.01),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.832 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.275 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        processParameters = cms.vstring('MSEL=0             ! User defined processes', 
            "MSUB(141)   = 1    ! ff -> gamma/Z0/Z\'", 
            'MSTP(44)    = 3    ! no Zprime/Z/gamma interference', 
            "PMAS(32,1)  = 500 ! Z\' mass (GeV)", 
            'CKIN(1)     = -1  ! lower invariant mass cutoff (GeV)', 
            'CKIN(2)     = -1   ! no upper invariant mass cutoff', 
            'MDME(289,1) = 0    ! d dbar', 
            'MDME(290,1) = 0    ! u ubar', 
            'MDME(291,1) = 0    ! s sbar', 
            'MDME(292,1) = 0    ! c cbar', 
            'MDME(293,1) = 0    ! b bar', 
            'MDME(294,1) = 0    ! t tbar', 
            'MDME(295,1) = -1   ! 4th gen q qbar', 
            'MDME(296,1) = -1   ! 4th gen q qbar', 
            'MDME(297,1) = 0    ! e-     e+', 
            'MDME(298,1) = 0    ! nu_e   nu_ebar', 
            'MDME(299,1) = 1    ! mu-    mu+', 
            'MDME(300,1) = 0    ! nu_mu  nu_mubar', 
            'MDME(301,1) = 0    ! tau    tau', 
            'MDME(302,1) = 0    ! nu_tau nu_taubar', 
            'MDME(303,1) = -1   ! 4th gen l- l+', 
            'MDME(304,1) = -1   ! 4th gen nu nubar', 
            'MDME(305,1) = -1   ! W+ W-', 
            'MDME(306,1) = -1   ! H+ H-', 
            'MDME(307,1) = -1   ! Z0 gamma', 
            'MDME(308,1) = -1   ! Z0 h0', 
            'MDME(309,1) = -1   ! h0 A0', 
            'MDME(310,1) = -1   ! H0 A0'),
        parameterSets = cms.vstring('pythiaUESettings', 'processParameters')
    )
)

process.ProductionFilterSequence = cms.Sequence(process.generator)

process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step])
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

process.XMLFromDBSource.label='Extended'
process.GlobalTag.toGet = cms.VPSet(
	cms.PSet(record = cms.string("GeometryFileRcd"),
		 tag = cms.string("XMLFILE_Geometry_380V3_Extended_mc"),
		 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_GEOMETRY"),
		 label = cms.untracked.string("Extended")
		 )
	)
process.newSeedFromTriplets.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = cms.uint32(100000)
process.newSeedFromPairs.OrderedHitsFactoryPSet.maxElement = cms.uint32(100000)
process.secTriplets.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = cms.uint32(100000)
process.thTripletsA.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = cms.uint32(100000)
process.thTripletsB.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = cms.uint32(100000)
process.fourthPLSeeds.OrderedHitsFactoryPSet.maxElement = cms.uint32(100000)
process.fifthSeeds.OrderedHitsFactoryPSet.maxElement = cms.uint32(100000)
process.pixelVertices.useBeamConstraint = False

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = None
pset = %(pset_fn)s
get_edm_output = 1
number_of_jobs = 100
events_per_job = 550
first_lumi = 1

[USER]
ui_working_dir = crab/crab_gensimreco_PU_%(name)s
copy_data = 1
storage_element = T3_US_FNALLPC
check_user_remote_dir = 0
publish_data = 1
publish_data_name = %(name)s-RECO-413p2-START311V2_PU
dbs_url_for_publication = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
'''

    pythia = {
        'dy': '''[
            'MSEL=0            !User defined processes', 
            'MSUB(1)=1         !Incl Z0/gamma* production', 
            'MSTP(43)=3        !Both Z0 and gamma*', 
            'MDME(174,1)=0     !Z decay into d dbar', 
            'MDME(175,1)=0     !Z decay into u ubar', 
            'MDME(176,1)=0     !Z decay into s sbar', 
            'MDME(177,1)=0     !Z decay into c cbar', 
            'MDME(178,1)=0     !Z decay into b bbar', 
            'MDME(179,1)=0     !Z decay into t tbar', 
            'MDME(182,1)=0     !Z decay into e- e+', 
            'MDME(183,1)=0     !Z decay into nu_e nu_ebar', 
            'MDME(184,1)=1     !Z decay into mu- mu+', 
            'MDME(185,1)=0     !Z decay into nu_mu nu_mubar', 
            'MDME(186,1)=0     !Z decay into tau- tau+', 
            'MDME(187,1)=0     !Z decay into nu_tau nu_taubar', 
            'CKIN(2)=-1        !Maximum sqrt(s_hat) value (=Z mass)'
            ]''',
        'zssm': '''[
            'MSEL=0             ! User defined processes',
            'MSUB(141)   = 1    ! ff -> gamma/Z0/Zprime',
            'MSTP(44)    = 3    ! no Zprime/Z/gamma interference',
            'CKIN(2)     = -1   ! no upper invariant mass cutoff',
            'MDME(289,1) = 0    ! d dbar',
            'MDME(290,1) = 0    ! u ubar',
            'MDME(291,1) = 0    ! s sbar',
            'MDME(292,1) = 0    ! c cbar',
            'MDME(293,1) = 0    ! b bar',
            'MDME(294,1) = 0    ! t tbar',
            'MDME(295,1) = -1   ! 4th gen q qbar',
            'MDME(296,1) = -1   ! 4th gen q qbar',
            'MDME(297,1) = 0    ! e-     e+',
            'MDME(298,1) = 0    ! nu_e   nu_ebar',
            'MDME(299,1) = 1    ! mu-    mu+',
            'MDME(300,1) = 0    ! nu_mu  nu_mubar',
            'MDME(301,1) = 0    ! tau    tau',
            'MDME(302,1) = 0    ! nu_tau nu_taubar',
            'MDME(303,1) = -1   ! 4th gen l- l+',
            'MDME(304,1) = -1   ! 4th gen nu nubar',
            'MDME(305,1) = -1   ! W+ W-',
            'MDME(306,1) = -1   ! H+ H-',
            'MDME(307,1) = -1   ! Z0 gamma',
            'MDME(308,1) = -1   ! Z0 h0',
            'MDME(309,1) = -1   ! h0 A0',
            'MDME(310,1) = -1   ! H0 A0',
            ]''',
        'zpsi': '''[
            # Couplings: first generation
            "PARU(121)=  0.        ! vd",
            "PARU(122)= -0.506809  ! ad",
            "PARU(123)=  0.        ! vu",
            "PARU(124)= -0.506809  ! au",
            "PARU(125)=  0.        ! ve",
            "PARU(126)= -0.506809  ! ae",
            "PARU(127)= -0.253405  ! vnu",
            "PARU(128)= -0.253405  ! anu",
            # Couplings: second generation
            "PARJ(180)=  0.        ! vd",
            "PARJ(181)= -0.506809  ! ad",
            "PARJ(182)=  0.        ! vu",
            "PARJ(183)= -0.506809  ! au",
            "PARJ(184)=  0.        ! ve",
            "PARJ(185)= -0.506809  ! ae",
            "PARJ(186)= -0.253405  ! vnu",
            "PARJ(187)= -0.253405  ! anu",
            # Couplings: third generation
            "PARJ(188)=  0.        ! vd",
            "PARJ(189)= -0.506809  ! ad",
            "PARJ(190)=  0.        ! vu",
            "PARJ(191)= -0.506809  ! au",
            "PARJ(192)=  0.        ! ve",
            "PARJ(193)= -0.506809  ! ae",
            "PARJ(194)= -0.253405  ! vnu",
            "PARJ(195)= -0.253405  ! anu",
            ]'''
        }


    x = [
#       ('dy', None, 120, '7.9'),
        ('dy', None, 200, '0.965'),
        ('dy', None, 500, '0.0269'),
        ('dy', None, 800, '0.0031'),
#       ('dy', None, 1200, None),
#       ('dy', None, 1500, None),
#       ('dy', None, 1800, None),
#       ('zssm',  250, None, None),
        ('zssm',  500, None, '2.01'),
        ('zssm',  750, None, '0.355'),
        ('zssm', 1000, None, '0.0923'),
        ('zssm', 1250, None, '0.028'),
        ('zssm', 1500, None, '0.0099'),
        ('zssm', 1750, None, '0.0037'),
#       ('zssm', 2000, None, None),
#       ('zpsi',  250, None, None),
        ('zpsi',  500, None, None),
        ('zpsi',  750, None, None),
        ('zpsi', 1000, None, None),
        ('zpsi', 1250, None, None),
        ('zpsi', 1500, None, None),
        ('zpsi', 1750, None, None),
#       ('zpsi', 2000, None, None),
        ]
        
    os.system('mkdir -p psets')
    os.system('mkdir -p crab')

    for kind, mass, minmass, sigma in x:
        name = '%s%i' % (kind, mass if mass is not None else minmass)
        pset_fn = 'psets/gensimreco_%(name)s.py' % locals()
        pset = open('gensimreco.py').read()
        pars = pythia[kind] if kind != 'zpsi' else pythia['zssm'] + ' + ' + pythia['zpsi']
        pset += '\n\nprocess.generator.PythiaParameters.processParameters = ' + pars + '\n'
        if mass:
            pset += 'process.generator.PythiaParameters.processParameters.append("PMAS(32,1)  = %(mass)i ! Z\' mass (GeV)")\n' % locals()
        if minmass:
            pset += 'process.generator.PythiaParameters.processParameters.append("CKIN(1)     = %(minmass)i  ! lower invariant mass cutoff (GeV)")\n' % locals()
        if sigma:
            pset += 'process.generator.crossSection = %(sigma)s\n' % locals()
        open(pset_fn, 'wt').write(pset)
        print pset_fn
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not 'testing' in sys.argv:
            os.system('crab -create -submit')
