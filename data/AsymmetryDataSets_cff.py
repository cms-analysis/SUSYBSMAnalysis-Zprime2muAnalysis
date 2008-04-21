import FWCore.ParameterSet.Config as cms

dataSets = cms.PSet(
    
Zssm1000 = cms.PSet(
    outputFileBase = cms.untracked.string('Zssm1000'),
    # EDM ROOT files with generated data in them to obtain
    # mistag parameterizations
    genSampleFiles = cms.vstring('rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_0/Zssm1000_fm_GEN_SIM.forParam.root'),
    peakMass = cms.double(1000.0),
    # massDistType: 1 = falling exponential, 2 = lorentzian peak, 3 = 1+2
    massDistType = cms.int32(3),
    maxPt = cms.double(2000.0),
    maxRapidity = cms.double(2.5),
    genWindow = cms.vdouble(400.0, 2000.0),
    fitWindowOnPeak = cms.vdouble(900.0, 1400.0),
    fitWindowOffPeak = cms.vdouble(500.0, 900.0),
    # recSigma are the sigma of the gaussians used in the smear to simulate
    # detector resolution
    #   order: cos_cs, rapidity, pT, phi, mass, phi_cs
    recSigma = cms.vdouble(0.00032, 0.017, 16.0, 0.2, 33.0, 0.21)
),

Zssm3000 = cms.PSet(
    outputFileBase = cms.untracked.string('Zssm3000'),
    genSampleFiles = cms.vstring('rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_0/Zssm3000_fm_GEN_SIM.forParam.root'),
    peakMass = cms.double(3000.0),
    massDistType = cms.int32(3),
    maxPt = cms.double(3000.0),
    maxRapidity = cms.double(2.5),
    genWindow = cms.vdouble(1500.0, 5000.0),
    fitWindowOnPeak = cms.vdouble(2700.0, 4000.0),
    fitWindowOffPeak = cms.vdouble(1800.0, 2700.0),
    recSigma = cms.vdouble(0.00031, 0.024, 59.0, 0.32, 129.0, 0.36)
),

dy_above400 = cms.PSet(
    genSampleFiles = cms.vstring('rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_0/dy_above400_GEN_SIM.forParam.root'),
    peakMass = cms.double(450.0), # no peak, but set this for completeness
    massDistType = cms.int32(1),
    maxPt = cms.double(1000.0),
    maxRapidity = cms.double(3.5),
    outputFileBase = cms.untracked.string('dy_above400'),
    genWindow = cms.vdouble(400.0, 500.0),
    fitWindowOnPeak = cms.vdouble(400.0, 500.0),
    fitWindowOffPeak = cms.vdouble(400.0, 500.0),
    recSigma = cms.vdouble(0.00024, 0.017, 25.0, 0.29, 46.0, 0.35)
)

)


