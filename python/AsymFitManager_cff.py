import FWCore.ParameterSet.Config as cms

# Below are defaults for DY->mumu with 200 < M_mumu < 400.
AsymFitManager = cms.PSet(
    # which type of leptons (and thus which mass to use)
    doingElectrons = cms.bool(False),

    # whether to use the mistag correction (i.e. to include
    # omega(y) in the pdf)
    correctMistags = cms.bool(True),

    # if correctMistags, and not calculateMistag, useMistagHist
    # decides between using the values of the 2D mistag
    # probability histogram
    useMistagHist = cms.bool(False),

    # whether to use the calculation of the mistag prob instead of
    # the parameterization
    calculateMistag = cms.bool(True),

    # specifies the pdf for the dilepton mass: 1 = falling
    # exponential, 2 = lorentzian, 3 = exp + lorentzian. if type
    # is 2 or 3 then also need to specify peakMass below.
    massDistType = cms.int32(1),

    # where the lorentzian peak should be found; only used if
    # massDistType is 2 or 3. This also gets used as M_ll in the
    # mistag calculation for the 2D fits.
    peakMass = cms.double(250.0),

    # max pT and rapidity expected
    maxPt = cms.double(1000.0),
    maxRapidity = cms.double(3.5),

    # the mass window to use for the measurement (counting/fits)
    fitWindow = cms.vdouble(200.0, 400.0),

    # the resolution used in the Gaussian smearing; order is
    # cos_cs, rapidity, pT, phi, mass, phi_cs
    recSigma = cms.vdouble(0.00023, 0.005, 2.3, 0.13, 4.2, 0.13),
)
