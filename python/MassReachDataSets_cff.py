import FWCore.ParameterSet.Config as cms

dataSets = cms.PSet(

# 50pb-1 Early paper excercise
Zssm1200_EPE = cms.PSet(
    # Model name and its id (true peak mass).
    resModel = cms.string('Zssm'),
    resMassId = cms.uint32(1200),
    # Number of bins and mass window for plots.
    nBins = cms.uint32(24),
    massWin = cms.vdouble(600.0, 1800.0),
    # Lower and upper mass limits, used to pre-select events based on
    # true dimuon masses to avoid double-counting.
    # Order of GenMass and other arrays: signal, bckg #1, bckg #2
    # (see doc/README). 
    lowerGenMass = cms.vdouble(600.0, 600.0, 200.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 600.0),
    # Production cross-section times branching ratio, in fb (PYTHIA 6.227S,
    # CTEQ6L)
    XSec = cms.vdouble(111, 54, 1620),
    # Event weights for 50pb-1 (0.2775, 0.135, 4.050)
    
    # K-factors.  Currently use the mass-average value of QCD NNLO/LO K-factor,
    # 1.35 +/- 0.05, calculated by M. Schmitt (SUSY/BSM meeting, Dec.16, 2005).
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    # Filter efficiency.
    effFilter = cms.vdouble(1.0, 1.0, 1.0),
    # Number of *generated* events in each sample.  Used to switch between
    # samples, and for normalization of histograms.
    nGenEvents = cms.vuint32(20000, 20000, 20000)
),


# Physics TDR parameters.
Zssm1000_PTDR = cms.PSet(
    # Model name and its id (true peak mass).
    resModel = cms.string('Zssm'),
    resMassId = cms.uint32(1000),
    # Number of bins and mass window for plots.
    nBins = cms.uint32(24),
    massWin = cms.vdouble(400.0, 1600.0),
    # Lower and upper mass limits, used to pre-select events based on
    # true dimuon masses to avoid double-counting.
    # Order of GenMass and other arrays: signal, bckg #1, bckg #2
    # (see doc/README). 
    lowerGenMass = cms.vdouble(400.0, 400.0, 200.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 400.0),
    # Production cross-section times branching ratio, in fb (PYTHIA 6.227S,
    # CTEQ6L)
    XSec = cms.vdouble(611.0, 221.0, 2479.0),
    # K-factors.  Currently use the mass-average value of QCD NNLO/LO K-factor,
    # 1.35 +/- 0.05, calculated by M. Schmitt (SUSY/BSM meeting, Dec.16, 2005).
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    # Filter efficiency.
    effFilter = cms.vdouble(1.0, 1.0, 1.0),
    # Number of *generated* events in each sample.  Used to switch between
    # samples, and for normalization of histograms.
    nGenEvents = cms.vuint32(1000, 1000, 1000)
),

# CMSSW_1_6_X parameters.
Zssm1000 = cms.PSet(
    resModel = cms.string('Zssm'),
    resMassId = cms.uint32(1000),
    nBins = cms.uint32(24),
    massWin = cms.vdouble(400.0, 1600.0),
    lowerGenMass = cms.vdouble(400.0, 400.0, 200.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 400.0),
    XSec = cms.vdouble(624.0, 100.0, 2558.0),
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    effFilter = cms.vdouble(1.0, 0.777, 0.587),
    nGenEvents = cms.vuint32(1000, 1000, 1000),
),

Zpsi1000 = cms.PSet(
    resModel = cms.string('Zpsi'),
    resMassId = cms.uint32(1000),
    nBins = cms.uint32(24),
    massWin = cms.vdouble(400.0, 1600.0),
    lowerGenMass = cms.vdouble(400.0, 400.0, 200.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 400.0),
    XSec = cms.vdouble(362.0, 100.0, 2558.0),
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    effFilter = cms.vdouble(1.0, 0.777, 0.587),
    nGenEvents = cms.vuint32(10000, 1000, 10000)
),

Zssm1500 = cms.PSet(
    resModel = cms.string('Zssm'),
    resMassId = cms.uint32(1500),
    nBins = cms.uint32(32),
    massWin = cms.vdouble(600.0, 2200.0),
    lowerGenMass = cms.vdouble(600.0, 600.0, 400.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 600.0),
    XSec = cms.vdouble(120.3, 50.6, 220.0),
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    effFilter = cms.vdouble(1.0, 1.0, 0.777),
    nGenEvents = cms.vuint32(2000, 2000, 2000)
),

Zpsi1500 = cms.PSet(
    resModel = cms.string('Zpsi'),
    resMassId = cms.uint32(1500),
    nBins = cms.uint32(32),
    massWin = cms.vdouble(600.0, 2200.0),
    lowerGenMass = cms.vdouble(600.0, 600.0, 400.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 600.0),
    XSec = cms.vdouble(75.4, 50.6, 220.0),
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    effFilter = cms.vdouble(1.0, 1.0, 0.777),
    nGenEvents = cms.vuint32(2000, 2000, 2000)
),

Zssm2000 = cms.PSet(
    resModel = cms.string('Zssm'),
    resMassId = cms.uint32(2000),
    nBins = cms.uint32(20),
    massWin = cms.vdouble(1000.0, 3000.0),
    lowerGenMass = cms.vdouble(1000.0, 1000.0, 500.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 1000.0),
    XSec = cms.vdouble(25.0, 6.6, 100.0),
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    effFilter = cms.vdouble(1.0, 1.0, 0.777),
    nGenEvents = cms.vuint32(10000, 1000, 10000)
),

Zpsi2000 = cms.PSet(
    resModel = cms.string('Zpsi'),
    resMassId = cms.uint32(2000),
    nBins = cms.uint32(20),
    massWin = cms.vdouble(1000.0, 3000.0),
    lowerGenMass = cms.vdouble(1000.0, 1000.0, 500.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 1000.0),
    XSec = cms.vdouble(12.9, 6.6, 100.0),
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    effFilter = cms.vdouble(1.0, 1.0, 0.777),
    nGenEvents = cms.vuint32(10000, 1000, 10000)
),

Zssm3000 = cms.PSet(
    resModel = cms.string('Zssm'),
    resMassId = cms.uint32(3000),
    nBins = cms.uint32(30),
    massWin = cms.vdouble(1500.0, 4500.0),
    lowerGenMass = cms.vdouble(1500.0, 1500.0, 1000.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 1500.0),
    XSec = cms.vdouble(2.8, 1.1, 6.6),
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    nGenEvents = cms.vuint32(1000, 1000, 1000)
),

Zssm5000 = cms.PSet(
    resModel = cms.string('Zssm'),
    resMassId = cms.uint32(5000),
    nBins = cms.uint32(40),
    massWin = cms.vdouble(3000.0, 7000.0),
    lowerGenMass = cms.vdouble(3000.0, 3000.0, 2000.0),
    upperGenMass = cms.vdouble(99999.0, 99999.0, 3000.0),
    XSec = cms.vdouble(0.05, 0.02, 0.24),
    KFactor = cms.vdouble(1.35, 1.35, 1.35),
    nGenEvents = cms.vuint32(1000, 1000, 1000)
)

)
