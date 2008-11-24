import FWCore.ParameterSet.Config as cms

from AsymmetryDataSets_cff import dataSets

def attachAsymmetry(process):
    process.Zprime2muAsymmetry = cms.EDAnalyzer(
        'Zprime2muAsymmetry',
        process.Zprime2muAnalysisCommon,
        process.plainAnalysisPSet,
        dataSets,
        dataSet = cms.string('dy_above400'),
        
        # if noFit is true, only make the histograms -- useful for
        # getting the recSigma information for above
        noFit = cms.bool(False),
        
        # if onlyEvalLLR is true, only evaluate the log-likelihood
        # ratio useful in spin discrimination studies
        onlyEvalLLR = cms.bool(False),

        # fitType determines the form of the cos_cs pdf used in the
        # fits it is a magic number here, but is defined in the
        # FITTYPE enum in Zprime2muAsymmetry.h
        fitType = cms.int32(0),
        
        # only do so many of the fits (up to 6; from gen+gen to
        # rec+rec)
        numFits = cms.int32(6),
        
        # max number of events to read from the parameterization
        # sample
        maxParamEvents = cms.int32(-1),
        
        # whether to use the cached parameterization
        useCachedParams = cms.bool(False),
        
        # the name of the root file that holds the cached parameters
        paramCacheFile = cms.string('cached.root'),
        
        # whether only to calculate the parameterization
        calcParamsOnly = cms.bool(False),
        
        # whether to use the on-peak fit window or the off-peak one
        onPeak = cms.bool(True),
        
        # whether bremsstrahlung was turned on for the generated
        # events
        internalBremOn = cms.bool(True),
        
        # whether to fix b for the simple 1-D fits
        fixbIn1DFit = cms.bool(False),
        
        # whether to use cos_true in 2/6-D fits
        useCosTrueInFit = cms.bool(False),
        
        # whether to correct the cos_cs values using MC truth
        artificialCosCS = cms.bool(False),
        
        # whether to use the mistag correction (i.e. to include
        # omega(y) in the pdf)
        correctMistags = cms.bool(True),
        
        # whether to use the calculation of the mistag prob instead of
        # the parameterization
        calculateMistag = cms.bool(True)
        )

    process.analysis = cms.Path(process.Zprime2muAsymmetry)
