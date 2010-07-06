import FWCore.ParameterSet.Config as cms

from AsymFitManager_cff import AsymFitManager

def attachAsymmetry(process, **kwargs):
    module = process.Zprime2muAsymmetry = cms.EDAnalyzer(
        'Zprime2muAsymmetry',
        process.Zprime2muAnalysisCommon,
        AsymFitManager,

        # if noFit is true, only make the histograms -- useful for
        # getting the recSigma information for above
        noFit = cms.bool(False),
        
        # only do so many of the fits (up to 6; from gen+gen to
        # rec+rec)
        numFits = cms.int32(6),
        
        # the name of the root file that holds the cached parameters
        paramCacheFile = cms.string('AsymmetryParametrizer.root'),
        
        # whether bremsstrahlung was turned on for the generated
        # events
        internalBremOn = cms.bool(True),
        
        # whether to fix b for the simple 1-D fits
        fixbIn1DFit = cms.bool(False),
        
        # whether to use cos_true in 2/6-D fits
        useCosTrueInFit = cms.bool(False),
        
        # whether to correct the cos_cs values using MC truth
        artificialCosCS = cms.bool(False),
        )

    process.analysis = cms.Path(module)

    for key, val in kwargs.items():
        setattr(module, key, val)

__all__ = ['attachAsymmetry']
