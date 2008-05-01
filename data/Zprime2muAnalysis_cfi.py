import FWCore.ParameterSet.Config as cms

def attachAnalysis(process):
    process.Zprime2muAnalysis = cms.EDAnalyzer(
        'Zprime2muAnalysis',
        process.Zprime2muAnalysisCommon,
        process.recLevelHelperPSet,
        verbosity = cms.untracked.int32(2)
        )

    process.analysis = cms.Path(process.Zprime2muAnalysis)
