import FWCore.ParameterSet.Config as cms

def makeAnalysis(process):
    process.Zprime2muAnalysis = cms.EDAnalyzer(
        'Zprime2muAnalysis',
        Zprime2muAnalysisCommon,
        verbosity = cms.untracked.int32(2)
        )

    process.analysis = cms.Path(process.Zprime2muAnalysis)
