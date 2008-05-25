import FWCore.ParameterSet.Config as cms

def attachTrigComparison(process):
    process.Zprime2muTrigComparison = cms.EDAnalyzer(
        'Zprime2muTrigComparison',
        process.Zprime2muAnalysisCommon,
        process.plainAnalysisPSet,
        process.recLevelHelperPSet,
        )

    process.analysis = cms.Path(process.Zprime2muTrigComparison)
