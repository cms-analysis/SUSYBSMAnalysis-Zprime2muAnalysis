import FWCore.ParameterSet.Config as cms

from ResolutionDataSets_cff import dataSets

# Idea: pass in the parameters below as arguments to the function? 
# Stick with the replace paradigm for now.
def attachResolution(process):
    process.Zprime2muResolution = cms.EDAnalyzer(
        'Zprime2muResolution',
        process.Zprime2muAnalysisCommon,
        dataSets,
        dataSet = cms.string('Zp1000'),
        )
    process.analysisResolution = cms.Path(process.Zprime2muResolution)
