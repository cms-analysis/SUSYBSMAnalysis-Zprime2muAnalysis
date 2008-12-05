import FWCore.ParameterSet.Config as cms

from ResolutionDataSets_cff import dataSets

# Idea: pass in the parameters below as arguments to the function? 
# Stick with the replace paradigm for now.
def attachData(process):
    process.Zprime2muHistos = cms.EDAnalyzer(
        'Zprime2muHistos',
        process.Zprime2muAnalysisCommon,
        dataSets,
        dataSet = cms.string('Zp1000'),
        )
    process.analysis = cms.Path(process.Zprime2muHistos)
