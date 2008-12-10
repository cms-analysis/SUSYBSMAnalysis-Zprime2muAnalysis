import FWCore.ParameterSet.Config as cms

from ResolutionDataSets_cff import dataSets

def attachResolution(process, **kwargs):
    module = process.Zprime2muResolution = cms.EDAnalyzer(
        'Zprime2muResolution',
        process.Zprime2muAnalysisCommon,
        dataSets,
        dataSet = cms.string('Zp1000'),
        )

    process.analysisResolution = cms.Path(module)

    for key, val in kwargs.items():
        setattr(module, key, val)

__all__ = ['attachResolution']
