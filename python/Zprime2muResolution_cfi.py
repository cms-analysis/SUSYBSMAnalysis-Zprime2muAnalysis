import FWCore.ParameterSet.Config as cms

from ResolutionDataSets_cff import dataSets

def attachResolution(process, **kwargs):
    module = process.Zprime2muResolution = cms.EDAnalyzer(
        'Zprime2muResolution',
        process.Zprime2muAnalysisCommon,
        dataSets,
        dataSet              = cms.string('Zp1000'),     # Which set of parameters from dataSets to use.
        leptonsFromDileptons = cms.bool(False),          # Whether only to fill the lepton plots from leptons that make it into dileptons.
        doQoverP             = cms.bool(False)           # Whether the inverse momenta plots will be e.g. q/pt.
        )

    process.analysisResolution = cms.Path(module)

    for key, val in kwargs.items():
        setattr(module, key, val)

__all__ = ['attachResolution']
