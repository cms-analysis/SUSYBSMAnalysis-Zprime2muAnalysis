import FWCore.ParameterSet.Config as cms

#from SUSYBSMAnalysis.Zprime2muAnalysis.data.Zprime2muAnalysisCommon_cff import process, Zprime2muAnalysisCommon
from Zprime2muAnalysisCommon_cff import process, Zprime2muAnalysisCommon
from ResolutionDataSets_cff import dataSets

process.Zprime2muResolution = cms.EDAnalyzer('Zprime2muResolution',
                                             Zprime2muAnalysisCommon,
                                             dataSets,
                                             verbosity         = cms.untracked.int32(0),
                                             outputFile        = cms.untracked.string('muon_resolution.ps'),
                                             histoFile         = cms.untracked.string('resolution_histos.root'),
                                             useHistosFromFile = cms.untracked.bool(False),
                                             dataSet           = cms.string('Zp1000'),
                                             )

process.analysis = cms.Path(process.Zprime2muResolution)
