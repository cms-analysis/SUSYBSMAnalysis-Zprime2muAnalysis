import FWCore.ParameterSet.Config as cms

#from SUSYBSMAnalysis.Zprime2muAnalysis.data.Zprime2muAnalysisCommon_cff import process, Zprime2muAnalysisCommon
from Zprime2muAnalysisCommon_cff import process, Zprime2muAnalysisCommon

process.Zprime2muAnalysis = cms.EDAnalyzer('Zprime2muAnalysis',
                                           Zprime2muAnalysisCommon,
                                           verbosity = cms.untracked.int32(0)
)

process.analysis = cms.Path(process.Zprime2muAnalysis)
