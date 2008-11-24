import FWCore.ParameterSet.Config as cms

from ResolutionDataSets_cff import dataSets

# Idea: pass in the parameters below as arguments to the function? 
# Stick with the replace paradigm for now.
def attachResolution(process):
    process.Zprime2muResolution = cms.EDAnalyzer(
        'Zprime2muResolution',
        process.Zprime2muAnalysisCommon,
        process.plainAnalysisPSet,
        process.recLevelHelperPSet,
        dataSets,
        dataSet           = cms.string('Zp1000'),

        # The postscript file for ~87 pages of plots.
        outputFile        = cms.untracked.string('muon_resolution.ps'),

        # If useHistosFromFile is true, pull the histograms from the
        # ROOT file specified (made earlier by TFileService and
        # perhaps merged from many jobs) instead of creating them from
        # the data in the event, draw the postscript file of plots and
        # quit.
        useHistosFromFile = cms.untracked.bool(False),
        histoFile         = cms.untracked.string('resolution_histos.root'),
        )

    process.analysis = cms.Path(process.Zprime2muResolution)
