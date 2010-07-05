import FWCore.ParameterSet.Config as cms

CocktailAnalyzer = cms.EDAnalyzer('CocktailAnalyzer',
                                  muon_src = cms.InputTag('leptons', 'muons')
                                  )
