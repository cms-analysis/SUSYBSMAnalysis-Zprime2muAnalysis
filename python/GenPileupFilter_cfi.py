import FWCore.ParameterSet.Config as cms

GenPileupFilter = cms.EDFilter('GenPileupFilter',
                               min_intime = cms.int32(0),
                               max_intime = cms.int32(50),
                               min_late = cms.int32(0),
                               max_late = cms.int32(50),
                               )
