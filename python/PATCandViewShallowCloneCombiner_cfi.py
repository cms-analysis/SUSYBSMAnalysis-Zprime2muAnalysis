import FWCore.ParameterSet.Config as cms

allDimuons = cms.EDProducer('PATCandViewShallowCloneCombiner',
                            decay = cms.string('leptons:muons@+ leptons:muons@-'),
                            cut = cms.string('')
                            )
