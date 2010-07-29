import FWCore.ParameterSet.Config as cms

# The below is mostly copied from
# https://twiki.cern.ch/twiki/bin/view/CMS/Collisions2010Recipes ; the
# user is responsible for checking that what's used is up-to-date and
# appropriate to their analysis.
#
# Special reminder about MC: bit 0 in L1T1 is not emulated as it would
# just fire all the time, so for running on MC the
# L1T1.L1SeedsLogicalExpression must be changed.

hltPhysicsDeclared = cms.EDFilter('HLTPhysicsDeclared',
                                  invert = cms.bool(False),
                                  L1GtReadoutRecordTag = cms.InputTag('gtDigis')
                                  )

primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
                                   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                   minimumNDOF = cms.uint32(4),
                                   maxAbsZ = cms.double(15),
                                   maxd0 = cms.double(2)
                                   )

noscraping = cms.EDFilter('FilterOutScraping',
                          applyfilter = cms.untracked.bool(True),
                          debugOn = cms.untracked.bool(False),
                          numtrack = cms.untracked.uint32(10),
                          thresh = cms.untracked.double(0.25)
                          )

from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
L1T1 = hltLevel1GTSeed.clone()
L1T1.L1TechTriggerSeeding = cms.bool(True)
L1T1.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

goodData = cms.Sequence(hltPhysicsDeclared * primaryVertexFilter * noscraping * L1T1)
