import FWCore.ParameterSet.Config as cms

# The below is mostly copied from
# https://twiki.cern.ch/twiki/bin/view/CMS/Collisions2010Recipes ; the
# user is responsible for checking that what's used is up-to-date and
# appropriate to their analysis.

hltPhysicsDeclared = cms.EDFilter('HLTPhysicsDeclared',
                                  invert = cms.bool(False),
                                  L1GtReadoutRecordTag = cms.InputTag('gtDigis')
                                  )

primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
                                   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                   minimumNDOF = cms.uint32(4),
                                   maxAbsZ = cms.double(24),
                                   maxd0 = cms.double(2)
                                   )

# Obsolete; not included in goodDataAll
noscraping = cms.EDFilter('FilterOutScraping',
                          applyfilter = cms.untracked.bool(True),
                          debugOn = cms.untracked.bool(False),
                          numtrack = cms.untracked.uint32(10),
                          thresh = cms.untracked.double(0.25)
                          )

# Some more filters are defined in RecoMET/METFilters/python/metFilters_cff.py
