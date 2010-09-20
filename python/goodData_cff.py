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

noscraping = cms.EDFilter('FilterOutScraping',
                          applyfilter = cms.untracked.bool(True),
                          debugOn = cms.untracked.bool(False),
                          numtrack = cms.untracked.uint32(10),
                          thresh = cms.untracked.double(0.25)
                          )

# Instead of filtering out events at PAT-tupling time, schedule
# separate paths for all the "good data" filters so that the results
# of them get stored in a small TriggerResults::PAT object that can be
# read and used to filter events in the analyzer process using
# e.g. the filter in Zprime2muAnalysis_cff.py. (Useful so we don't
# have to keep around the entire generalTracks collection for
# noscraping, for example.)
#
# Make one for each so they can be accessed separately in the
# TriggerResults object; the "All" path isn't necessary because it
# could be emulated using the AND of all of the separate ones, but
# it's nice for convenience.
goodDataHLTPhysicsDeclared = cms.Path(hltPhysicsDeclared)
goodDataPrimaryVertexFilter = cms.Path(primaryVertexFilter)
goodDataNoScraping = cms.Path(noscraping)
goodDataAll = cms.Path(hltPhysicsDeclared * primaryVertexFilter * noscraping)
