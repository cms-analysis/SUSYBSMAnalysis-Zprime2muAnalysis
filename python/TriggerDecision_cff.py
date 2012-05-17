import FWCore.ParameterSet.Config as cms

triggerDecision = cms.PSet(doingElectrons = cms.bool(False),
                           useTrigger = cms.bool(True),
                           l1GtObjectMap = cms.InputTag('hltL1GtObjectMap'),
#                           hltResults = cms.InputTag('TriggerResults', '', 'HLT'),
                           hltResults = cms.InputTag('TriggerResults', '', 'PAT'),
                           # Note: these next two paths should not be
                           # changed just because the trigger menu for
                           # data-taking changes. This module is meant
                           # for MC studies, so if the pT threshold to
                           # be studied is e.g. 40 GeV, that has to be
                           # taken care of in code as there is no such
                           # path in the Spring11/Summer11 MC.
                           l1Paths = cms.vstring('L1_SingleMu12'),
                           hltPaths = cms.vstring('HLT_Mu40_eta2p1_v9')
                           )
