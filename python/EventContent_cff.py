import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import *
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

zPrimeEventContent = [
#    'keep recoGenParticles_genParticles_*_*',
    'keep GenEventInfoProduct_*_*_*',
    'keep GenRunInfoProduct_*_*_*',
    'keep recoTrackExtras_standAloneMuons_*_*',
    'keep recoTracks_standAloneMuons__*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep edmTriggerResults_TriggerResults__HLT*',
    'keep L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT*',
    'keep *_hltTriggerSummaryAOD__HLT*',
#    'keep *_genSimParticles_*_*'
    ]

