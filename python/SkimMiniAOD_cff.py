import FWCore.ParameterSet.Config as cms

from math import fabs
##process.load('PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi')
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import*
selectedPatMuons.src = cms.InputTag('slimmedMuons')
selectedPatMuons.cut = cms.string('isGlobalMuon && pt>20 && abs(eta)<2.4')

def addUserData(patMuonProducer, tag=cms.InputTag('selectedPatMuons')):
    patMuonProducer.userData.userCands.src.append(tag)
