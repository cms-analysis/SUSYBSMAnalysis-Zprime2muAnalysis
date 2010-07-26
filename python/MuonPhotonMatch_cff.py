import FWCore.ParameterSet.Config as cms

muonPhotonMatch = cms.EDProducer('TrivialDeltaRViewMatcher',
                                 src     = cms.InputTag('cleanPatMuons'),
                                 matched = cms.InputTag('cleanPatPhotons'),
                                 distMin = cms.double(0.1)
                                 )

## Helper function to add this info into a pat::Muon
def addUserData(patMuonProducer, label="muonPhotonMatch"):
    patMuonProducer.userData.userCands.src += [
        cms.InputTag(label,"")
        ]
