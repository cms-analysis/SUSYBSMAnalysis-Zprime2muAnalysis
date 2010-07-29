import FWCore.ParameterSet.Config as cms

# JMTBAD muon-photon matching is done here using the default muon
# momentum and not whichever refit momentum will be eventually used in
# the analysis. For most muons this shouldn't make much difference as
# the matching is in eta-phi, but if the default is really bad and the
# selected refit manages to recover, the photon match may be
# wrong. The photon brem recovery needs to be studied more anyway.
muonPhotonMatch = cms.EDProducer('TrivialDeltaRViewMatcher',
                                 src     = cms.InputTag('cleanPatMuons'),
                                 matched = cms.InputTag('cleanPatPhotons'),
                                 distMin = cms.double(0.1)
                                 )

def addUserData(patMuonProducer, tag=cms.InputTag('muonPhotonMatch')):
    patMuonProducer.userData.userCands.src.append(tag)
