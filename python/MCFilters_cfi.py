import FWCore.ParameterSet.Config as cms

# If not running on mini aod, need to switch src input tag to prunedMCLeptons
DYPtZskim = cms.EDFilter('DyPt_ZSkim',
                   src = cms.InputTag('prunedGenParticles'),
                   min_mass = cms.double(0),
                   max_mass = cms.double(100), 
                   )
TTbarGenMassFilter = cms.EDFilter('TTbarSelection',
                   src = cms.InputTag('prunedGenParticles'),
                   min_mass = cms.double(50),
                   max_mass = cms.double(500), 
                   )
DibosonGenMassFilter= cms.EDFilter('DibosonGenMass',
                       src = cms.InputTag('prunedGenParticles'),
                       min_mass = cms.double(50),
                       max_mass = cms.double(200), 
                       )
TauTauFilter= cms.EDFilter('TauTauSelection',
                       src = cms.InputTag('prunedGenParticles'),
                       )

