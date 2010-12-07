import FWCore.ParameterSet.Config as cms

dy_gen_mass_cut = 'status == 3 && (pdgId == 23 || pdgId == 32 || pdgId == 39 || pdgId == 5000039) && mass > %(lo)i && mass < %(hi)i'

DYGenMassFilter = cms.EDFilter('CandViewSelector',
                               src = cms.InputTag('prunedGenSimLeptons'),
                               cut = cms.string(dy_gen_mass_cut % {'lo': 0, 'hi': 10000}),
                               filter = cms.bool(True),
                               )
