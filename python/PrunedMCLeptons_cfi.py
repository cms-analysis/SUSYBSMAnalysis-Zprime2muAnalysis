import FWCore.ParameterSet.Config as cms

# Do some further pruning of the genSimLeptons as
# detailed in the comments below.
prunedMCLeptons = cms.EDProducer('GenParticlePruner',
                                 src = cms.InputTag('genParticles'),
                                 select = cms.vstring(
                                     'drop *',
                                     # Keep all doc-line particles (e.g. the Z0/Z', its parents, daughters, etc.).
                                     'keep status == 3',
                                     # Keep any final-state muon (either from the generator with status 1,
                                     # or from GEANT with status 8 as set above) that has potential to
                                     # reach the muon system (should reevaluate for muons produced far from
                                     # the IP), and all of their ancestors.
                                     '++keep abs(pdgId) == 13 && (status == 1 || status == 8)',
                                     # Keep any final-state non-GEANT electron and all of their ancestors.
                                     '++keep abs(pdgId) == 11 && status == 1',
                                     # Just to be sure, keep all the descendants of the Z0/Z'/G*.
                                     'keep++ pdgId == 23 || pdgId == 32 || pdgId == 39 || pdgId == 5000039',
                                     )
                                 )
