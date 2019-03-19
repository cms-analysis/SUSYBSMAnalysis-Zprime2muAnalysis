import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction,hardInteraction_MiniAOD

ResolutionUsingMC = cms.EDAnalyzer('ResolutionUsingMC',
                                  nbins = cms.uint32(500),
                                  min_mass = cms.double(0),
                                  max_mass = cms.double(10000),
                                   hardInteraction = hardInteraction,
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   use_vertex_mass = cms.bool(False),
                                   )

ResolutionUsingMC_MiniAOD = cms.EDAnalyzer('ResolutionUsingMC',
                                  nbins = cms.uint32(500),
                                  min_mass = cms.double(0),
                                  max_mass = cms.double(10000),
                                   hardInteraction = hardInteraction_MiniAOD,
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   use_vertex_mass = cms.bool(False),
                                   )
