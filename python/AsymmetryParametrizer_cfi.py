import FWCore.ParameterSet.Config as cms

from AsymFitManager_cff import AsymFitManager

AsymmetryParametrizer = cms.EDFilter('AsymmetryParametrizer',
                                     AsymFitManager,
                                     gen_particle_src = cms.InputTag('genParticles'),
                                     internal_brem_on = cms.bool(True),
                                     assemble_only = cms.bool(False),
                                     histos_fn = cms.string('Parametrizer.root'),
                                     postscript_fn = cms.string('Parametrizer.ps'),
                                     )
