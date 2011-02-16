import FWCore.ParameterSet.Config as cms

from AsymFitManager_cff import AsymFitManager

Zprime2muAsymmetry = cms.EDAnalyzer('Zprime2muAsymmetry',
                                    AsymFitManager,
                                    gen_particle_src = cms.InputTag('genParticles'),
                                    gen_dilepton_src = cms.InputTag('genDimuons'),
                                    dilepton_src = cms.InputTag('dimuons'),
                                    useGen = cms.bool(True),
                                    noFit = cms.bool(False),
                                    numFits = cms.int32(6),
                                    internalBremOn = cms.bool(True),
                                    fixbIn1DFit = cms.bool(False),
                                    useCosTrueInFit = cms.bool(False),
                                    artificialCosCS = cms.bool(False),
                                    paramCacheFile = cms.string('Parametrizer.root'),
                                    )
