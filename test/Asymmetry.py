import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAsymmetry_cfi')

process.genDimuons = cms.EDProducer('GenDilCandViewShallowCloneCombiner',
                                    decay = cms.string('genParticles@+ genParticles@-'),
                                    cut = cms.string('')
                                    )

process.p = cms.Path(process.hltFilter * process.Zprime2muAnalysisSequence * process.genDimuons * process.Zprime2muAsymmetry)

process.Zprime2muAsymmetry.numFits = 4

#process.Zprime2muAsymmetry.peakMass = cms.double(91.2)
#process.Zprime2muAsymmetry.massDistType = cms.int32(3),
#process.Zprime2muAsymmetry.maxPt = cms.double(200.0),
#process.Zprime2muAsymmetry.maxRapidity = cms.double(3.5),
#process.Zprime2muAsymmetry.fitWindow = cms.vdouble(80.0, 100.0),
##process.Zprime2muAsymmetry.fitWindow = cms.vdouble(200.0, 500.0),
##   order: cos_cs, rapidity, pT, phi, mass, phi_cs
#process.Zprime2muAsymmetry.recSigma = cms.vdouble(5.8e-4, 8.6e-3, 3.9, 8.6e-2, 6.14, 9.5e-2)


