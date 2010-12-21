import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
from SUSYBSMAnalysis.Zprime2muAnalysis.TriggerDecision_cff import triggerDecision

EfficiencyFromMC = cms.EDAnalyzer('EfficiencyFromMC',
                                  hardInteraction = hardInteraction,
                                  triggerDecision = triggerDecision,
                                  nbins = cms.uint32(2000),
                                  min_mass = cms.double(0),
                                  max_mass = cms.double(2000),
                                  use_resonance_mass = cms.bool(False),
                                  use_resonance_mass_denom = cms.bool(False),
                                  dimuon_src = cms.InputTag('dimuons'),
                                  hlt_obj_src = cms.InputTag(''),
                                  hlt_single_min_pt = cms.double(-1),
                                  acceptance_max_eta = cms.double(2.4),
                                  acceptance_min_pt = cms.double(20),
                                  )
