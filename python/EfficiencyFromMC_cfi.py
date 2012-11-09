import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
from SUSYBSMAnalysis.Zprime2muAnalysis.TriggerDecision_cff import triggerDecision

EfficiencyFromMC = cms.EDAnalyzer('EfficiencyFromMC',
                                  hardInteraction = hardInteraction,
                                  triggerDecision = triggerDecision,
                                  check_l1 = cms.bool(True),
                                  trigger_summary_src = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
                                  nbins = cms.uint32(3200),
                                  min_mass = cms.double(0),
                                  max_mass = cms.double(3200),
                                  use_resonance_mass = cms.bool(False),
                                  use_resonance_mass_denom = cms.bool(False),
                                  dimuon_src = cms.InputTag('dimuons'),
                                  hlt_obj_src = cms.InputTag(''),
                                  hlt_single_min_pt = cms.double(40),
                                  hlt_single_max_eta = cms.double(2.1),
                                  checking_prescaled_path = cms.bool(False),
                                  acceptance_max_eta_1 = cms.double(2.1),
                                  acceptance_max_eta_2 = cms.double(2.4),
                                  acceptance_min_pt = cms.double(45),
                                  )
