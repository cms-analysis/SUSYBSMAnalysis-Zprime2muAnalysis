import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction

HardInteractionFilter = cms.EDFilter('HardInteractionFilter',
                                     hardInteraction = hardInteraction,
                                     min_mass = cms.double(0),
                                     max_mass = cms.double(1e99),
                                     use_resonance_mass = cms.bool(False),
                                     max_muon_eta = cms.double(1e99),
                                     min_muon_pt = cms.double(0),
                                     )
