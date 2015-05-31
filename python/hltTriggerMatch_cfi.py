import FWCore.ParameterSet.Config as cms


### ==== Unpack trigger, and match ====
#from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
#patTriggerFull = patTrigger.clone(onlyStandAlone = False)
#patTriggerMuon =  patTrigger.clone(src = cms.InputTag("patTriggerFull"),
 #                                  collections = cms.vstring("hltL3MuonCandidates")
  #                                ) 


patTrigger = cms.EDProducer( "PATTriggerProducer", 
                             onlyStandAlone = cms.bool( False ), 
                             processName    = cms.string( "HLT" )             
                             )

muonTriggerMatchHLTMuons = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                          src = cms.InputTag( 'cleanPatMuons' ),
                                          matched = cms.InputTag( 'patTrigger' ),
                                          matchedCuts = cms.string('type("TriggerMuon") && (path("HLT_Mu9*",1,0) || path("HLT_Mu15*",1,0) || path("HLT_Mu24_v*",1,0)|| path("HLT_Mu24*",1,0) || path("HLT_Mu30*",1,0) || path("HLT_Mu40*",1,0))'),
                                          maxDPtRel   = cms.double( 1. ), 
                                          maxDeltaR   = cms.double( 0.2 ),
                                          maxDeltaEta = cms.double( 0.2 ),
                                          resolveAmbiguities    = cms.bool( True ),
                                          resolveByMatchQuality = cms.bool( True )
                                          )

## ==== Embed ====
#from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import patMuonsWithTrigger
cleanPatMuonsTriggerMatch =  cms.EDProducer("PATTriggerMatchMuonEmbedder",                                           
                                            src = cms.InputTag("cleanPatMuons"),                                         
                                            matches = cms.VInputTag("muonTriggerMatchHLTMuons")  
                                            )

## ==== Trigger Sequence ====
patTriggerMatching = cms.Sequence(patTrigger *
                                  muonTriggerMatchHLTMuons *
                                  cleanPatMuonsTriggerMatch)










# JMTBAD to drop
#muonTriggerMatchHLTMuons = cms.EDProducer('PATTriggerMatcherDRLessByR',
 #                                         src = cms.InputTag('cleanPatMuons'),
  #                                        matched = cms.InputTag('patTrigger'),
   #                                       matchedCuts           = cms.string('type("TriggerMuon") && (path("HLT_Mu9*",1,0) || path("HLT_Mu15*",1,0) || path("HLT_Mu24_v*",1,0)|| path("HLT_Mu24*",1,0) || path("HLT_Mu30*",1,0) || path("HLT_Mu40*",1,0))'),
    #maxDPtRel             = cms.double(1),
    #maxDeltaR             = cms.double(0.2),
   # resolveAmbiguities    = cms.bool(True),
   # resolveByMatchQuality = cms.bool(True)
#)

trigger_pt_threshold = 40
offline_pt_threshold = 45
trigger_paths = ['HLT_Mu40_eta2p1_v%i' % i for i in (9,10,11)]
trigger_match = 'userFloat("TriggerMatchPt") > %(trigger_pt_threshold)i && abs(userFloat("TriggerMatchEta")) < 2.1' % locals()

overall_prescale = 300
prescaled_trigger_pt_threshold = 24
prescaled_offline_pt_threshold = 27
prescaled_trigger_paths = ['HLT_Mu24_eta2p1_v%i' % i for i in (3,4,5)]
prescaled_trigger_match = trigger_match.replace('Trigger', 'prescaledTrigger').replace('%i' % trigger_pt_threshold, '%i' % prescaled_trigger_pt_threshold)
