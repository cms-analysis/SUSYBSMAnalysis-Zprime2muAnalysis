import FWCore.ParameterSet.Config as cms


muonTriggerMatchHLTMuons = cms.EDProducer('PATTriggerMatcherDRLessByR',
                                          src = cms.InputTag( 'cleanPatMuons' ),
                                          matched = cms.InputTag( 'patTrigger' ),
                                          matchedCuts = cms.string('type("TriggerMuon") && (path("HLT_Mu9*",1,0) || path("HLT_Mu15*",1,0) || path("HLT_Mu24_v*",1,0)|| path("HLT_Mu24*",1,0) || path("HLT_Mu30*",1,0) || path("HLT_Mu40*",1,0) || path("HLT_Mu45*",1,0) || path("HLT_Mu50*",1,0))'),
                                          maxDPtRel   = cms.double( 1. ), 
                                          maxDeltaR   = cms.double( 0.2 ),
                                          resolveAmbiguities    = cms.bool( True ),
                                          resolveByMatchQuality = cms.bool( True )
                                          )

muonTriggerMatchHLTMuonsMiniAOD = cms.EDProducer('PATTriggerMatcherDRLessByR',
                                          src = cms.InputTag( 'slimmedMuons' ),
                                          matched = cms.InputTag( 'patTrigger' ),
                                          matchedCuts = cms.string('type("TriggerMuon") && (path("HLT_Mu9*",1,0) || path("HLT_Mu15*",1,0) || path("HLT_Mu24_v*",1,0)|| path("HLT_Mu24*",1,0) || path("HLT_Mu30*",1,0) || path("HLT_Mu40*",1,0) || path("HLT_Mu45*",1,0) || path("HLT_Mu50*",1,0))'),
                                          maxDPtRel   = cms.double( 1. ), 
                                          maxDeltaR   = cms.double( 0.2 ),
                                          resolveAmbiguities    = cms.bool( True ),
                                          resolveByMatchQuality = cms.bool( True )
)

#trigger_pt_threshold = 45
#offline_pt_threshold = 48 #?
trigger_pt_threshold = 50
offline_pt_threshold = 53 #?
trigger_paths = ['HLT_Mu50_v%i' % i for i in (6, 7, 8, 9, 10, 11)]
#trigger_paths = ['HLT_Mu45_eta2p1_v%i' % i for i in (1,2)]
#trigger_paths = ['HLT_Mu45_eta2p1_v1']
#trigger_paths = ['HLT_Mu50_v1']
#trigger_match = 'userFloat("TriggerMatchPt") > %(trigger_pt_threshold)i && abs(userFloat("TriggerMatchEta")) < 2.1' % locals()
trigger_match = 'userFloat("TriggerMatchPt") > %(trigger_pt_threshold)i ' % locals()
#trigger_match = '1>0'

#overall_prescale = 1
overall_prescale = 350
prescaled_trigger_pt_threshold = 24
prescaled_offline_pt_threshold = 27
#prescaled_trigger_paths = ['HLT_L1SingleMuOpen_v1']# for first collisions
#prescaled_trigger_paths = ['HLT_L1SingleMuOpen_v1','HLT_L1SingleMu3p5_v1']
#prescaled_trigger_paths = ['HLT_Mu27_v1']
#HLT_Mu24_eta2p1_v
prescaled_trigger_paths = ['HLT_Mu27_v%i' % i for i in (6, 7, 8, 9, 10, 11)]
prescaled_trigger_match = trigger_match.replace('Trigger', 'prescaledTrigger').replace('%i' % trigger_pt_threshold, '%i' % prescaled_trigger_pt_threshold)


# -- for updated plugins/Zprime2muLeptonProducer_miniAOD.cc
# Mu50:     hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q::HLT
# OldMu100: hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q::HLT
# TkMu100:  hltL3fL1sMu25f0TkFiltered100Q::HLT

# Mu27:     hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::HLT
# Mu27      hltL3fL1sMu25L1f0L2f10QL3Filtered27Q::HLT (???)
# TkMu27:   hltL3fL1sMu22Or25f0TkFiltered27Q::HLT

def make_string_cut_for_trigger_matching( list_filters, list_filters_pt ):
  cut = ''
  if len(list_filters) != len(list_filters_pt):
    print 'len(list_filters) != len(list_filters_pt) -> return ', cut
    return cut
  for i, f in enumerate(list_filters):
    if f != list_filters[-1]:
      cut += 'userFloat("%sTriggerMatchPt")>=%i || ' % (list_filters[i], list_filters_pt[i])
    else:
      cut += 'userFloat("%sTriggerMatchPt")>=%i ' % (list_filters[i], list_filters_pt[i])
  return cut

trigger_filters = [
                    'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
                    'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q',
                    'hltL3fL1sMu25f0TkFiltered100Q'
                  ]
trigger_filters_pt = [
                    50,
                    100,
                    100
                  ]
prescaled_trigger_filters = [
                    'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q'
                  ]
prescaled_trigger_filters_pt = [
                    27
                  ]

trigger_match_2018 = make_string_cut_for_trigger_matching( trigger_filters, trigger_filters_pt )
prescaled_trigger_match_2018 = make_string_cut_for_trigger_matching( prescaled_trigger_filters, prescaled_trigger_filters_pt )

#trigger_match_2018 = 'userFloat("%sTriggerMatchPt") >= %i || ' \
#                     'userFloat("%sTriggerMatchPt") >= %i || ' \
#                     'userFloat("%sTriggerMatchPt") >= %i ' % tuple([i for pair in zip(trigger_filters, trigger_filters_pt) for i in pair])



### ==== Unpack trigger, and match ==== Not needed

#from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
#patTrigger = cms.EDProducer( "PATTriggerProducer", 
                           #  onlyStandAlone = cms.bool( False ), 
                           #  processName    = cms.string( "HLT" )             
                           #  )
#patTriggerFull = patTrigger.clone(onlyStandAlone = False)
#patTriggerMuon =  patTrigger.clone(src = cms.InputTag("patTriggerFull"),
 #                                  collections = cms.vstring("hltL3MuonCandidates")
  #                                ) 


## ==== Embed ==== Already done in the package, not needed

#cleanPatMuonsTriggerMatchEmbedded =  cms.EDProducer("PATTriggerMatchMuonEmbedder",                                           
                                          #  src = cms.InputTag("cleanPatMuons"),                                         
                                          #  matches = cms.VInputTag("muonTriggerMatchHLTMuons")  
                                          #  )

## ==== Trigger Sequence ====
#patTriggerMatching = cms.Sequence(patTrigger *
                                  #muonTriggerMatchHLTMuons *
                                 # cleanPatMuonsTriggerMatchEmbedded)












