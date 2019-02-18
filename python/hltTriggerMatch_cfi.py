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
#trigger_paths = ['HLT_Mu50_v%i' % i for i in (6, 7, 8, 9, 10, 11)]
#trigger_paths = ['HLT_Mu45_eta2p1_v%i' % i for i in (1,2)]
#trigger_paths = ['HLT_Mu45_eta2p1_v1']
#trigger_paths = ['HLT_Mu50_v1']
#trigger_match = 'userFloat("TriggerMatchPt") > %(trigger_pt_threshold)i && abs(userFloat("TriggerMatchEta")) < 2.1' % locals()
trigger_match = 'userFloat("TriggerMatchPt") > %(trigger_pt_threshold)i ' % locals()

#overall_prescale = 1
prescaled_trigger_pt_threshold = 27
prescaled_offline_pt_threshold = 27

# http://fwyzard.web.cern.ch/fwyzard/hlt/2015/summary
#Path HLT_Mu27: 2015
# - first seen online on run 248036 (/cdaq/physics/Run2015/5e33/v1.0/HLT/V1)
# - last  seen online on run 260627 (/cdaq/physics/Run2015/25ns14e33/v4.4.5/HLT/V1)
# - V1: (runs 248036 - 252126)
# - V2: (runs 254227 - 260627)
#Path HLT_Mu50:
# - first seen online on run 248036 (/cdaq/physics/Run2015/5e33/v1.0/HLT/V1)
# - last  seen online on run 260627 (/cdaq/physics/Run2015/25ns14e33/v4.4.5/HLT/V1)
# - V1: (runs 248036 - 252126)
# - V2: (runs 254227 - 260627)
#prescaled_trigger_paths = ['HLT_Mu27_v%i' % i for i in (1,2)]
#trigger_paths = ['HLT_Mu50_v%i' % i for i in (1,2)]

# http://fwyzard.web.cern.ch/fwyzard/hlt/2016/summary
#Path HLT_Mu27: 2016
# - first seen online on run 272760 (/cdaq/physics/Run2016/25ns10e33/v1.0/HLT/V16)
# - last  seen online on run 284044 (/cdaq/physics/Run2016/25ns15e33/v4.2.3/HLT/V2)
# - V2: (runs 272760 - 274443)
# - V3: (runs 274954 - 276244)
# - V4: (runs 276282 - 280385)
# - V5: (runs 281613 - 284044)
#prescaled_trigger_paths = ['HLT_Mu27_v%i' % i for i in (2,3,4,5)]
#Path HLT_TkMu27:
# - first seen online on run 272760 (/cdaq/physics/Run2016/25ns10e33/v1.0/HLT/V16)
# - last  seen online on run 284044 (/cdaq/physics/Run2016/25ns15e33/v4.2.3/HLT/V2)
# - V2: (runs 272760 - 274443)
# - V3: (runs 274954 - 275376)
# - V4: (runs 275656 - 276244)
# - V5: (runs 276282 - 284044)
#Path HLT_Mu50:
# - first seen online on run 272760 (/cdaq/physics/Run2016/25ns10e33/v1.0/HLT/V16)
# - last  seen online on run 284044 (/cdaq/physics/Run2016/25ns15e33/v4.2.3/HLT/V2)
# - V2: (runs 272760 - 274443)
# - V3: (runs 274954 - 276244)
# - V4: (runs 276282 - 280385)
# - V5: (runs 281613 - 284044)
#Path HLT_TkMu50:
# - first seen online on run 274954 (/cdaq/physics/Run2016/25ns10e33/v2.1.0/HLT/V14)
# - last  seen online on run 284044 (/cdaq/physics/Run2016/25ns15e33/v4.2.3/HLT/V2)
# - V1: (runs 274954 - 275376)
# - V2: (runs 275656 - 276244)
# - V3: (runs 276282 - 284044)
#trigger_paths += ['HLT_Mu50_v%i' % i for i in (2,3,4,5)]
#overall_prescale = 320 # 196 pb-1
#overall_prescale = 290  # 1828 pb-1

# http://fwyzard.web.cern.ch/fwyzard/hlt/2017/summary
#Path HLT_Mu27: 2017
# - first seen online on run 296070 (/cdaq/physics/Run2017/2e34/v1.0.0/HLT/V2)
# - last  seen online on run 306460 (/cdaq/physics/Run2017/2e34/v4.2.1/HLT/V2)
# - V6: (runs 296070 - 297057)
# - V7: (runs 297099 - 297505)
# - V8: (runs 297557 - 299329)
# - V9: (runs 299368 - 299649)
# - V10: (runs 300079 - 302019)
# - V11: (runs 302026 - 306171)
# - V12: (runs 306416 - 306460)
#Path HLT_Mu50:
# - first seen online on run 296070 (/cdaq/physics/Run2017/2e34/v1.0.0/HLT/V2)
# - last  seen online on run 306460 (/cdaq/physics/Run2017/2e34/v4.2.1/HLT/V2)
# - V6: (runs 296070 - 297057)
# - V7: (runs 297099 - 297505)
# - V8: (runs 297557 - 299329)
# - V9: (runs 299368 - 299649)
# - V10: (runs 300079 - 302019)
# - V11: (runs 302026 - 306171)
# - V12: (runs 306416 - 306460)
#prescaled_trigger_paths = ['HLT_Mu27_v%i' % i for i in (6, 7, 8, 9, 10, 11)]
#trigger_paths = ['HLT_Mu50_v%i' % i for i in (6, 7, 8, 9, 10, 11)]
#overall_prescale = ?

# http://fwyzard.web.cern.ch/fwyzard/hlt/2018/summary
#Path HLT_Mu27: 2018
# - first seen online on run 315252 (/cdaq/physics/Run2018/2e34/v1.1.0/HLT/V4)
# - last  seen online on run 325175 (/cdaq/physics/Run2018/2e34/v3.6.1/HLT/V2)
# - V12: (runs 315252 - 316271)
# - V13: (runs 316361 - 325175)
#Path HLT_Mu50:
# - first seen online on run 315252 (/cdaq/physics/Run2018/2e34/v1.1.0/HLT/V4)
# - last  seen online on run 325175 (/cdaq/physics/Run2018/2e34/v3.6.1/HLT/V2)
# - V12: (runs 315252 - 316271)
# - V13: (runs 316361 - 325175)
#Path HLT_OldMu100:
# - first seen online on run 315252 (/cdaq/physics/Run2018/2e34/v1.1.0/HLT/V4)
# - last  seen online on run 325175 (/cdaq/physics/Run2018/2e34/v3.6.1/HLT/V2)
# - V3: (runs 315252 - 325175)
#Path HLT_TkMu100:
# - first seen online on run 315252 (/cdaq/physics/Run2018/2e34/v1.1.0/HLT/V4)
# - last  seen online on run 325175 (/cdaq/physics/Run2018/2e34/v3.6.1/HLT/V2)
# - V2: (runs 315252 - 325175)
prescaled_trigger_paths = ['HLT_Mu27_v%i' % i for i in (12,13)]
trigger_paths = ['HLT_Mu50_v%i' % i for i in (12,13)]
#trigger_paths += ['HLT_OldMu100_v%i' % i for i in (3)]
#trigger_paths += ['HLT_TkMu100_v%i' % i for i in (2)]
overall_prescale = 500 

prescaled_trigger_match = trigger_match.replace('Trigger', 'prescaledTrigger').replace('%i' % trigger_pt_threshold, '%i' % prescaled_trigger_pt_threshold)


# -- for updated plugins/Zprime2muLeptonProducer_miniAOD.cc
# Mu50:     hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q::HLT
# OldMu100: hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q::HLT
# TkMu100:  hltL3fL1sMu25f0TkFiltered100Q::HLT

# Mu27:     hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::HLT
# Mu27      hltL3fL1sMu25L1f0L2f10QL3Filtered27Q::HLT (???)
# TkMu27:   hltL3fL1sMu22Or25f0TkFiltered27Q::HLT

def make_string_cut_for_trigger_matching( list_path_names, list_filters_pt, extra=''):
  cut = ''
  if len(list_path_names) != len(list_filters_pt):
    print 'len(list_path_names) != len(list_filters_pt) -> return ', cut
    return cut
  for i, f in enumerate(list_path_names):
    if f != list_path_names[-1]:
      cut += 'userFloat("%s%s_TriggerMatchPt")>=%i || ' % (extra,list_path_names[i], list_filters_pt[i])
    else:                  
      cut += 'userFloat("%s%s_TriggerMatchPt")>=%i ' % (extra,list_path_names[i], list_filters_pt[i])
  return cut

trigger_filters = [
                    'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
                    'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q',
                    'hltL3fL1sMu25f0TkFiltered100Q'
                  ]
trigger_path_names = [
        'Mu50',
        'OldMu100',
        'TkMu100'
        ]
trigger_filters_pt = [
                    50,
                    100,
                    100
                  ]
prescaled_trigger_filters = [
                    'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q'
                  ]
prescaled_trigger_path_names = [
        'Mu27'
        ]
prescaled_trigger_filters_pt = [
                    27
                  ]

trigger_match_2018 = make_string_cut_for_trigger_matching( trigger_path_names, trigger_filters_pt )
prescaled_trigger_match_2018 = make_string_cut_for_trigger_matching( prescaled_trigger_path_names, prescaled_trigger_filters_pt, extra='prescaled')

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












