#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

def addGenSimLeptons(process):
    # For muon and electron MC matching, want to be able to match to
    # decays-in-flight produced by GEANT, so make some GenParticles out of
    # the simTracks.
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.GenPlusSim_cfi')
    
    for x in (process.muonMatch, process.electronMatch):
        # Use the new gen + sim MC particles created above.
        x.matched = cms.InputTag('prunedGenSimLeptons')
        
        # PAT muon/electron-MC matching requires, in addition to deltaR < 0.5,
        # the MC and reconstructed leptons to have the same charge, and (reco
        # pt - gen pt)/gen pt < 0.5. Disable these two cuts. If using the other 
        x.checkCharge = False
        x.maxDPtRel = 1e6

    process.patDefaultSequence = cms.Sequence(process.genSimLeptons * process.prunedGenSimLeptons * process.patDefaultSequence._seq)

def addMuonMCClassification(process):
    # Run and embed in the patMuons the classification of muons by
    # their GEANT hits.
    process.load('MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi')
    from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
    addClassByHits(process.patMuons, extraInfo=True)
    process.patDefaultSequence = cms.Sequence(process.muonClassificationByHits * process.patDefaultSequence._seq)

def addMuonStations(process):
    # Embed the muon station counts. This should be obsolete in 37X
    # with the addition to reco::HitPattern -- needs checking.
    process.load('MuonAnalysis.Examples.muonStations_cfi')
    from MuonAnalysis.Examples.muonStations_cfi import addUserData as addStations
    addStations(process.patMuons)
    process.patDefaultSequence.replace(process.patCandidates, process.muonStations * process.patCandidates)

def addMuonHitCount(process):
    # Embed the muon hit counts. Redundant with reco::HitPattern?
    process.load('UserCode.Examples.muonHitCount_cfi')
    from UserCode.Examples.muonHitCount_cfi import addUserData as addHitCount
    addHitCount(process.patMuons)
    process.patDefaultSequence.replace(process.patCandidates, process.muonHitCounts * process.patCandidates)

def removeMuonMCClassification(process):
    process.patDefaultSequence.remove(process.muonClassificationByHits)

    # Remove the InputTags that were added to the userData of the
    # patMuons for the muonClassification.
    def filter(v, s):
        v2 = []
        for x in v:
            if type(x) == cms.InputTag and x.moduleLabel != s:
                v2.append(x)
        return v2
    
    i = process.patMuons.userData.userInts.src.value()
    f = process.patMuons.userData.userFloats.src.value()
    for s in ['classByHitsGlb', 'classByHitsTM', 'classByHitsTMLSAT', 'classByHitsSta']:
        i = filter(i, s)
        f = filter(f, s)
    process.patMuons.userData.userInts.src = i
    process.patMuons.userData.userFloats.src = f

def removeMCUse(process):
    # Remove anything that requires MC truth.
    from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
    removeMCMatching(process, ['All'])
    removeMuonMCClassification(process)
    process.patDefaultSequence.remove(process.genSimLeptons)
    process.patDefaultSequence.remove(process.prunedGenSimLeptons)
    
def changeMuonHLTMatch(process):
    # Configure the PAT trigger matcher as we want it.
    process.load('PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff')
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi')
    process.patTriggerMatcher += process.muonTriggerMatchHLTMuons
    # Removing these three really isn't necessary as the trigger match
    # embedding is by name anyway...
    process.patTriggerMatcher.remove(process.patTriggerMatcherElectron)
    process.patTriggerMatcher.remove(process.patTriggerMatcherMuon)
    process.patTriggerMatcher.remove(process.patTriggerMatcherTau)
    process.patTriggerEvent.patTriggerMatches = ['muonTriggerMatchHLTMuons']
    process.cleanPatMuonsTriggerMatch.matches = [cms.InputTag('muonTriggerMatchHLTMuons')]
    for x in [process.cleanPatPhotonsTriggerMatch,
              process.cleanPatElectronsTriggerMatch,
              process.cleanPatTausTriggerMatch,
              process.cleanPatJetsTriggerMatch,
              process.patMETsTriggerMatch]:
        process.patTriggerMatchEmbedder.remove(x)

def switchHLTProcessName(process, name):
    process.patTrigger.processName = name
    process.patTriggerEvent.processName = name

def addHEEPId(process):
    # Run the HEEP electron id. This must be done at PAT tuple making time
    # and cannot be done later unless some modifications are done the
    # GsfElectron/GsfElectronCore classes.
    from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import heepBarrelCuts, heepEndcapCuts
    process.HEEPId = cms.EDProducer('HEEPIdValueMapProducer',
                                    eleLabel = cms.InputTag('gsfElectrons'),
                                    barrelCuts = heepBarrelCuts,
                                    endcapCuts = heepEndcapCuts
                                    )

    # Embed the HEEP cut bitwords into the userData of the patElectrons.
    process.patElectrons.userData.userInts.src.append('HEEPId')
    
    process.patDefaultSequence.replace(process.patCandidates, process.HEEPId * process.patCandidates)

def AODOnly(process):
    raise NotImplementedError('JMTBAD')

    from PhysicsTools.PatAlgos.tools.coreTools import restrictInputToAOD
    restrictInputToAOD(process)

    removeMuonMCClassification(process) # throw the baby out with the bathwater...

    for x in (process.muonMatch, process.electronMatch):
        x.matched = cms.InputTag('prunedGenLeptons')
    process.patDefaultSequence.replace(process.genSimLeptons * process.prunedGenSimLeptons, process.prunedGenLeptons)


# Some scraps to aid in debugging that can be put in your top-level
# config (could be turned into functions a la the above):
'''
process.patDefaultSequence.remove(process.selectedPatCandidateSummary)
process.selectedPatCandidateSummary.perEvent = cms.untracked.bool(True)
process.selectedPatCandidateSummary.dumpItems = cms.untracked.bool(True)
process.patDefaultSequence *= process.selectedPatCandidateSummary


process.Tracer = cms.Service('Tracer')
process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck')


process.eca = cms.EDAnalyzer('EventContentAnalyzer')
process.peca = cms.Path(process.eca)


process.load('L1Trigger.L1ExtraFromDigis.l1extratest_cfi')
process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
#process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD', '', hltProcessName)
process.ptrigAnalyzer = cms.Path(process.l1extratest*process.triggerSummaryAnalyzerAOD)


process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.printTree = cms.EDAnalyzer(
    'ParticleListDrawer',
    maxEventsToPrint = cms.untracked.int32(-1),
    src = cms.InputTag('genParticles'),
    printOnlyHardInteraction = cms.untracked.bool(True),
    useMessageLogger = cms.untracked.bool(True)
    )
process.MessageLogger.categories.append('ParticleListDrawer')
process.ptree = cms.Path(process.printTree)


process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


process.MessageLogger.cerr.FwkReport.reportEvery = 1
'''
