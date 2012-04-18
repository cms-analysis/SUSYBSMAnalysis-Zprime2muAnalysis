#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

def pruneMCLeptons(process, use_sim=False):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
    obj = process.prunedMCLeptons
    
    if use_sim:
        # For muon and electron MC matching, want to be able to match to
        # decays-in-flight produced in SIM (whether by GEANT or FastSim),
        # so make some genParticles out of the simTracks.
        process.load('SUSYBSMAnalysis.Zprime2muAnalysis.GenPlusSim_cfi')
        process.prunedMCLeptons.src = 'genSimLeptons'
        obj = process.genSimLeptons * process.prunedMCLeptons

    for x in (process.muonMatch, process.electronMatch):
        # Switch to Use the new GEN+SIM particles created above.
        x.matched = cms.InputTag('prunedMCLeptons')
        
        # Default PAT muon/electron-MC matching requires, in addition
        # to deltaR < 0.5, the MC and reconstructed leptons to have
        # the same charge, and (reco pt - gen pt)/gen pt <
        # 0.5. Disable these two cuts.
        x.checkCharge = False
        x.maxDPtRel = 1e6

    process.patDefaultSequence = cms.Sequence(obj * process.patDefaultSequence._seq)

def addMuonMCClassification(process):
    # Run and embed in the patMuons the classification of muons by
    # their GEANT hits.
    process.load('MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi')
    from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
    addClassByHits(process.patMuons, extraInfo=True)
    process.patDefaultSequence = cms.Sequence(process.muonClassificationByHits * process.patDefaultSequence._seq)

def removeMuonMCClassification(process):
    # Some magic to undo the use of MC added above in
    # addMuonMCClassification.
    
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

def removeSimLeptons(process):
    if hasattr(process, 'genSimLeptons'):
        process.patDefaultSequence.remove(process.genSimLeptons)
        if hasattr(process, 'prunedMCLeptons'):
            process.prunedMCLeptons.src = 'genParticles'

def removePrunedMCLeptons(process):
    if hasattr(process, 'prunedMCLeptons'):
        process.patDefaultSequence.remove(process.prunedMCLeptons)

def removeMCUse(process):
    # Remove anything that requires MC truth.
    from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
    removeMCMatching(process, ['All'])
    removeMCMatching(process, ['METs'], postfix='TC')
    removeMCMatching(process, ['METs'], postfix='PF')
    removeMuonMCClassification(process)
    removeSimLeptons(process)
    removePrunedMCLeptons(process)
    
def switchHLTProcessName(process, name):
    # As the correct trigger process name is different from the
    # default "HLT" for some MC samples, this is a simple tool to
    # switch this in all places that it's needed.
    process.patTrigger.processName = name
    process.patTriggerEvent.processName = name

def addHEEPId(process):
    # Run the HEEP electron id. This must be done at PAT tuple making
    # time and cannot be done later unless some modifications are done
    # the GsfElectron/GsfElectronCore classes.
    from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import heepBarrelCuts, heepEndcapCuts
    from SHarper.HEEPAnalyzer.HEEPEventParameters_cfi import heepEventPara
    process.HEEPId = cms.EDProducer('HEEPIdValueMapProducer',
                                    eleLabel = cms.InputTag('gsfElectrons'),
                                    barrelCuts = heepBarrelCuts,
                                    endcapCuts = heepEndcapCuts,
                                    eleIsolEffectiveAreas = heepEventPara.eleIsolEffectiveAreas,
                                    applyRhoCorrToEleIsol = heepEventPara.applyRhoCorrToEleIsol,
                                    eleRhoCorrLabel = heepEventPara.eleRhoCorrTag,
                                    writeIdAsInt = cms.bool(True),
                                    )

    # Embed the HEEP cut bitwords into the userData of the
    # patElectrons so we can use it to cut on at dilepton-making time
    # instead of necessarily dropping them in the selectedPatElectrons
    # or cleanPatElectrons steps.
    process.patElectrons.userData.userInts.src.append('HEEPId')
    
    process.patDefaultSequence.replace(process.patCandidates, process.HEEPId * process.patCandidates)

def AODOnly(process):
    from PhysicsTools.PatAlgos.tools.coreTools import restrictInputToAOD
    restrictInputToAOD(process)
    removeMuonMCClassification(process) # throw the baby out with the bathwater...
    removeSimLeptons(process)

# Some scraps to aid in debugging that can be put in your top-level
# config (could be turned into functions a la the above):
'''
# At the end of the job, print out a table summarizing the PAT
# candidates seen/made.
process.patDefaultSequence.remove(process.selectedPatCandidateSummary)
process.selectedPatCandidateSummary.perEvent = cms.untracked.bool(True)
process.selectedPatCandidateSummary.dumpItems = cms.untracked.bool(True)
process.patDefaultSequence *= process.selectedPatCandidateSummary

# Print messages tracing through the execution of the
# analyzers/producers (e.g. beginJob/beginRun/analyze/endRun/endJob).
process.Tracer = cms.Service('Tracer')
process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck')

# To print every event the content (the branch names) of the
# edm::Event.
process.eca = cms.EDAnalyzer('EventContentAnalyzer')
process.peca = cms.Path(process.eca)

# Dump extensive L1 and HLT trigger info (objects, path results).
process.load('L1Trigger.L1ExtraFromDigis.l1extratest_cfi')
process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
#process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD', '', hltProcessName)
process.ptrigAnalyzer = cms.Path(process.l1extratest*process.triggerSummaryAnalyzerAOD)

# Dump the list of genParticles in a format similar to that from
# turning on PYTHIA's verbosity.
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

# Print tables of the results of module/path execution and timing info.
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Make MessageLogger print a message every event with (run, lumi,
# event) numbers.
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Extra options for controlling how CMSSW works.
process.source.noEventSort = cms.untracked.bool(True)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

'''
