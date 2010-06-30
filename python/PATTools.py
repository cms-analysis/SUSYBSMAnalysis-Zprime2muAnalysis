#!/usr/bin/env python

def configure(process, useAODOnly=False, muonsOnly=False, useMonteCarlo=False):
    from PhysicsTools.PatAlgos.tools.coreTools import restrictInputToAOD, removeAllPATObjectsBut, removeMCMatching
    if useAODOnly:
        restrictInputToAOD(process)
    if muonsOnly:
        removeAllPATObjectsBut(process, ['Muons'])
    if not useMonteCarlo:
        x = 'All'
        if muonsOnly:
            x = 'Muons'
        # Take the simParticles and genSimParticles modules out of the
        # sequence that we added.
        if hasattr(process, 'patDefaultSequenceBeforeZp2mu'):
            process.patDefaultSequence = process.patDefaultSequenceBeforeZp2mu
        removeMCMatching(process, [x])
        

__all__ = [
    'configure'
    ]

# Some scraps to aid in debugging that can be put in your top-level
# config:
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
process.ptrigAnalyzer = cms.Path(process.l1extratest*process.trigAnalyzer)


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
