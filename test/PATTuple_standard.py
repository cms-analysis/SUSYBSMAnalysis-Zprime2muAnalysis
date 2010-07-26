## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

## ------------------------------------------------------
#  NOTE: you can make some typical event selection cuts
#  such as good vertex, beam scraping, and L1 technical bits
#  by uncommenting the lines below and including it in
#  the path.
## ------------------------------------------------------
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
## Modify L1 trigger bits if MC (bit 0 not defined for MC)
process.L1T1.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

## ------------------------------------------------------
#  NOTE: you can use a bunch of core tools of PAT
#  and some Exotica tools to
#  taylor your PAT configuration; for a few examples
#  uncomment the lines below
## ------------------------------------------------------
from PhysicsTools.PatAlgos.tools.coreTools import *
## restrict the input file to run off of AOD
#restrictInputToAOD(process)


from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime_PATTools import *
## configure the PAT to be run with or without MC information
runMC(process,True)

## extend the patMuon to have the embedding that we want
extendMuon(process)

## remove certain objects from the default sequence
#removeAllPATObjectsBut(process, ['Muons','Electrons'])
#removeSpecificPATObjects(process, ['Taus'])

## ------------------------------------------------------
#  Add collections that we want to keep on top of the
#  standard event content.
## ------------------------------------------------------
from SUSYBSMAnalysis.Zprime2muAnalysis.EventContent_cff import zPrimeEventContent
process.out.outputCommands += zPrimeEventContent



## ------------------------------------------------------
#  NOTE: you can still run PAT in the 36X version on
#  input files produced within the 35X series. This
#  implies some reconfigurations, example are given
#  below.
## ------------------------------------------------------
#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *

## uncomment this line to run on an 35X input sample
#run36xOn35xInput(process)

## let it run
process.p = cms.Path(
    process.goodData *
    process.patDefaultSequence
    )

## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
process.source.fileNames = [
    '/store/mc/Summer10/Zmumu/GEN-SIM-RECO/START36_V9_S09-v1/0045/FC750941-147B-DF11-BA8C-00261834B569.root',
    ]
#
#process.GlobalTag.globaltag = 'START36_V10::All'
process.GlobalTag.globaltag = 'MC_37Y_V5::All'
#
process.maxEvents.input = 1000        ##  (e.g. -1 to run on all events)
#
#   process.out.outputCommands = [ ... ]  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#                                         ##
#   process.out.fileName = ...            ##  (e.g. 'myTuple.root')
#                                         ##
#   process.options.wantSummary = True    ##  (to suppress the long output at the end of the job)    
#
## switch on PAT trigger
changeMuonHLTMatch( process )
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
process.out.outputCommands += patTriggerEventContent
