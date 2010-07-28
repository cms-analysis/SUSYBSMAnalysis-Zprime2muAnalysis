from FWCore.GuiBrowsers.ConfigToolBase import *

from PhysicsTools.PatAlgos.tools.helpers import *

class AddGenSimTracks(ConfigToolBase):
    """ create GenParticles for all of the particles in a given
    simulator collection
    """
    _label_='addGenSimTracks'
    _defaultParameters=dicttypes.SortedKeysDict()
    def __init__(self):
        ConfigToolBase.__init__(self)
        self.addParameter(self._defaultParameters,'theSrc',cms.InputTag("g4SimHits"),"collection of SimTracks and SimVetices")
        self.addParameter(self._defaultParameters,'theFilter',"(abs(pdgId) == 13 || abs(pdgId) == 11) && pt > 2","default filter to use elctrons and muons")
        self.addParameter(self._defaultParameters,'theGenParticles',cms.InputTag("genParticles"),"original genParticle list")
        self.addParameter(self._defaultParameters,'theSetStatus',8,"set status = 8 for GEANT GPs")
        self._parameters=copy.deepcopy(self._defaultParameters)
        self._comment = ""

    def getDefaultParameters(self):
        return self._defaultParameters

    def __call__(self,process,
                 theSrc = None,
                 theFilter = None,
                 theSetStatus = None,
                 theGenParticles = None) :
        if  theSrc is None:
            theSrc=self._defaultParameters['theSrc'].value
        if  theFilter is None:
            theFilter=self._defaultParameters['theFilter'].value
        if  theGenParticles  is None:
            theGenParticles=self._defaultParameters['theGenParticles'].value
        if  theSetStatus is None:
            theSetStatus=self._defaultParameters['theSetStatus'].value
        self.setParameter('theSrc',theSrc)
        self.setParameter('theGenParticles',theGenParticles)
        self.setParameter('theFilter',theFilter)
        self.setParameter('theSetStatus',theSetStatus)
        
        self.apply(process) 

    def toolCode(self, process):        
        theSrc=self._parameters['theSrc'].value
        theGenParticles=self._parameters['theGenParticles'].value
        theFilter=self._parameters['theFilter'].value
        theSetStatus=self._parameters['theSetStatus'].value

        genSimLabel = "genParticles" + 'SimLep'
        print "Making new GenPlusSimParticles with label " + genSimLabel

        process.load("SUSYBSMAnalysis.Zprime2muAnalysis.GenPlusSim_cfi")
        setattr(process, genSimLabel, process.genSimLeptons.clone(src=theSrc))
        process.prunedGenParticles = cms.EDProducer(
            "GenParticlePruner",
            src = cms.InputTag("genParticlesSimLep"),
            select = cms.vstring(
            "drop * ",
            "keep status <= 3",
            "++keep++ pdgId = {mu-} && pt > 2",
            "++keep++ pdgId = {e-} && pt > 2",
            #"+drop pdgId = {mu-}",
            #"+drop pdgId = {e-}"
            )
            )
        process.patDefaultSequence.replace(process.patCandidates, (getattr(process,genSimLabel) + process.prunedGenParticles) + process.patCandidates)

addGenSimTracks=AddGenSimTracks()

def addGenSimParticles(process):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.GenPlusSim_cfi')
    process.electronMatch.matched = "genSimParticles"
    process.muonMatch.matched = "genSimParticles"
    process.tauMatch.matched = "genSimParticles"
    process.tauGenJets.GenParticles = "genSimParticles"
    process.photonMatch.matched = "genSimParticles"
    process.patJetPartonMatch.matched = "genSimParticles"
    ##process.patJetPartons.src = "genSimParticles"

def addMuonMCTruth(process):
    process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
    process.load("SimGeneral.MixingModule.mixNoPU_cfi")
    process.load("SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi")    # On RECO
    process.load("SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi")  # On RECO
    process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
    from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
    addClassByHits(process.patMuons, extraInfo=True)
    return(process)

def addMuonStations(process):
    process.load("MuonAnalysis.Examples.muonStations_cfi")
    from MuonAnalysis.Examples.muonStations_cfi import addUserData as addStations
    addStations(process.patMuons) 
    return(process)

def addMuonHitCount(process):
    process.load("UserCode.Examples.muonHitCount_cfi")
    from UserCode.Examples.muonHitCount_cfi import addUserData as addHitCount
    addHitCount(process.patMuons)
    
def runMC(process,useMonteCarlo=False):
    if useMonteCarlo:
        addGenSimParticles(process)
        addMuonMCTruth(process)
        process.patDefaultSequence.replace(process.patCandidates, getattr(process,'genSimMergerSequence') +  getattr(process,'muonClassificationByHits') + process.patCandidates)
        return(process)
    else :
        from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
        removeMCMatching(process, ['All'])
        return(process)

def extendMuon(process):
    process.patMuons.embedTrack = True
    process.selectedPatMuons.cut = 'isGlobalMuon || isTrackerMuon'
    process.countPatMuons.minNumber = 1
    addMuonStations(process)
    addMuonHitCount(process)
    process.patDefaultSequence.replace(process.patCandidates,getattr(process,'muonStations') + getattr(process,'muonHitCounts')  + process.patCandidates)
    return(process)

def changeMuonHLTMatch(process):
    process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
    process.load("SUSYBSMAnalysis.Zprime2muAnalysis.Zprime_hltTriggerMatch_cfi")
    process.patTriggerMatcher += process.muonTriggerMatchHLTMuons
    process.patTriggerMatcher.remove( process.patTriggerMatcherElectron )
    process.patTriggerMatcher.remove( process.patTriggerMatcherMuon )
    process.patTriggerMatcher.remove( process.patTriggerMatcherTau )
    process.patTriggerEvent.patTriggerMatches  = [ "muonTriggerMatchHLTMuons" ]
