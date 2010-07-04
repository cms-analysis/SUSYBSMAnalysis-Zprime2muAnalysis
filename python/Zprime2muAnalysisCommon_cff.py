import os
import FWCore.ParameterSet.Config as cms

from tools import parse_enum, rec_levels

__debug = False

# Info about the defined rec levels.
recLevelDict, recLevels = rec_levels()
numRecLevels = len(recLevels)
# Make constants for the reclevels to avoid magic numbers for them below.
for enum, value in recLevelDict.items():
    exec '%s = %i' % (enum, value)

# The InputTag names for the muon collections.
muonCollections = ['muCand%s' % rec for rec in recLevels]

# The InputTag names for the electron collections. (None specifies
# that that rec level is to be skipped, especially during dilepton
# construction since the CandCombiner modules throw exceptions if the
# collection is not found.)
electronCollections = [c.replace('mu','el') for c in muonCollections]
for rec in xrange(lGR, len(electronCollections)):
    electronCollections[rec] = None
electronCollections[lOP] = 'heepElectrons'

# The InputTags for the L1 collections (needs to be configurable to
# support FastSim).
l1InputLabel = 'l1extraParticles' # 'l1extraParticles' for FastSim
l1MapLabel   = 'hltL1GtObjectMap'    # 'gtDigis'           "     "

# Dilepton construction specifiers (see docstring in the main method
# make...() below).
oppSign = ('+','-')
oppSignMP = ('-','+') # not really necessary
likeSignPos = ('+','+')
likeSignNeg = ('-','-')
diMuons = ('muons', 'muons')
diElectrons = ('electrons', 'electrons')
muonElectron = ('muons', 'electrons')
electronMuon = ('electrons', 'muons')

# Lepton/dilepton cut specifiers, for inclusion in a bitmask.  Magic
# happens here to get the cut bit values directly from CutHelper.h so
# we don't have to worry about keeping them consistent between the two
# files.
cuts = parse_enum('src/SUSYBSMAnalysis/Zprime2muAnalysis/src/CutHelper.h', 'CutResult')

# Now use the cuts so defined and make some sets of them:
# AN 2007/038 cuts (pT > 20 and isolation dR sumPt < 10)
cuts['TeVmu'] = cuts['PT'] | cuts['ISO'];
# AN 2008/044 e mu bkg study cuts (pT > 80).
cuts['HEEP']  = cuts['PTOTHER'];
# Various combinations of AN 2008/015 cuts ("top group").
cuts['TopSkim']      = cuts['PT'] # | cuts['TOPTRIGGER']
cuts['TopEl']        = cuts['ELTIGHT'] | cuts['D0EL'] | cuts['COLLEMU']
cuts['TopMu']        = cuts['D0'] | cuts['NSIHITS'] | cuts['CHI2DOF']
cuts['TopLeptonId']  = cuts['TopEl'] | cuts['TopMu']
cuts['TopIsolation'] = cuts['ISOS']
cuts['Top']          = cuts['TopSkim'] | cuts['TopLeptonId'] | cuts['TopIsolation']

# We want to be able to add modules to the path in loops, with their
# names depending on the loop variables.  Keep track of the modules we
# make so we can set them all into a real CMSSW path in the process
# later. The order of the modules in this list will be obeyed by
# CMSSW.
__path = None
def pathAppend(prodobj):
    global __path
    if __path is None: __path  = prodobj
    else:              __path *= prodobj
    
def addToPath(process, name, prodobj):
    pathAppend(prodobj)
    setattr(process, name, prodobj)

def finalizePath(process, pathName):
    global __path
    if __path is not None:
        p = cms.Path(__path)
        __path = None
        setattr(process, pathName, p)

def makeZprime2muAnalysisProcess(fileNames=[],
                                 maxEvents=-1,
                                 doingElectrons=False,
                                 useGen=True,
                                 useSim=True,
                                 useRaw=True,
                                 useReco=True,
                                 usingAODOnly=False,
                                 useTrigger=True,
                                 useOtherMuonRecos=True,
                                 recoverBrem=True,
                                 disableElectrons=True,
                                 conditionsGlobalTag='MC_3XY_V27::All',
                                 dumpHardInteraction=False,
                                 strictGenDilepton=True,
                                 dumpTriggerSummary=False,
                                 flavorsForDileptons=diMuons,
                                 chargesForDileptons=oppSign,
                                 maxDileptons=1,
                                 skipPAT=True,
                                 runOnPATTuple=False,
                                 useHLTDEBUG=False,
                                 minGenPt=3.0,
                                 cutMask=cuts['TeVmu'],
                                 bestRecLevel=lOP,
                                 inputIsFastSim=False,
                                 muons=muonCollections,
                                 electrons=electronCollections,
                                 photons='photons',
                                 defaultMuons='muons',
                                 tevMuons='tevMuons',
                                 hltProcessName='HLT',
                                 processName='Zprime2muAnalysis',
                                 fromNtuple=False,
                                 processObj=None):
    '''Return a CMSSW process for running Zprime2muAnalysis-derived
    code. See e.g. testZprime2muResolution_cfg.py for example use.

    Arguments:
    
    fileNames: list of strings containing filenames for PoolSource. If
    left as [] (i.e. nothing passed in), we assume the user is going
    to set up process.source explicitly.
    
    maxEvents: number of events to run (same convention as in the
    original configuration files). Default is all events in all files.

      Comment: even though we set up process.source,
      process.maxEvents, process.MessageLogger, etc., once the caller
      has the process object, any of these can be replaced. The only
      thing to be careful of when replacing modules is that the order
      with which paths are put into the process matters, i.e. one
      cannot simply replace any of the paths below after they are
      already made and get the modules to run in the right
      order. Other than paths everything should be replaceable.
    
    doingElectrons: whether to run on electron or muon collections
    (i.e. which collections are in allLeptons in the code). Default is
    False, so we run on muons.
    
    useGen/Sim/Trigger/Reco: whether to expect information from GEN,
    SIM, HLT, or RECO branches. (The expectation from RECO is modified
    by the usingAODOnly parameter below.)
    
    usingAODOnly: if True, do not expect to be able to get information
    that is not in AOD, such as rechits. Also turns off looking at
    GEN/SIM collections.

    useOtherMuonRecos: whether to expect the extra TeV muon tracks and
    the track-to-track maps for them in the input files. If True, then
    the OP rec level is a copy of GR, and the TR level is a copy of
    TK.

      Comment: the use* flags listed above are merely shortcuts used
      to disable the appropriate rec levels. They are used in the
      "sanity check" section below to set in the muons/electrons list
      the appropriate rec levels to None if they are to be
      disabled. After that, they are not checked directly, but rather
      the individual rec levels themselves are checked for None-ness.
      
    recoverBrem: whether to expect a collection of reconstructed
    photons in the input files in order to construct "resonances",
    i.e. add closest photons to the dimuons made.

    disableElectrons: whether to completely kill electrons from being
    accessed.
            
    conditionsGlobalTag: this sets the globalTag for database
    conditions (e.g. alignment, calibration, etc.) needed for trigger
    redigitization/l1extra production, or for the PAT.

    dumpHardInteraction: if True, use ParticleListDrawer to dump the
    generator info on the hard interaction particles.

    strictGenDilepton: if True, require the dilepton produced at gen
    level to have leptons that have the same generated mother as
    specified in the MC record (or, in the case of leptons that have
    radiated photons, try to find the first non-lepton ancestor of
    each and compare those). Set this to False for CompHEP or other
    samples which do not have the resonance listed in the MC record;
    the dilepton will then be constructed using final state muons in
    the same way as for all the other rec levels.

    dumpTriggerSummary: if True, use TriggerSummaryAnalyzerAOD to dump
    the information in the hltTriggerSummaryAOD object (which is the
    input the L2 and L3 collections when not running on HLTDEBUG
    files).
    
    flavorsForDileptons: the pair of collections to be used for the
    "official" dileptons (i.e. those accessed by allDileptons in the
    code); one of diMuons, diElectrons, etc. above.
    
    chargesForDileptons: the pair of charges (in corresponding order
    to the flavors) to be used for the "official" dileptons (see
    previous); one of oppSign, etc. above.

      Example: if flavorsForDileptons == muonElectron and
      chargesForDileptons == oppSign, then the official dileptons
      will be all pairs of mu+ e-.

    maxDileptons: the maximum number of dileptons that are kept in the
    collections after cuts, sorting, and overlap removal.

    skipPAT: whether to skip running layers 0 and 1 of the PAT. Useful
    if running on files that already have PAT collections in them.

    useHLTDEBUG: can be set True if running on HLTDEBUG files, in
    which case the L2 and L3 collections are taken directly from the
    appropriate collections in the event. Otherwise, L2 and L3
    collections are produced from the hltTriggerSummaryAOD object.

    minGenPt: a simple cut on generator-level pT to cut out a lot of
    leptons that would never have a chance of being reconstructed. We
    want to keep high-pT leptons from in-flight decays, so we add in
    the GEANT leptons, but without a pT cut lots of junk leptons would
    be added.

    cutMask: an unsigned value representing the chosen cuts on
    leptons/dileptons. The values must correspond to those set in
    src/CutHelper.h. The "cuts" dictionary defined at the top of this
    file gives some examples of common ones. The default is the
    "TeVmu" selection: lepton pT > 20 GeV, and tracker isolation sumPt
    < 10 GeV.

    inputIsFastSim: whether the input file has been produced by
    FastSim or not. If so, then some collection names need to be
    changed (the L1/HLT collections and trigger bitmaps), and TeV
    muons should be turned off since they are not produced (this
    latter part will be fixed in 3_1_X).
    
    muons/electrons/photons/defaultMuons/tevMuons: for experts who
    wish to directly specify the collections to be used for each rec
    level.

    hltProcessName: now that HLT is typically run twice for 1E31 and
    8E29 menus and the results for both are kept, this specifies which
    is to be used (e.g. "HLT" or "HLT8E29").
    '''

    ####################################################################
    ## Sanity checks.
    ####################################################################

    if len(muons) != numRecLevels or len(electrons) != numRecLevels:
        raise RuntimeError, 'at least one of the muon and electron collections is not the right length'

    if runOnPATTuple:
        skipPAT = True
        disableElectrons = True
        useSim = False
        useOtherMuonRecos = False
        photons = 'selectedLayer1Photons'
        defaultMuons = 'selectedLayer1Muons'

    if inputIsFastSim:
        useOtherMuonRecos = False
        global l1InputLabel
        global l1MapLabel
        l1InputLabel = 'l1extraParticles'
        l1MapLabel = 'gtDigis'
        
    if disableElectrons:
        for i in xrange(numRecLevels):
            electrons[i] = None
        
    if usingAODOnly:
        useGen = useSim = useRaw = False
        
    if not useGen:
        useSim = False
        muons[lGN] = electrons[lGN] = None

    if not useReco:
        useTrigger = False
        useOtherMuonRecos = False
        for i in xrange(lL1, lOP+1):
            muons[i] = None
            electrons[i] = None

    if not useOtherMuonRecos:
        for i in xrange(lTK+1, lPR+1):
            muons[i] = None
        muons[lOP] = muons[lGR]
        muons[lTR] = muons[lTK]

    if not useTrigger:
        for i in xrange(lL1, lL3+1):
            muons[i] = electrons[i] = None

    print 'After sanity checks, using these muon collections:', muons
    print 'And these electron collections:', electrons
    
    if flavorsForDileptons not in [diMuons, diElectrons, muonElectron, electronMuon]:
        raise RuntimeError, 'bad input for flavorsForDileptons'
    
    if chargesForDileptons not in [oppSign, oppSignMP, likeSignPos, likeSignNeg]:
        raise RuntimeError, 'bad input for chargesForDileptons'        

    need_global_tag = False
    need_l1extra = False
    
    if 'muCandL1' in muons or 'elCandL1' in electrons:
        need_global_tag = True
        need_l1extra = True

    if processObj is None:
        ####################################################################
        ## Set up the CMSSW process and useful services.
        ####################################################################

        if fromNtuple and processName == 'Zprime2muAnalysis':
            processName += '2'

        process = cms.Process(processName)

        if fileNames:
            process.source = cms.Source(
                'PoolSource',
                fileNames = cms.untracked.vstring(*fileNames)
                )
            print 'Setting noDuplicateCheck mode...'
            process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

        process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))

        process.options = cms.untracked.PSet(
            #IgnoreCompletely = cms.untracked.vstring('ProductNotFound')
            )

        process.load('FWCore.MessageLogger.MessageLogger_cfi')
        process.MessageLogger.cerr.FwkReport.reportEvery = 200

        if __debug:
            process.options.wantSummary = cms.untracked.bool(True)
            process.Tracer = cms.Service('Tracer')
            process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck')

        process.TFileService = cms.Service('TFileService', fileName = cms.string('zp2mu_histos.root'))

        # If doing reconstruction or the PAT, need some common things
        # like geometry, magnetic field, and conditions.
        if not skipPAT:
            process.load("Configuration.StandardSequences.Geometry_cff")
            process.load("Configuration.StandardSequences.MagneticField_cff")
            need_global_tag = True

        if need_global_tag:
            process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
            process.GlobalTag.globaltag = cms.string(conditionsGlobalTag)

        if useGen and dumpHardInteraction:
            process.include("SimGeneral/HepPDTESSource/data/pythiapdt.cfi")
            process.printTree = cms.EDAnalyzer(
                'ParticleListDrawer',
                maxEventsToPrint = cms.untracked.int32(-1),
                src = cms.InputTag('genParticles'),
                printOnlyHardInteraction = cms.untracked.bool(True),
                useMessageLogger = cms.untracked.bool(True)
                )
            process.ptree = cms.Path(process.printTree)

        if useTrigger and dumpTriggerSummary:
            process.load('L1Trigger.L1ExtraFromDigis.l1extratest_cfi')
            process.trigAnalyzer = cms.EDAnalyzer(
                'TriggerSummaryAnalyzerAOD',
                inputTag = cms.InputTag('hltTriggerSummaryAOD', '', hltProcessName)
                )
            process.ptrigAnalyzer = cms.Path(process.l1extratest*process.trigAnalyzer)

        ####################################################################
        ## Run PAT layers 0 and 1 (unless instructed to skip them).
        ####################################################################

        if not skipPAT:
            process.load("PhysicsTools.PatAlgos.patLayer0_cff")
            process.load("PhysicsTools.PatAlgos.patLayer1_cff")
            process.pPAT = cms.Path(process.patLayer0 + process.patLayer1)
    else:
        process = processObj

    # Since some version of CMSSW 3, l1extraParticles are either not
    # made or are transient. Make them.
    if False and useRaw and need_l1extra:
        # For the samples on CASTOR in directory ~tucker/CMSSW_3_1_2,
        # the digitization step in the RECO sequence (the digis that
        # are extant in these files) was done using the 8E29 menu. If
        # we're not using that menu, redo digitization (raw for both
        # HLT and HLT8E29 are in these files).
        if 'STARTUP' not in conditionsGlobalTag:
            process.load("Configuration.StandardSequences.Geometry_cff")
            process.load('Configuration.StandardSequences.RawToDigi_cff')
            digitag = 'rawDataCollector::HLT'
            process.scalersRawToDigi.scalersInputTag = digitag
            process.csctfDigis.producer = digitag
            process.dttfDigis.DTTF_FED_Source = digitag
            process.gctDigis.inputLabel = digitag
            process.gtDigis.DaqGtInputTag = digitag
            process.siPixelDigis.InputLabel = digitag
            # ecal digi producer can't take an InputTag, and using the
            # string in the ::HLT form above just makes it throw an
            # exception (it doesn't see that HLT is supposed to be the
            # process name, instead taking ::HLT as part of the module
            # label.)
            #process.ecalDigis.InputLabel = digitag 
            process.ecalPreshowerDigis.sourceTag = digitag
            process.hcalDigis.InputLabel = digitag
            process.muonCSCDigis.InputObjects = digitag
            process.muonDTDigis.inputLabel = digitag
            process.muonRPCDigis.InputLabel = digitag
            process.gtEvmDigis.EvmGtInputTag = digitag
            process.praw2digi = cms.Path(process.RawToDigi)

        process.load('L1Trigger.Configuration.L1Extra_cff')
        process.pl1extra = cms.Path(process.L1Extra)


    ####################################################################
    ## Make a CandidateCollection out of the GEANT tracks.
    ####################################################################

    if useSim:
        process.simParticleCandidates = cms.EDProducer(
            'PATGenCandsFromSimTracksProducer',
            src = cms.InputTag('g4SimHits'),
            setStatus = cms.int32(1)
            )
        
        process.psimParticleCandidates = cms.Path(process.simParticleCandidates)

    ####################################################################
    ## Below we want to make the EDProducer object and add it to a
    ## path only if the collection is declared to be in use by being
    ## in the muons or electrons list (these collections define the
    ## rec levels). To simplify the logic (remove a lot of "if 'X' in
    ## muons" calls and indentation), make a helper function to
    ## encapsulate this functionality.
    ####################################################################

    EDProducerType = type(cms.EDProducer('dummy'))
    
    def appendIfUsing(name, *args, **kwargs):
        if name in muons or name in electrons:
            if len(args) == 1 and type(args[0]) == EDProducerType:
                # If the first argument after name is an EDProducer
                # object, go ahead and use it as passed in. This
                # enables using cms.EDProducer.clone().
                prodobj = args[0]
            else:
                # Otherwise, make a new EDProducer object, passing the
                # rest of the positional and keyword arguments to it.
                # *args, **kwargs is a python idiom.
                prodobj = cms.EDProducer(*args, **kwargs)
            addToPath(process, name, prodobj)
    
    ####################################################################
    ## Produce/clone/copy/do whatever to the needed muon collections.
    ####################################################################

    if 'muCandGN' in muons:
        process.genMuons = cms.EDFilter(
            'CandViewSelector',
            src = cms.InputTag('genParticles'),
            cut = cms.string('abs(pdgId) = 13 & status = 1 & pt > %f' % minGenPt)
            )
        pathAppend(process.genMuons)
        GNtags = cms.VInputTag(cms.InputTag('genMuons'))
        
        if useSim:
            process.simMuons = cms.EDFilter(
                'CandViewSelector',
                src = cms.InputTag('simParticleCandidates'),
                cut = cms.string('abs(pdgId) = 13 & status = 1 & pt > %f' % minGenPt)
                )
            pathAppend(process.simMuons)
            GNtags.append(cms.InputTag('simMuons'))

        process.muCandGN = cms.EDProducer(
            'CandViewMerger',
            src = GNtags
            )
        pathAppend(process.muCandGN)

    
    # 'Sanitize' the L1 muons, i.e. shift their phi values from
    # the bin edge to the bin center (an offset of 0.0218 rad).
    appendIfUsing('muCandL1', 'L1MuonSanitizer',
                  src = cms.InputTag(l1InputLabel)
                  )

    if useHLTDEBUG:
        # 'Sanitize' the L2 muons, i.e. in case there is no
        # hltL2MuonCandidates collection in the event, put an empty
        # one in (otherwise just copy the existing one).
        appendIfUsing('muCandL2', 'L2MuonSanitizer',
                      src = cms.InputTag('hltL2MuonCandidates')
                      )

        # 'Sanitize' the L3 muons, i.e. make up reco::Muons from
        # reco::MuonTrackLinks (the hltL3MuonCandidates drop some of
        # the extra track information, while the
        # MuonTrackLinksCollection hltL3Muons appropriately has links
        # to all 3 tracks).
        appendIfUsing('muCandL3', 'L3MuonSanitizer',
                      src = cms.InputTag('hltL3Muons')
                      )
    else:
        # Extract the L2 and L3 muons from the hltTriggerSummaryAOD
        # object.
        appendIfUsing('muCandL2', 'HLTLeptonsFromTriggerEvent',
                      summary = cms.InputTag('hltTriggerSummaryAOD', '', hltProcessName),
                      leptons = cms.VInputTag(cms.InputTag('hltL2MuonCandidates', '', hltProcessName))
                      )

        appendIfUsing('muCandL3', 'HLTLeptonsFromTriggerEvent',
                      summary = cms.InputTag('hltTriggerSummaryAOD', '', hltProcessName),
                      leptons = cms.VInputTag(cms.InputTag('hltL3MuonCandidates', '', hltProcessName))
                      )

    from RecoMuon.MuonIdentification.refitMuons_cfi import refitMuons

    # Copy only the GlobalMuons from the default muon collection,
    # ignoring the TrackerMuons and CaloMuons for now. Since 3_1_0,
    # this is now what was called the sigma switch, which uses the
    # tracker-only track if both the tracker-only and global pT were
    # below 200 GeV, or if the q/p of the two tracks differ by more
    # than twice the tracker-only error on q/p.
    appendIfUsing('muCandGR', refitMuons.clone(
        src           = defaultMuons,
        fromCocktail  = False,
        tevMuonTracks = 'none'
        ))

    # Make tracker-only reco::Muons out of the tracker tracks in
    # the muons collection.
    appendIfUsing('muCandTK', refitMuons.clone(
        src              = defaultMuons,
        fromCocktail     = False,
        tevMuonTracks    = 'none',
        fromTrackerTrack = cms.bool(True)
        ))

    # Make first-muon-station (FMS) reco::Muons using the supplied
    # TeV refit tracks.
    appendIfUsing('muCandFS', refitMuons.clone(
        src           = defaultMuons,
        fromCocktail  = False,
        tevMuonTracks = tevMuons + ':firstHit'
        ))

    # Make picky-muon-reconstructor (PMR) reco::Muons using the
    # supplied TeV refit tracks.
    appendIfUsing('muCandPR', refitMuons.clone(
        src           = defaultMuons,
        fromCocktail  = False,
        tevMuonTracks = tevMuons + ':picky'
        ))

    # Use the official TeV muon cocktail code to pick the best
    # muons using the supplied TeV refit tracks.
    appendIfUsing('muCandOP', refitMuons.clone(
        src           = defaultMuons,
        fromCocktail  = True,
        tevMuonTracks = tevMuons
        ))

    # Use the TMR cocktail code to pick muons using the
    # supplied TeV refit tracks.
    appendIfUsing('muCandTR', refitMuons.clone(
        src           = defaultMuons,
        fromCocktail  = False,
        fromTMR       = True,
        tevMuonTracks = tevMuons,
        ))

    finalizePath(process, 'pMuons')

    ####################################################################
    ## Same for the electrons
    ####################################################################

    if 'elCandGN' in electrons:
        process.genElectrons = cms.EDFilter(
            'CandViewSelector',
            src = cms.InputTag('genParticles'),
            cut = cms.string('abs(pdgId) = 11 & status = 1 & pt > %f' % minGenPt)
            )
        pathAppend(process.genElectrons)
        GNtags = cms.VInputTag(cms.InputTag('genElectrons'))
        
        if useSim:
            process.simElectrons = cms.EDFilter(
                'CandViewSelector',
                src = cms.InputTag('simParticleCandidates'),
                cut = cms.string('abs(pdgId) = 11 & status = 1 & pt > %f' % minGenPt)
                )
            pathAppend(process.simElectrons)
            GNtags.append(cms.InputTag('simElectrons'))
    
        process.elCandGN = cms.EDProducer(
            'CandViewMerger',
            src = GNtags
            )
        pathAppend(process.elCandGN)

    appendIfUsing('elCandL1', 'CandViewMerger',
                  src = cms.VInputTag(cms.InputTag(l1InputLabel,'NonIsolated'), cms.InputTag(l1InputLabel,'Isolated'))
                  )

    if useHLTDEBUG:
        # For electrons, the HLT makes its decision based off of the
        # l1(Non)IsoRecoEcalCandidate and the
        # pixelMatchElectronsL1(Non)IsoForHLT collections. I choose to
        # call them L2 and L3 electrons here, respectively.

        # 'Sanitize' the L2 & L3 electrons, i.e. in case there is no
        # collection in the event, put an empty one in (otherwise just
        # copy the existing one).

        if 'elCandL2' in electrons:
            process.elCandL2NonIso = cms.EDProducer(
                'L2ElectronSanitizer',
                src = cms.InputTag('hltL1NonIsoRecoEcalCandidate')
                )
            
            process.elCandL2Iso = cms.EDProducer(
                'L2ElectronSanitizer',
                src = cms.InputTag('hltL1IsoRecoEcalCandidate')
                )

            process.elCandL2 = cms.EDProducer(
                'CandViewMerger',
                src = cms.VInputTag(cms.InputTag('elCandL2NonIso'), cms.InputTag('elCandL2Iso'))
                )
            
            pathAppend(process.elCandL2NonIso)
            pathAppend(process.elCandL2Iso)
            pathAppend(process.elCandL2)

        if 'elCandL3' in electrons:
            process.elCandL3NonIso = cms.EDProducer(
                'L3ElectronSanitizer',
                src = cms.InputTag('hltPixelMatchElectronsL1NonIso')
                )

            process.elCandL3Iso = cms.EDProducer(
                'L3ElectronSanitizer',
                src = cms.InputTag('hltPixelMatchElectronsL1Iso')
                )

            process.elCandL3 = cms.EDProducer(
                'CandViewMerger',
                src = cms.VInputTag(cms.InputTag('elCandL3NonIso'), cms.InputTag('elCandL3Iso'))
                )

            pathAppend(process.elCandL2NonIso)
            pathAppend(process.elCandL2Iso)
            pathAppend(process.elCandL3)
    else:
        appendIfUsing('elCandL2', 'HLTLeptonsFromTriggerEvent',
                      summary = cms.InputTag('hltTriggerSummaryAOD', '', hltProcessName),
                      leptons = cms.VInputTag(cms.InputTag('hltL1NonIsoRecoEcalCandidate', '', hltProcessName), cms.InputTag('hltL1IsoRecoEcalCandidate', '', hltProcessName))
                      )

        appendIfUsing('elCandL3', 'HLTLeptonsFromTriggerEvent',
                      summary = cms.InputTag('hltTriggerSummaryAOD', '', hltProcessName),
                      leptons = cms.VInputTag(cms.InputTag('hltPixelMatchElectronsL1NonIso', '', hltProcessName), cms.InputTag('hltPixelMatchElectronsL1Iso', '', hltProcessName))
                      )

    try:
        from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import heepBarrelCuts, heepEndcapCuts
        appendIfUsing('heepElectrons', 'SimpleHEEPSelector',
                      src = cms.InputTag('gsfElectrons'),
                      barrelCuts = heepBarrelCuts,
                      endcapCuts = heepEndcapCuts
                      )
    except ImportError:
        # If the user hasn't checked out HEEPAnalyzer, skip it.
        pass

    finalizePath(process, 'pElectrons')
    
    ####################################################################
    ## Common module configuration parameters.
    ####################################################################

    process.Zprime2muAnalysisCommon = cms.PSet(
        ################################################################
        ## Analysis configuration 
        ################################################################
        verbosity         = cms.untracked.int32(0),
        maxDileptons      = cms.uint32(maxDileptons),
        doingElectrons    = cms.bool(doingElectrons),
        useGen            = cms.bool(useGen),
        useSim            = cms.bool(useSim),
        useReco           = cms.bool(useReco),
        usingAODOnly      = cms.bool(usingAODOnly),
        useTrigger        = cms.bool(useTrigger),
        useOtherMuonRecos = cms.bool(useOtherMuonRecos),
        cutMask           = cms.uint32(cutMask),
        bestRecLevel      = cms.int32(bestRecLevel),
        dateHistograms    = cms.untracked.bool(True),

        ################################################################
        ## Input tags for trigger paths and particles.
        ################################################################
        l1GtObjectMap = cms.InputTag(l1MapLabel),
        # Every process puts a TriggerResults product into the event;
        # pick the HLT one.
        hltResults = cms.InputTag('TriggerResults', '', hltProcessName),
        # Trigger paths we use: the result is the OR of these. These
        # are the highest pT muon triggers in the 8E29 menu. The 1E31
        # menu has in addition L1_SingleMu10 which seeds HLT_Mu15.
        l1Paths = cms.vstring('L1_SingleMu7', 'L1_DoubleMu3'),
        hltPaths = cms.vstring('HLT_Mu9', 'HLT_DoubleMu3'),
        )

    if doingElectrons:
        # Need to find out what triggers HEEP is using -- was
        # L1_SingleEG15, HLT_EM80, HLT_EM200 in 2E30 menu in 2_X_Y.
        process.Zprime2muAnalysisCommon.l1Paths = cms.vstring('L1_SingleEG8')
        process.Zprime2muAnalysisCommon.hltPaths = cms.vstring('HLT_Ele15_SW_L1R')

    ####################################################################
    ## Input tags for leptons at the different rec levels, in order;
    ## filled below using the collections defined at the top of the
    ## file, and according to whether the default leptons are
    ## electrons or muons (determined by doingElectrons).
    ####################################################################
    if doingElectrons:
        lepcoll = electrons
    else:
        lepcoll = muons

    for i in xrange(numRecLevels):
        if lepcoll[i] is not None:
            tag = cms.InputTag(lepcoll[i])    
            setattr(process.Zprime2muAnalysisCommon, 'leptons' + recLevels[i], tag)

    ####################################################################
    ## Do all the closest-in-deltaR matching -- between each rec level
    ## and MC truth, and between each offline rec level and
    ## reconstructed photons.
    ####################################################################

    # Closest (deltaR) matching to MC truth.
    for irec in xrange(lL1, numRecLevels):
        jrec = lGN

        fromColl = lepcoll[irec]
        toColl   = lepcoll[jrec]
        if fromColl is None or toColl is None:
            continue

        prod = cms.EDProducer(
            'TrivialDeltaRViewMatcher',
            src     = cms.InputTag(fromColl),
            matched = cms.InputTag(toColl),
            distMin = cms.double(0.5)
            )
        addToPath(process, 'genMatch' + recLevels[irec], prod)

    finalizePath(process, 'genMatchPath')

    # Match all muons to their single closest photons (within dR <
    # 0.1) to try and recover energy lost to brem later. (Don't bother
    # doing this for electrons, since the GSF algorithm already takes
    # brem losses into account.)
    if recoverBrem:
        for rec in xrange(lGR, numRecLevels):
            if muons[rec] is None: continue
            prod = cms.EDProducer(
                'TrivialDeltaRViewMatcher',
                src     = cms.InputTag(muons[rec]),
                matched = cms.InputTag(photons),
                distMin = cms.double(0.1)
                )
            addToPath(process, 'photonMatch' + recLevels[rec], prod)

        finalizePath(process, 'photonMatchPath')

    ####################################################################
    ## Dilepton construction. Make the requested dileptons, by default
    ## opposite-sign dimuons.
    ####################################################################

    for rec in xrange(numRecLevels):
        # Do not make dileptons for trigger levels.
        if rec >= lL1 and rec <= lL3:
            continue

        # Set up the flavor and charge pair for this rec level.
        collections = []
        for j in xrange(2):
            if flavorsForDileptons[j] == 'muons':
                collection_list = muons
            elif flavorsForDileptons[j] == 'electrons':
                collection_list = electrons
            collections.append(collection_list[rec])

        # We finally come to the point of the Nones in the electron
        # collection above: to skip rec levels that have no
        # collections here when producing dileptons. Otherwise,
        # CandCombiner throws an exception when it cannot find the
        # input collection.
        if None in collections:
            continue

        # Use CandViewShallowCloneCombiner for the ShallowClone bit;
        # this enables us to get the reference to the original lepton
        # that makes up the dilepton we're looking at in the code.
        combiner = 'CandViewShallowCloneCombiner'
        if rec == 0 and strictGenDilepton:
            # At generator level, look only at the leptons which came
            # from the resonance.
            combiner = 'GenDil' + combiner

        # Put the decay string in the format expected by CandCombiner:
        # e.g. 'muons@+ muons@-'.
        decay = '%s@%s %s@%s' % (collections[0], chargesForDileptons[0],
                                 collections[1], chargesForDileptons[1])

        # These are the "raw" unsorted, uncut dileptons...
        prod = cms.EDProducer(combiner,
                              decay = cms.string(decay),
                              cut = cms.string('')
                              )
        # ...except at generator level.
        if rec != lGN: name = 'rawDileptons'
        else:          name = 'dileptons'
        name += recLevels[rec]
        addToPath(process, name, prod)

        if rec != lGN:
            # Produce the dileptons after cuts and overlap removal.
            prod = cms.EDProducer('DileptonPicker',
                                  src          = cms.InputTag(name),
                                  doingGen     = cms.bool(rec == 0 and strictGenDilepton),
                                  maxDileptons = cms.uint32(maxDileptons),
                                  cutMask      = cms.uint32(cutMask)
                                  )
            addToPath(process, 'dileptons' + recLevels[rec], prod)

    finalizePath(process, 'dileptonPath')

    # Set up the "main" dileptons -- the ones the analysis module is
    # supposed to use by default.
    for rec in xrange(numRecLevels):
        # If we skipped making this dilepton collection in the logic
        # above, skip using it.
        tag = 'dileptons' + recLevels[rec]
        if hasattr(process, tag):
            setattr(process.Zprime2muAnalysisCommon, tag, cms.InputTag(tag))
      
    ####################################################################
    ## Done!
    ####################################################################

    if __debug:
        process.EventContentAnalyzer = cms.EDAnalyzer('EventContentAnalyzer')
        process.pECA = cms.Path(process.EventContentAnalyzer)        

    if fromNtuple:
        # Should probably just never make the modules and paths in the
        # above, but this works for now.
        paths_to_del = ['praw2digi','pl1extra','psimParticleCandidates','pMuons','genMatchPath','photonMatchPath','dileptonPath']
        for x in paths_to_del:
            if hasattr(process, x):
                delattr(process, x)

    print 'Zprime2muAnalysis base process built!'
    return process

########################################################################
## Utility functions.
########################################################################

# Function to attach simple EDAnalyzers that don't need any
# parameters. Allows skipping of making data/Zprime2mu*_cfi.py files.
def attachAnalysis(process, name, mod_name=None, **kwargs):
    if mod_name is None:
        mod_name = name
    analyzer = cms.EDAnalyzer(name, process.Zprime2muAnalysisCommon)
    setattr(process, mod_name, analyzer)
    setattr(process, mod_name + 'Path', cms.Path(analyzer))
    for key, val in kwargs.items():
        setattr(analyzer, key, val)

# Attach a PoolOutputModule to the process, by default dumping all
# branches (useful for inspecting the file).
def attachOutputModule(process, fileName='/scratchdisk2/tucker/zp2muout.root'):
    process.out = cms.OutputModule(
        'PoolOutputModule',
        fileName = cms.untracked.string(fileName)
        )
    process.endp = cms.EndPath(process.out)

def dumpNtuple(process, fileName):
    attachOutputModule(process, fileName)
    process.out.outputCommands = cms.untracked.vstring(
        'drop *',
        #'keep edmHepMCProduct_generator__HLT',
        'keep edmTriggerResults_TriggerResults__HLT',
        #'keep GenEventInfoProduct_generator__HLT',
        #'keep ints_genParticles__HLT',
        'keep L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT',
        'keep recoGenParticles_genParticles__HLT',
        #'keep SimTracks_g4SimHits__HLT',
        #'keep SimVertexs_g4SimHits__HLT',
        #'keep triggerTriggerEvent_hltTriggerSummaryAOD__HLT',
        'keep edmTriggerResults_TriggerResults__HLT8E29',
        #'keep triggerTriggerEvent_hltTriggerSummaryAOD__HLT8E29',
        'keep edmErrorSummaryEntrys_logErrorHarvester__*',
        'keep recoMuons_muons__*',
        'keep recoPhotonCores_photonCore__*',
        'keep recoPhotons_photons__*',
        'keep recoTrackExtras_generalTracks__*',
        'keep recoTrackExtras_globalMuons__*',
        'keep recoTrackExtras_standAloneMuons__*',
        'keep recoTrackExtras_tevMuons_default_*',
        'keep recoTrackExtras_tevMuons_firstHit_*',
        'keep recoTrackExtras_tevMuons_picky_*',
        'keep recoTracks_generalTracks__*',
        'keep recoTracks_globalMuons__*',
        'keep recoTracks_standAloneMuons_UpdatedAtVtx_*',
        'keep recoTracks_tevMuons_default_*',
        'keep recoTracks_tevMuons_firstHit_*',
        'keep recoTracks_tevMuons_picky_*',
        'keep recoTracksToOnerecoTracksAssociation_standAloneMuons__*',
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_default_*',
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_firstHit_*',
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_picky_*',
        'keep TrackingRecHitsOwned_generalTracks__*',
        'keep TrackingRecHitsOwned_globalMuons__*',
        'keep TrackingRecHitsOwned_standAloneMuons__*',
        'keep TrackingRecHitsOwned_tevMuons_default_*',
        'keep TrackingRecHitsOwned_tevMuons_firstHit_*',
        'keep TrackingRecHitsOwned_tevMuons_picky_*',
        #'keep edmTriggerResults_TriggerResults__Zprime2muAnalysis',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchFS__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchGR__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchL1__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchL2__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchL3__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchOP__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchPR__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchTK__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_genMatchTR__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_photonMatchFS__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_photonMatchGR__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_photonMatchOP__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_photonMatchPR__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_photonMatchTK__*',
        'keep recoCandidateedmViewrecoCandidateedmViewuintrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaseProdrecoCandidateedmRefToBaserecoCandidateedmRefToBaseedmOneToOneGenericedmAssociationMap_photonMatchTR__*',
        #'keep recoCandidatesOwned_genMuons__*',
        #'keep recoGenParticles_simParticleCandidates__*',
        #'keep recoCandidatesOwned_simMuons__*',
        'keep recoCompositeCandidates_dileptonsFS__*',
        'keep recoCompositeCandidates_dileptonsGN__*',
        'keep recoCompositeCandidates_dileptonsGR__*',
        'keep recoCompositeCandidates_dileptonsOP__*',
        'keep recoCompositeCandidates_dileptonsPR__*',
        'keep recoCompositeCandidates_dileptonsTK__*',
        'keep recoCompositeCandidates_dileptonsTR__*',
        #'keep recoCompositeCandidates_rawDileptonsFS__*',
        #'keep recoCompositeCandidates_rawDileptonsGR__*',
        #'keep recoCompositeCandidates_rawDileptonsOP__*',
        #'keep recoCompositeCandidates_rawDileptonsPR__*',
        #'keep recoCompositeCandidates_rawDileptonsTK__*',
        #'keep recoCompositeCandidates_rawDileptonsTR__*',
        'keep recoCandidatesOwned_muCandGN__*',
        'keep l1extraL1MuonParticles_muCandL1__*',
        'keep recoRecoChargedCandidates_muCandL2__*',
        'keep recoRecoChargedCandidates_muCandL3__*',
        'keep recoMuons_muCandFS__*',
        'keep recoMuons_muCandGR__*',
        'keep recoMuons_muCandOP__*',
        'keep recoMuons_muCandPR__*',
        'keep recoMuons_muCandTK__*',
        'keep recoMuons_muCandTR__*',
        )
    
# Return a list of files suitable for passing into a PoolSource,
# obtained by globbing using the pattern passed in (which includes the
# path to them).
def poolAllFiles(pattern):
    import glob
    return ['file:%s' % x for x in glob.glob(pattern)]

# Function to select a set of alignment constants. Useful when doing
# track re-reconstruction.
def setAlignment(process, name, connect, **kwargs):
    if len(kwargs) == 0:
        raise ValueError, 'must specify at least one record parameter'
    
    from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
    alignment_source = cms.ESSource('PoolDBESSource',
        CondDBSetup,
        connect = cms.string(connect),
        toGet = cms.VPSet()
    )

    for record, tag in kwargs.iteritems():
        alignment_source.toGet.append(cms.PSet(record = cms.string(record), tag = cms.string(tag)))

    setattr(process, name, alignment_source)
    process.prefer(name)

    #code = "cms.ESPrefer('PoolDBESSource', name, %s)" % ', '.join(['%s=cms.vstring("%s")' % x for x in kwargs.iteritems()])
    #setattr(process, name + '_es_prefer', eval(code))

def setTrackerAlignment(process, connectString, tagTrackerAlignmentRcd, tagTrackerAlignmentErrorRcd):
    setAlignment(process, 'trackerAlignment', connectString,
                 TrackerAlignmentRcd      = tagTrackerAlignmentRcd,
                 TrackerAlignmentErrorRcd = tagTrackerAlignmentErrorRcd)
                        
def setMuonAlignment(process, connectString, tagDTAlignmentRcd, tagDTAlignmentErrorRcd, tagCSCAlignmentRcd, tagCSCAlignmentErrorRcd):
    setAlignment(process, 'muonAlignment', connectString,
                 DTAlignmentRcd  = tagDTAlignmentRcd,  DTAlignmentErrorRcd  = tagDTAlignmentErrorRcd,
                 CSCAlignmentRcd = tagCSCAlignmentRcd, CSCAlignmentErrorRcd = tagCSCAlignmentErrorRcd)
    
################################################################################

# so 'from module import *' doesn't clutter the user's namespace
__all__ = [
    'muonCollections', 'electronCollections',
    'oppSign', 'oppSignMP', 'likeSignPos', 'likeSignNeg', 'diMuons',
    'diElectrons', 'muonElectron', 'electronMuon', 'cuts',
    'makeZprime2muAnalysisProcess', 'attachAnalysis',
    'attachOutputModule', 'dumpNtuple', 'poolAllFiles',
    'setAlignment', 'setTrackerAlignment', 'setMuonAlignment',
    'lGN','lL1','lL2','lL3','lGR','lTK','lFS','lPR','lOP','lTR'
    ]

# To debug this config, you can do from the shell:
#   python Zprime2muAnalysisCommon_cff.py
if __name__ == '__main__':
    import sys
    process = makeZprime2muAnalysisProcess(['file:/dev/null'], 42, skipPAT=True)
    if 'silent' not in sys.argv:
        print process.dumpConfig()
