import os
import FWCore.ParameterSet.Config as cms

__debug = False

# Info about the defined rec levels.
recLevels = ['GN','L1','L2','L3','GR','TK','FS','PR','OP','TR']
numRecLevels = len(recLevels)
# Make enums for the reclevels to avoid magic numbers below.
for i, rec in enumerate(recLevels):
    locals()['l' + rec] = i

# The InputTag names for the muon collections.
muonCollections = [
    'muCandGN',
    'muCandL1', # hltL1extraParticles
    'muCandL2', # hltL2MuonCandidates
    'muCandL3', # hltL3MuonCandidates
    'muCandGR',
    'muCandTK',
    'muCandFS',
    'muCandPR',
    'muCandOP',
    'muCandTR'
    ]

# The InputTag names for the electron collections.
#
# (None specifies that that rec level is to be skipped, especially
# during dilepton construction since the CandCombiner modules throw
# exceptions if the collection is not found.)
electronCollections = [
    'elCandGN',
    'elCandL1',
    'elCandL2',
    'elCandL3',
    'pixelMatchGsfElectrons',
    None,
    None,
    None,
    'selectedLayer1Electrons',
    None
    ]

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
cutHelperText = file(os.path.join(os.getenv('CMSSW_BASE'),
                                  'src/SUSYBSMAnalysis/Zprime2muAnalysis/src/CutHelper.h')).read().split('CutCutCut')[1]
cutHelperText = ''.join([line.split('//')[0].strip() for line in cutHelperText.split('\n')])
cutHelperText = [enum.strip().split('=') for enum in cutHelperText.split(',') if enum]
cutHelperText = [(x, eval(y)) for x, y in cutHelperText]
cuts = dict(cutHelperText)

# Now that that nastiness is over, use the cuts so defined and make
# some sets of cuts:
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
                                 useReco=True,
                                 usingAODOnly=False,
                                 useTrigger=True,
                                 useOtherMuonRecos=True,
                                 recoverBrem=True,
                                 disableElectrons=False,
                                 performTrackReReco=False,
                                 conditionsGlobalTag='IDEAL_V9::All',
                                 dumpHardInteraction=False,
                                 dumpTriggerSummary=False,
                                 flavorsForDileptons=diMuons,
                                 chargesForDileptons=oppSign,
                                 maxDileptons=1,
                                 skipPAT=False,
                                 runOnPATTuple=False,
                                 useHLTDEBUG=False,
                                 minGenPt=3.0,
                                 cutMask=cuts['TeVmu'],
                                 bestRecLevel=lOP,
                                 muons=muonCollections,
                                 electrons=electronCollections,
                                 photons='photons',
                                 defaultMuons='muons',
                                 tevMuons='tevMuons'):
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
    that is not in AOD, such as rechits.

    useOtherMuonRecos: whether to expect the extra TeV muon tracks and
    the track-to-track maps for them in the input files.

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

    performTrackReReco: if set, re-run track, muon, and photon
    reconstruction on-the-fly before running any analysis
    paths. Does this work in 2_1/2_2?
            
    conditionsGlobalTag: if performing re-reconstruction on the fly,
    this sets the globalTag for database conditions (e.g. alignment,
    calibration, etc.)

    dumpHardInteraction: if True, use ParticleListDrawer to dump the
    generator info on the hard interaction particles.

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

    muons/electrons/photons/defaultMuons/tevMuons: for experts who
    wish to directly specify the collections to be used for each rec
    level.
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

    if disableElectrons:
        for i in xrange(numRecLevels):
            electrons[i] = None
        
    if usingAODOnly:
        useGen = useSim = False
        
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

    if not useTrigger:
        for i in xrange(lL1, lL3+1):
            muons[i] = electrons[i] = None

    print 'After sanity checks, using these muon collections:', muons
    print 'And these electron collections:', electrons
    
    if flavorsForDileptons not in [diMuons, diElectrons, muonElectron, electronMuon]:
        raise RuntimeError, 'bad input for flavorsForDileptons'
    
    if chargesForDileptons not in [oppSign, oppSignMP, likeSignPos, likeSignNeg]:
        raise RuntimeError, 'bad input for chargesForDileptons'        
    
    ####################################################################
    ## Set up the CMSSW process and useful services.
    ####################################################################

    process = cms.Process('Zprime2muAnalysis')

    if fileNames:
        process.source = cms.Source(
            'PoolSource',
            fileNames = cms.untracked.vstring(*fileNames)
            )

    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))

    process.options = cms.untracked.PSet(
        #IgnoreCompletely = cms.untracked.vstring('ProductNotFound')
        )

    if __debug:
        process.options.wantSummary = cms.untracked.bool(True)
        outFile = 'cout'
    else:
        outFile = 'Zprime'
        
    process.MessageLogger = cms.Service(
        'MessageLogger',
        destinations = cms.untracked.vstring(outFile),
        categories = cms.untracked.vstring(
            #'FwkJob', 'FwkReport', 'Root_Warning',
            'Root_NoDictionary', 'RFIOFileDebug',
            'DDLParser', 'PixelGeom', 'EcalGeom', 'TIBGeom', 'TIDGeom',
            'TOBGeom', 'TECGeom', 'SFGeom', 'HCalGeom', 'TrackerGeom',
            'GeometryConfiguration', 'HcalHardcodeGeometry',
            'PoolDBESSource', 'TkDetLayers', 'TkNavigation',
            'Done', 'CSC', 'EcalTrivialConditionRetriever',
            'Geometry', 'GlobalMuonTrajectoryBuilder', 'HCAL', 'Muon',
            'RecoMuon', 'setEvent', 'Starting', 'TrackProducer',
            'trajectories', 'DetLayers', 'RadialStripTopology',
            'SiStripPedestalsFakeESSource', 'TrackingRegressionTest',
            'CaloExtractorByAssociator', 'CaloGeometryBuilder'
            'CompositeTrajectoryFilterESProducer', 'NavigationSetter',
            'SiStripPedestalsFakeESSource', 'TrajectoryFilterESProducer',
            'ZDC', 'ZdcHardcodeGeometry', 'RunLumiMerging',
            'TrackAssociator', #'RecoVertex/PrimaryVertexProducer',
            'ConversionTrackCandidateProducer','GsfTrackProducer',
            'PhotonProducer','TrackProducerWithSCAssociation',
            'PartonSelector', 'JetPartonMatcher', 'Alignments',
            'L1GtConfigProducers', 'SiStripQualityESProducer',
            'LikelihoodPdf', 'HemisphereAlgo', 'LikelihoodPdfProduct',
            'ObjectResolutionCalc'
            ),
        Zprime = cms.untracked.PSet(
            threshold    = cms.untracked.string('INFO'),
            lineLength   = cms.untracked.int32(132),
            noLineBreaks = cms.untracked.bool(True)
            )
        )

    if not __debug:
        process.MessageLogger.Zprime.extension = cms.untracked.string('.out')
        setattr(process.MessageLogger.Zprime, 'FwkReport', cms.untracked.PSet(reportEvery = cms.untracked.int32(500)))
        
        # Instead of line after line of limit psets in Zprime above, set
        # them all here.
        limitZero = cms.untracked.PSet(limit = cms.untracked.int32(0))
        for cat in process.MessageLogger.categories:
            setattr(process.MessageLogger.Zprime, cat, limitZero)

    process.TFileService = cms.Service(
        'TFileService',
        fileName = cms.string('zp2mu_histos.root')
        )
    
    if __debug:
        process.Tracer = cms.Service('Tracer')
        process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck')
    
    if performTrackReReco:
        process.include('SUSYBSMAnalysis/Zprime2muAnalysis/data/TrackReReco.cff')
        process.pReReco = cms.Path(process.trackReReco)

        process.include('Configuration/StandardSequences/data/FrontierConditions_GlobalTag.cff')
        process.GlobalTag.globaltag = conditionsGlobalTag

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
        process.trigAnalyzer = cms.EDAnalyzer(
            'TriggerSummaryAnalyzerAOD',
            inputTag = cms.InputTag('hltTriggerSummaryAOD')
            )
        process.ptrigAnalyzer = cms.Path(process.trigAnalyzer)

    ####################################################################
    ## Run PAT layers 0 and 1 (unless instructed to skip them).
    ####################################################################

    if not skipPAT:
        process.load("Configuration.StandardSequences.Geometry_cff")
        process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
        process.GlobalTag.globaltag = cms.string(conditionsGlobalTag)
        process.load("Configuration.StandardSequences.MagneticField_cff")

        process.load("PhysicsTools.PatAlgos.patLayer0_cff")
        process.load("PhysicsTools.PatAlgos.patLayer1_cff")

        process.pPAT = cms.Path(process.patLayer0 + process.patLayer1)

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
        process.genMuons = cms.EDProducer(
            'CandViewSelector',
            src = cms.InputTag('genParticles'),
            cut = cms.string('abs(pdgId) = 13 & status = 1 & pt > %f' % minGenPt)
            )
        pathAppend(process.genMuons)
        GNtags = cms.VInputTag(cms.InputTag('genMuons'))
        
        if useSim:
            process.simMuons = cms.EDProducer(
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
                  src = cms.InputTag('hltL1extraParticles')
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
                      src = cms.VInputTag(cms.InputTag('hltL2MuonCandidates::HLT'))
                      )

        appendIfUsing('muCandL3', 'HLTLeptonsFromTriggerEvent',
                      src = cms.VInputTag(cms.InputTag('hltL3MuonCandidates::HLT'))
                      )

    # Use the right collection of muons depending on whether we've
    # done track re-reconstruction.
    if performTrackReReco:
        defaultMuons += '2'
        tevMuons += '2'

    from RecoMuon.MuonIdentification.refitMuons_cfi import refitMuons

    # Copy only the GlobalMuons from the default muon collection,
    # ignoring the TrackerMuons and CaloMuons for now.
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
        process.genElectrons = cms.EDProducer(
            'CandViewSelector',
            src = cms.InputTag('genParticles'),
            cut = cms.string('abs(pdgId) = 11 & status = 1 & pt > %f' % minGenPt)
            )
        pathAppend(process.genElectrons)
        GNtags = cms.VInputTag(cms.InputTag('genElectrons'))
        
        if useSim:
            process.simElectrons = cms.EDProducer(
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
                  src = cms.VInputTag(cms.InputTag('hltL1extraParticles','NonIsolated'), cms.InputTag('hltL1extraParticles','Isolated'))
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
                      src = cms.VInputTag(cms.InputTag('hltL1NonIsoRecoEcalCandidate::HLT'), cms.InputTag('hltL1IsoRecoEcalCandidate::HLT'))
                      )

        appendIfUsing('elCandL3', 'HLTLeptonsFromTriggerEvent',
                      src = cms.VInputTag(cms.InputTag('hltPixelMatchElectronsL1NonIso::HLT'), cms.InputTag('hltPixelMatchElectronsL1Iso::HLT'))
                      )

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
        l1GtObjectMap = cms.InputTag('hltL1GtObjectMap'),
        # Every process puts a TriggerResults product into the event;
        # pick the HLT one.
        hltResults = cms.InputTag('TriggerResults','','HLT')
        )

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
        
        if doingElectrons:
            fromColl = electrons[irec]
            toColl   = electrons[jrec]
        else:
            fromColl = muons[irec]
            toColl   = muons[jrec]

        if fromColl is None or toColl is None:
            continue

        prod = cms.EDProducer(
            'TrivialDeltaRViewMatcher',
            src     = cms.InputTag(muons[irec]),
            matched = cms.InputTag(muons[jrec]),
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
        if rec == 0:
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

    print 'Zprime2muAnalysis base process built!'
    return process

########################################################################
## Utility functions.
########################################################################

# Function to attach simple EDAnalyzers that don't need any
# parameters. Allows skipping of making data/Zprime2mu*_cfi.py files.
def attachAnalysis(process, name, **kwargs):
    analyzer = cms.EDAnalyzer(name, process.Zprime2muAnalysisCommon)
    setattr(process, name, analyzer)
    setattr(process, name + 'Path', cms.Path(getattr(process, name)))
    for key, val in kwargs.items():
        setattr(analyzer, key, val)

# Attach a PoolOutputModule to the process, by default dumping all
# branches (useful for inspecting the file).
def attachOutputModule(process, fileName='/scratchdisk1/tucker/out.root'):
    process.out = cms.OutputModule(
        'PoolOutputModule',
        fileName = cms.untracked.string(fileName)
        )
    process.endp = cms.EndPath(process.out)

# Return a list of files suitable for passing into a PoolSource,
# obtained by globbing using the pattern passed in (which includes the
# path to them).
def poolAllFiles(pattern):
    import glob
    return ['file:%s' % x for x in glob.glob(pattern)]

# Function to select a set of alignment constants. Useful when doing
# track re-reconstruction.
def setTrackerAlignment(process, connectString, tagTrackerAlignmentRcd, tagTrackerAlignmentErrorRcd):
    process.include('CondCore/DBCommon/data/CondDBSetup.cfi')

    process.trackerAlignment = cms.ESSource('PoolDBESSource',
                                            process.CondDBSetup,
                                            connect = cms.string(connectString),
                                            toGet = cms.VPSet(
        cms.PSet(record = cms.string('TrackerAlignmentRcd'), tag = cms.string(tagTrackerAlignmentRcd)),
        cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'), tag = cms.string(tagTrackerAlignmentErrorRcd)),
        ))

    process.prefer("trackerAlignment")

def setMuonAlignment(process, connectString, tagDTAlignmentRcd, tagDTAlignmentErrorRcd, tagCSCAlignmentRcd, tagCSCAlignmentErrorRcd):
    process.include('CondCore/DBCommon/data/CondDBSetup.cfi')

    process.muonAlignment = cms.ESSource('PoolDBESSource',
                                         process.CondDBSetup,
                                         connect = cms.string(connectString),
                                         toGet = cms.VPSet(
        cms.PSet(record = cms.string('DTAlignmentRcd'), tag = cms.string(tagDTAlignmentRcd)),
        cms.PSet(record = cms.string('DTAlignmentErrorRcd'), tag = cms.string(tagDTAlignmentErrorRcd)),
        cms.PSet(record = cms.string('CSCAlignmentRcd'), tag = cms.string(tagCSCAlignmentRcd)),
        cms.PSet(record = cms.string('CSCAlignmentErrorRcd'), tag = cms.string(tagCSCAlignmentErrorRcd)),
        ))

    process.prefer("muonAlignment")


################################################################################

# so 'from module import *' doesn't clutter the user's namespace
__all__ = [
    'muonCollections', 'electronCollections',
    'oppSign', 'oppSignMP', 'likeSignPos', 'likeSignNeg', 'diMuons',
    'diElectrons', 'muonElectron', 'electronMuon', 'cuts',
    'makeZprime2muAnalysisProcess', 'attachAnalysis',
    'attachOutputModule', 'poolAllFiles', 'setTrackerAlignment',
    'setMuonAlignment',
    'lGN','lL1','lL2','lL3','lGR','lTK','lFS','lPR','lOP','lTR'
    ]

# To debug this config, you can do from the shell:
#   python Zprime2muAnalysisCommon_cff.py
if __name__ == '__main__':
    process = makeZprime2muAnalysisProcess(['file:/dev/null'], 42, skipPAT=True)
    print process.dumpConfig()
