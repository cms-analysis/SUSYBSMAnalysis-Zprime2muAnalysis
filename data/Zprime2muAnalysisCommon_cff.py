import FWCore.ParameterSet.Config as cms

# Info about the defined rec levels.
recLevels = ['GN','L1','L2','L3','GR','TK','FS','PR','OP']
numRecLevels = len(recLevels)
# Make enums for the reclevels to avoid magic numbers below.
for i, rec in enumerate(recLevels):
    locals()['l' + rec] = i

# The InputTag names for the muon collections.
muonCollections = [
    'muCandGN',
    'muCandL1', # l1extraParticles
    'muCandL2', # hltL2MuonCandidates
    'muCandL3', # hltL3MuonCandidates
    'muCandGR',
    'muCandTK', # muonsTK
    'muCandFS',
    'muCandPR',
    'bestMuons'
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
    'heepSelector'
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

def makeZprime2muAnalysisProcess(fileNames=[], maxEvents=-1,
                                 doingElectrons=False,
                                 useGen=True,
                                 useSim=True,
                                 useReco=True,
                                 usingAODOnly=False,
                                 useTrigger=False,
                                 useOtherMuonRecos=False,
                                 useHEEPSelector=False,
                                 recoverBrem=True,
                                 disableElectrons=False,
                                 lumiForCSA07=0.,
                                 dumpHardInteraction=False,
                                 flavorsForDileptons=diMuons,
                                 chargesForDileptons=oppSign,
                                 maxDileptons=1,
                                 # These last two are passed in just
                                 # so they can be put standalone at
                                 # the top of the file so they stick
                                 # out.
                                 muons=muonCollections,
                                 electrons=electronCollections):
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
      order. Other than paths everything is (should be?) replaceable.
    
    doingElectrons: whether to run on electron or muon collections
    (i.e. which collections are in allLeptons in the code). Default is
    False, so we run on muons.
    
    useGen/Sim/Reco: whether to expect information from GEN, SIM or
    RECO branches. (The expectation from RECO is modified by the
    usingAODOnly parameter below.)
    
    usingAODOnly: if True, do not expect to be able to get information
    that is not in AOD, such as rechits.

    useOtherMuonRecos: whether to expect the extra TeV muon
    reconstructors in the input files. Default is True, but if running
    on "official" files that do not have these extra collections
    (i.e. most likely every file that was not produced with the
    scripts in test/makeScripts of this package), this needs to be set
    False.

    useHEEPSelector: whether to use D.L. Evans\'s HEEPSelector for
    "best" electrons. It requires quantities not in AOD, so it may
    need to be disabled via this flag.
    
    useTrigger: whether to expect to be able to get trigger
    collections (i.e. those destined for L1-L3 rec levels) from the
    file.

    lumiForCSA07: if > 0, set up the CSA07EventWeightProducer with
    weights corresponding to the given lumi (in pb^-1).

    dumpHardInteraction: if True, use the (modified)
    ParticleListDrawer to dump the generator info on the hard
    interaction particles.
    
    flavorsForDileptons: the pair of collections to be used for the
    "official" dileptons (i.e. those accessed by allDileptons in the
    code); one of diMuons, diElectrons, etc. above.
    
    chargesForDileptons: the pair of charges (in corresponding order
    to the flavors) to be used for the "official" dileptons (see
    previous); one of oppSign, etc. above.

      Example: if flavorsForDileptons == muonElectron and
      chargesForDileptons == oppSign, then the official dileptons
      will be all pairs of mu+ e-.
    
    '''

    ####################################################################
    ## Sanity checks.
    ####################################################################

    if len(muons) != numRecLevels or len(electrons) != numRecLevels:
        raise RuntimeError, 'at least one of the muon and electron collections is not the right length'

    if disableElectrons:
        for i in xrange(numRecLevels):
            electrons[i] = None
        useHEEPSelector = False
        
    if usingAODOnly:
        useGen = useSim = False
        useOtherMuonRecos = False
        useHEEPSelector = False
        
    if not useGen:
        useSim = False
        muons[lGN] = electrons[lGN] = None

    if not useReco:
        useTrigger = False
        useOtherMuonRecos = False
        useHEEPSelector = False
        for i in xrange(lL1, lOP+1):
            muons[i] = None
            electrons[i] = None

    if not useHEEPSelector:
        electrons[lOP] = electrons[lGR]

    if not useOtherMuonRecos:
        for i in xrange(lTK, lPR+1):
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

    #process.options = cms.untracked.PSet(IgnoreCompletely = cms.untracked.vstring('ProductNotFound'))
    
    process.MessageLogger = cms.Service(
        'MessageLogger',
        destinations = cms.untracked.vstring('Zprime'),
        categories = cms.untracked.vstring(
            #'FwkJob', 'FwkReport', 'Root_Warning',
            'Root_NoDictionary', 'RFIOFileDebug',
            # For some reason these next ones come up (without being
            # asked to in debugModules below) when doing electrons,
            # i.e. when the HEEPSelector.cfi is included.
            'DDLParser', 'PixelGeom', 'EcalGeom', 'TIBGeom', 'TIDGeom',
            'TOBGeom', 'TECGeom', 'SFGeom', 'HCalGeom', 'TrackerGeom',
            'GeometryConfiguration', 'HcalHardcodeGeometry'
            ),
        Zprime = cms.untracked.PSet(
            extension    = cms.untracked.string('.out'),
            threshold    = cms.untracked.string('DEBUG'),
            lineLength   = cms.untracked.int32(132),
            noLineBreaks = cms.untracked.bool(True)
            ),
        debugModules = cms.untracked.vstring(
            'bestMuons', 'Zprime2muAnalysis', 'Zprime2muResolution',
            'Zprime2muAsymmetry', 'Zprime2muMassReach',
            'Zprime2muBackgrounds', 'Zprime2muMatchStudy')
        )

    # Instead of line after line of limit psets in Zprime above, set
    # them all here.
    limitZero = cms.untracked.PSet(limit = cms.untracked.int32(0))
    for cat in process.MessageLogger.categories:
        setattr(process.MessageLogger.Zprime, cat, limitZero)
    setattr(process.MessageLogger.Zprime, 'FwkReport', cms.untracked.PSet(reportEvery = cms.untracked.int32(500)))

    process.TFileService = cms.Service(
        'TFileService',
        fileName = cms.string('zp2mu_histos.root')
        )
    
    #process.Tracer = cms.Service('Tracer')
    #process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck')

    if useGen and lumiForCSA07 > 0:
        process.csa07EventWeightProducer = cms.EDProducer(
            'CSA07EventWeightProducer',
            src         = cms.InputTag('source'),
            talkToMe    = cms.untracked.bool(False),
            overallLumi = cms.double(lumiForCSA07), # in pb^-1
            ttKfactor   = cms.double(1) # 1.85
            )

        process.csa07 = cms.Path(process.csa07EventWeightProducer)

    if useGen and dumpHardInteraction:
        process.include("SimGeneral/HepPDTESSource/data/pythiapdt.cfi")
        
        process.printTree = cms.EDAnalyzer(
            'ParticleListDrawer',
            maxEventsToPrint = cms.untracked.int32(-1),
            src = cms.InputTag('genParticles') #,
            #printOnlyHardInteraction = cms.untracked.bool(True),
            #useMessageLogger = cms.untracked.bool(True),
            )
        
        process.ptree = cms.Path(process.printTree)

    ####################################################################
    ## Make a CandidateCollection out of the GEANT tracks.
    ####################################################################

    if useSim:
        process.simParticleCandidates = cms.EDProducer(
            'GenCandsFromSimTracksProducer',
            src = cms.InputTag('g4SimHits')
            )
        
        process.psimParticleCandidates = cms.Path(process.simParticleCandidates)

    ####################################################################
    ## Produce/clone/copy/do whatever to the needed muon collections.
    ####################################################################

    if useGen:
        process.genMuons = cms.EDProducer(
            'CandViewSelector',
            src = cms.InputTag('genParticles'),
            cut = cms.string('abs(pdgId) = 13 & status = 1') # & pt > 0.001'
            )

        GNtags = cms.VInputTag(cms.InputTag('genMuons'))
        if useSim:
            process.simMuons = cms.EDProducer(
                'CandViewSelector',
                src = cms.InputTag('simParticleCandidates'),
                cut = cms.string('abs(pdgId) = 13 & status = 1') # & pt > 0.001'
                )
            
            process.psimMuons = cms.Path(process.simMuons)

            GNtags.append(cms.InputTag('simMuons'))
            
        process.muCandGN = cms.EDProducer(
            'CandViewMerger',
            src = GNtags
            )
        
        process.pmuCandGN = cms.Path(process.genMuons * process.muCandGN)

    if useTrigger:
        # 'Sanitize' the L1 muons, i.e. shift their phi values from
        # the bin edge to the bin center (an offset of 0.0218 rad).
        process.muCandL1 = cms.EDProducer(
            'L1MuonSanitizer',
            src = cms.InputTag('l1extraParticles')
            )
    
        # 'Sanitize' the L2 muons, i.e. in case there is no
        # hltL2MuonCandidates collection in the event, put an empty
        # one in (otherwise just copy the existing one).
        process.muCandL2 = cms.EDProducer(
            'L2MuonSanitizer',
            src = cms.InputTag('hltL2MuonCandidates')
            )
    
        # 'Sanitize' the L3 muons, i.e. make up reco::Muons from
        # reco::MuonTrackLinks (the hltL3MuonCandidates drop some of
        # the extra track information, while the
        # MuonTrackLinksCollection hltL3Muons appropriately has links
        # to all 3 tracks).
        process.muCandL3 = cms.EDProducer(
            'L3MuonSanitizer',
            src = cms.InputTag('hltL3Muons')
            )

        process.pmuTrig = cms.Path(process.muCandL1 * process.muCandL2 *
                                   process.muCandL3)
        
    if useReco:
        # Copy only the GlobalMuons from the default muon collection,
        # ignoring the TrackerMuons and CaloMuons for now.
        process.muCandGR = cms.EDProducer(
            'GlobalOnlyMuonProducer',
            src              = cms.InputTag('muons'),
            fromTrackerTrack = cms.bool(False)
            )

        # Make tracker-only reco::Muons out of the tracker tracks in
        # the muons collection.
        process.muCandTK = cms.EDProducer(
            'GlobalOnlyMuonProducer',
            src              = cms.InputTag('muons'),
            fromTrackerTrack = cms.bool(True)
            )
    
        process.pmuCandTK = cms.Path(process.muCandGR * process.muCandTK)

    ####################################################################
    ## Same for the electrons
    ####################################################################

    if useGen:
        process.genElectrons = cms.EDProducer(
            'CandViewSelector',
            src = cms.InputTag('genParticles'),
            cut = cms.string('abs(pdgId) = 11 & status = 1')
            )

        GNtags = cms.VInputTag(cms.InputTag('genElectrons'))
        if useSim:
            process.simElectrons = cms.EDProducer(
                'CandViewSelector',
                src = cms.InputTag('simParticleCandidates'),
                cut = cms.string('abs(pdgId) = 11 & status = 1')
                )
            
            process.psimElectrons = cms.Path(process.simElectrons)

            GNtags.append(cms.InputTag('simElectrons'))
    
        process.elCandGN = cms.EDProducer(
            'CandViewMerger',
            src = GNtags
            )

        process.pelCandGN = cms.Path(process.genElectrons * process.elCandGN)

    if useTrigger:
        process.elCandL1 = cms.EDProducer(
            'CandViewMerger',
            src = cms.VInputTag(
            cms.InputTag('l1extraParticles','NonIsolated'),
            cms.InputTag('l1extraParticles','Isolated'))
            )

        # For electrons, the HLT makes its decision based off of the
        # l1(Non)IsoRecoEcalCandidate and the
        # pixelMatchElectronsL1(Non)IsoForHLT collections. I choose to
        # call them L2 and L3 electrons here, respectively.
        
        # 'Sanitize' the L2 & L3 electrons, i.e. in case there is no
        # collection in the event, put an empty one in (otherwise just
        # copy the existing one).
        
        process.elCandL2NonIso = cms.EDProducer(
            'L2ElectronSanitizer',
            src = cms.InputTag('l1NonIsoRecoEcalCandidate')
            )
        
        process.elCandL2Iso = cms.EDProducer(
            'L2ElectronSanitizer',
            src = cms.InputTag('l1IsoRecoEcalCandidate')
            )
        
        process.elCandL3NonIso = cms.EDProducer(
            'L3ElectronSanitizer',
            src = cms.InputTag('pixelMatchElectronsL1NonIsoForHLT')
            )
        
        process.elCandL3Iso = cms.EDProducer(
            'L3ElectronSanitizer',
            src = cms.InputTag('pixelMatchElectronsL1IsoForHLT')
            )
        
        process.elCandL2 = cms.EDProducer(
            'CandViewMerger',
            src = cms.VInputTag(
            cms.InputTag('elCandL2NonIso'),
            cms.InputTag('elCandL2Iso'))
            )
        
        process.elCandL3 = cms.EDProducer(
            'CandViewMerger',
            src = cms.VInputTag(
            cms.InputTag('elCandL3NonIso'),
            cms.InputTag('elCandL3Iso'))
            )

        process.pelTrig = cms.Path(process.elCandL1 * process.elCandL2NonIso *
                                   process.elCandL2Iso * process.elCandL3NonIso *
                                   process.elCandL3Iso * process.elCandL2 *
                                   process.elCandL3)

    if useHEEPSelector:
        ## Hack for now to include old-style cfgs in python ones.
        #process.include('DLEvans/HEEPSelector/data/heepSelection_1_1.cfi')
        #process.bestEl = cms.Path(process.heepSelection)

        # Use the output of running dumpPython() on the above for
        # faster loading (will have to be redumped when there is a new
        # HEEPSelector!).
        import HEEPSelector11_cfi
        HEEPSelector11_cfi.attachHEEPSelector(process)

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

        dateHistograms    = cms.untracked.bool(True),

        #resonanceIds      = cms.vint32(32, 23, 39, 5000039),

        ################################################################
        ## Input tags for trigger paths and particles.
        ################################################################
        l1ParticleMap = cms.InputTag('l1extraParticleMap'),
        # Every process puts a TriggerResults product into the event;
        # pick the HLT one.
        hltResults = cms.InputTag('TriggerResults','','HLT')
        )

    process.recLevelHelperPSet = cms.PSet()
    process.plainAnalysisPSet  = cms.PSet()
    
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
        tagname = 'leptons' + recLevels[i]
        
        if lepcoll[i] != None:
            tag = cms.InputTag(lepcoll[i])    
        else:
            # The code (via RecLevelHelper) will know to ignore this
            # (we cannot have empty InputTags, apparently).
            tag = cms.InputTag('skip')
            
        setattr(process.recLevelHelperPSet, tagname, tag)

    ####################################################################
    ## After this point, we will add a *lot* of modules to the path in
    ## loops, with their names depending on the loop variables.  Keep
    ## track of the modules we make so we can set them all into a real
    ## CMSSW path in the process later. The order of the modules in
    ## this list will be obeyed by CMSSW.
    ####################################################################
    path = []
    def addToPath(name, prodobj):
        path.append('process.' + name)
        setattr(process, name, prodobj)
                    
    ####################################################################
    ## Do all the by-seed matching between each pair of globally-fit
    ## muon collections, two of which are needed by
    ## CocktailMuonProducer, so these must be run first.
    ####################################################################

    for irec in xrange(lL3, numRecLevels-1):
        for jrec in xrange(lL3, numRecLevels-1):
            # Don't store an identity map, and skip missing collections.
            if irec == jrec or muons[irec] is None or muons[jrec] is None:
                continue
            
            prod = cms.EDProducer(
                'MuonBySeedMatcher',
                # Stand-alone muon tracks are used as 'seeds'.
                seedTracks = cms.InputTag('standAloneMuons','UpdatedAtVtx'),
                src        = cms.InputTag(muons[irec]),
                matched    = cms.InputTag(muons[jrec])
                )
            addToPath('seedMatch' + recLevels[irec] + recLevels[jrec], prod)

    ####################################################################
    ## Pick the 'best' muons according to a cocktail (currently
    ## implemented are Piotr's old cocktail and TMR), and match them
    ## to the other rec levels.
    ####################################################################

    if useOtherMuonRecos:
        process.bestMuons = cms.EDProducer(
            'CocktailMuonProducer',
            useTMR           = cms.bool(False),
            trackerOnlyMuons = cms.InputTag(muons[lTK]),
            toFMSMap         = cms.InputTag('seedMatchTKFS'),
            toPMRMap         = cms.InputTag('seedMatchTKPR')
            )
        path.append('process.bestMuons')

    ####################################################################
    ## Do all the closest-in-deltaR matching -- between each pair of
    ## rec levels (including MC), and between each rec level and
    ## reconstructed photons. 
    ####################################################################

    # Closest (deltaR) matching between all pairs of rec levels.
    for irec in xrange(numRecLevels):
        for jrec in xrange(numRecLevels):
            # Again, don't store an identity map.
            if irec == jrec: continue

            if doingElectrons:
                fromColl = electrons[irec]
                toColl   = electrons[jrec]
            else:
                fromColl = muons[irec]
                toColl   = muons[jrec]

            if fromColl is None or toColl is None:
                continue
            
            matcher = 'DeltaRViewMatcher'
            if irec == lGN or jrec == lGN:
                # Could switch to MCTruth matcher, but needs studying.
                matcher = 'PtMin' + matcher
            else:
                matcher = 'Trivial' + matcher

            prod = cms.EDProducer(
                matcher,
                src     = cms.InputTag(muons[irec]),
                matched = cms.InputTag(muons[jrec]),
                distMin = cms.double(0.7072)
                )
            addToPath('closestMatch' + recLevels[irec] + recLevels[jrec], prod)

    # Match all muons to closest photons to try and recover energy
    # lost to brem later. (Don't bother doing this for electrons,
    # since the GSF algorithm already takes brem losses into account.)
    if recoverBrem:
        for rec in xrange(lL3, numRecLevels):
            if muons[rec] is None: continue
            prod = cms.EDProducer(
                'TrivialDeltaRViewMatcher',
                src     = cms.InputTag(muons[rec]),
                matched = cms.InputTag('correctedPhotons'),
                distMin = cms.double(999)
                )
            addToPath('photonMatch' + recLevels[rec], prod)

    # Also do the by-seed matching from the best level to the others,
    # since we skipped it before. (Again, not available for
    # electrons.)
    irec = lOP
    if muons[irec] is not None:
        for jrec in xrange(lL3, numRecLevels-1):
            if muons[jrec] is None: continue
            prod = cms.EDProducer(
                'MuonBySeedMatcher',
                # Stand-alone muon tracks are used as 'seeds'.
                seedTracks = cms.InputTag('standAloneMuons','UpdatedAtVtx'),
                src        = cms.InputTag(muons[irec]),
                matched    = cms.InputTag(muons[jrec])
                )
            addToPath('seedMatch' + recLevels[irec] + recLevels[jrec], prod)

    ####################################################################
    ## Dilepton construction. Make all possible dileptons:
    ## {mu, e} x {+, -}.
    ####################################################################

    def nameDilCollection(flavors, charges, reclevel):
        cdict = {'+': 'P', '-': 'M'}
        # I would put a _ in the label name but that is illegal in
        # CMSSW.
        return 'dileptons%s%s%s%s%s' % (flavors[0][:2], flavors[1][:2],
                                        cdict[charges[0]], cdict[charges[1]],
                                        reclevel)

    # Enumerate the valid combinations with no double counting.
    combos = [
        (diMuons,      oppSign),     # mu+ mu-
        (diMuons,      likeSignPos), # mu+ mu+
        (diMuons,      likeSignNeg), # mu- mu-
        (diElectrons,  oppSign),     # e+  e-
        (diElectrons,  likeSignPos), # e+  e+
        (diElectrons,  likeSignNeg), # e-  e-
        (muonElectron, likeSignPos), # mu+ e+
        (muonElectron, likeSignNeg), # mu- e-
        (muonElectron, oppSign),     # mu+ e-
        (muonElectron, oppSignMP)    # mu- e+
        ]

    # Make all the dileptons.
    for theseFlavors, theseCharges in combos:
        for i in xrange(numRecLevels):
            # Skip the other TeV muon collections if desired.
            if not useOtherMuonRecos and i >= lTK and i <= lPR:
                continue

            # Set up the flavor and charge pair for this rec level.
            collections = []
            charges = []
            for j in xrange(2):
                collections.append(eval(theseFlavors[j])[i])
                # Electrons at L1 and L2 do not have charges (they are
                # just ECAL superclusters); try to use them to make
                # dileptons anyway.
                if theseFlavors[j] == 'electrons' and (i == lL1 or i == lL2):
                    charges.append('')
                else:
                    charges.append('@%s' % theseCharges[j])

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
            if i == 0:
                # At generator level, look only at the leptons which came
                # from the resonance.
                combiner = 'GenDil' + combiner

            # Put the decay string in the format expected by CandCombiner:
            # e.g. 'muons@+ muons@-'.
            decay = '%s%s %s%s' % (collections[0], charges[0],
                                   collections[1], charges[1])

            # A dummy cut, otherwise CandCombiner crashes.
            cut = 'mass > 0'

            # Encode the name as e.g. 'dileptonselmuMPOP' meaning
            # e-mu+ at OP rec level.
            dilname = nameDilCollection(theseFlavors, theseCharges,
                                        recLevels[i])
            rawname = dilname + 'Raw'
            resname = dilname + 'Res'
            
            # These are the "raw" dileptons, i.e. before overlap
            # removal and selection cuts.
            prod = cms.EDProducer(combiner,
                                  decay = cms.string(decay),
                                  cut = cms.string(cut))
            addToPath(rawname, prod)

            # Remove dilepton overlap and apply selection cuts, also
            # making sure we get the requested combination by PDG id.
            # This doesn't work for L1 or L2.
            if i != lL1 and i != lL2:
                pdgIds = []
                for flavor, charge in zip(theseFlavors, theseCharges):
                    pid =  {'muons': 13, 'electrons': 11}[flavor]
                    if charge == '+': pid *= -1
                    pdgIds.append(pid)
            else:
                pdgIds = [0,0]
            
            prod = cms.EDProducer(
                'DileptonPicker',
                src          = cms.InputTag(rawname),
                doingGen     = cms.bool(i == 0),
                maxDileptons = cms.uint32(maxDileptons),
                pdgIds       = cms.vint32(*pdgIds)
                )
            addToPath(dilname, prod)
            
            # For offline dimuons, attempt to recover some of the
            # energy lost to bremsstrahlung by using the daughter
            # muons' closest photons (those within dRmax).
            if i >= lL3 and theseFlavors == diMuons:
                prod = cms.EDProducer(
                    'DimuonBremRecoverer',
                    dimuons        = cms.InputTag(dilname),
                    photonMatchMap = cms.InputTag('photonMatch' + recLevels[i]),
                    dRmax          = cms.double(0.1)
                    )
                addToPath(resname, prod)

            # For generator-level dileptons, we don't need to attempt
            # brem recovery, we can just take the generated resonance
            # from the event.
            if i == lGN:
                if doingElectrons: lf = 11
                else:              lf = 13
                prod = cms.EDProducer(
                    'GenResonanceProducer',
                    leptonFlavor = cms.int32(lf)
                    )
                addToPath(resname, prod)

    # Set up the "main" dileptons -- the ones the analysis module is
    # supposed to use by default.

    # For RecLevelAnalyzers:
    for i in xrange(numRecLevels):
        # If we skipped making this dilepton collection in the logic
        # above, skip using it.
        name = nameDilCollection(flavorsForDileptons, chargesForDileptons,
                                 recLevels[i])
        if hasattr(process, name): tag = name
        else:                      tag = 'skip'
            
        tagname = 'dileptons' + recLevels[i]
        setattr(process.recLevelHelperPSet, tagname, cms.string(tag))

    # For plain analyzers, point them at gen, hlt, default offline,
    # and "best" offline dileptons and resonances.
    names = [nameDilCollection(flavorsForDileptons, chargesForDileptons, recLevels[l]) for l in [lGN, lL3, lGR, lOP]]
    kinds = ['gen', 'hlt', 'rec', 'best']
    for kind, name in zip(kinds, names):
        setattr(process.plainAnalysisPSet, '%sDileptons' % kind,
                cms.InputTag(name))
        setattr(process.plainAnalysisPSet, '%sResonances' % kind,
                cms.InputTag(name))

    ####################################################################
    ## Make a path that runs all the producers we just made.
    ####################################################################

    # When the config is converted to the old style configuration
    # language, and then reparsed (as happens if you run this module
    # with CRAB), Python complains about too much recursion if there
    # are too many modules in one path. So, split up the path into
    # reasonable chunks.
    count = 0
    chunkSize = 5
    while path:
        pathCodeStr = ' * '.join(path[:chunkSize])
        path = path[chunkSize:]
        setattr(process, 'zp2muPath%i' % count, cms.Path(eval(pathCodeStr)))
        count += 1

    ####################################################################
    ## Done!
    ####################################################################

    print 'Zprime2muAnalysis base process built!'
    #print process.dumpConfig()
    return process

# Function to attach simple EDAnalyzers that don't need any
# parameters. Allows skipping of making data/Zprime2mu*_cfi.py files.
def attachAnalysis(process, name, isRecLevelAnalysis=False):
    args = (name, process.Zprime2muAnalysisCommon, process.plainAnalysisPSet)
    if isRecLevelAnalysis:
        args = args + (process.recLevelHelperPSet,)
        
    analyzer = cms.EDAnalyzer(*args)
    setattr(process, name, analyzer)
    setattr(process, name + 'Path', cms.Path(getattr(process, name)))
