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
    'muons',
    'muCandTK', # muonsTK
    'muonsFMS',
    'muonsPMR',
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
                                 useTrigger=True,
                                 useOtherMuonRecos=True,
                                 useHEEPSelector=True,
                                 lumiForCSA07=0.,
                                 dumpHardInteraction=False,
                                 flavorsForDileptons=diMuons,
                                 chargesForDileptons=oppSign,
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

    if useHEEPSelector:
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
            'leptonMatches', 'bestMuons', 'bestMatches',
            'Zprime2muAnalysis', 'Zprime2muResolution',
            'Zprime2muAsymmetry', 'Zprime2muMassReach',
            'Zprime2muBackgrounds')
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
            src = cms.InputTag('genParticleCandidates') #,
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
            'CandSelector',
            src = cms.InputTag('genParticleCandidates'),
            cut = cms.string('abs(pdgId) = 13 & status = 1') # & pt > 0.001'
            )

        GNtags = cms.VInputTag(cms.InputTag('genMuons'))
        if useSim:
            process.simMuons = cms.EDProducer(
                'CandSelector',
                src = cms.InputTag('simParticleCandidates'),
                cut = cms.string('abs(pdgId) = 13 & status = 1') # & pt > 0.001'
                )
            
            process.psimMuons = cms.Path(process.simMuons)

            GNtags.append(cms.InputTag('simMuons'))
            
        process.muCandGN = cms.EDProducer(
            'CandMerger',
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
        # Make tracker-only reco::Muons out of the tracker tracks in
        # the muons collection.
        process.muCandTK = cms.EDProducer(
            'TrackerOnlyMuonProducer',
            src = cms.InputTag('muons')
            )
    
        process.pmuCandTK = cms.Path(process.muCandTK)

    ####################################################################
    ## Same for the electrons
    ####################################################################

    if useGen:
        process.genElectrons = cms.EDProducer(
            'CandSelector',
            src = cms.InputTag('genParticleCandidates'),
            cut = cms.string('abs(pdgId) = 11 & status = 1')
            )

        GNtags = cms.VInputTag(cms.InputTag('genElectrons'))
        if useSim:
            process.simElectrons = cms.EDProducer(
                'CandSelector',
                src = cms.InputTag('simParticleCandidates'),
                cut = cms.string('abs(pdgId) = 11 & status = 1')
                )
            
            process.psimElectrons = cms.Path(process.simElectrons)

            GNtags.append(cms.InputTag('simElectrons'))
    
        process.elCandGN = cms.EDProducer(
            'CandMerger',
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
        # Hack for now to include old-style cfgs in python ones.
        process.include('DLEvans/HEEPSelector/data/heepSelection_1_1.cfi')
        process.bestEl = cms.Path(process.heepSelection)

    ####################################################################
    ## Common module configuration parameters.
    ####################################################################

    process.Zprime2muAnalysisCommon = cms.PSet(
        ################################################################
        ## Analysis configuration 
        ################################################################
        doingElectrons    = cms.bool(doingElectrons),
        useGen            = cms.bool(useGen),
        useSim            = cms.bool(useSim),
        useReco           = cms.bool(useReco),
        usingAODOnly      = cms.bool(usingAODOnly),
        useTrigger        = cms.bool(useTrigger),
        useOtherMuonRecos = cms.bool(useOtherMuonRecos),

        dateHistograms    = cms.untracked.bool(True),
        doingHiggs        = cms.bool(False),
        constructGenDil   = cms.bool(False),

        ################################################################
        ## Input tags for trigger paths and particles.
        ################################################################
        l1ParticleMap = cms.InputTag('l1extraParticleMap'),
        # Every process puts a TriggerResults product into the event;
        # pick the HLT one.
        hltResults = cms.InputTag('TriggerResults','','HLT')
        )

    process.recLevelHelperPSet = cms.PSet()
    
    ################################################################
    ## Input tags for leptons at the different rec levels, in
    ## order; filled below using the collections defined at the
    ## top of the file, and according to whether the default leptons
    ## are electrons or muons (determined by doingElectrons).
    ################################################################
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
    ## Do all the matching -- between each pair of rec levels
    ## (including MC), to seed leptons (for global fits), and to
    ## photons
    ####################################################################

    leptonMatchPSet = cms.PSet(
        verbosity      = cms.untracked.int32(0),
        doingElectrons = cms.bool(doingElectrons),
        # Stand-alone muons are used as 'seeds' to try to match the
        # various levels of global reconstruction.
        standAloneMuons = cms.InputTag('standAloneMuons','UpdatedAtVtx'),
        # Reconstructed photons to try to recover brem losses when
        # running on muons.
        photons = cms.InputTag('correctedPhotons')
        )
        
    process.leptonMatches = cms.EDProducer(
        'LeptonAssociator',
        leptonMatchPSet,
        process.recLevelHelperPSet,
        fromBest  = cms.bool(False)
        )

    process.matches = cms.Path(process.leptonMatches)

    ####################################################################
    ## Pick the 'best' muons according to a cocktail (currently
    ## implemented are Piotr's old cocktail and TMR), and match them
    ## to the other rec levels.
    ####################################################################

    process.bestMuons = cms.EDProducer(
        'CocktailMuonProducer',
        process.recLevelHelperPSet,
        useOtherMuonRecos = cms.bool(useOtherMuonRecos),
        verbosity = cms.untracked.int32(0),
        useTMR    = cms.bool(False),
        )
    
    # Need to match the "best" leptons to all the others. (And,
    # CocktailMuonProducer needs the results of the by-seed
    # matching to do its thing, so the matching must be done in
    # two steps.)
    process.bestMatches = cms.EDProducer(
        'LeptonAssociator',
        leptonMatchPSet,
        process.recLevelHelperPSet,
        fromBest  = cms.bool(True)
        )
    
    process.bestMu = cms.Path(process.bestMuons * process.bestMatches)

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

    # Keep track of the modules we make so we can set them all into a
    # path in the process later.
    path = []
    
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

            dilProd = cms.EDProducer(combiner,
                                     decay = cms.string(decay),
                                     cut = cms.string(cut))


            # Attach this producer to the process, encoding the
            # name as e.g. 'dileptonselmuMPOP' meaning e-mu+ at OP rec
            # level.
            name = nameDilCollection(theseFlavors, theseCharges, recLevels[i])
            setattr(process, name, dilProd)
            path.append('process.' + name)

    # Make a path that runs all the producers we just made.
    #
    # When the config is converted to the old style configuration
    # language, and then reparsed (as happens if you run this module
    # with CRAB), Python complains about too much recursion if there
    # are too many modules in one path. So, split up the path.
    count = 0
    while path:
        sub = path[:10]
        path = path[10:]
        pathCodeStr = ' * '.join(sub)
        setattr(process, 'dileptons%i' % count, cms.Path(eval(pathCodeStr)))
        count += 1

    # Set up the "main" dileptons -- the ones the analysis module is
    # supposed to use by default.
    for i in xrange(numRecLevels):
        # If we skipped making this dilepton collection in the logic
        # above, skip using it.
        name = nameDilCollection(flavorsForDileptons, chargesForDileptons, recLevels[i])
        if hasattr(process, name):
            tag = cms.InputTag(name)
        else:
            tag = cms.InputTag('skip')
            
        tagname = 'dileptons' + recLevels[i]
        setattr(process.recLevelHelperPSet, tagname, tag)

    ####################################################################
    ## Done!
    ####################################################################

    return process
