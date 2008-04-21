import FWCore.ParameterSet.Config as cms

# Info about the defined rec levels.
recLevels = ['GN','L1','L2','L3','GR','TK','FS','PR','OP']
numRecLevels = len(recLevels)

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
    None,
    None,
    None,
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
                                 useOtherMuonRecos=True,
                                 useTriggerInfo=True,
                                 flavorsForDileptons=diMuons,
                                 chargesForDileptons=oppSign,
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
    
    doingElectrons: whether to run on electron or muon
    collections. Default is False, so we run on muons.
    
    useOtherMuonRecos: whether to expect the extra TeV muon
    reconstructors in the input files. Default is True, but if running
    on "official" files that do not have these extra collections
    (i.e. most likely every file that was not produced with the
    scripts in test/makeScripts of this package), this needs to be set
    False.

    useTriggerInfo: whether to expect to be able to get trigger
    collections (i.e. those destined for L1-L3 rec levels) from the
    file. Default is True, but they are not yet implemented for
    electrons.
    
    flavorsForDileptons: the pair of collections to be used in making
    dileptons; one of diMuons, diElectrons, etc. above. 
    
    chargesForDileptons: the pair of charges (in corresponding order
    to the flavors) to be used in making dileptons; one of oppSign,
    etc. above.

      Example: if flavorsForDileptons == muonElectron and
      chargesForDileptons == oppSign, then the dilepton construction
      will produce all pairs of mu+ e-.
    
    '''

    ####################################################################
    ## Sanity checks.
    ####################################################################

    if len(muons) != numRecLevels or len(electrons) != numRecLevels:
        raise RuntimeError, 'at least one of the muon and electron collections is not the right length'

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

    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

    process.MessageLogger = cms.Service(
        'MessageLogger',
        destinations = cms.untracked.vstring('Zprime'),
        categories = cms.untracked.vstring(
            'FwkJob', 'FwkReport', 'Root_Warning',
            'Root_NoDictionary', 'RFIOFileDebug'),
        Zprime = cms.untracked.PSet(
            extension    = cms.untracked.string('.out'),
            threshold    = cms.untracked.string('DEBUG'),
            lineLength   = cms.untracked.int32(132),
            noLineBreaks = cms.untracked.bool(True),
            FwkReport = cms.untracked.PSet(reportEvery = cms.untracked.int32(500)),
            RFIOFileDebug = cms.untracked.PSet(limit = cms.untracked.int32(0)),
            Root_NoDictionary = cms.untracked.PSet(limit = cms.untracked.int32(0))
            ),
        debugModules = cms.untracked.vstring(
            'leptonMatches', 'bestMuons', 'bestMatches',
            'Zprime2muAnalysis', 'Zprime2muResolution',
            'Zprime2muAsymmetry', 'Zprime2muMassReach')
        )

    process.TFileService = cms.Service(
        'TFileService',
        fileName = cms.string('zp2mu_histos.root')
        )
    
    #process.Tracer = cms.Service('Tracer')
    #process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck')

    ####################################################################
    ## Make a CandidateCollection out of the GEANT tracks.
    ####################################################################

    process.simParticleCandidates = cms.EDProducer(
        'GenCandsFromSimTracksProducer',
        src = cms.InputTag('g4SimHits')
        )

    process.makeSimCands = cms.Path(process.simParticleCandidates)

    ####################################################################
    ## Produce/clone/copy/do whatever to the needed muon collections.
    ####################################################################

    if not doingElectrons:
        process.genMuons = cms.EDFilter(
            'CandSelector',
            src = cms.InputTag('genParticleCandidates'),
            cut = cms.string('abs(pdgId) = 13 & status = 1') # & pt > 0.001'
            )
        
        process.simMuons = cms.EDFilter(
            'CandSelector',
            src = cms.InputTag('simParticleCandidates'),
            cut = cms.string('abs(pdgId) = 13 & status = 1') # & pt > 0.001'
            )
        
        process.muCandGN = cms.EDFilter(
            'CandMerger',
            src = cms.VInputTag(cms.InputTag('genMuons'),
                                cms.InputTag('simMuons'))
            )
        
        # 'Sanitize' the L1 muons, i.e. shift their phi values from
        # the bin edge to the bin center (an offset of 0.0218 rad).
        process.muCandL1 = cms.EDFilter(
            'L1MuonSanitizer',
            src = cms.InputTag('l1extraParticles')
            )
        
        # 'Sanitize' the L2 muons, i.e. in case there is no
        # hltL2MuonCandidates collection in the event, put an empty
        # one in (otherwise just copy the existing one).
        process.muCandL2 = cms.EDFilter(
            'L2MuonSanitizer',
            src = cms.InputTag('hltL2MuonCandidates')
            )
        
        # 'Sanitize' the L3 muons, i.e. make up reco::Muons from
        # reco::MuonTrackLinks (the hltL3MuonCandidates drop some of
        # the extra track information, while the
        # MuonTrackLinksCollection hltL3Muons appropriately has links
        # to all 3 tracks).
        process.muCandL3 = cms.EDFilter(
            'L3MuonSanitizer',
            src = cms.InputTag('hltL3Muons')
            )
        
        # Make tracker-only reco::Muons out of the tracker tracks in
        # the muons collection.
        process.muCandTK = cms.EDProducer(
            'TrackerOnlyMuonProducer',
            src = cms.InputTag('muons')
            )
        
        process.makeMuonCandidates = cms.Sequence(        
            process.genMuons * process.simMuons * process.muCandGN *
            process.muCandL1 * process.muCandL2 * process.muCandL3 *
            process.muCandTK
            )

    ####################################################################
    ## Same for the electrons
    ####################################################################

    if doingElectrons:
        process.genElectrons = cms.EDFilter(
            'CandSelector',
            src = cms.InputTag('genParticleCandidates'),
            cut = cms.string('abs(pdgId) = 11 & status = 1')
            )
        
        process.simElectrons = cms.EDFilter(
            'CandSelector',
            src = cms.InputTag('simParticleCandidates'),
            cut = cms.string('abs(pdgId) = 11 & status = 1')
            )
        
        process.elCandGN = cms.EDFilter(
            'CandMerger',
            src = cms.VInputTag(cms.InputTag('genElectrons'),
                                cms.InputTag('simElectrons'))
            )
        
        # Hack for now to include old-style cfgs in python ones. (Doing
        # this takes ~30 seconds on my machine! Hopefully a python config
        # version will appear soon.)
        process.include('DLEvans/HEEPSelector/data/heepSelection_1_1.cfi')
        
        process.makeElectronCandidates = cms.Sequence(
            process.genElectrons * process.simElectrons * process.elCandGN
            * process.heepSelection
            )

    ####################################################################
    ## Make the appropriate candidates defined in the sequences above.
    ####################################################################

    if doingElectrons:
        process.makeCandidates = cms.Path(process.makeElectronCandidates)
    else:
        process.makeCandidates = cms.Path(process.makeMuonCandidates)

    ####################################################################
    ## Common module configuration parameters.
    ####################################################################

    process.Zprime2muAnalysisCommon = cms.PSet(
        ################################################################
        ## Analysis configuration 
        ################################################################
        doingElectrons    = cms.bool(doingElectrons),
        useOtherMuonRecos = cms.bool(useOtherMuonRecos),
        useTriggerInfo    = cms.bool(useTriggerInfo),
        dateHistograms    = cms.untracked.bool(True),
        doingHiggs        = cms.bool(False),
        constructGenDil   = cms.bool(False),
        generatedOnly     = cms.bool(False),
        doingGeant4       = cms.bool(True),
        reconstructedOnly = cms.bool(False),
        usingAODOnly      = cms.bool(False),

        ################################################################
        ## Input tags for particles/trigger paths
        ################################################################
        l1ParticleMap = cms.InputTag('l1extraParticleMap'),
        # Every process puts a TriggerResults product into the event;
        # pick the HLT one.
        hltResults = cms.InputTag('TriggerResults','','HLT'),
        # Stand-alone muons are used as 'seeds' to try to match the
        # various levels of global reconstruction.
        standAloneMuons = cms.InputTag('standAloneMuons','UpdatedAtVtx'),
        # Reconstructed photons to try to recover brem losses when
        # running on muons.
        photons = cms.InputTag('correctedPhotons'),
        
        ################################################################
        ## Input tags for leptons at the different rec levels, in
        ## order; filled below using the collections defined at the
        ## top of the file.
        ################################################################
        muInputs = cms.VInputTag(),
        elInputs = cms.VInputTag()
        )

    for name in muons:
        if name != None:
            process.Zprime2muAnalysisCommon.muInputs.append(cms.InputTag(name))
    for name in electrons:
        if name != None:
            process.Zprime2muAnalysisCommon.elInputs.append(cms.InputTag(name))

    ####################################################################
    ## Do all the matching -- between each pair of rec levels
    ## (including MC), to seed leptons (for global fits), and to
    ## photons
    ####################################################################
            
    process.leptonMatches = cms.EDFilter(
        'LeptonAssociator',
        process.Zprime2muAnalysisCommon,
        verbosity    = cms.untracked.int32(0),
        fromCocktail = cms.bool(False)
        )

    process.matches = cms.Path(process.leptonMatches)

    ####################################################################
    ## Pick the 'best' muons according to a cocktail (currently
    ## implemented are Piotr's old cocktail and TMR), and match them
    ## to the other rec levels.
    ####################################################################

    process.bestMuons = cms.EDProducer(
        'CocktailMuonProducer',
        process.Zprime2muAnalysisCommon,
        verbosity = cms.untracked.int32(0),
        useTMR    = cms.bool(False),
        )

    process.bestMatches = cms.EDFilter(
        'LeptonAssociator',
        process.Zprime2muAnalysisCommon,
        verbosity    = cms.untracked.int32(0),
        fromCocktail = cms.bool(True)
        )

    process.best = cms.Path(process.bestMuons * process.bestMatches)

    ####################################################################
    ## Dilepton construction.
    ####################################################################

    path = []
    for i in xrange(numRecLevels):
        # Skip the other TeV muon collections if desired.
        if not useOtherMuonRecos and i >= 5 and i <= 7:
            continue

        # We finally come to the point of the Nones in the electron
        # collection above: to skip rec levels that have no
        # collections here when producing dileptons. Otherwise,
        # CandCombiner throws an exception when it cannot find the
        # input collection.
        coll0 = eval(flavorsForDileptons[0])[i]
        coll1 = eval(flavorsForDileptons[1])[i]
        if coll0 == None or coll1 == None:
            continue

        charge0, charge1 = chargesForDileptons
        
        # Use CandViewShallowCloneCombiner for the ShallowClone bit;
        # this enables us to get the reference to the original lepton
        # that makes up the dilepton we're looking at in the code.
        combiner = 'CandViewShallowCloneCombiner'

        # Put the decay string in the format expected by CandCombiner:
        # e.g. 'muons@+ muons@-'.
        decay = '%s@%s %s@%s' % (coll0, charge0, coll1, charge1)

        # A dummy cut, otherwise CandCombiner crashes.
        cut = 'mass > 0'

        dilProd = cms.EDProducer(combiner,
                                 decay = cms.string(decay),
                                 cut = cms.string(cut))

        # Attach this producer to the process.
        name = 'dileptons%s' % recLevels[i]
        setattr(process, name, dilProd)
        path.append('process.' + name)

    # Make a path that runs all the producers we just made.
    pathCodeStr = ' * '.join(path)
    process.dileptons = cms.Path(eval(pathCodeStr))

    ####################################################################
    ## Done!
    ####################################################################

    return process
