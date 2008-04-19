import FWCore.ParameterSet.Config as cms

'''
# Set up convenient aliases by modifying the locals dict, allowing one
# to write e.g. 'cu_string' instead of 'cms.untracked.string', or
# 'c_PSet' for 'cms.PSet'.
cmsTypes = ['bool','int32','uint32','double','string','PSet','InputTag']
for t in cmsTypes:
    for vectorized in ['v','']:
        for tracked in ['u','']:
            newName = 'c%s_%s%s' % (tracked, vectorized, t)

ubool = cms.untracked.
uint32 = cms.untracked.int32
ustring = cms.untracked.string
uvstring = cms.untracked.vstring
uPSet = cms.untracked.PSet
'''

####################################################################
## Set up the CMSSW process and services.
####################################################################

process = cms.Process('Zprime2muAnalysis')

process.MessageLogger = cms.Service('MessageLogger',
                                    destinations = cms.untracked.vstring('Zprime'),
                                    categories   = cms.untracked.vstring('FwkJob', 'FwkReport', 'Root_Warning', 'Root_NoDictionary', 'RFIOFileDebug'),
                                    Zprime       = cms.untracked.PSet(extension    = cms.untracked.string('.out'),
                                                                      threshold    = cms.untracked.string('DEBUG'),
                                                                      lineLength   = cms.untracked.int32(132),
                                                                      noLineBreaks = cms.untracked.bool(True)
                                                                      ),
                                    FwkReport    = cms.untracked.PSet(reportEvery = cms.untracked.int32(500)),
                                    debugModules = cms.untracked.vstring('leptonMatches', 'bestMuons', 'bestMatches',
                                                                         'Zprime2muAnalysis', 'Zprime2muResolution',
                                                                         'Zprime2muAsymmetry', 'Zprime2muMassReach')
                                    )

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('zp2mu_histos.root')
                                   )

#process.Tracer = cms.Service('Tracer')
#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck')
    
####################################################################
## Make a CandidateCollection out of the GEANT tracks.
####################################################################

process.simParticleCandidates = cms.EDProducer("GenCandsFromSimTracksProducer",
                                               src = cms.InputTag("g4SimHits")
                                               )

process.makeSimCands = cms.Sequence(process.simParticleCandidates)

####################################################################
## Produce/clone/copy/do whatever to the needed muon collections
####################################################################

process.genMuons = cms.EDFilter("CandSelector",
                                src = cms.InputTag("genParticleCandidates"),
                                cut = cms.string('abs(pdgId) = 13 & status = 1') # & pt > 0.001"
                                )

process.simMuons = cms.EDFilter("CandSelector",
                                src = cms.InputTag("simParticleCandidates"),
                                cut = cms.string('abs(pdgId) = 13 & status = 1') # & pt > 0.001"
                                )

process.muCandGN = cms.EDFilter("CandMerger",
                                src = cms.VInputTag(cms.InputTag("genMuons"),
                                                    cms.InputTag("simMuons")
                                                    )
                                )

# "Sanitize" the L1 muons, i.e. shift their phi values from the bin
# edge to the bin center (an offset of 0.0218 rad).
process.muCandL1 = cms.EDFilter("L1MuonSanitizer",
                                src = cms.InputTag("l1extraParticles")
                                )

# "Sanitize" the L2 muons, i.e. in case there is no
# hltL2MuonCandidates collection in the event, put an empty one in
# (otherwise just copy the existing one).
process.muCandL2 = cms.EDFilter("L2MuonSanitizer",
                                src = cms.InputTag("hltL2MuonCandidates")
                                )

# "Sanitize" the L3 muons, i.e. make up reco::Muons from
# reco::MuonTrackLinks (the hltL3MuonCandidates drop some of the extra
# track information, while the MuonTrackLinksCollection hltL3Muons
# appropriately has links to all 3 tracks)
process.muCandL3 = cms.EDFilter("L3MuonSanitizer",
                                src = cms.InputTag("hltL3Muons")
                                )

# Make tracker-only reco::Muons out of the tracker tracks in the muons
# collection.
process.muCandTK = cms.EDProducer("TrackerOnlyMuonProducer",
                                  src = cms.InputTag("muons")
                                  )

process.makeMuonCandidates = cms.Sequence(process.genMuons * process.simMuons * process.muCandGN * process.muCandL1 * process.muCandL2 * process.muCandL3 * process.muCandTK)

####################################################################
## Same for the electrons
####################################################################

process.genElectrons = cms.EDFilter("CandSelector",
                                    src = cms.InputTag("genParticleCandidates"),
                                    cut = cms.string('abs(pdgId) = 11 & status = 1')
                                    )

process.simElectrons = cms.EDFilter("CandSelector",
                                    src = cms.InputTag("simParticleCandidates"),
                                    cut = cms.string('abs(pdgId) = 11 & status = 1')
                                    )

process.elCandGN = cms.EDFilter("CandMerger",
                                src = cms.VInputTag(cms.InputTag("genElectrons"), cms.InputTag("simElectrons"))
                                )

process.makeElectronCandidates = cms.Sequence(process.genElectrons * process.simElectrons * process.elCandGN)

####################################################################
## Make all the candidates defined in the sequences above.
####################################################################

process.makeCandidates = cms.Path(process.makeSimCands * process.makeMuonCandidates * process.makeElectronCandidates)

####################################################################
## Module configuration parameters
####################################################################

Zprime2muAnalysisCommon = cms.PSet(
    ################################################################
    ## Analysis configuration 
    ################################################################
    dateHistograms    = cms.untracked.bool(True),
    doingElectrons    = cms.bool(False),
    doingHiggs        = cms.bool(False),
    constructGenDil   = cms.bool(False),
    generatedOnly     = cms.bool(False),
    doingGeant4       = cms.bool(True),
    reconstructedOnly = cms.bool(False),
    useOtherMuonRecos = cms.bool(True),
    usingAODOnly      = cms.bool(False),
    useTriggerInfo    = cms.bool(True),

    ################################################################
    ## Input tags for particles/trigger paths
    ################################################################
    l1ParticleMap = cms.InputTag("l1extraParticleMap"),
    # Every process puts a TriggerResults product into the event;
    # pick the HLT one.
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    # Stand-alone muons are used as "seeds" to try to match the
    # various levels of global reconstruction.
    standAloneMuons = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    # Reconstructed photons to try to recover brem losses.
    photons = cms.InputTag("correctedPhotons"),

    ################################################################
    ## Input tags for leptons at the different rec levels, in order;
    ## filled below.
    ################################################################
    muInputs = cms.VInputTag(),
    elInputs = cms.VInputTag()
)

recLevels = ['GN','L1','L2','L3','GR','TK','FS','PR','OP']
numRecLevels = len(recLevels)

muons = ['muCandGN',
         'muCandL1', # l1extraParticles
         'muCandL2', # hltL2MuonCandidates
         'muCandL3', # hltL3MuonCandidates
         'muons',
         'muCandTK', # muonsTK
         'muonsFMS',
         'muonsPMR',
         'bestMuons']

electrons = ['elCandGN',
             None,
             None,
             None,
             'pixelMatchGsfElectrons',
             None,
             None,
             None,
             'pixelMatchGsfElectrons'
             ]

if len(muons) != numRecLevels or len(electrons) != numRecLevels:
    raise RuntimeError, 'at least one of the muon and electron collections is not the right length'

for name in muons:
    if name != None:
        Zprime2muAnalysisCommon.muInputs.append(cms.InputTag(name))
for name in electrons:
    if name != None:
        Zprime2muAnalysisCommon.elInputs.append(cms.InputTag(name))
    
####################################################################
## Do all the matching -- between each pair of rec levels (including
## MC), to seed leptons (for global fits), and to photons
####################################################################

process.leptonMatches = cms.EDFilter("LeptonAssociator",
                                     Zprime2muAnalysisCommon,
                                     verbosity    = cms.untracked.int32(0),
                                     fromCocktail = cms.bool(False)
                                     )

process.matches = cms.Path(process.leptonMatches)

####################################################################
## Pick the "best" muons according to a cocktail (currently ##
## implemented are Piotr's old cocktail, and TMR), and match them to
## the other rec levels.
####################################################################

process.bestMuons = cms.EDProducer("CocktailMuonProducer",
                                   Zprime2muAnalysisCommon,
                                   verbosity = cms.untracked.int32(0),
                                   useTMR = cms.bool(False),
                                   )

process.bestMatches = cms.EDFilter("LeptonAssociator",
                                   Zprime2muAnalysisCommon,
                                   verbosity = cms.untracked.int32(0),
                                   fromCocktail = cms.bool(True)
                                   )

process.best = cms.Path(process.bestMuons * process.bestMatches)

####################################################################
## Dilepton construction.
####################################################################

oppSign = ('+','-')
likeSignPos = ('+','+')
likeSignNeg = ('-','-')
diMuons = ('muons', 'muons')
diElectrons = ('electrons', 'electrons')
muonElectron = ('muon', 'electron')
electronMuon = ('electron', 'muon')

collections = diMuons
charges = oppSign

def setDileptonPath(process, collections, charges):
#    global process
    path = []
    for i in xrange(numRecLevels):
        coll0 = eval(collections[0])[i]
        coll1 = eval(collections[1])[i]
        if coll0 == None or coll1 == None:
            continue
        
        combiner = ''
        decay = '%s@%s %s@%s' % (coll0, charges[0], coll1, charges[1])
        cut = 'mass > 0'
        
        dilProd = cms.EDProducer('CandViewShallowCloneCombiner',
                                 decay = cms.string(decay),
                             cut = cms.string(cut))
        name = 'dileptons%s' % recLevels[i]
        setattr(process, name, dilProd)
        path.append('process.' + name)
    
    pathCodeStr = ' * '.join(path)
    process.dileptons = cms.Path(eval(pathCodeStr))

#setDileptonPath(process, diMuons, oppSign)

####################################################################
## Done!
####################################################################
