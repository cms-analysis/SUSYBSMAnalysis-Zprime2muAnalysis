import FWCore.ParameterSet.Config as cms

def attachHEEPSelector(process):
 process.cscSegments = cms.EDProducer("CSCSegmentProducer",
    inputObjects = cms.InputTag("csc2DRecHits"),
    algo_type = cms.int32(3),
    algo_psets = cms.VPSet(cms.PSet(
        parameters_per_chamber_type = cms.vint32(2, 1, 1, 1, 1, 
            1, 1, 1, 1),
        algo_name = cms.string('CSCSegAlgoSK'),
        algo_psets = cms.VPSet(cms.PSet(
            dPhiFineMax = cms.double(0.025),
            verboseInfo = cms.untracked.bool(True),
            chi2Max = cms.double(99999.0),
            dPhiMax = cms.double(0.003),
            wideSeg = cms.double(3.0),
            minLayersApart = cms.int32(2),
            dRPhiFineMax = cms.double(8.0),
            dRPhiMax = cms.double(8.0)
        ), 
            cms.PSet(
                dPhiFineMax = cms.double(0.025),
                verboseInfo = cms.untracked.bool(True),
                chi2Max = cms.double(99999.0),
                dPhiMax = cms.double(0.025),
                wideSeg = cms.double(3.0),
                minLayersApart = cms.int32(2),
                dRPhiFineMax = cms.double(3.0),
                dRPhiMax = cms.double(8.0)
            ))
    ), 
        cms.PSet(
            parameters_per_chamber_type = cms.vint32(2, 1, 1, 1, 1, 
                1, 1, 1, 1),
            algo_name = cms.string('CSCSegAlgoTC'),
            algo_psets = cms.VPSet(cms.PSet(
                dPhiFineMax = cms.double(0.02),
                verboseInfo = cms.untracked.bool(True),
                SegmentSorting = cms.int32(1),
                chi2Max = cms.double(6000.0),
                dPhiMax = cms.double(0.003),
                chi2ndfProbMin = cms.double(0.0001),
                minLayersApart = cms.int32(2),
                dRPhiFineMax = cms.double(6.0),
                dRPhiMax = cms.double(1.2)
            ), 
                cms.PSet(
                    dPhiFineMax = cms.double(0.013),
                    verboseInfo = cms.untracked.bool(True),
                    SegmentSorting = cms.int32(1),
                    chi2Max = cms.double(6000.0),
                    dPhiMax = cms.double(0.00198),
                    chi2ndfProbMin = cms.double(0.0001),
                    minLayersApart = cms.int32(2),
                    dRPhiFineMax = cms.double(3.0),
                    dRPhiMax = cms.double(0.6)
                ))
        ), 
        cms.PSet(
            parameters_per_chamber_type = cms.vint32(3, 1, 2, 2, 1, 
                2, 1, 2, 1),
            algo_name = cms.string('CSCSegAlgoDF'),
            algo_psets = cms.VPSet(cms.PSet(
                minHitsPerSegment = cms.untracked.int32(3),
                dPhiFineMax = cms.untracked.double(0.025),
                dXclusBoxMax = cms.untracked.double(4.0),
                tanThetaMax = cms.untracked.double(1.2),
                BrutePruning = cms.untracked.bool(False),
                preClustering = cms.untracked.bool(False),
                tanPhiMax = cms.untracked.double(0.5),
                nSigmaFromSegment = cms.untracked.double(5.0),
                minLayersApart = cms.untracked.int32(2),
                dRPhiFineMax = cms.untracked.double(8.0),
                CSCSegmentDebug = cms.untracked.bool(False),
                Pruning = cms.untracked.bool(False),
                dYclusBoxMax = cms.untracked.double(8.0)
            ), 
                cms.PSet(
                    minHitsPerSegment = cms.untracked.int32(3),
                    dPhiFineMax = cms.untracked.double(0.025),
                    dXclusBoxMax = cms.untracked.double(4.0),
                    tanThetaMax = cms.untracked.double(2.0),
                    BrutePruning = cms.untracked.bool(False),
                    preClustering = cms.untracked.bool(False),
                    tanPhiMax = cms.untracked.double(0.8),
                    nSigmaFromSegment = cms.untracked.double(5.0),
                    minLayersApart = cms.untracked.int32(2),
                    dRPhiFineMax = cms.untracked.double(12.0),
                    CSCSegmentDebug = cms.untracked.bool(False),
                    Pruning = cms.untracked.bool(False),
                    dYclusBoxMax = cms.untracked.double(8.0)
                ), 
                cms.PSet(
                    minHitsPerSegment = cms.untracked.int32(3),
                    dPhiFineMax = cms.untracked.double(0.025),
                    dXclusBoxMax = cms.untracked.double(4.0),
                    tanThetaMax = cms.untracked.double(1.2),
                    BrutePruning = cms.untracked.bool(False),
                    preClustering = cms.untracked.bool(False),
                    tanPhiMax = cms.untracked.double(0.5),
                    nSigmaFromSegment = cms.untracked.double(5.0),
                    minLayersApart = cms.untracked.int32(2),
                    dRPhiFineMax = cms.untracked.double(8.0),
                    CSCSegmentDebug = cms.untracked.bool(False),
                    Pruning = cms.untracked.bool(False),
                    dYclusBoxMax = cms.untracked.double(8.0)
                ))
        ), 
        cms.PSet(
            parameters_per_chamber_type = cms.vint32(2, 1, 1, 1, 1, 
                1, 1, 1, 1),
            algo_name = cms.string('CSCSegAlgoST'),
            algo_psets = cms.VPSet(cms.PSet(
                minHitsPerSegment = cms.untracked.int32(3),
                dXclusBoxMax = cms.untracked.double(4.0),
                BrutePruning = cms.untracked.bool(False),
                preClustering = cms.untracked.bool(True),
                maxRecHitsInCluster = cms.untracked.int32(20),
                onlyBestSegment = cms.untracked.bool(False),
                CSCDebug = cms.untracked.bool(False),
                Pruning = cms.untracked.bool(False),
                dYclusBoxMax = cms.untracked.double(8.0)
            ), 
                cms.PSet(
                    minHitsPerSegment = cms.untracked.int32(3),
                    dXclusBoxMax = cms.untracked.double(4.0),
                    BrutePruning = cms.untracked.bool(False),
                    preClustering = cms.untracked.bool(True),
                    maxRecHitsInCluster = cms.untracked.int32(20),
                    onlyBestSegment = cms.untracked.bool(False),
                    CSCDebug = cms.untracked.bool(False),
                    Pruning = cms.untracked.bool(False),
                    dYclusBoxMax = cms.untracked.double(8.0)
                ))
        ))
)


 process.heepElectronTkIsolation = cms.EDProducer("EgammaElectronTkIsolationProducer",
    absolut = cms.bool(True),
    trackProducer = cms.InputTag("ctfWithMaterialTracks"),
    intRadius = cms.double(0.1),
    electronProducer = cms.InputTag("pixelMatchGsfElectrons"),
    extRadius = cms.double(0.2),
    ptMin = cms.double(1.5),
    maxVtxDist = cms.double(0.1)
)


 process.heepHcalIsolation = cms.EDProducer("EgammaHcalIsolationProducer",
    absolut = cms.bool(True),
    intRadius = cms.double(0.15),
    emObjectProducer = cms.InputTag("pixelMatchGsfElectrons"),
    extRadius = cms.double(0.3),
    etMin = cms.double(0.0),
    hcalRecHitProducer = cms.InputTag("hbhereco")
)


 process.ecalRegionalEgammaWeightUncalibRecHit = cms.EDProducer("EcalWeightUncalibRecHitProducer",
    EBdigiCollection = cms.InputTag("ecalRegionalEgammaDigis","ebDigis"),
    EEhitCollection = cms.string('EcalUncalibRecHitsEE'),
    EEdigiCollection = cms.InputTag("ecalRegionalEgammaDigis","eeDigis"),
    EBhitCollection = cms.string('EcalUncalibRecHitsEB')
)


 process.siPixelClusters = cms.EDProducer("SiPixelClusterProducer",
    PedestalValue = cms.double(28.2),
    src = cms.InputTag("siPixelDigis"),
    ChannelThreshold = cms.int32(2500),
    VerbosityLevel = cms.untracked.int32(1),
    MissCalibrate = cms.untracked.bool(True),
    MaxGain = cms.double(3.0),
    UseCalibDataFromDB = cms.bool(False),
    MaxPed = cms.double(35.0),
    MinPed = cms.double(25.0),
    MinGain = cms.double(2.5),
    SeedThreshold = cms.int32(3000),
    GainValue = cms.double(2.8),
    ClusterThreshold = cms.double(5050.0)
)


 process.pixelTracks = cms.EDProducer("PixelTrackProducer",
    FitterPSet = cms.PSet(
        ComponentName = cms.string('PixelFitterByHelixProjections'),
        TTRHBuilder = cms.string('PixelTTRHBuilderWithoutAngle')
    ),
    FilterPSet = cms.PSet(
        chi2 = cms.double(1000.0),
        ComponentName = cms.string('PixelTrackFilterByKinematics'),
        ptMin = cms.double(0.0),
        tipMax = cms.double(1.0)
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            precise = cms.bool(True),
            originZPos = cms.double(0.0),
            originRadius = cms.double(0.2),
            ptMin = cms.double(0.9),
            originHalfLength = cms.double(15.9)
        )
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitTripletGenerator'),
        SeedingLayers = cms.string('PixelLayerTriplets'),
        GeneratorPSet = cms.PSet(
            useBending = cms.bool(True),
            useFixedPreFiltering = cms.bool(False),
            ComponentName = cms.string('PixelTripletHLTGenerator'),
            extraHitRPhitolerance = cms.double(0.06),
            useMultScattering = cms.bool(True),
            phiPreFiltering = cms.double(0.3),
            extraHitRZtolerance = cms.double(0.06)
        )
    ),
    CleanerPSet = cms.PSet(
        ComponentName = cms.string('PixelTrackCleanerBySharedHits')
    )
)


 process.ecalWeightUncalibRecHit = cms.EDProducer("EcalWeightUncalibRecHitProducer",
    EBdigiCollection = cms.InputTag("ecalDigis","ebDigis"),
    EEhitCollection = cms.string('EcalUncalibRecHitsEE'),
    EEdigiCollection = cms.InputTag("ecalDigis","eeDigis"),
    EBhitCollection = cms.string('EcalUncalibRecHitsEB')
)


 process.ecalRegionalEgammaFEDs = cms.EDProducer("EcalListOfFEDSProducer",
    Muon = cms.untracked.bool(False),
    Ptmin_iso = cms.untracked.double(5.0),
    EM_l1TagIsolated = cms.untracked.InputTag("l1extraParticles","Isolated"),
    OutputLabel = cms.untracked.string(''),
    EM_regionEtaMargin = cms.untracked.double(0.25),
    EM_doNonIsolated = cms.untracked.bool(True),
    EM_l1TagNonIsolated = cms.untracked.InputTag("l1extraParticles","NonIsolated"),
    debug = cms.untracked.bool(False),
    EM_regionPhiMargin = cms.untracked.double(0.4),
    EGamma = cms.untracked.bool(True),
    Ptmin_noniso = cms.untracked.double(5.0),
    EM_doIsolated = cms.untracked.bool(True)
)


 process.ecalRegionalMuonsWeightUncalibRecHit = cms.EDProducer("EcalWeightUncalibRecHitProducer",
    EBdigiCollection = cms.InputTag("ecalRegionalMuonsDigis","ebDigis"),
    EEhitCollection = cms.string('EcalUncalibRecHitsEE'),
    EEdigiCollection = cms.InputTag("ecalRegionalMuonsDigis","eeDigis"),
    EBhitCollection = cms.string('EcalUncalibRecHitsEB')
)


 process.heepElectronTkNumIsolation = cms.EDProducer("EgammaElectronTkNumIsolationProducer",
    trackProducer = cms.InputTag("ctfWithMaterialTracks"),
    intRadius = cms.double(0.1),
    electronProducer = cms.InputTag("pixelMatchGsfElectrons"),
    extRadius = cms.double(0.2),
    ptMin = cms.double(1.5),
    maxVtxDist = cms.double(0.1)
)


 process.ecalRegionalMuonsRecHit = cms.EDProducer("EcalRecHitProducer",
    EErechitCollection = cms.string('EcalRecHitsEE'),
    EEuncalibRecHitCollection = cms.InputTag("ecalRegionalMuonsWeightUncalibRecHit","EcalUncalibRecHitsEE"),
    EBuncalibRecHitCollection = cms.InputTag("ecalRegionalMuonsWeightUncalibRecHit","EcalUncalibRecHitsEB"),
    EBrechitCollection = cms.string('EcalRecHitsEB')
)


 process.heepEcalIsolation = cms.EDProducer("EgammaEcalIsolationProducer",
    absolut = cms.bool(True),
    basicClusterProducer = cms.InputTag("egammaBasicClusterMerger"),
    superClusterProducer = cms.InputTag("egammaSuperClusterMerger"),
    extRadius = cms.double(0.3),
    etMin = cms.double(0.0),
    emObjectProducer = cms.InputTag("pixelMatchGsfElectrons")
)


 process.ctfWithMaterialTracks = cms.EDProducer("TrackProducer",
    src = cms.string('ckfTrackCandidates'),
    producer = cms.string(''),
    Fitter = cms.string('KFFittingSmoother'),
    TrajectoryInEvent = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    Propagator = cms.string('PropagatorWithMaterial')
)


 process.rsWithMaterialTracks = cms.EDProducer("TrackProducer",
    src = cms.string('rsTrackCandidates'),
    producer = cms.string(''),
    Fitter = cms.string('KFFittingSmoother'),
    TrajectoryInEvent = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    Propagator = cms.string('PropagatorWithMaterial')
)


 process.globalMixedSeeds = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            precise = cms.bool(True),
            originZPos = cms.double(0.0),
            originRadius = cms.double(0.2),
            ptMin = cms.double(0.9),
            originHalfLength = cms.double(15.9)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitPairGenerator'),
        SeedingLayers = cms.string('MixedLayerPairs')
    ),
    TTRHBuilder = cms.string('WithTrackAngle')
)


 process.ecalTriggerPrimitiveDigis = cms.EDProducer("EcalTrigPrimProducer",
    BarrelOnly = cms.bool(False),
    TTFHighEnergyEB = cms.double(1.0),
    InstanceEB = cms.string(''),
    InstanceEE = cms.string(''),
    TTFHighEnergyEE = cms.double(1.0),
    Famos = cms.bool(False),
    Debug = cms.bool(False),
    TcpOutput = cms.bool(False),
    TTFLowEnergyEE = cms.double(1.0),
    Label = cms.string('ecalUnsuppressedDigis'),
    TTFLowEnergyEB = cms.double(1.0)
)


 process.dt4DSegments = cms.EDProducer("DTRecSegment4DProducer",
    debug = cms.untracked.bool(False),
    Reco4DAlgoName = cms.string('DTCombinatorialPatternReco4D'),
    Reco4DAlgoConfig = cms.PSet(
        segmCleanerMode = cms.int32(1),
        nSharedHitsMax = cms.int32(2),
        debug = cms.untracked.bool(False),
        nUnSharedHitsMin = cms.int32(2),
        AllDTRecHits = cms.bool(True),
        Reco2DAlgoConfig = cms.PSet(
            segmCleanerMode = cms.int32(1),
            AlphaMaxPhi = cms.double(1.0),
            MaxAllowedHits = cms.uint32(50),
            nSharedHitsMax = cms.int32(2),
            AlphaMaxTheta = cms.double(0.1),
            debug = cms.untracked.bool(False),
            nUnSharedHitsMin = cms.int32(2),
            recAlgoConfig = cms.PSet(
                tTrigMode = cms.string('DTTTrigSyncFromDB'),
                minTime = cms.double(-3.0),
                interpolate = cms.bool(True),
                debug = cms.untracked.bool(False),
                tTrigModeConfig = cms.PSet(
                    vPropWire = cms.double(24.4),
                    doTOFCorrection = cms.bool(True),
                    tofCorrType = cms.int32(1),
                    kFactor = cms.double(-2.0),
                    wirePropCorrType = cms.int32(1),
                    doWirePropCorrection = cms.bool(True),
                    doT0Correction = cms.bool(True),
                    debug = cms.untracked.bool(False)
                ),
                maxTime = cms.double(415.0)
            ),
            recAlgo = cms.string('DTParametrizedDriftAlgo')
        ),
        Reco2DAlgoName = cms.string('DTCombinatorialPatternReco'),
        recAlgoConfig = cms.PSet(
            tTrigMode = cms.string('DTTTrigSyncFromDB'),
            minTime = cms.double(-3.0),
            interpolate = cms.bool(True),
            debug = cms.untracked.bool(False),
            tTrigModeConfig = cms.PSet(
                vPropWire = cms.double(24.4),
                doTOFCorrection = cms.bool(True),
                tofCorrType = cms.int32(1),
                kFactor = cms.double(-2.0),
                wirePropCorrType = cms.int32(1),
                doWirePropCorrection = cms.bool(True),
                doT0Correction = cms.bool(True),
                debug = cms.untracked.bool(False)
            ),
            maxTime = cms.double(415.0)
        ),
        recAlgo = cms.string('DTParametrizedDriftAlgo')
    ),
    recHits1DLabel = cms.string('dt1DRecHits'),
    recHits2DLabel = cms.string('dt2DSegments')
)


 process.pixelVertices = cms.EDProducer("PixelVertexProducer",
    WtAverage = cms.bool(True),
    ZOffset = cms.double(5.0),
    Verbosity = cms.int32(0),
    UseError = cms.bool(True),
    TrackCollection = cms.string('pixelTracks'),
    ZSeparation = cms.double(0.05),
    NTrkMin = cms.int32(2),
    Method2 = cms.bool(True),
    Finder = cms.string('DivisiveVertexFinder'),
    PtMin = cms.double(1.0)
)


 process.ecalPreshowerRecHit = cms.EDProducer("ESRecHitProducer",
    ESrechitCollection = cms.string('EcalRecHitsES'),
    ESdigiCollection = cms.InputTag("ecalPreshowerDigis")
)


 process.ecalRegionalMuonsFEDs = cms.EDProducer("EcalListOfFEDSProducer",
    Muon = cms.untracked.bool(True),
    MuonSource = cms.untracked.InputTag("l1extraParticles"),
    MU_regionPhiMargin = cms.untracked.double(1.0),
    OutputLabel = cms.untracked.string(''),
    Ptmin_muon = cms.untracked.double(0.0),
    debug = cms.untracked.bool(False),
    EGamma = cms.untracked.bool(False),
    MU_regionEtaMargin = cms.untracked.double(1.0)
)


 process.dt2DSegments = cms.EDProducer("DTRecSegment2DProducer",
    debug = cms.untracked.bool(False),
    Reco2DAlgoConfig = cms.PSet(
        segmCleanerMode = cms.int32(1),
        AlphaMaxPhi = cms.double(1.0),
        MaxAllowedHits = cms.uint32(50),
        nSharedHitsMax = cms.int32(2),
        AlphaMaxTheta = cms.double(0.1),
        debug = cms.untracked.bool(False),
        nUnSharedHitsMin = cms.int32(2),
        recAlgoConfig = cms.PSet(
            tTrigMode = cms.string('DTTTrigSyncFromDB'),
            minTime = cms.double(-3.0),
            interpolate = cms.bool(True),
            debug = cms.untracked.bool(False),
            tTrigModeConfig = cms.PSet(
                vPropWire = cms.double(24.4),
                doTOFCorrection = cms.bool(True),
                tofCorrType = cms.int32(1),
                kFactor = cms.double(-2.0),
                wirePropCorrType = cms.int32(1),
                doWirePropCorrection = cms.bool(True),
                doT0Correction = cms.bool(True),
                debug = cms.untracked.bool(False)
            ),
            maxTime = cms.double(415.0)
        ),
        recAlgo = cms.string('DTParametrizedDriftAlgo')
    ),
    recHits1DLabel = cms.string('dt1DRecHits'),
    Reco2DAlgoName = cms.string('DTCombinatorialPatternReco')
)


 process.heepRobustElectronId = cms.EDProducer("ElectronIDProducer",
    electronLabel = cms.string(''),
    doNeuralNet = cms.bool(False),
    doLikelihood = cms.bool(False),
    endcapClusterShapeAssociation = cms.InputTag("islandBasicClusters","islandEndcapShapeAssoc"),
    doPtdrId = cms.bool(False),
    doCutBased = cms.bool(True),
    electronIDAssociationLabel = cms.string(''),
    electronProducer = cms.string('pixelMatchGsfElectrons'),
    electronIDLabel = cms.string(''),
    barrelClusterShapeAssociation = cms.InputTag("hybridSuperClusters","hybridShapeAssoc"),
    algo_psets = cms.VPSet(cms.PSet(
        useBremFraction = cms.vint32(0, 0, 0),
        tightEleIDCuts = cms.PSet(
            invEMinusInvP = cms.vdouble(0.02, 0.02, 0.02, 0.02, 0.02, 
                0.02, 0.02, 0.02),
            EoverPInMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0),
            EoverPOutMin = cms.vdouble(0.6, 0.75, 0.75, 0.75, 0.5, 
                0.8, 0.5, 0.8),
            sigmaEtaEtaMin = cms.vdouble(0.005, 0.005, 0.005, 0.005, 0.008, 
                0.008, 0.008, 0.008),
            EoverPOutMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
                999.0, 999.0, 999.0),
            EoverPInMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
                999.0, 999.0, 999.0),
            deltaPhiOut = cms.vdouble(0.02, 999.0, 0.02, 999.0, 0.02, 
                999.0, 0.02, 999.0),
            sigmaEtaEtaMax = cms.vdouble(0.011, 0.011, 0.011, 0.011, 0.03, 
                0.03, 0.03, 0.022),
            deltaPhiIn = cms.vdouble(0.02, 0.03, 0.02, 0.04, 0.04, 
                0.04, 0.04, 0.05),
            HoverE = cms.vdouble(0.05, 0.05, 0.05, 0.05, 0.07, 
                0.07, 0.07, 0.07),
            sigmaPhiPhiMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0),
            bremFraction = cms.vdouble(0.0, 0.1, 0.1, 0.1, 0.0, 
                0.2, 0.2, 0.2),
            deltaEtaIn = cms.vdouble(0.004, 0.004, 0.004, 0.005, 0.005, 
                0.005, 0.005, 0.005),
            E9overE25 = cms.vdouble(0.8, 0.65, 0.75, 0.65, 0.8, 
                0.7, 0.7, 0.65),
            sigmaPhiPhiMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
                999.0, 999.0, 999.0)
        ),
        useSigmaPhiPhi = cms.vint32(0, 1, 0),
        useEoverPIn = cms.vint32(0, 1, 0),
        useHoverE = cms.vint32(1, 1, 1),
        useDeltaPhiOut = cms.vint32(0, 1, 1),
        useInvEMinusInvP = cms.vint32(0, 0, 0),
        useDeltaEtaIn = cms.vint32(1, 1, 1),
        useSigmaEtaEta = cms.vint32(0, 1, 1),
        useDeltaPhiIn = cms.vint32(1, 1, 1),
        looseEleIDCuts = cms.PSet(
            invEMinusInvP = cms.vdouble(0.02, 0.02, 0.02, 0.02, 0.02, 
                0.02, 0.02, 0.02),
            EoverPInMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0),
            EoverPOutMin = cms.vdouble(0.6, 1.7, 0.9, 0.5, 0.6, 
                1.7, 0.9, 0.5),
            sigmaEtaEtaMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0),
            EoverPOutMax = cms.vdouble(2.5, 999.0, 2.2, 999.0, 2.5, 
                999.0, 2.2, 999.0),
            EoverPInMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
                999.0, 999.0, 999.0),
            deltaPhiOut = cms.vdouble(0.011, 999.0, 999.0, 999.0, 0.02, 
                999.0, 999.0, 999.0),
            sigmaEtaEtaMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
                999.0, 999.0, 999.0),
            deltaPhiIn = cms.vdouble(0.06, 0.06, 0.06, 0.08, 0.06, 
                0.06, 0.06, 0.09),
            HoverE = cms.vdouble(0.09, 0.06, 0.07, 0.12, 0.09, 
                0.06, 0.07, 0.12),
            sigmaPhiPhiMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0),
            bremFraction = cms.vdouble(0.0, 0.1, 0.1, 0.1, 0.0, 
                0.2, 0.2, 0.2),
            deltaEtaIn = cms.vdouble(0.008, 0.008, 0.008, 0.009, 0.008, 
                0.008, 0.008, 0.009),
            E9overE25 = cms.vdouble(0.7, 0.7, 0.7, 0.5, 0.8, 
                0.8, 0.8, 0.5),
            sigmaPhiPhiMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
                999.0, 999.0, 999.0)
        ),
        useE9overE25 = cms.vint32(1, 1, 1),
        useEoverPOut = cms.vint32(1, 1, 1),
        mediumEleIDCuts = cms.PSet(
            invEMinusInvP = cms.vdouble(0.02, 0.02, 0.02, 0.02, 0.02, 
                0.02, 0.02, 0.02),
            EoverPInMin = cms.vdouble(0.9, 0.9, 0.9, 0.6, 0.9, 
                0.9, 0.9, 0.7),
            EoverPOutMin = cms.vdouble(0.6, 1.8, 1.0, 0.75, 0.6, 
                1.5, 1.0, 0.8),
            sigmaEtaEtaMin = cms.vdouble(0.005, 0.005, 0.005, 0.005, 0.008, 
                0.008, 0.008, 0.0),
            EoverPOutMax = cms.vdouble(2.5, 999.0, 999.0, 999.0, 2.0, 
                999.0, 999.0, 999.0),
            EoverPInMax = cms.vdouble(1.3, 1.2, 1.3, 999.0, 999.0, 
                999.0, 999.0, 999.0),
            deltaPhiOut = cms.vdouble(0.011, 999.0, 999.0, 999.0, 0.02, 
                999.0, 999.0, 999.0),
            sigmaEtaEtaMax = cms.vdouble(0.011, 0.011, 0.011, 0.011, 0.022, 
                0.022, 0.022, 0.3),
            deltaPhiIn = cms.vdouble(0.04, 0.07, 0.04, 0.08, 0.06, 
                0.07, 0.06, 0.07),
            HoverE = cms.vdouble(0.06, 0.05, 0.06, 0.14, 0.1, 
                0.1, 0.1, 0.12),
            sigmaPhiPhiMin = cms.vdouble(0.005, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0),
            bremFraction = cms.vdouble(0.0, 0.1, 0.1, 0.1, 0.0, 
                0.2, 0.2, 0.2),
            deltaEtaIn = cms.vdouble(0.004, 0.006, 0.005, 0.007, 0.007, 
                0.008, 0.007, 0.008),
            E9overE25 = cms.vdouble(0.7, 0.75, 0.8, 0.0, 0.85, 
                0.75, 0.8, 0.0),
            sigmaPhiPhiMax = cms.vdouble(0.015, 999.0, 999.0, 999.0, 0.02, 
                999.0, 999.0, 999.0)
        ),
        electronQuality = cms.string('loose')
    ), 
        cms.PSet(
            robustEleIDCuts = cms.PSet(
                barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005),
                endcap = cms.vdouble(0.08, 0.0275, 0.092, 0.007)
            ),
            tightEleIDCuts = cms.PSet(
                eSeedOverPinMax = cms.vdouble(99999.0, 99999.0, 99999.0, 99999.0, 99999.0, 
                    99999.0, 99999.0, 99999.0),
                eSeedOverPinMin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
                    0.83, 0.0, 0.0),
                deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
                    0.035, 0.065, 0.092),
                hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
                    0.037, 0.05, 0.0),
                sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
                    0.0252, 0.026, 0.0),
                deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
                    0.0055, 0.0075, 0.0)
            ),
            electronQuality = cms.string('robust'),
            looseEleIDCuts = cms.PSet(
                deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
                    0.03, 0.092, 0.092),
                hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
                    0.12, 0.15, 0.0),
                sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
                    0.0265, 0.0265, 0.0),
                deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
                    0.0068, 0.01, 0.0),
                eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
                    0.85, 0.0, 0.0)
            )
        ))
)


 process.rpcRecHits = cms.EDProducer("RPCRecHitProducer",
    recAlgo = cms.string('RPCRecHitStandardAlgo'),
    recAlgoConfig = cms.PSet(

    ),
    rpcDigiLabel = cms.string('muonRPCDigis')
)


 process.ecalRegionalEgammaRecHit = cms.EDProducer("EcalRecHitProducer",
    EErechitCollection = cms.string('EcalRecHitsEE'),
    EEuncalibRecHitCollection = cms.InputTag("ecalRegionalEgammaWeightUncalibRecHit","EcalUncalibRecHitsEE"),
    EBuncalibRecHitCollection = cms.InputTag("ecalRegionalEgammaWeightUncalibRecHit","EcalUncalibRecHitsEB"),
    EBrechitCollection = cms.string('EcalRecHitsEB')
)


 process.globalPixelSeeds = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            precise = cms.bool(True),
            originZPos = cms.double(0.0),
            originRadius = cms.double(0.2),
            ptMin = cms.double(0.9),
            originHalfLength = cms.double(15.9)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitPairGenerator'),
        SeedingLayers = cms.string('PixelLayerPairs')
    ),
    TTRHBuilder = cms.string('WithTrackAngle')
)


 process.ecalRecHit = cms.EDProducer("EcalRecHitProducer",
    EErechitCollection = cms.string('EcalRecHitsEE'),
    EEuncalibRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit","EcalUncalibRecHitsEE"),
    EBuncalibRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit","EcalUncalibRecHitsEB"),
    EBrechitCollection = cms.string('EcalRecHitsEB')
)


 process.dt1DRecHits = cms.EDProducer("DTRecHitProducer",
    debug = cms.untracked.bool(False),
    recAlgo = cms.string('DTParametrizedDriftAlgo'),
    dtDigiLabel = cms.string('muonDTDigis'),
    recAlgoConfig = cms.PSet(
        tTrigMode = cms.string('DTTTrigSyncFromDB'),
        minTime = cms.double(-3.0),
        interpolate = cms.bool(True),
        debug = cms.untracked.bool(False),
        tTrigModeConfig = cms.PSet(
            vPropWire = cms.double(24.4),
            doTOFCorrection = cms.bool(True),
            tofCorrType = cms.int32(1),
            kFactor = cms.double(-2.0),
            wirePropCorrType = cms.int32(1),
            doWirePropCorrection = cms.bool(True),
            doT0Correction = cms.bool(True),
            debug = cms.untracked.bool(False)
        ),
        maxTime = cms.double(415.0)
    )
)


 process.csc2DRecHits = cms.EDProducer("CSCRecHitBProducer",
    CSCStripMaxDistance = cms.untracked.double(5.0),
    CSCSegmentPerChamberMax = cms.untracked.int32(3),
    CSCWireMaxDistance = cms.untracked.double(2.0),
    CSCStripxtalksOffset = cms.untracked.double(0.0),
    CSCUseCalibrations = cms.untracked.bool(True),
    CSCStripSegmentDeltaT = cms.untracked.int32(1),
    CSCminGattiStepSize = cms.untracked.double(0.0001),
    CSCIsRunningOnData = cms.untracked.bool(False),
    CSCUseGattiFit = cms.untracked.bool(True),
    CSCStripClusterChargeCut = cms.untracked.double(30.0),
    CSCuseCleanWireCollection = cms.untracked.bool(False),
    CSCchamberIdPrefix = cms.untracked.int32(0),
    theMappingFile = cms.FileInPath('CondFormats/CSCObjects/data/csc_slice_test_map.txt'),
    CSCstripWireDeltaTime = cms.untracked.int32(1),
    CSCWireminLayersApart = cms.untracked.int32(2),
    CSCminGattiError = cms.untracked.double(0.015),
    CSCMaxGattiChi2 = cms.untracked.double(200.0),
    CSCStripClusterSize = cms.untracked.int32(3),
    CSCStripPeakThreshold = cms.untracked.double(10.0),
    CSCminStripHitsPerSegment = cms.untracked.int32(3),
    CSCStripCloseToSegment = cms.untracked.double(8.0),
    CSCproduce1DHits = cms.untracked.bool(False),
    CSCuseLeftOverStripHits = cms.untracked.bool(False),
    CSCuseCleanStripCollection = cms.untracked.bool(False),
    CSCuseLeftOverWireHits = cms.untracked.bool(False),
    CSCDebug = cms.untracked.bool(False),
    CSCWireSegmentDeltaT = cms.untracked.int32(1),
    CSCuseStripHitsFromFits = cms.untracked.bool(False),
    CSCStripxtalksSystematics = cms.untracked.double(0.0),
    CSCStripminLayersApart = cms.untracked.int32(2),
    CSCuseWireHitsFromFits = cms.untracked.bool(False),
    CSCWireClusterDeltaT = cms.untracked.int32(1),
    CSCminWireHitsPerSegment = cms.untracked.int32(3),
    CSCWireClusterMaxSize = cms.untracked.int32(999),
    CSCStripDigiProducer = cms.string('muonCSCDigis'),
    CSCWireDigiProducer = cms.string('muonCSCDigis'),
    CSCCalibrationSystematics = cms.untracked.double(0.0)
)


 process.horeco = cms.EDFilter("HcalSimpleReconstructor",
    correctionPhaseNS = cms.double(13.0),
    digiLabel = cms.InputTag("hcalZeroSuppressedDigis"),
    samplesToAdd = cms.int32(4),
    Subdetector = cms.string('HO'),
    correctForTimeslew = cms.bool(True),
    correctForPhaseContainment = cms.bool(True),
    firstSample = cms.int32(4)
)


 process.ecalPreshowerDigis = cms.EDFilter("ESRawToDigi",
    debugMode = cms.untracked.bool(False),
    InstanceES = cms.string(''),
    ESdigiCollection = cms.string(''),
    Label = cms.string('rawDataCollector')
)


 process.l1SeedBegin = cms.EDFilter("HLTLevel1Seed",
    L1ExtraParticleMap = cms.InputTag("l1extraParticleMap"),
    L1GTReadoutRecord = cms.InputTag("l1extraParticleMap"),
    L1ExtraCollections = cms.InputTag("l1extraParticles"),
    L1Seeds = cms.vstring(),
    andOr = cms.bool(True)
)


 process.hbhereco = cms.EDFilter("HcalSimpleReconstructor",
    correctionPhaseNS = cms.double(13.0),
    digiLabel = cms.InputTag("hcalZeroSuppressedDigis"),
    samplesToAdd = cms.int32(4),
    Subdetector = cms.string('HBHE'),
    correctForTimeslew = cms.bool(True),
    correctForPhaseContainment = cms.bool(True),
    firstSample = cms.int32(4)
)


 process.ecalRegionalEgammaDigis = cms.EDFilter("EcalRawToDigiDev",
    tccUnpacking = cms.untracked.bool(True),
    FedLabel = cms.untracked.string('ecalRegionalEgammaFEDs'),
    srpUnpacking = cms.untracked.bool(False),
    syncCheck = cms.untracked.bool(False),
    headerUnpacking = cms.untracked.bool(False),
    DCCMapFile = cms.untracked.string('EventFilter/EcalRawToDigiDev/data/DCCMap.txt'),
    eventPut = cms.untracked.bool(True),
    feUnpacking = cms.untracked.bool(True),
    InputLabel = cms.untracked.string('rawDataCollector'),
    DoRegional = cms.untracked.bool(True),
    memUnpacking = cms.untracked.bool(False)
)


 process.caloTowersForMuons = cms.EDFilter("CaloTowerCandidateCreator",
    src = cms.InputTag("towerMakerForMuons"),
    minimumEt = cms.double(-1.0),
    minimumE = cms.double(-1.0)
)


 process.egammaSuperClusterMerger = cms.EDFilter("SuperClusterMerger",
    src = cms.VInputTag(cms.InputTag("islandSuperClusters","islandBarrelSuperClusters"), cms.InputTag("correctedEndcapSuperClustersWithPreshower"))
)


 process.siPixelDigis = cms.EDFilter("SiPixelRawToDigi",
    Timing = cms.untracked.bool(False),
    InputLabel = cms.untracked.string('rawDataCollector'),
    IncludeErrors = cms.untracked.bool(False)
)


 process.muonCSCDigis = cms.EDFilter("CSCDCCUnpacker",
    PrintEventNumber = cms.untracked.bool(False),
    theMappingFile = cms.FileInPath('CondFormats/CSCObjects/data/csc_slice_test_map.txt'),
    UseExaminer = cms.untracked.bool(False),
    ErrorMask = cms.untracked.uint32(3754946559),
    InputObjects = cms.InputTag("rawDataCollector"),
    ExaminerMask = cms.untracked.uint32(133921782),
    UnpackStatusDigis = cms.untracked.bool(False),
    isMTCCData = cms.untracked.bool(False),
    Debug = cms.untracked.bool(False)
)


 process.hlt2GetRaw = cms.EDFilter("HLTGetRaw",
    RawDataCollection = cms.InputTag("rawDataCollector")
)


 process.hcalZeroSuppressedDigis = cms.EDFilter("HcalRawToDigi",
    FilterDataQuality = cms.bool(True),
    lastSample = cms.int32(9),
    InputLabel = cms.InputTag("rawDataCollector"),
    ComplainEmptyData = cms.untracked.bool(False),
    UnpackCalib = cms.untracked.bool(False),
    ExceptionEmptyData = cms.untracked.bool(False),
    firstSample = cms.int32(0)
)


 process.hltMakeSummaryObjects = cms.EDFilter("HLTMakeSummaryObjects")


 process.SiStripRawToClustersFacility = cms.EDFilter("SiStripRawToClusters",
    ProductLabel = cms.untracked.string('rawDataCollector'),
    ChannelThreshold = cms.untracked.double(2.0),
    MaxHolesInCluster = cms.untracked.uint32(0),
    ClusterizerAlgorithm = cms.untracked.string('ThreeThreshold'),
    SeedThreshold = cms.untracked.double(3.0),
    ProductInstance = cms.untracked.string(''),
    ClusterThreshold = cms.untracked.double(5.0)
)


 process.towerMaker = cms.EDFilter("CaloTowersCreator",
    EBSumThreshold = cms.double(0.2),
    EBWeight = cms.double(1.0),
    hfInput = cms.InputTag("hfreco"),
    AllowMissingInputs = cms.untracked.bool(False),
    EESumThreshold = cms.double(0.45),
    HOThreshold = cms.double(1.1),
    HBThreshold = cms.double(0.9),
    EcutTower = cms.double(-1000.0),
    HcalThreshold = cms.double(-1000.0),
    HEDWeight = cms.double(1.0),
    EEWeight = cms.double(1.0),
    UseHO = cms.bool(True),
    HF1Weight = cms.double(1.0),
    HOWeight = cms.double(1.0),
    HESWeight = cms.double(1.0),
    hbheInput = cms.InputTag("hbhereco"),
    HF2Weight = cms.double(1.0),
    HF2Threshold = cms.double(1.8),
    EEThreshold = cms.double(0.45),
    HESThreshold = cms.double(1.4),
    hoInput = cms.InputTag("horeco"),
    HF1Threshold = cms.double(1.2),
    HEDThreshold = cms.double(1.4),
    EBThreshold = cms.double(0.09),
    ecalInputs = cms.VInputTag(cms.InputTag("ecalRecHit","EcalRecHitsEB"), cms.InputTag("ecalRecHit","EcalRecHitsEE")),
    HBWeight = cms.double(1.0)
)


 process.l1tTrigReport = cms.EDFilter("L1TrigReport",
    L1ExtraParticleMap = cms.InputTag("l1extraParticleMap"),
    L1GTReadoutRecord = cms.InputTag("l1extraParticleMap")
)


 process.egammaBasicClusterMerger = cms.EDFilter("BasicClusterMerger",
    src = cms.VInputTag(cms.InputTag("islandBasicClusters","islandBarrelBasicClusters"), cms.InputTag("islandBasicClusters","islandEndcapBasicClusters"))
)


 process.caloTowers = cms.EDFilter("CaloTowerCandidateCreator",
    src = cms.InputTag("towerMaker"),
    minimumEt = cms.double(-1.0),
    minimumE = cms.double(-1.0)
)


 process.roadSearchClouds = cms.EDFilter("RoadSearchCloudMaker",
    MinimalFractionOfUsedLayersPerCloud = cms.double(0.5),
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    StraightLineNoBeamSpotCloud = cms.bool(False),
    UsePixelsinRS = cms.bool(True),
    SeedProducer = cms.InputTag("roadSearchSeeds"),
    DoCloudCleaning = cms.bool(True),
    IncreaseMaxNumberOfConsecutiveMissedLayersPerCloud = cms.uint32(4),
    rphiStripRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    UseStereoRecHits = cms.bool(False),
    ZPhiRoadSize = cms.double(0.06),
    MaximalFractionOfConsecutiveMissedLayersPerCloud = cms.double(0.15),
    stereoStripRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    MaximalFractionOfMissedLayersPerCloud = cms.double(0.3),
    scalefactorRoadSeedWindow = cms.double(1.5),
    MaxDetHitsInCloudPerDetId = cms.uint32(32),
    IncreaseMaxNumberOfMissedLayersPerCloud = cms.uint32(3),
    RoadsLabel = cms.string(''),
    MaxRecHitsInCloud = cms.int32(100),
    UseRphiRecHits = cms.bool(False),
    MergingFraction = cms.double(0.8),
    RPhiRoadSize = cms.double(0.02),
    matchedStripRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    MinimumHalfRoad = cms.double(0.55)
)


 process.roadSearchSeeds = cms.EDFilter("RoadSearchSeedFinder",
    OuterSeedRecHitAccessMode = cms.string('RPHI'),
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    MaximalEndcapImpactParameter = cms.double(1.2),
    MergeSeedsCenterCut_C = cms.double(0.4),
    MergeSeedsCenterCut_B = cms.double(0.25),
    MergeSeedsCenterCut_A = cms.double(0.05),
    MergeSeedsDifferentHitsCut = cms.uint32(1),
    rphiStripRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    MaximalBarrelImpactParameter = cms.double(0.2),
    stereoStripRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    RoadsLabel = cms.string(''),
    OuterSeedRecHitAccessUseStereo = cms.bool(False),
    MinimalReconstructedTransverseMomentum = cms.double(1.5),
    PhiRangeForDetIdLookupInRings = cms.double(0.5),
    Mode = cms.string('STANDARD'),
    MergeSeedsRadiusCut_A = cms.double(0.05),
    InnerSeedRecHitAccessMode = cms.string('RPHI'),
    InnerSeedRecHitAccessUseStereo = cms.bool(False),
    OuterSeedRecHitAccessUseRPhi = cms.bool(False),
    MergeSeedsRadiusCut_B = cms.double(0.25),
    MergeSeedsRadiusCut_C = cms.double(0.4),
    matchedStripRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    InnerSeedRecHitAccessUseRPhi = cms.bool(False)
)


 process.rsTrackCandidates = cms.EDFilter("RoadSearchTrackCandidateMaker",
    NumHitCut = cms.int32(5),
    HitChi2Cut = cms.double(30.0),
    StraightLineNoBeamSpotCloud = cms.bool(False),
    MeasurementTrackerName = cms.string(''),
    MinimumChunkLength = cms.int32(7),
    TTRHBuilder = cms.string('WithTrackAngle'),
    nFoundMin = cms.int32(4),
    CloudProducer = cms.InputTag("roadSearchClouds")
)


 process.ecalDigis = cms.EDFilter("EcalRawToDigiDev",
    tccUnpacking = cms.untracked.bool(True),
    FedLabel = cms.untracked.string('ecalRegionalFEDs'),
    srpUnpacking = cms.untracked.bool(False),
    syncCheck = cms.untracked.bool(False),
    headerUnpacking = cms.untracked.bool(False),
    DCCMapFile = cms.untracked.string('EventFilter/EcalRawToDigiDev/data/DCCMap.txt'),
    eventPut = cms.untracked.bool(True),
    feUnpacking = cms.untracked.bool(True),
    InputLabel = cms.untracked.string('rawDataCollector'),
    DoRegional = cms.untracked.bool(False),
    memUnpacking = cms.untracked.bool(False)
)


 process.heepSelector = cms.EDFilter("HEEPSelector",
    trackIsolation = cms.bool(True),
    classification = cms.bool(True),
    ecalIsolation = cms.bool(True),
    kinematicsEt = cms.bool(True),
    cleanDuplicates = cms.int32(1),
    hcalIsolation = cms.bool(True),
    inputProducer = cms.InputTag("pixelMatchGsfElectrons"),
    trackNumIsolation = cms.bool(True),
    debug = cms.bool(False),
    kinematicsEta = cms.bool(True),
    HEEPHelper = cms.PSet(
        ebMax = cms.double(1.4442),
        trackIsolation_p0 = cms.double(0.2),
        hcalIsolation_p1 = cms.double(0.005),
        hcalIsolation_p0 = cms.double(4.0),
        trackNumIsolation = cms.bool(True),
        endcapClusterShapeAssociation = cms.InputTag("islandBasicClusters","islandEndcapShapeAssoc"),
        electronIDProducer = cms.InputTag("heepRobustElectronId"),
        hltName = cms.InputTag("hltL1NonIsoSingleElectronTrackIsolFilter"),
        hcalIsolationProducer = cms.InputTag("heepHcalIsolation"),
        trackNumIsolationProducer = cms.InputTag("heepElectronTkNumIsolation"),
        ecalIsolationProducer = cms.InputTag("heepEcalIsolation"),
        trackIsolationProducer = cms.InputTag("heepElectronTkIsolation"),
        ecalIsolation_p0 = cms.double(6.0),
        eeMax = cms.double(2.5),
        barrelClusterShapeAssociation = cms.InputTag("hybridSuperClusters","hybridShapeAssoc"),
        etCut = cms.double(20.0),
        badClassifications = cms.vint32(40),
        eeMin = cms.double(1.56),
        trackNumIsolation_p0 = cms.double(4.0),
        ecalIsolation_p1 = cms.double(0.01)
    ),
    electronID = cms.bool(True)
)


 process.hltTrigReport = cms.EDFilter("HLTrigReport",
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
)


 process.ecalRegionalMuonsDigis = cms.EDFilter("EcalRawToDigiDev",
    tccUnpacking = cms.untracked.bool(True),
    FedLabel = cms.untracked.string('ecalRegionalMuonsFEDs'),
    srpUnpacking = cms.untracked.bool(False),
    syncCheck = cms.untracked.bool(False),
    headerUnpacking = cms.untracked.bool(False),
    DCCMapFile = cms.untracked.string('EventFilter/EcalRawToDigiDev/data/DCCMap.txt'),
    eventPut = cms.untracked.bool(True),
    feUnpacking = cms.untracked.bool(True),
    InputLabel = cms.untracked.string('rawDataCollector'),
    DoRegional = cms.untracked.bool(True),
    memUnpacking = cms.untracked.bool(False)
)


 process.towerMakerForMuons = cms.EDFilter("CaloTowersCreator",
    EBSumThreshold = cms.double(0.2),
    EBWeight = cms.double(1.0),
    hfInput = cms.InputTag("hfreco"),
    AllowMissingInputs = cms.untracked.bool(False),
    EESumThreshold = cms.double(0.45),
    HOThreshold = cms.double(1.1),
    HBThreshold = cms.double(0.9),
    EcutTower = cms.double(-1000.0),
    HcalThreshold = cms.double(-1000.0),
    HEDWeight = cms.double(1.0),
    EEWeight = cms.double(1.0),
    UseHO = cms.bool(True),
    HF1Weight = cms.double(1.0),
    HOWeight = cms.double(1.0),
    HESWeight = cms.double(1.0),
    hbheInput = cms.InputTag("hbhereco"),
    HF2Weight = cms.double(1.0),
    HF2Threshold = cms.double(1.8),
    EEThreshold = cms.double(0.45),
    HESThreshold = cms.double(1.4),
    hoInput = cms.InputTag("horeco"),
    HF1Threshold = cms.double(1.2),
    HEDThreshold = cms.double(1.4),
    EBThreshold = cms.double(0.09),
    ecalInputs = cms.VInputTag(cms.InputTag("ecalRegionalMuonsRecHit","EcalRecHitsEB"), cms.InputTag("ecalRegionalMuonsRecHit","EcalRecHitsEE")),
    HBWeight = cms.double(1.0)
)


 process.ckfTrackCandidates = cms.EDFilter("CkfTrackCandidateMaker",
    NavigationPSet = cms.PSet(
        ComponentName = cms.string('SimpleNavigationSchool')
    ),
    TransientInitialStateEstimatorParameters = cms.PSet(
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SeedLabel = cms.string(''),
    SeedProducer = cms.string('globalMixedSeeds'),
    TrajectoryBuilder = cms.string('CkfTrajectoryBuilder')
)


 process.siPixelRecHits = cms.EDFilter("SiPixelRecHitConverter",
    src = cms.InputTag("siPixelClusters"),
    CPE = cms.string('Initial'),
    VerboseLevel = cms.untracked.int32(0),
    TanLorentzAnglePerTesla = cms.double(0.106),
    Alpha2Order = cms.bool(True),
    speed = cms.int32(0)
)


 process.siStripClusters = cms.EDFilter("SiStripRawToClustersRoI",
    All = cms.untracked.bool(True),
    Random = cms.untracked.bool(False),
    DeltaR = cms.untracked.double(0.5),
    InputModuleLabel = cms.untracked.string('SiStripRawToClustersFacility'),
    Electron = cms.untracked.bool(False)
)


 process.muonRPCDigis = cms.EDFilter("RPCUnpackingModule",
    InputLabel = cms.untracked.InputTag("rawDataCollector")
)


 process.hfreco = cms.EDFilter("HcalSimpleReconstructor",
    correctionPhaseNS = cms.double(0.0),
    digiLabel = cms.InputTag("hcalZeroSuppressedDigis"),
    samplesToAdd = cms.int32(1),
    Subdetector = cms.string('HF'),
    correctForTimeslew = cms.bool(False),
    correctForPhaseContainment = cms.bool(False),
    firstSample = cms.int32(3)
)


 process.heepEcalIsolationSequenceSpecial = cms.Sequence( process.heepEcalIsolation)


 process.recopixelvertexing = cms.Sequence( process.pixelTracks* process.pixelVertices)


 process.dtlocalreco_with_2DSegments = cms.Sequence( process.dt1DRecHits* process.dt2DSegments* process.dt4DSegments)


 process.hcalLocalRecoSequence = cms.Sequence( process.hbhereco+ process.hfreco+ process.horeco)


 process.hcalLocalRecoWithoutHO = cms.Sequence( process.hbhereco+ process.hfreco)


 process.ecalLocalRecoSequence = cms.Sequence( process.ecalWeightUncalibRecHit* process.ecalRecHit+ process.ecalPreshowerRecHit)


 process.calolocalreco = cms.Sequence( process.ecalLocalRecoSequence+ process.hcalLocalRecoSequence)


 process.csclocalreco = cms.Sequence( process.csc2DRecHits* process.cscSegments)


 process.rawtoclusters = cms.Sequence( process.SiStripRawToClustersFacility* process.siStripClusters)


 process.caloTowersRec = cms.Sequence( process.towerMaker* process.caloTowers)


 process.muonlocalreco_with_2DSegments = cms.Sequence( process.dtlocalreco_with_2DSegments+ process.csclocalreco+ process.rpcRecHits)


 process.doLocalRPC = cms.Sequence( process.muonRPCDigis* process.rpcRecHits)


 process.SiStripRawToClusters = cms.Sequence( process.rawtoclusters)


 process.rstracks = cms.Sequence( process.roadSearchSeeds* process.roadSearchClouds* process.rsTrackCandidates* process.rsWithMaterialTracks)


 process.hltBegin = cms.Sequence( process.hlt2GetRaw+ process.l1SeedBegin)


 process.dtlocalreco = cms.Sequence( process.dt1DRecHits* process.dt4DSegments)


 process.doLocalStrip = cms.Sequence( process.SiStripRawToClustersFacility* process.siStripClusters)


 process.ckftracks = cms.Sequence( process.globalMixedSeeds* process.globalPixelSeeds* process.ckfTrackCandidates* process.ctfWithMaterialTracks)


 process.EcalESRawToDigi = cms.Sequence( process.ecalPreshowerDigis)


 process.doLocalHcalWithoutHO = cms.Sequence( process.hcalZeroSuppressedDigis* process.hcalLocalRecoWithoutHO)


 process.ecalRegionalEgammaRecoSequence = cms.Sequence( process.ecalRegionalEgammaFEDs* process.ecalRegionalEgammaDigis* process.ecalRegionalEgammaWeightUncalibRecHit* process.ecalRegionalEgammaRecHit+ process.ecalPreshowerRecHit)


 process.doLocalCSC = cms.Sequence( process.muonCSCDigis* process.csclocalreco)


 process.CSCRawToDigi = cms.Sequence( process.muonCSCDigis)


 process.HcalRawToDigi = cms.Sequence( process.hcalZeroSuppressedDigis)


 process.ecalRegionalMuonsRecoSequence = cms.Sequence( process.ecalRegionalMuonsFEDs* process.ecalRegionalMuonsDigis* process.ecalRegionalMuonsWeightUncalibRecHit* process.ecalRegionalMuonsRecHit+ process.ecalPreshowerRecHit)


 process.doLocalEcal = cms.Sequence( process.ecalDigis+ process.ecalPreshowerDigis+ process.ecalLocalRecoSequence)


 process.SiPixelRawToDigi = cms.Sequence( process.siPixelDigis)


 process.doLocalCaloWithoutHO = cms.Sequence( process.doLocalEcal+ process.doLocalHcalWithoutHO)


 process.EcalRawToDigi = cms.Sequence( process.ecalDigis)


 process.RPCRawToDigi = cms.Sequence( process.muonRPCDigis)


 process.doLocalPixel = cms.Sequence( process.siPixelDigis* process.siPixelClusters* process.siPixelRecHits)


 process.heepEcalIsolationSequence = cms.Sequence( process.egammaSuperClusterMerger* process.egammaBasicClusterMerger* process.heepEcalIsolation)


 process.doRegionalMuonsEcal = cms.Sequence( process.ecalPreshowerDigis+ process.ecalRegionalMuonsRecoSequence)


 process.doLocalHcal = cms.Sequence( process.hcalZeroSuppressedDigis* process.hcalLocalRecoSequence)


 process.doLocalDT = cms.Sequence( process.dtlocalreco)


 process.heepSelectionProducersSpecial = cms.Sequence( process.heepRobustElectronId* process.heepEcalIsolationSequenceSpecial* process.heepHcalIsolation* process.heepElectronTkNumIsolation* process.heepElectronTkIsolation)


 process.doRegionalCaloForMuons = cms.Sequence(( process.doRegionalMuonsEcal+ process.doLocalHcal)* process.towerMakerForMuons* process.caloTowersForMuons)


 process.heepSelection = cms.Sequence( process.heepRobustElectronId* process.heepEcalIsolationSequence* process.heepHcalIsolation* process.heepElectronTkNumIsolation* process.heepElectronTkIsolation* process.heepSelector)


 process.heepSelectionProducersOnly = cms.Sequence( process.heepRobustElectronId* process.heepEcalIsolationSequence* process.heepHcalIsolation* process.heepElectronTkNumIsolation* process.heepElectronTkIsolation)


 process.doLocalCalo = cms.Sequence( process.doLocalEcal+ process.doLocalHcal)


 process.muonlocalreco = cms.Sequence( process.dtlocalreco+ process.csclocalreco+ process.rpcRecHits)


 process.doRegionalEgammaEcal = cms.Sequence( process.ecalPreshowerDigis+ process.ecalRegionalEgammaRecoSequence)


 process.RawToDigi = cms.Sequence( process.SiPixelRawToDigi+ process.SiStripRawToClusters+ process.EcalRawToDigi+ process.EcalESRawToDigi+ process.HcalRawToDigi+ process.CSCRawToDigi+ process.RPCRawToDigi)


 process.doLocalMuon = cms.Sequence( process.doLocalDT+ process.doLocalCSC+ process.doLocalRPC)


 process.doCalo = cms.Sequence( process.doLocalCalo* process.caloTowersRec)


 process.doLocalTracker = cms.Sequence( process.doLocalPixel+ process.doLocalStrip)


 process.HLTEndpath = cms.EndPath( process.hltMakeSummaryObjects+ process.l1tTrigReport+ process.hltTrigReport)


 process.Chi2MeasurementEstimator = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('Chi2'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(30.0)
)


 process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder")


 process.myTTRHBuilderWithoutAngle4MixedPairs = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    StripCPE = cms.string('Fake'),
    ComponentName = cms.string('TTRHBuilderWithoutAngle4MixedPairs'),
    PixelCPE = cms.string('PixelCPEfromTrackAngle'),
    Matcher = cms.string('StandardMatcher')
)


 process.pixellayertriplets = cms.ESProducer("PixelLayerTripletsESProducer",
    ComponentName = cms.string('PixelLayerTriplets'),
    layerList = cms.vstring('BPix1+BPix2+BPix3', 
        'BPix1+BPix2+FPix1_pos', 
        'BPix1+BPix2+FPix1_neg', 
        'BPix1+FPix1_pos+FPix2_pos', 
        'BPix1+FPix1_neg+FPix2_neg'),
    BPix = cms.PSet(
        useErrorsFromParam = cms.untracked.bool(True),
        hitErrorRPhi = cms.double(0.0027),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits'),
        hitErrorRZ = cms.double(0.006)
    ),
    FPix = cms.PSet(
        useErrorsFromParam = cms.untracked.bool(True),
        hitErrorRPhi = cms.double(0.0051),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits'),
        hitErrorRZ = cms.double(0.0036)
    )
)


 process.CaloTowerHardcodeGeometryEP = cms.ESProducer("CaloTowerHardcodeGeometryEP")


 process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducer",
    useParametrizedTrackerField = cms.bool(False),
    findVolumeTolerance = cms.double(0.0),
    timerOn = cms.untracked.bool(False),
    debugBuilder = cms.untracked.bool(False),
    cacheLastVolume = cms.untracked.bool(True)
)


 process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEfromTrackAngleESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle')
)


 process.KFFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    Fitter = cms.string('KFFitter'),
    ComponentName = cms.string('KFFittingSmoother'),
    Smoother = cms.string('KFSmoother'),
    EstimateCut = cms.double(-1.0),
    MinNumberOfHits = cms.int32(5)
)


 process.myTTRHBuilderWithoutAngle4PixelPairs = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    StripCPE = cms.string('Fake'),
    ComponentName = cms.string('TTRHBuilderWithoutAngle4PixelPairs'),
    PixelCPE = cms.string('PixelCPEfromTrackAngle'),
    Matcher = cms.string('StandardMatcher')
)


 process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


 process.roads = cms.ESProducer("RoadMapMakerESProducer",
    GeometryStructure = cms.string('FullDetector'),
    ComponentName = cms.string(''),
    RingsLabel = cms.string(''),
    WriteOutRoadMapToAsciiFile = cms.untracked.bool(False),
    SeedingType = cms.string('FourRingSeeds'),
    RoadMapAsciiFile = cms.untracked.string('roads.dat')
)


 process.KFUpdatorESProducer = cms.ESProducer("KFUpdatorESProducer",
    ComponentName = cms.string('KFUpdator')
)


 process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0)
)


 process.pixellayerpairs = cms.ESProducer("PixelLayerPairsESProducer",
    ComponentName = cms.string('PixelLayerPairs'),
    layerList = cms.vstring('BPix1+BPix2', 
        'BPix1+BPix3', 
        'BPix2+BPix3', 
        'BPix1+FPix1_pos', 
        'BPix1+FPix1_neg', 
        'BPix1+FPix2_pos', 
        'BPix1+FPix2_neg', 
        'BPix2+FPix1_pos', 
        'BPix2+FPix1_neg', 
        'BPix2+FPix2_pos', 
        'BPix2+FPix2_neg', 
        'FPix1_pos+FPix2_pos', 
        'FPix1_neg+FPix2_neg'),
    BPix = cms.PSet(
        useErrorsFromParam = cms.untracked.bool(True),
        hitErrorRPhi = cms.double(0.0027),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelPairs'),
        HitProducer = cms.string('siPixelRecHits'),
        hitErrorRZ = cms.double(0.006)
    ),
    FPix = cms.PSet(
        useErrorsFromParam = cms.untracked.bool(True),
        hitErrorRPhi = cms.double(0.0051),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelPairs'),
        HitProducer = cms.string('siPixelRecHits'),
        hitErrorRZ = cms.double(0.0036)
    )
)


 process.EcalTrigPrimESProducer = cms.ESProducer("EcalTrigPrimESProducer",
    DatabaseFileEE = cms.untracked.string('TPG_EE.txt'),
    DatabaseFileEB = cms.untracked.string('TPG.txt')
)


 process.MeasurementTracker = cms.ESProducer("MeasurementTrackerESProducer",
    StripCPE = cms.string('SimpleStripCPE'),
    UseStripNoiseDB = cms.bool(False),
    ComponentName = cms.string(''),
    stripClusterProducer = cms.string('siStripClusters'),
    Regional = cms.bool(True),
    UseStripCablingDB = cms.bool(False),
    pixelClusterProducer = cms.string('siPixelClusters'),
    HitMatcher = cms.string('StandardMatcher'),
    PixelCPE = cms.string('PixelCPEfromTrackAngle')
)


 process.rings = cms.ESProducer("RingMakerESProducer",
    DumpDetIds = cms.untracked.bool(False),
    ComponentName = cms.string(''),
    RingAsciiFileName = cms.untracked.string('rings.dat'),
    DetIdsDumpFileName = cms.untracked.string('tracker_detids.dat'),
    WriteOutRingsToAsciiFile = cms.untracked.bool(False),
    Configuration = cms.untracked.string('FULL')
)


 process.EcalPreshowerGeometryEP = cms.ESProducer("EcalPreshowerGeometryEP")


 process.CkfTrajectoryBuilder = cms.ESProducer("CkfTrajectoryBuilderESProducer",
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    maxConsecLostHit = cms.int32(1),
    maxLostHit = cms.int32(1),
    minimumNumberOfHits = cms.int32(5),
    maxCand = cms.int32(5),
    ComponentName = cms.string('CkfTrajectoryBuilder'),
    intermediateCleaning = cms.bool(True),
    MeasurementTrackerName = cms.string(''),
    maxNumberOfHits = cms.int32(-1),
    TTRHBuilder = cms.string('WithTrackAngle'),
    updator = cms.string('KFUpdator'),
    alwaysUseInvalidHits = cms.bool(True),
    ptCut = cms.double(0.9),
    estimator = cms.string('Chi2'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    lostHitPenalty = cms.double(30.0)
)


 process.myTTRHBuilderWithoutAngle4PixelTriplets = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    StripCPE = cms.string('Fake'),
    ComponentName = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
    PixelCPE = cms.string('PixelCPEfromTrackAngle'),
    Matcher = cms.string('StandardMatcher')
)


 process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True),
    useGangedStripsInME1a = cms.bool(True),
    useCentreTIOffsets = cms.bool(False),
    applyAlignment = cms.untracked.bool(False)
)


 process.EcalEndcapGeometryEP = cms.ESProducer("EcalEndcapGeometryEP")


 process.SiStripRegionConnectivity = cms.ESProducer("SiStripRegionConnectivity",
    TIBLayers = cms.vuint32(1, 2, 3, 4),
    TECWheels = cms.vuint32(1, 2, 3, 4, 5, 
        6, 7, 8, 9),
    TOBLayers = cms.vuint32(1, 2, 3, 4, 5, 
        6),
    EtaDivisions = cms.untracked.uint32(50),
    TIDWheels = cms.vuint32(1, 2, 3),
    PhiDivisions = cms.untracked.uint32(50),
    EtaMax = cms.untracked.double(2.4)
)


 process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


 process.MaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    MaxDPhi = cms.double(1.6),
    ComponentName = cms.string('PropagatorWithMaterial'),
    Mass = cms.double(0.105),
    PropagationDirection = cms.string('alongMomentum'),
    useRungeKutta = cms.bool(False)
)


 process.GroupedCkfTrajectoryBuilder = cms.ESProducer("GroupedCkfTrajectoryBuilderESProducer",
    bestHitOnly = cms.bool(True),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    maxConsecLostHit = cms.int32(1),
    minimumNumberOfHits = cms.int32(5),
    maxCand = cms.int32(5),
    maxLostHit = cms.int32(1),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    MeasurementTrackerName = cms.string(''),
    lockHits = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    foundHitBonus = cms.double(5.0),
    updator = cms.string('KFUpdator'),
    alwaysUseInvalidHits = cms.bool(True),
    ptCut = cms.double(0.9),
    requireSeedHitsInRebuild = cms.bool(True),
    ComponentName = cms.string('GroupedCkfTrajectoryBuilder'),
    estimator = cms.string('Chi2'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    minNrOfHitsForRebuild = cms.int32(5),
    maxNumberOfHits = cms.int32(-1)
)


 process.HcalHardcodeGeometryEP = cms.ESProducer("HcalHardcodeGeometryEP")


 process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True)
)


 process.ttrhbwr = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    ComponentName = cms.string('WithTrackAngle'),
    PixelCPE = cms.string('PixelCPEfromTrackAngle'),
    Matcher = cms.string('StandardMatcher')
)


 process.KFTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('KFFitter'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('PropagatorWithMaterial'),
    Updator = cms.string('KFUpdator')
)


 process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EEMap.txt')
)


 process.HcalTopologyIdealEP = cms.ESProducer("HcalTopologyIdealEP")


 process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)


 process.mixedlayerpairs = cms.ESProducer("MixedLayerPairsESProducer",
    ComponentName = cms.string('MixedLayerPairs'),
    layerList = cms.vstring('BPix1+BPix2', 
        'BPix1+BPix3', 
        'BPix2+BPix3', 
        'BPix1+FPix1_pos', 
        'BPix1+FPix1_neg', 
        'BPix1+FPix2_pos', 
        'BPix1+FPix2_neg', 
        'BPix2+FPix1_pos', 
        'BPix2+FPix1_neg', 
        'BPix2+FPix2_pos', 
        'BPix2+FPix2_neg', 
        'FPix1_pos+FPix2_pos', 
        'FPix1_neg+FPix2_neg', 
        'FPix2_pos+TEC1_pos', 
        'FPix2_pos+TEC2_pos', 
        'TEC2_pos+TEC3_pos', 
        'FPix2_neg+TEC1_neg', 
        'FPix2_neg+TEC2_neg', 
        'TEC2_neg+TEC3_neg'),
    TEC = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        useRingSlector = cms.untracked.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'),
        minRing = cms.int32(1),
        maxRing = cms.int32(2)
    ),
    BPix = cms.PSet(
        useErrorsFromParam = cms.untracked.bool(True),
        hitErrorRPhi = cms.double(0.0027),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4MixedPairs'),
        HitProducer = cms.string('siPixelRecHits'),
        hitErrorRZ = cms.double(0.006)
    ),
    FPix = cms.PSet(
        useErrorsFromParam = cms.untracked.bool(True),
        hitErrorRPhi = cms.double(0.0051),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4MixedPairs'),
        HitProducer = cms.string('siPixelRecHits'),
        hitErrorRZ = cms.double(0.0036)
    )
)


 process.StripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('SimpleStripCPE')
)


 process.TrackerDigiGeometryESModule = cms.ESProducer("TrackerDigiGeometryESModule",
    applyAlignment = cms.untracked.bool(False)
)


 process.EcalBarrelGeometryEP = cms.ESProducer("EcalBarrelGeometryEP")


 process.TrackerGeometricDetESModule = cms.ESProducer("TrackerGeometricDetESModule")


 process.l1CaloGeometry = cms.ESProducer("L1CaloGeometryProd")


 process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    applyAlignment = cms.untracked.bool(False)
)


 process.PixelCPEParmErrorESProducer = cms.ESProducer("PixelCPEParmErrorESProducer",
    ComponentName = cms.string('PixelCPEfromTrackAngle'),
    UseNewParametrization = cms.bool(True),
    TanLorentzAnglePerTesla = cms.double(0.106),
    UseSigma = cms.bool(True),
    Alpha2Order = cms.bool(True),
    PixelErrorParametrization = cms.string('NOTcmsim')
)


 process.SiStripConnectivity = cms.ESProducer("SiStripConnectivity")


 process.OppositeMaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    MaxDPhi = cms.double(1.6),
    ComponentName = cms.string('PropagatorWithMaterialOpposite'),
    Mass = cms.double(0.105),
    PropagationDirection = cms.string('oppositeToMomentum'),
    useRungeKutta = cms.bool(False)
)


 process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


 process.MuonNumberingInitialization = cms.ESProducer("MuonNumberingInitialization")


 process.myTTRHBuilderWithoutAngle = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    StripCPE = cms.string('Fake'),
    ComponentName = cms.string('PixelTTRHBuilderWithoutAngle'),
    PixelCPE = cms.string('PixelCPEfromTrackAngle'),
    Matcher = cms.string('StandardMatcher')
)


 process.KFTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('KFSmoother'),
    errorRescaling = cms.double(100.0),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('PropagatorWithMaterial'),
    Updator = cms.string('KFUpdator')
)


 process.magfield = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/cms.xml', 
        'Geometry/CMSCommonData/data/cmsMagneticField.xml', 
        'Geometry/CMSCommonData/data/MagneticFieldVolumes.xml'),
    rootNodeName = cms.string('MagneticFieldVolumes:MAGF')
)


 process.eegeom = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalMappingRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)


 process.XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/materials.xml', 
        'Geometry/CMSCommonData/data/rotations.xml', 
        'Geometry/CMSCommonData/data/totem_rotations.xml', 
        'Geometry/CMSCommonData/data/cms.xml', 
        'Geometry/CMSCommonData/data/cmsMother.xml', 
        'Geometry/CMSCommonData/data/cmsTracker.xml', 
        'Geometry/CMSCommonData/data/caloBase.xml', 
        'Geometry/CMSCommonData/data/cmsCalo.xml', 
        'Geometry/CMSCommonData/data/muonBase.xml', 
        'Geometry/CMSCommonData/data/cmsMuon.xml', 
        'Geometry/CMSCommonData/data/mgnt.xml', 
        'Geometry/CMSCommonData/data/beampipe.xml', 
        'Geometry/CMSCommonData/data/cmsBeam.xml', 
        'Geometry/CMSCommonData/data/muonMB.xml', 
        'Geometry/CMSCommonData/data/muonMagnet.xml', 
        'Geometry/TrackerCommonData/data/pixfwdMaterials.xml', 
        'Geometry/TrackerCommonData/data/pixfwdCommon.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq1x2.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq1x5.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq2x3.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq2x4.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq2x5.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPanelBase.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPanel.xml', 
        'Geometry/TrackerCommonData/data/pixfwdBlade.xml', 
        'Geometry/TrackerCommonData/data/pixfwdNipple.xml', 
        'Geometry/TrackerCommonData/data/pixfwdDisk.xml', 
        'Geometry/TrackerCommonData/data/pixfwdCylinder.xml', 
        'Geometry/TrackerCommonData/data/pixfwd.xml', 
        'Geometry/TrackerCommonData/data/pixbarmaterial.xml', 
        'Geometry/TrackerCommonData/data/pixbarladder.xml', 
        'Geometry/TrackerCommonData/data/pixbarladderfull.xml', 
        'Geometry/TrackerCommonData/data/pixbarladderhalf.xml', 
        'Geometry/TrackerCommonData/data/pixbarlayer.xml', 
        'Geometry/TrackerCommonData/data/pixbarlayer0.xml', 
        'Geometry/TrackerCommonData/data/pixbarlayer1.xml', 
        'Geometry/TrackerCommonData/data/pixbarlayer2.xml', 
        'Geometry/TrackerCommonData/data/pixbar.xml', 
        'Geometry/TrackerCommonData/data/tibtidcommonmaterial.xml', 
        'Geometry/TrackerCommonData/data/tibmaterial.xml', 
        'Geometry/TrackerCommonData/data/tibmodpar.xml', 
        'Geometry/TrackerCommonData/data/tibmodule0.xml', 
        'Geometry/TrackerCommonData/data/tibmodule0a.xml', 
        'Geometry/TrackerCommonData/data/tibmodule0b.xml', 
        'Geometry/TrackerCommonData/data/tibmodule2.xml', 
        'Geometry/TrackerCommonData/data/tibstringpar.xml', 
        'Geometry/TrackerCommonData/data/tibstringds.xml', 
        'Geometry/TrackerCommonData/data/tibstring0c.xml', 
        'Geometry/TrackerCommonData/data/tibstring0ll.xml', 
        'Geometry/TrackerCommonData/data/tibstring0lr.xml', 
        'Geometry/TrackerCommonData/data/tibstring0ul.xml', 
        'Geometry/TrackerCommonData/data/tibstring0ur.xml', 
        'Geometry/TrackerCommonData/data/tibstring0.xml', 
        'Geometry/TrackerCommonData/data/tibstring1c.xml', 
        'Geometry/TrackerCommonData/data/tibstring1ll.xml', 
        'Geometry/TrackerCommonData/data/tibstring1lr.xml', 
        'Geometry/TrackerCommonData/data/tibstring1ul.xml', 
        'Geometry/TrackerCommonData/data/tibstring1ur.xml', 
        'Geometry/TrackerCommonData/data/tibstring1.xml', 
        'Geometry/TrackerCommonData/data/tibstringss.xml', 
        'Geometry/TrackerCommonData/data/tibstring2c.xml', 
        'Geometry/TrackerCommonData/data/tibstring2ll.xml', 
        'Geometry/TrackerCommonData/data/tibstring2lr.xml', 
        'Geometry/TrackerCommonData/data/tibstring2ul.xml', 
        'Geometry/TrackerCommonData/data/tibstring2ur.xml', 
        'Geometry/TrackerCommonData/data/tibstring2.xml', 
        'Geometry/TrackerCommonData/data/tibstring3c.xml', 
        'Geometry/TrackerCommonData/data/tibstring3ll.xml', 
        'Geometry/TrackerCommonData/data/tibstring3lr.xml', 
        'Geometry/TrackerCommonData/data/tibstring3ul.xml', 
        'Geometry/TrackerCommonData/data/tibstring3ur.xml', 
        'Geometry/TrackerCommonData/data/tibstring3.xml', 
        'Geometry/TrackerCommonData/data/tiblayerpar.xml', 
        'Geometry/TrackerCommonData/data/tiblayer0.xml', 
        'Geometry/TrackerCommonData/data/tiblayer1.xml', 
        'Geometry/TrackerCommonData/data/tiblayer2.xml', 
        'Geometry/TrackerCommonData/data/tiblayer3.xml', 
        'Geometry/TrackerCommonData/data/tib.xml', 
        'Geometry/TrackerCommonData/data/tidmaterial.xml', 
        'Geometry/TrackerCommonData/data/tidmodpar.xml', 
        'Geometry/TrackerCommonData/data/tidmodule0.xml', 
        'Geometry/TrackerCommonData/data/tidmodule0r.xml', 
        'Geometry/TrackerCommonData/data/tidmodule0l.xml', 
        'Geometry/TrackerCommonData/data/tidmodule1.xml', 
        'Geometry/TrackerCommonData/data/tidmodule1r.xml', 
        'Geometry/TrackerCommonData/data/tidmodule1l.xml', 
        'Geometry/TrackerCommonData/data/tidmodule2.xml', 
        'Geometry/TrackerCommonData/data/tidringpar.xml', 
        'Geometry/TrackerCommonData/data/tidring0.xml', 
        'Geometry/TrackerCommonData/data/tidring0f.xml', 
        'Geometry/TrackerCommonData/data/tidring0b.xml', 
        'Geometry/TrackerCommonData/data/tidring1.xml', 
        'Geometry/TrackerCommonData/data/tidring1f.xml', 
        'Geometry/TrackerCommonData/data/tidring1b.xml', 
        'Geometry/TrackerCommonData/data/tidring2.xml', 
        'Geometry/TrackerCommonData/data/tid.xml', 
        'Geometry/TrackerCommonData/data/tidf.xml', 
        'Geometry/TrackerCommonData/data/tidb.xml', 
        'Geometry/TrackerCommonData/data/tibtidservices.xml', 
        'Geometry/TrackerCommonData/data/tibtidservicesf.xml', 
        'Geometry/TrackerCommonData/data/tibtidservicesb.xml', 
        'Geometry/TrackerCommonData/data/tobmaterial.xml', 
        'Geometry/TrackerCommonData/data/tobmodpar.xml', 
        'Geometry/TrackerCommonData/data/tobmodule0.xml', 
        'Geometry/TrackerCommonData/data/tobmodule2.xml', 
        'Geometry/TrackerCommonData/data/tobmodule4.xml', 
        'Geometry/TrackerCommonData/data/tobrodpar.xml', 
        'Geometry/TrackerCommonData/data/tobrod0c.xml', 
        'Geometry/TrackerCommonData/data/tobrod0l.xml', 
        'Geometry/TrackerCommonData/data/tobrod0h.xml', 
        'Geometry/TrackerCommonData/data/tobrod0.xml', 
        'Geometry/TrackerCommonData/data/tobrod1l.xml', 
        'Geometry/TrackerCommonData/data/tobrod1h.xml', 
        'Geometry/TrackerCommonData/data/tobrod1.xml', 
        'Geometry/TrackerCommonData/data/tobrod2c.xml', 
        'Geometry/TrackerCommonData/data/tobrod2l.xml', 
        'Geometry/TrackerCommonData/data/tobrod2h.xml', 
        'Geometry/TrackerCommonData/data/tobrod2.xml', 
        'Geometry/TrackerCommonData/data/tobrod3l.xml', 
        'Geometry/TrackerCommonData/data/tobrod3h.xml', 
        'Geometry/TrackerCommonData/data/tobrod3.xml', 
        'Geometry/TrackerCommonData/data/tobrod4c.xml', 
        'Geometry/TrackerCommonData/data/tobrod4l.xml', 
        'Geometry/TrackerCommonData/data/tobrod4h.xml', 
        'Geometry/TrackerCommonData/data/tobrod4.xml', 
        'Geometry/TrackerCommonData/data/tobrod5l.xml', 
        'Geometry/TrackerCommonData/data/tobrod5h.xml', 
        'Geometry/TrackerCommonData/data/tobrod5.xml', 
        'Geometry/TrackerCommonData/data/tob.xml', 
        'Geometry/TrackerCommonData/data/tecmaterial.xml', 
        'Geometry/TrackerCommonData/data/tecmodpar.xml', 
        'Geometry/TrackerCommonData/data/tecmodule0.xml', 
        'Geometry/TrackerCommonData/data/tecmodule0r.xml', 
        'Geometry/TrackerCommonData/data/tecmodule0s.xml', 
        'Geometry/TrackerCommonData/data/tecmodule1.xml', 
        'Geometry/TrackerCommonData/data/tecmodule1r.xml', 
        'Geometry/TrackerCommonData/data/tecmodule1s.xml', 
        'Geometry/TrackerCommonData/data/tecmodule2.xml', 
        'Geometry/TrackerCommonData/data/tecmodule3.xml', 
        'Geometry/TrackerCommonData/data/tecmodule4.xml', 
        'Geometry/TrackerCommonData/data/tecmodule4r.xml', 
        'Geometry/TrackerCommonData/data/tecmodule4s.xml', 
        'Geometry/TrackerCommonData/data/tecmodule5.xml', 
        'Geometry/TrackerCommonData/data/tecmodule6.xml', 
        'Geometry/TrackerCommonData/data/tecpetpar.xml', 
        'Geometry/TrackerCommonData/data/tecring0.xml', 
        'Geometry/TrackerCommonData/data/tecring1.xml', 
        'Geometry/TrackerCommonData/data/tecring2.xml', 
        'Geometry/TrackerCommonData/data/tecring3.xml', 
        'Geometry/TrackerCommonData/data/tecring4.xml', 
        'Geometry/TrackerCommonData/data/tecring5.xml', 
        'Geometry/TrackerCommonData/data/tecring6.xml', 
        'Geometry/TrackerCommonData/data/tecring0f.xml', 
        'Geometry/TrackerCommonData/data/tecring1f.xml', 
        'Geometry/TrackerCommonData/data/tecring2f.xml', 
        'Geometry/TrackerCommonData/data/tecring3f.xml', 
        'Geometry/TrackerCommonData/data/tecring4f.xml', 
        'Geometry/TrackerCommonData/data/tecring5f.xml', 
        'Geometry/TrackerCommonData/data/tecring6f.xml', 
        'Geometry/TrackerCommonData/data/tecring0b.xml', 
        'Geometry/TrackerCommonData/data/tecring1b.xml', 
        'Geometry/TrackerCommonData/data/tecring2b.xml', 
        'Geometry/TrackerCommonData/data/tecring3b.xml', 
        'Geometry/TrackerCommonData/data/tecring4b.xml', 
        'Geometry/TrackerCommonData/data/tecring5b.xml', 
        'Geometry/TrackerCommonData/data/tecring6b.xml', 
        'Geometry/TrackerCommonData/data/tecpetalf.xml', 
        'Geometry/TrackerCommonData/data/tecpetalb.xml', 
        'Geometry/TrackerCommonData/data/tecpetal0.xml', 
        'Geometry/TrackerCommonData/data/tecpetal0f.xml', 
        'Geometry/TrackerCommonData/data/tecpetal0b.xml', 
        'Geometry/TrackerCommonData/data/tecpetal3.xml', 
        'Geometry/TrackerCommonData/data/tecpetal3f.xml', 
        'Geometry/TrackerCommonData/data/tecpetal3b.xml', 
        'Geometry/TrackerCommonData/data/tecpetal6f.xml', 
        'Geometry/TrackerCommonData/data/tecpetal6b.xml', 
        'Geometry/TrackerCommonData/data/tecpetal8f.xml', 
        'Geometry/TrackerCommonData/data/tecpetal8b.xml', 
        'Geometry/TrackerCommonData/data/tecwheel.xml', 
        'Geometry/TrackerCommonData/data/tecwheela.xml', 
        'Geometry/TrackerCommonData/data/tecwheelb.xml', 
        'Geometry/TrackerCommonData/data/tecwheelc.xml', 
        'Geometry/TrackerCommonData/data/tecwheeld.xml', 
        'Geometry/TrackerCommonData/data/tec.xml', 
        'Geometry/TrackerCommonData/data/trackermaterial.xml', 
        'Geometry/TrackerCommonData/data/tracker.xml', 
        'Geometry/TrackerCommonData/data/trackerpixbar.xml', 
        'Geometry/TrackerCommonData/data/trackerpixfwd.xml', 
        'Geometry/TrackerCommonData/data/trackertibtidservices.xml', 
        'Geometry/TrackerCommonData/data/trackertib.xml', 
        'Geometry/TrackerCommonData/data/trackertid.xml', 
        'Geometry/TrackerCommonData/data/trackertob.xml', 
        'Geometry/TrackerCommonData/data/trackertec.xml', 
        'Geometry/TrackerCommonData/data/trackerother.xml', 
        'Geometry/EcalCommonData/data/eregalgo.xml', 
        'Geometry/EcalCommonData/data/ebalgo.xml', 
        'Geometry/EcalCommonData/data/ebcon.xml', 
        'Geometry/EcalCommonData/data/ebrot.xml', 
        'Geometry/EcalCommonData/data/eealgo.xml', 
        'Geometry/EcalCommonData/data/esalgo.xml', 
        'Geometry/HcalCommonData/data/hcalrotations.xml', 
        'Geometry/HcalCommonData/data/hcalalgo.xml', 
        'Geometry/HcalCommonData/data/hcalbarrelalgo.xml', 
        'Geometry/HcalCommonData/data/hcalendcapalgo.xml', 
        'Geometry/HcalCommonData/data/hcalouteralgo.xml', 
        'Geometry/HcalCommonData/data/hcalforwardalgo.xml', 
        'Geometry/HcalCommonData/data/hcalforwardmaterial.xml', 
        'Geometry/MuonCommonData/data/mbCommon.xml', 
        'Geometry/MuonCommonData/data/mb1.xml', 
        'Geometry/MuonCommonData/data/mb2.xml', 
        'Geometry/MuonCommonData/data/mb3.xml', 
        'Geometry/MuonCommonData/data/mb4.xml', 
        'Geometry/MuonCommonData/data/muonYoke.xml', 
        'Geometry/MuonCommonData/data/mf.xml', 
        'Geometry/MuonCommonData/data/muonNumbering.xml', 
        'Geometry/TrackerCommonData/data/trackerStructureTopology.xml', 
        'Geometry/TrackerSimData/data/trackersens.xml', 
        'Geometry/TrackerRecoData/data/trackerRecoMaterial.xml', 
        'Geometry/EcalSimData/data/ecalsens.xml', 
        'Geometry/HcalCommonData/data/hcalsens.xml', 
        'Geometry/HcalSimData/data/CaloUtil.xml', 
        'Geometry/MuonSimData/data/muonSens.xml', 
        'Geometry/DTGeometryBuilder/data/dtSpecsFilter.xml', 
        'Geometry/CSCGeometryBuilder/data/cscSpecsFilter.xml', 
        'Geometry/CSCGeometryBuilder/data/cscSpecs.xml', 
        'Geometry/RPCGeometryBuilder/data/RPCSpecs.xml', 
        'Geometry/CMSCommonData/data/cmsMagneticField.xml', 
        'Geometry/CMSCommonData/data/MagneticFieldVolumes.xml', 
        'Geometry/HcalSimData/data/HcalProdCuts.xml', 
        'Geometry/EcalSimData/data/EcalProdCuts.xml', 
        'Geometry/TrackerSimData/data/trackerProdCuts.xml', 
        'Geometry/TrackerSimData/data/trackerProdCutsBEAM.xml', 
        'Geometry/MuonSimData/data/muonProdCuts.xml', 
        'Geometry/CMSCommonData/data/FieldParameters.xml'),
    rootNodeName = cms.string('cms:OCMS')
)


 process.tpparams = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPParametersRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)


 process.DTCabling = cms.ESSource("PoolDBESSource",
    siteLocalConfig = cms.untracked.bool(True),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('DTReadOutMappingRcd'),
        tag = cms.string('DTROMap')
    )),
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(30),
        loadBlobStreamer = cms.untracked.bool(False),
        messageLevel = cms.untracked.int32(0),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableConnectionSharing = cms.untracked.bool(False),
        connectionRetrialTimeOut = cms.untracked.int32(180),
        connectionTimeOut = cms.untracked.int32(600),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False)
    ),
    catalog = cms.untracked.string('relationalcatalog_frontier://cms_conditions_data/CMS_COND_CSA07_FRONTIER'),
    timetype = cms.string('runnumber'),
    connect = cms.string('frontier://FrontierCSA07/CMS_COND_CSA07_DT')
)


 process.l1CaloGeomRecordSource = cms.ESSource("EmptyESSource",
    recordName = cms.string('L1CaloGeometryRecord'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)


 process.SiStripDBCabling = cms.ESSource("PoolDBESSource",
    siteLocalConfig = cms.untracked.bool(True),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('SiStripFedCablingRcd'),
        tag = cms.string('CSA07_SiStripFedCabling')
    )),
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(30),
        loadBlobStreamer = cms.untracked.bool(False),
        messageLevel = cms.untracked.int32(0),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableConnectionSharing = cms.untracked.bool(False),
        connectionRetrialTimeOut = cms.untracked.int32(180),
        connectionTimeOut = cms.untracked.int32(600),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False)
    ),
    catalog = cms.untracked.string('relationalcatalog_frontier://cms_conditions_data/CMS_COND_CSA07_FRONTIER'),
    timetype = cms.string('runnumber'),
    connect = cms.string('frontier://FrontierCSA07/CMS_COND_CSA07_STRIP')
)


 process.siPixelCabling = cms.ESSource("PoolDBESSource",
    siteLocalConfig = cms.untracked.bool(True),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelFedCablingMapRcd'),
        tag = cms.string('SiPixelFedCablingMap_v1')
    )),
    catalog = cms.untracked.string('relationalcatalog_frontier://cms_conditions_data/CMS_COND_CSA07_FRONTIER'),
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(30),
        loadBlobStreamer = cms.untracked.bool(False),
        messageLevel = cms.untracked.int32(0),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableConnectionSharing = cms.untracked.bool(False),
        connectionRetrialTimeOut = cms.untracked.int32(180),
        connectionTimeOut = cms.untracked.int32(600),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False)
    ),
    timetype = cms.string('runnumber'),
    connect = cms.string('frontier://FrontierCSA07/CMS_COND_CSA07_PIXEL'),
    authenticationMethod = cms.untracked.uint32(0)
)


 process.RPCCabling = cms.ESSource("PoolDBESSource",
    siteLocalConfig = cms.untracked.bool(True),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('RPCReadOutMappingRcd'),
        tag = cms.string('RPCReadOutMapping_v1')
    )),
    catalog = cms.untracked.string('relationalcatalog_frontier://cms_conditions_data/CMS_COND_CSA07_FRONTIER'),
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(30),
        loadBlobStreamer = cms.untracked.bool(False),
        messageLevel = cms.untracked.int32(0),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableConnectionSharing = cms.untracked.bool(False),
        connectionRetrialTimeOut = cms.untracked.int32(180),
        connectionTimeOut = cms.untracked.int32(600),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False)
    ),
    timetype = cms.string('runnumber'),
    connect = cms.string('frontier://FrontierCSA07/CMS_COND_CSA07_RPC'),
    authenticationMethod = cms.untracked.uint32(0)
)


 process.DTCombinatorialPatternReco2DAlgo_ParamDrift_CSA07 = cms.PSet(
    Reco2DAlgoConfig = cms.PSet(
        segmCleanerMode = cms.int32(1),
        AlphaMaxPhi = cms.double(1.0),
        MaxAllowedHits = cms.uint32(50),
        nSharedHitsMax = cms.int32(2),
        AlphaMaxTheta = cms.double(0.1),
        debug = cms.untracked.bool(False),
        nUnSharedHitsMin = cms.int32(2),
        recAlgoConfig = cms.PSet(
            tTrigMode = cms.string('DTTTrigSyncFromDB'),
            minTime = cms.double(-3.0),
            interpolate = cms.bool(True),
            debug = cms.untracked.bool(False),
            tTrigModeConfig = cms.PSet(
                vPropWire = cms.double(24.4),
                doTOFCorrection = cms.bool(True),
                tofCorrType = cms.int32(1),
                kFactor = cms.double(-2.0),
                wirePropCorrType = cms.int32(1),
                doWirePropCorrection = cms.bool(True),
                doT0Correction = cms.bool(True),
                debug = cms.untracked.bool(False)
            ),
            maxTime = cms.double(415.0)
        ),
        recAlgo = cms.string('DTParametrizedDriftAlgo')
    ),
    Reco2DAlgoName = cms.string('DTCombinatorialPatternReco')
)

 process.CondDBSetup = cms.PSet(
    catalog = cms.untracked.string(''),
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(30),
        loadBlobStreamer = cms.untracked.bool(False),
        messageLevel = cms.untracked.int32(0),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableConnectionSharing = cms.untracked.bool(False),
        connectionRetrialTimeOut = cms.untracked.int32(180),
        connectionTimeOut = cms.untracked.int32(600),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False)
    )
)

 process.TC_ME1234 = cms.PSet(
    dPhiFineMax = cms.double(0.02),
    verboseInfo = cms.untracked.bool(True),
    SegmentSorting = cms.int32(1),
    chi2Max = cms.double(6000.0),
    dPhiMax = cms.double(0.003),
    chi2ndfProbMin = cms.double(0.0001),
    minLayersApart = cms.int32(2),
    dRPhiFineMax = cms.double(6.0),
    dRPhiMax = cms.double(1.2)
)

 process.CSCSegAlgoDF = cms.PSet(
    algo_name = cms.string('CSCSegAlgoDF'),
    algo_psets = cms.VPSet(cms.PSet(
        minHitsPerSegment = cms.untracked.int32(3),
        dPhiFineMax = cms.untracked.double(0.025),
        dXclusBoxMax = cms.untracked.double(4.0),
        tanThetaMax = cms.untracked.double(1.2),
        BrutePruning = cms.untracked.bool(False),
        preClustering = cms.untracked.bool(False),
        tanPhiMax = cms.untracked.double(0.5),
        nSigmaFromSegment = cms.untracked.double(5.0),
        minLayersApart = cms.untracked.int32(2),
        dRPhiFineMax = cms.untracked.double(8.0),
        CSCSegmentDebug = cms.untracked.bool(False),
        Pruning = cms.untracked.bool(False),
        dYclusBoxMax = cms.untracked.double(8.0)
    ), 
        cms.PSet(
            minHitsPerSegment = cms.untracked.int32(3),
            dPhiFineMax = cms.untracked.double(0.025),
            dXclusBoxMax = cms.untracked.double(4.0),
            tanThetaMax = cms.untracked.double(2.0),
            BrutePruning = cms.untracked.bool(False),
            preClustering = cms.untracked.bool(False),
            tanPhiMax = cms.untracked.double(0.8),
            nSigmaFromSegment = cms.untracked.double(5.0),
            minLayersApart = cms.untracked.int32(2),
            dRPhiFineMax = cms.untracked.double(12.0),
            CSCSegmentDebug = cms.untracked.bool(False),
            Pruning = cms.untracked.bool(False),
            dYclusBoxMax = cms.untracked.double(8.0)
        ), 
        cms.PSet(
            minHitsPerSegment = cms.untracked.int32(3),
            dPhiFineMax = cms.untracked.double(0.025),
            dXclusBoxMax = cms.untracked.double(4.0),
            tanThetaMax = cms.untracked.double(1.2),
            BrutePruning = cms.untracked.bool(False),
            preClustering = cms.untracked.bool(False),
            tanPhiMax = cms.untracked.double(0.5),
            nSigmaFromSegment = cms.untracked.double(5.0),
            minLayersApart = cms.untracked.int32(2),
            dRPhiFineMax = cms.untracked.double(8.0),
            CSCSegmentDebug = cms.untracked.bool(False),
            Pruning = cms.untracked.bool(False),
            dYclusBoxMax = cms.untracked.double(8.0)
        )),
    parameters_per_chamber_type = cms.vint32(3, 1, 2, 2, 1, 
        2, 1, 2, 1)
)

 process.SiStripClusterization = cms.PSet(
    ClusterizerAlgorithm = cms.untracked.string('ThreeThreshold'),
    ChannelThreshold = cms.untracked.double(2.0),
    SeedThreshold = cms.untracked.double(3.0),
    MaxHolesInCluster = cms.untracked.uint32(0),
    ClusterThreshold = cms.untracked.double(5.0)
)

 process.HEEPCutBased_ID = cms.PSet(
    looseEleIDCuts = cms.PSet(
        deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
            0.03, 0.092, 0.092),
        hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
            0.12, 0.15, 0.0),
        sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
            0.0265, 0.0265, 0.0),
        deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
            0.0068, 0.01, 0.0),
        eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
            0.85, 0.0, 0.0)
    ),
    robustEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005),
        endcap = cms.vdouble(0.08, 0.0275, 0.092, 0.007)
    ),
    tightEleIDCuts = cms.PSet(
        eSeedOverPinMax = cms.vdouble(99999.0, 99999.0, 99999.0, 99999.0, 99999.0, 
            99999.0, 99999.0, 99999.0),
        eSeedOverPinMin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
            0.83, 0.0, 0.0),
        deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
            0.035, 0.065, 0.092),
        hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
            0.037, 0.05, 0.0),
        sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
            0.0252, 0.026, 0.0),
        deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
            0.0055, 0.0075, 0.0)
    ),
    electronQuality = cms.string('robust')
)

 process.ST_ME1234 = cms.PSet(
    minHitsPerSegment = cms.untracked.int32(3),
    dXclusBoxMax = cms.untracked.double(4.0),
    BrutePruning = cms.untracked.bool(False),
    preClustering = cms.untracked.bool(True),
    maxRecHitsInCluster = cms.untracked.int32(20),
    onlyBestSegment = cms.untracked.bool(False),
    CSCDebug = cms.untracked.bool(False),
    Pruning = cms.untracked.bool(False),
    dYclusBoxMax = cms.untracked.double(8.0)
)

 process.DTCombinatorialPatternReco4DAlgo_ParamDrift_CSA07 = cms.PSet(
    Reco4DAlgoName = cms.string('DTCombinatorialPatternReco4D'),
    Reco4DAlgoConfig = cms.PSet(
        segmCleanerMode = cms.int32(1),
        nSharedHitsMax = cms.int32(2),
        debug = cms.untracked.bool(False),
        nUnSharedHitsMin = cms.int32(2),
        AllDTRecHits = cms.bool(True),
        Reco2DAlgoConfig = cms.PSet(
            segmCleanerMode = cms.int32(1),
            AlphaMaxPhi = cms.double(1.0),
            MaxAllowedHits = cms.uint32(50),
            nSharedHitsMax = cms.int32(2),
            AlphaMaxTheta = cms.double(0.1),
            debug = cms.untracked.bool(False),
            nUnSharedHitsMin = cms.int32(2),
            recAlgoConfig = cms.PSet(
                tTrigMode = cms.string('DTTTrigSyncFromDB'),
                minTime = cms.double(-3.0),
                interpolate = cms.bool(True),
                debug = cms.untracked.bool(False),
                tTrigModeConfig = cms.PSet(
                    vPropWire = cms.double(24.4),
                    doTOFCorrection = cms.bool(True),
                    tofCorrType = cms.int32(1),
                    kFactor = cms.double(-2.0),
                    wirePropCorrType = cms.int32(1),
                    doWirePropCorrection = cms.bool(True),
                    doT0Correction = cms.bool(True),
                    debug = cms.untracked.bool(False)
                ),
                maxTime = cms.double(415.0)
            ),
            recAlgo = cms.string('DTParametrizedDriftAlgo')
        ),
        Reco2DAlgoName = cms.string('DTCombinatorialPatternReco'),
        recAlgoConfig = cms.PSet(
            tTrigMode = cms.string('DTTTrigSyncFromDB'),
            minTime = cms.double(-3.0),
            interpolate = cms.bool(True),
            debug = cms.untracked.bool(False),
            tTrigModeConfig = cms.PSet(
                vPropWire = cms.double(24.4),
                doTOFCorrection = cms.bool(True),
                tofCorrType = cms.int32(1),
                kFactor = cms.double(-2.0),
                wirePropCorrType = cms.int32(1),
                doWirePropCorrection = cms.bool(True),
                doT0Correction = cms.bool(True),
                debug = cms.untracked.bool(False)
            ),
            maxTime = cms.double(415.0)
        ),
        recAlgo = cms.string('DTParametrizedDriftAlgo')
    )
)

 process.SiStripFull = cms.PSet(
    TOBLayers = cms.vuint32(1, 2, 3, 4, 5, 
        6),
    TIBLayers = cms.vuint32(1, 2, 3, 4),
    TECWheels = cms.vuint32(1, 2, 3, 4, 5, 
        6, 7, 8, 9),
    TIDWheels = cms.vuint32(1, 2, 3)
)

 process.SK_ME1234 = cms.PSet(
    dPhiFineMax = cms.double(0.025),
    verboseInfo = cms.untracked.bool(True),
    chi2Max = cms.double(99999.0),
    dPhiMax = cms.double(0.003),
    wideSeg = cms.double(3.0),
    minLayersApart = cms.int32(2),
    dRPhiFineMax = cms.double(8.0),
    dRPhiMax = cms.double(8.0)
)

 process.ST_ME1A = cms.PSet(
    minHitsPerSegment = cms.untracked.int32(3),
    dXclusBoxMax = cms.untracked.double(4.0),
    BrutePruning = cms.untracked.bool(False),
    preClustering = cms.untracked.bool(True),
    maxRecHitsInCluster = cms.untracked.int32(20),
    onlyBestSegment = cms.untracked.bool(False),
    CSCDebug = cms.untracked.bool(False),
    Pruning = cms.untracked.bool(False),
    dYclusBoxMax = cms.untracked.double(8.0)
)

 process.SiPixelGainCalibrationServiceParameters = cms.PSet(
    PedestalValue = cms.double(28.2),
    MaxGain = cms.double(3.0),
    UseCalibDataFromDB = cms.bool(False),
    MaxPed = cms.double(35.0),
    MinPed = cms.double(25.0),
    MinGain = cms.double(2.5),
    GainValue = cms.double(2.8)
)

 process.SK_ME1A = cms.PSet(
    dPhiFineMax = cms.double(0.025),
    verboseInfo = cms.untracked.bool(True),
    chi2Max = cms.double(99999.0),
    dPhiMax = cms.double(0.025),
    wideSeg = cms.double(3.0),
    minLayersApart = cms.int32(2),
    dRPhiFineMax = cms.double(3.0),
    dRPhiMax = cms.double(8.0)
)

 process.CSCSegAlgoST = cms.PSet(
    algo_name = cms.string('CSCSegAlgoST'),
    algo_psets = cms.VPSet(cms.PSet(
        minHitsPerSegment = cms.untracked.int32(3),
        dXclusBoxMax = cms.untracked.double(4.0),
        BrutePruning = cms.untracked.bool(False),
        preClustering = cms.untracked.bool(True),
        maxRecHitsInCluster = cms.untracked.int32(20),
        onlyBestSegment = cms.untracked.bool(False),
        CSCDebug = cms.untracked.bool(False),
        Pruning = cms.untracked.bool(False),
        dYclusBoxMax = cms.untracked.double(8.0)
    ), 
        cms.PSet(
            minHitsPerSegment = cms.untracked.int32(3),
            dXclusBoxMax = cms.untracked.double(4.0),
            BrutePruning = cms.untracked.bool(False),
            preClustering = cms.untracked.bool(True),
            maxRecHitsInCluster = cms.untracked.int32(20),
            onlyBestSegment = cms.untracked.bool(False),
            CSCDebug = cms.untracked.bool(False),
            Pruning = cms.untracked.bool(False),
            dYclusBoxMax = cms.untracked.double(8.0)
        )),
    parameters_per_chamber_type = cms.vint32(2, 1, 1, 1, 1, 
        1, 1, 1, 1)
)

 process.PTDR_ID = cms.PSet(
    useEoverPOut = cms.vint32(1, 1, 1),
    useDeltaPhiIn = cms.vint32(1, 1, 1),
    tightEleIDCuts = cms.PSet(
        invEMinusInvP = cms.vdouble(0.02, 0.02, 0.02, 0.02, 0.02, 
            0.02, 0.02, 0.02),
        EoverPInMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0),
        EoverPOutMin = cms.vdouble(0.6, 0.75, 0.75, 0.75, 0.5, 
            0.8, 0.5, 0.8),
        sigmaEtaEtaMin = cms.vdouble(0.005, 0.005, 0.005, 0.005, 0.008, 
            0.008, 0.008, 0.008),
        EoverPOutMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
            999.0, 999.0, 999.0),
        EoverPInMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
            999.0, 999.0, 999.0),
        deltaPhiOut = cms.vdouble(0.02, 999.0, 0.02, 999.0, 0.02, 
            999.0, 0.02, 999.0),
        sigmaEtaEtaMax = cms.vdouble(0.011, 0.011, 0.011, 0.011, 0.03, 
            0.03, 0.03, 0.022),
        deltaPhiIn = cms.vdouble(0.02, 0.03, 0.02, 0.04, 0.04, 
            0.04, 0.04, 0.05),
        HoverE = cms.vdouble(0.05, 0.05, 0.05, 0.05, 0.07, 
            0.07, 0.07, 0.07),
        sigmaPhiPhiMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0),
        bremFraction = cms.vdouble(0.0, 0.1, 0.1, 0.1, 0.0, 
            0.2, 0.2, 0.2),
        deltaEtaIn = cms.vdouble(0.004, 0.004, 0.004, 0.005, 0.005, 
            0.005, 0.005, 0.005),
        E9overE25 = cms.vdouble(0.8, 0.65, 0.75, 0.65, 0.8, 
            0.7, 0.7, 0.65),
        sigmaPhiPhiMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
            999.0, 999.0, 999.0)
    ),
    useBremFraction = cms.vint32(0, 0, 0),
    useSigmaPhiPhi = cms.vint32(0, 1, 0),
    useEoverPIn = cms.vint32(0, 1, 0),
    useHoverE = cms.vint32(1, 1, 1),
    looseEleIDCuts = cms.PSet(
        invEMinusInvP = cms.vdouble(0.02, 0.02, 0.02, 0.02, 0.02, 
            0.02, 0.02, 0.02),
        EoverPInMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0),
        EoverPOutMin = cms.vdouble(0.6, 1.7, 0.9, 0.5, 0.6, 
            1.7, 0.9, 0.5),
        sigmaEtaEtaMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0),
        EoverPOutMax = cms.vdouble(2.5, 999.0, 2.2, 999.0, 2.5, 
            999.0, 2.2, 999.0),
        EoverPInMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
            999.0, 999.0, 999.0),
        deltaPhiOut = cms.vdouble(0.011, 999.0, 999.0, 999.0, 0.02, 
            999.0, 999.0, 999.0),
        sigmaEtaEtaMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
            999.0, 999.0, 999.0),
        deltaPhiIn = cms.vdouble(0.06, 0.06, 0.06, 0.08, 0.06, 
            0.06, 0.06, 0.09),
        HoverE = cms.vdouble(0.09, 0.06, 0.07, 0.12, 0.09, 
            0.06, 0.07, 0.12),
        sigmaPhiPhiMin = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0),
        bremFraction = cms.vdouble(0.0, 0.1, 0.1, 0.1, 0.0, 
            0.2, 0.2, 0.2),
        deltaEtaIn = cms.vdouble(0.008, 0.008, 0.008, 0.009, 0.008, 
            0.008, 0.008, 0.009),
        E9overE25 = cms.vdouble(0.7, 0.7, 0.7, 0.5, 0.8, 
            0.8, 0.8, 0.5),
        sigmaPhiPhiMax = cms.vdouble(999.0, 999.0, 999.0, 999.0, 999.0, 
            999.0, 999.0, 999.0)
    ),
    useDeltaEtaIn = cms.vint32(1, 1, 1),
    useSigmaEtaEta = cms.vint32(0, 1, 1),
    useInvEMinusInvP = cms.vint32(0, 0, 0),
    useDeltaPhiOut = cms.vint32(0, 1, 1),
    useE9overE25 = cms.vint32(1, 1, 1),
    electronQuality = cms.string('loose'),
    mediumEleIDCuts = cms.PSet(
        invEMinusInvP = cms.vdouble(0.02, 0.02, 0.02, 0.02, 0.02, 
            0.02, 0.02, 0.02),
        EoverPInMin = cms.vdouble(0.9, 0.9, 0.9, 0.6, 0.9, 
            0.9, 0.9, 0.7),
        EoverPOutMin = cms.vdouble(0.6, 1.8, 1.0, 0.75, 0.6, 
            1.5, 1.0, 0.8),
        sigmaEtaEtaMin = cms.vdouble(0.005, 0.005, 0.005, 0.005, 0.008, 
            0.008, 0.008, 0.0),
        EoverPOutMax = cms.vdouble(2.5, 999.0, 999.0, 999.0, 2.0, 
            999.0, 999.0, 999.0),
        EoverPInMax = cms.vdouble(1.3, 1.2, 1.3, 999.0, 999.0, 
            999.0, 999.0, 999.0),
        deltaPhiOut = cms.vdouble(0.011, 999.0, 999.0, 999.0, 0.02, 
            999.0, 999.0, 999.0),
        sigmaEtaEtaMax = cms.vdouble(0.011, 0.011, 0.011, 0.011, 0.022, 
            0.022, 0.022, 0.3),
        deltaPhiIn = cms.vdouble(0.04, 0.07, 0.04, 0.08, 0.06, 
            0.07, 0.06, 0.07),
        HoverE = cms.vdouble(0.06, 0.05, 0.06, 0.14, 0.1, 
            0.1, 0.1, 0.12),
        sigmaPhiPhiMin = cms.vdouble(0.005, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0),
        bremFraction = cms.vdouble(0.0, 0.1, 0.1, 0.1, 0.0, 
            0.2, 0.2, 0.2),
        deltaEtaIn = cms.vdouble(0.004, 0.006, 0.005, 0.007, 0.007, 
            0.008, 0.007, 0.008),
        E9overE25 = cms.vdouble(0.7, 0.75, 0.8, 0.0, 0.85, 
            0.75, 0.8, 0.0),
        sigmaPhiPhiMax = cms.vdouble(0.015, 999.0, 999.0, 999.0, 0.02, 
            999.0, 999.0, 999.0)
    )
)

 process.TC_ME1A = cms.PSet(
    dPhiFineMax = cms.double(0.013),
    verboseInfo = cms.untracked.bool(True),
    SegmentSorting = cms.int32(1),
    chi2Max = cms.double(6000.0),
    dPhiMax = cms.double(0.00198),
    chi2ndfProbMin = cms.double(0.0001),
    minLayersApart = cms.int32(2),
    dRPhiFineMax = cms.double(3.0),
    dRPhiMax = cms.double(0.6)
)

 process.DF_ME1A = cms.PSet(
    preClustering = cms.untracked.bool(False),
    minHitsPerSegment = cms.untracked.int32(3),
    dPhiFineMax = cms.untracked.double(0.025),
    dXclusBoxMax = cms.untracked.double(4.0),
    BrutePruning = cms.untracked.bool(False),
    tanThetaMax = cms.untracked.double(1.2),
    tanPhiMax = cms.untracked.double(0.5),
    nSigmaFromSegment = cms.untracked.double(5.0),
    minLayersApart = cms.untracked.int32(2),
    dRPhiFineMax = cms.untracked.double(8.0),
    CSCSegmentDebug = cms.untracked.bool(False),
    Pruning = cms.untracked.bool(False),
    dYclusBoxMax = cms.untracked.double(8.0)
)

 process.CSCSegAlgoSK = cms.PSet(
    algo_name = cms.string('CSCSegAlgoSK'),
    algo_psets = cms.VPSet(cms.PSet(
        dPhiFineMax = cms.double(0.025),
        verboseInfo = cms.untracked.bool(True),
        chi2Max = cms.double(99999.0),
        dPhiMax = cms.double(0.003),
        wideSeg = cms.double(3.0),
        minLayersApart = cms.int32(2),
        dRPhiFineMax = cms.double(8.0),
        dRPhiMax = cms.double(8.0)
    ), 
        cms.PSet(
            dPhiFineMax = cms.double(0.025),
            verboseInfo = cms.untracked.bool(True),
            chi2Max = cms.double(99999.0),
            dPhiMax = cms.double(0.025),
            wideSeg = cms.double(3.0),
            minLayersApart = cms.int32(2),
            dRPhiFineMax = cms.double(3.0),
            dRPhiMax = cms.double(8.0)
        )),
    parameters_per_chamber_type = cms.vint32(2, 1, 1, 1, 1, 
        1, 1, 1, 1)
)

 process.DF_ME1234_1 = cms.PSet(
    preClustering = cms.untracked.bool(False),
    minHitsPerSegment = cms.untracked.int32(3),
    dPhiFineMax = cms.untracked.double(0.025),
    dXclusBoxMax = cms.untracked.double(4.0),
    BrutePruning = cms.untracked.bool(False),
    tanThetaMax = cms.untracked.double(1.2),
    tanPhiMax = cms.untracked.double(0.5),
    nSigmaFromSegment = cms.untracked.double(5.0),
    minLayersApart = cms.untracked.int32(2),
    dRPhiFineMax = cms.untracked.double(8.0),
    CSCSegmentDebug = cms.untracked.bool(False),
    Pruning = cms.untracked.bool(False),
    dYclusBoxMax = cms.untracked.double(8.0)
)

 process.DF_ME1234_2 = cms.PSet(
    preClustering = cms.untracked.bool(False),
    minHitsPerSegment = cms.untracked.int32(3),
    dPhiFineMax = cms.untracked.double(0.025),
    dXclusBoxMax = cms.untracked.double(4.0),
    BrutePruning = cms.untracked.bool(False),
    tanThetaMax = cms.untracked.double(2.0),
    tanPhiMax = cms.untracked.double(0.8),
    nSigmaFromSegment = cms.untracked.double(5.0),
    minLayersApart = cms.untracked.int32(2),
    dRPhiFineMax = cms.untracked.double(12.0),
    CSCSegmentDebug = cms.untracked.bool(False),
    Pruning = cms.untracked.bool(False),
    dYclusBoxMax = cms.untracked.double(8.0)
)

 process.DTParametrizedDriftAlgo_CSA07 = cms.PSet(
    recAlgoConfig = cms.PSet(
        tTrigMode = cms.string('DTTTrigSyncFromDB'),
        minTime = cms.double(-3.0),
        interpolate = cms.bool(True),
        debug = cms.untracked.bool(False),
        tTrigModeConfig = cms.PSet(
            vPropWire = cms.double(24.4),
            doTOFCorrection = cms.bool(True),
            tofCorrType = cms.int32(1),
            kFactor = cms.double(-2.0),
            wirePropCorrType = cms.int32(1),
            doWirePropCorrection = cms.bool(True),
            doT0Correction = cms.bool(True),
            debug = cms.untracked.bool(False)
        ),
        maxTime = cms.double(415.0)
    ),
    recAlgo = cms.string('DTParametrizedDriftAlgo')
)

 process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

 process.CSCSegAlgoTC = cms.PSet(
    algo_name = cms.string('CSCSegAlgoTC'),
    algo_psets = cms.VPSet(cms.PSet(
        dPhiFineMax = cms.double(0.02),
        verboseInfo = cms.untracked.bool(True),
        SegmentSorting = cms.int32(1),
        chi2Max = cms.double(6000.0),
        dPhiMax = cms.double(0.003),
        chi2ndfProbMin = cms.double(0.0001),
        minLayersApart = cms.int32(2),
        dRPhiFineMax = cms.double(6.0),
        dRPhiMax = cms.double(1.2)
    ), 
        cms.PSet(
            dPhiFineMax = cms.double(0.013),
            verboseInfo = cms.untracked.bool(True),
            SegmentSorting = cms.int32(1),
            chi2Max = cms.double(6000.0),
            dPhiMax = cms.double(0.00198),
            chi2ndfProbMin = cms.double(0.0001),
            minLayersApart = cms.int32(2),
            dRPhiFineMax = cms.double(3.0),
            dRPhiMax = cms.double(0.6)
        )),
    parameters_per_chamber_type = cms.vint32(2, 1, 1, 1, 1, 
        1, 1, 1, 1)
)

 process.HEEPSelectionPath = cms.Path(process.heepSelection)

