import FWCore.ParameterSet.Config as cms

## Define some utilities to declare bins easily
def _nBins(n,min,max): return cms.vdouble(*[min + (max-min)/n*i for i in range(n+1)])
def _evenBins(min,max,delta): 
    ret = cms.vdouble(min)
    x = min
    while x < max - 1e-4: # need a small hint otherwise for some numbers it will overstep due to numerical resolution
        x += delta
        ret.append(x)
    return ret 

def makeInclusiveDiMuonPlots(rebinFactor=1):
    return cms.PSet(
        # ---- Kinematics ----
        massBins = _nBins( 100, 0., 100.),
        ptBins = _evenBins( 0, 100, 2 * rebinFactor),
        pBins  = _evenBins( 0, 500, 2  * rebinFactor),
        etaBins = _evenBins( -2.6, 2.6, 0.2 * rebinFactor),
        phiBins = _evenBins(-3.2,  3.2, 0.2 * rebinFactor),
        chargeBins = cms.vdouble(-2,0,2),
        # ---- Vertex ----
        dxyFineBins = _evenBins(-0.2, 0.2, 0.005), #  50um
        dzFineBins  = _evenBins(-0.5, 0.5, 0.010), # 100um
        dxyCoarseBins = _evenBins( -4,  4, 0.1), # 1mm
        dzCoarseBins  = _evenBins(-10, 10, 0.1), # 1mm
        # ---- Tracks ----
        pixelHitsBins       = _nBins(8,0,8),
        pixelLayersBins     = _nBins(5,0,5),
        trackerHitsBins     = _nBins(33,0,33),
        trackerLostHitsBins = _nBins(10,0,10),
        muonHitsBins        = _nBins(50,0,50),
        muonStationHitsBins = _nBins(20,0,20),
        muonBadHitsBins     = _nBins(20,0,20),
        globalHitsBins      = _nBins(80,0,80),
        trackerChi2nBins = _evenBins(0, 20, 0.2 * rebinFactor),
        muonChi2nBins    = _evenBins(0, 20, 0.2 * rebinFactor),
        globalChi2nBins  = _evenBins(0, 20, 0.2 * rebinFactor),
        trackerChi2RelBins  = _evenBins(0, 4, 0.2 * rebinFactor),
        muonChi2RelBins  = _evenBins(0, 20, 0.2 * rebinFactor),
        # ---- Isolation ----
        isolationBins = _evenBins(0,  5, .25  * rebinFactor),
        relIsoBins    = _evenBins(0, .5, .025 * rebinFactor),
        # ---- Muon ID ----
        muonStationsBins    = _nBins(5,0,5), 
        segmentMatchesBins = _nBins(12,0,12),
        segmentCompatBins  = _evenBins(0, 1 + 0.1*rebinFactor, 0.1 * rebinFactor), # need one bin for ">= 1.0"
        caloCompatBins     = _evenBins(0, 1 + 0.1*rebinFactor, 0.1 * rebinFactor), # need one bin for ">= 1.0"
        # ---- Production Vertex ----
        zBins = _nBins(100,-500.,500.),
        rBins = _nBins(100,0.,500.),
        rzXBins = cms.uint32(1000),
        rzXRange = cms.vdouble(-500.,500.),
        rzYBins = cms.uint32(500),
        rzYRange = cms.vdouble(0.,500.),
        # ---- ----
        boolBins = _nBins(2,-0.5,1.5),
        deltaPtBins = _evenBins( -50., 50., 2 * rebinFactor),
        deltaPtnBins = _evenBins(-5.,5.,0.1 * rebinFactor),
        muonHitCountsratioBins = _evenBins(0.,1.2,0.1 * rebinFactor),
        muonHitCountsrpcratioBins = _evenBins(0.,1.2,0.1 * rebinFactor),
        ratioBins = _evenBins(0.,1.2,0.1 * rebinFactor),
    )

inclusiveDiMuonPlots = cms.EDAnalyzer("Zprime2muAnalysisPlots",
    makeInclusiveDiMuonPlots(),
    dilepton_src = cms.InputTag('zToMuMuGG'),
    selection = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    #normalization   = cms.InputTag("countCollisionEvents"), ## read normalization from output of cms.EDProducer("EventCountProducer") 
)
