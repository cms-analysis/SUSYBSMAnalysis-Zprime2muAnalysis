import FWCore.ParameterSet.Config as cms

from MassReachDataSets_cff import dataSets

def attachMassReach(process):
    process.Zprime2muMassReach = cms.EDAnalyzer(
        'Zprime2muMassReach',
        process.Zprime2muAnalysisCommon,
        dataSets,
        dataSet        = cms.string('Zssm1000'),
        
        psFile         = cms.untracked.string('mass_fits.ps'),

        DYEvents       = cms.bool(False),  # set to use Drell-Yan only     
        BackgroundFit  = cms.bool(False),  # set to fit background slope   
        GenuineEvents  = cms.bool(True),   # set to use genuine sim. events
        GMR            = cms.bool(False),  # set to use GMR instead of "bes
        FixedMass      = cms.bool(True),   # set to fix signal mean mass   
        FixedFWHM      = cms.bool(True),   # set to fix signal FWHM        
        SmoothedSample = cms.bool(False),  # currently only for DY to 5 TeV
        RandomSeed     = cms.bool(False),  # set to randomly modify gen. se
        ExpPlots       = cms.bool(False),  # set to plot individual expts  
        BinnedFit      = cms.bool(False),  # set to perform binned fit     
        Nexp           = cms.uint32(1000), # total number of experiments   
        intLumi        = cms.double(0.1),  # integrated luminosity, in fb-1

        # in the unbinned fit, fit either the generated masses, or the
        # reconstructed masses
        fitGenMass     = cms.bool(False),
        fitRecMass     = cms.bool(True),

        nGenEvents     = cms.vuint32(1000,1000,1000)
        )

    process.analysis = cms.Path(process.Zprime2muMassReach)
