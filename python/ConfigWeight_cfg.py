import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

                                      #DY50
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0033A97B-8707-E511-9D3B-008CFA1980B8.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/003826EE-8807-E511-9628-008CFA06470C.root',

                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/F416988B-832A-E511-8F6D-002590D9D9DA.root'
    )
)

process.maxEvents.input = -1


process.demo = cms.EDAnalyzer("WeightAnalyzer",

	#GenEventInfo = cms.InputTag("GenEventInfoProduct_generator","","SIM")
         genEventInfo = cms.InputTag('generator'),
        
	
)


process.p = cms.Path(process.demo)
