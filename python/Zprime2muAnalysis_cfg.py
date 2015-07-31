import FWCore.ParameterSet.Config as cms

process = cms.Process('Zprime2muAnalysis')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring(#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/788396C0-9D6F-E411-97DF-002590494E34.root',
	#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3023C1A6-D56F-E411-B210-002590AC4C08.root',
	#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/90575C8B-D56F-E411-A39B-0025904B1420.root',
	#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/FCEF82CA-9D6F-E411-A2FC-002481E0D5CE.root',
	#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/44B8EEE7-9A6F-E411-8858-002590DB0640.root',
	#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/F03F87D0-9A6F-E411-8BC4-00266CFFA5E0.root'))
#file:///afs/cern.ch/user/k/klarson/private/miniWork/CMSSW_7_2_0/src/SUSYBSMAnalysis/Zprime2muAnalysis/MyOutputFile.root', 
'/store/relval/CMSSW_7_2_0/RelValZMM_13/MINIAODSIM/PU25ns_PHYS14_25_V1_Phys14-v2/00000/EA3D8F7C-A059-E411-8858-0025905A6136.root',
'/store/relval/CMSSW_7_2_0/RelValZMM_13/MINIAODSIM/PU25ns_PHYS14_25_V1_Phys14-v2/00000/0A3BFB7A-A059-E411-A9C1-0025905A48D6.root'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service('TFileService', fileName=cms.string('zp2mu_histos.root'))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag =  'MCRUN2_74_V7::All'

flag = cms.string('AOD')
if flag == 'miniAOD':
	process.load('SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff_miniAOD')
if flag == 'AOD':
	process.load('SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff')
