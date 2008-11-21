#!/usr/bin/env python

import os, sys

crabcfg = '''
[CRAB]
jobtype = cmssw
scheduler = glitecoll

[CMSSW]
datasetpath = %(dataset)s
pset = tmp_dilskim_cmssw.cfg
total_number_of_events = -1
number_of_jobs = %(njobs)i
#output_file = out.root, heeptree.root
output_file = zp2mu_histos.root

[USER]
ui_working_dir = crab_%(shortname)s
return_data = 1
#copy_data = 1
#storage_element = srm-cms.cern.ch
#storage_path = /srm/managerv2?SFN=/castor/cern.ch/user/t/tucker/crab/%(shortname)s

use_central_bossDB = 0
use_boss_rt = 0

[EDG]
CE_black_list = fnal.gov, cern.ch, in2p3.fr, rl.ac.uk, pic.es, fzk.de, infn.it, asgc.tw
CE_white_list = %(whitelist)s

rb = CERN 
proxy_server = myproxy.cern.ch 
virtual_organization = cms
retry_count = 0
lcg_catalog_type = lfc
lfc_host = lfc-cms-test.cern.ch
lfc_home = /grid/cms
'''

cmsswweights = '''
  # Make the weights for 100/pb.

  module csa07EventWeightProducer = CSA07EventWeightProducer {
    InputTag src = source
    untracked bool talkToMe = false
    double overallLumi = 100.
    # K-factor for ttbar (= 1. for ALPGEN LO cross sections)
    # from MCFM NLO vs LO study, the K-factor is found to be 1.85
    double ttKfactor = 1 # 1.85
  }
  
  path weights = { csa07EventWeightProducer }

  # Ignore Z+jets.
  module csaids = CSA07ProcessIdFilter {
    vint32 csa07Ids = { 0,1,2,3,4,5,6,7,8,9,10,22,23,24,25,26 }
    double overallLumi = 100.
    string csa07EventWeightProducerLabel = "csa07EventWeightProducer"
  }
'''

cmsswcfg_intermediate = '''
process Skim = {
  include "FWCore/MessageLogger/data/MessageLogger.cfi"
  replace MessageLogger.cerr.threshold = "INFO"
#  replace MessageLogger.categories += "PATLayer0Summary"
  replace MessageLogger.categories += "FwkReport"
#  replace MessageLogger.cerr = {
#      untracked PSet PATLayer0Summary = { untracked int32 limit = -1 }
#  }
  replace MessageLogger.cerr.FwkReport.reportEvery = 100
  untracked PSet options = { untracked bool wantSummary = true }
  untracked PSet maxEvents = { untracked int32 input = -1 }

  source = PoolSource {
    untracked vstring fileNames = {
      "file:/scratchdisk3/tucker/CSA07-CSA07AllEvents-Tier0-A3-Chowder/002198E8-CC9E-DC11-BDD6-00304855D4B8.root"
    }
  }

  %(weights)s
  
  module myHLTfilter = hltHighLevel from "HLTrigger/HLTfilters/data/hltHighLevel.cfi"
  replace myHLTfilter.HLTPaths = {
    # us + topDiLeptonMuonX
    "HLT1MuonNonIso",
    "HLT2MuonNonIso",
    "HLTXElectronMuonRelaxed",
    # topDiLepton2Electron
    "HLT1ElectronRelaxed",
    "HLT2ElectronRelaxed",
    # HEEP
    "HLT1EMHighEt",
    "HLT1EMVeryHighEt"
  }

  module CSA07MuonFilter = hltHighLevel from "HLTrigger/HLTfilters/data/hltHighLevel.cfi"
  replace CSA07MuonFilter.HLTPaths = {
    "HLT1MuonIso",
    "HLT1MuonNonIso",
    "HLT2MuonNonIso",
    "HLT2MuonJPsi",
    "HLT2MuonUpsilon",
    "HLT2MuonZ",
    "HLTNMuonNonIso",
    "HLT2MuonSameSign",
    "HLTBJPsiMuMu",
    "HLTXMuonJets",
    "HLTXElectronMuon",
    "HLTXElectronMuonRelaxed",
    "HLTXMuonTau"
  }

  # Kinematic pre-selection: require there be any two leptons with pT
  # > 20 GeV.

  module allLeptons = CandViewMerger {
    VInputTag src = { muons, pixelMatchGsfElectrons }
  }

  module ptLeptons = PtMinCandSelector {
    InputTag src = allLeptons
    double ptMin = 20
  }

  module numberFilter = CandCountFilter {
    InputTag src = ptLeptons
    uint32 minNumber = 2
  }

  sequence kine = { allLeptons, ptLeptons, numberFilter }

  # Perform TeV re-reconstruction.
  include "SUSYBSMAnalysis/Zprime2muAnalysis/data/TeVMuonReReco.cff"

  service = TFileService {
    string fileName = "heeptree.root"
  }

#  include "PhysicsTools/PatAlgos/data/patLayer0.cff"
#  include "PhysicsTools/PatAlgos/data/patLayer1.cff"
#
#  # Switch the PAT-default "robust" electron id to "tight".
#  replace electronId.doPtdrId   = false
#  replace electronId.doCutBased = true
#  replace CutBased_ID.electronQuality = "tight"

  # Use HEEPSelector while we have the full RECO.
  include "DLEvans/HEEPSelector/data/heepSelection_1_2.cfi"

  # Go ahead and make simParticles before dropping SIM.
  module simParticleCandidates = GenCandsFromSimTracksProducer {
    InputTag src = g4SimHits
  }

#  sequence stuff = { patLayer0, patLayer1, heepSelection, simParticleCandidates }
  sequence stuff = { heepSelection, simParticleCandidates }

  # The filter that will determine the event output.
  path myfilter = { myHLTfilter & %(nocsa07muon)s %(weights2)s kine & TeVMuonReReco & stuff }

  block myEventContent = {
    untracked vstring outputCommands = {
      "drop *",

      "keep edmHepMCProduct_*_*_*",
      "keep recoCandidatesOwned_genParticleCandidates_*_*",
      "keep recoCandidatesOwned_simParticleCandidates_*_*",


      "keep *_ctfWithMaterialTracks_*_Rec1",
      "keep *_pixelMatchGsfFitBarrel_*_*",
      "keep *_pixelMatchGsfFitEndcap_*_*",

      "keep *_hybridSuperClusters_*_*",
      "keep *_islandBasicClusters_*_*",
      "keep *_islandSuperClusters_*_*",
      "keep recoSuperClusters_*_*_*",
      
      "keep *_globalMuons_*_*",
      "keep *_globalMuonsFMS_*_*",
      "keep *_globalMuonsPMR_*_*",
      "keep *_globalMuonsTK_*_*",
      "keep *_standAloneMuons_*_*",

      "keep recoTracksrecoMuIsoDeposituintedmOneToValueedmAssociationMap_*_*_*",
      "keep recoTracksrecoSuperClustersuintedmOneToOneedmAssociationMap_*_*_*",

      "keep recoCaloJets_*_*_*",
      "keep recoCaloMETs_*_*_*",
      "keep recoJetTags_*_*_*",
      
      "keep recoMuons_*_*_*",
      "keep recoPixelMatchGsfElectrons_*_*_*",
      "keep recoPhotons_correctedPhotons_*_*",
      
      "keep double_genEventScale_*_*",
      "keep double_genEventWeight_*_*",
      "keep double_genEventRunInfo_*_*",
      "keep int_genEventProcID_*_*",

      "keep *_*_*_Rec1",
      
      "drop *_*_*_HLT",
      "keep edmTriggerResults_TriggerResults__HLT"
    }
  }

#  # Keep all the PAT Layer-0 and Layer-1 stuff too.
#  include "PhysicsTools/PatAlgos/data/patLayer0_EventContent.cff"
#  replace myEventContent.outputCommands += patLayer0EventContent.outputCommands
#  include "PhysicsTools/PatAlgos/data/patLayer1_EventContent.cff"
#  replace myEventContent.outputCommands += patLayer1EventContent.outputCommands

  module out = PoolOutputModule {
    untracked string fileName = "out.root"
    untracked PSet SelectEvents = {
      vstring SelectEvents = { "myfilter" }
    }
    using myEventContent
  }

  endpath ep = { out }
}
'''

cmsswcfg = '''
process Skim = {
  include "FWCore/MessageLogger/data/MessageLogger.cfi"
  replace MessageLogger.cerr.threshold = "INFO"
  replace MessageLogger.categories += "PATLayer0Summary"
  replace MessageLogger.categories += "FwkReport"
#  replace MessageLogger.cerr = {
#      untracked PSet PATLayer0Summary = { untracked int32 limit = -1 }
#  }
  replace MessageLogger.cerr.FwkReport.reportEvery = 100
  untracked PSet options = { untracked bool wantSummary = true }
  untracked PSet maxEvents = { untracked int32 input = 100 }

  source = PoolSource {
    untracked vstring fileNames = {
      "file:/scratchdisk1/tucker/E07F5B8E-14E1-DC11-82AF-000423D990CC.root"
    }
  }

  module csa07EventWeightProducer = CSA07EventWeightProducer {
    InputTag src = source
    untracked bool talkToMe = false
    double overallLumi = 100.
    # K-factor for ttbar (= 1. for ALPGEN LO cross sections)
    # from MCFM NLO vs LO study, the K-factor is found to be 1.85
    double ttKfactor = 1 # 1.85
  }

  # Ignore Z+jets.
  module csaids = CSA07ProcessIdFilter {
    vint32 csa07Ids = { 0,1,2,3,4,5,6,7,8,9,10,22,23,24,25,26 }
    double overallLumi = 100.
    string csa07EventWeightProducerLabel = "csa07EventWeightProducer"
  }
  
  module myHLTfilter = hltHighLevel from "HLTrigger/HLTfilters/data/hltHighLevel.cfi"
  replace myHLTfilter.HLTPaths = {
    # us + topDiLeptonMuonX
    "HLT1MuonNonIso",
    "HLT2MuonNonIso",
    "HLTXElectronMuonRelaxed",
    # topDiLepton2Electron
    "HLT1ElectronRelaxed",
    "HLT2ElectronRelaxed",
    # HEEP
    "HLT1EMHighEt",
    "HLT1EMVeryHighEt"
  }

  module CSA07MuonFilter = hltHighLevel from "HLTrigger/HLTfilters/data/hltHighLevel.cfi"
  replace CSA07MuonFilter.HLTPaths = {
    "HLT1MuonIso",
    "HLT1MuonNonIso",
    "HLT2MuonNonIso",
    "HLT2MuonJPsi",
    "HLT2MuonUpsilon",
    "HLT2MuonZ",
    "HLTNMuonNonIso",
    "HLT2MuonSameSign",
    "HLTBJPsiMuMu",
    "HLTXMuonJets",
    "HLTXElectronMuon",
    "HLTXElectronMuonRelaxed",
    "HLTXMuonTau"
  }

  # Kinematic pre-selection: require there be any two leptons with pT
  # > 20 GeV.

  module allLeptons = CandViewMerger {
    VInputTag src = { muons, pixelMatchGsfElectrons }
  }

  module ptLeptons = PtMinCandSelector {
    InputTag src = allLeptons
    double ptMin = 20
  }

  module numberFilter = CandCountFilter {
    InputTag src = ptLeptons
    uint32 minNumber = 2
  }

  sequence kine = { allLeptons, ptLeptons, numberFilter }

  # Perform TeV re-reconstruction.
  include "SUSYBSMAnalysis/Zprime2muAnalysis/data/TeVMuonReReco.cff"

  service = TFileService {
    string fileName = "zp2mu_histos.root"
  }

  include "PhysicsTools/PatAlgos/data/patLayer0.cff"
  include "PhysicsTools/PatAlgos/data/patLayer1.cff"

  # Switch the PAT-default "robust" electron id to "tight".
  replace electronId.doPtdrId   = false
  replace electronId.doCutBased = true
  replace CutBased_ID.electronQuality = "tight"

  # Use HEEPSelector while we have the full RECO.
  include "DLEvans/HEEPSelector/data/heepSelection_1_2.cfi"

  # Go ahead and make simParticles before dropping SIM.
  module simParticleCandidates = GenCandsFromSimTracksProducer {
    InputTag src = g4SimHits
  }

  # Make cocktail muons. Need TK, and seed matches (muonsFMS/PMR are made in the TeVReReco above).
  module muCandTK = TrackerOnlyMuonProducer { 
    InputTag src = muons
  }

  module seedMatchTKFS = MuonBySeedMatcher { 
    InputTag seedTracks = standAloneMuons:UpdatedAtVtx
    InputTag src = muCandTK
    InputTag matched = muonsFMS
  }

  module seedMatchTKPR = MuonBySeedMatcher { 
    InputTag seedTracks = standAloneMuons:UpdatedAtVtx
    InputTag src = muCandTK
    InputTag matched = muonsPMR
  }

  module bestMuons = CocktailMuonProducer { 
    InputTag trackerOnlyMuons = muCandTK
    InputTag toFMSMap = seedMatchTKFS
    InputTag toPMRMap = seedMatchTKPR
    bool useTMR = False
  }

#  module printTree = ParticleListDrawer { 
#    InputTag src = genParticleCandidates
#    untracked int32 maxEventsToPrint = -1
#    untracked bool printOnlyHardInteraction = True
#  }
#  path pt = {printTree}

  module ntupledump = EMuBackgroundsNtupleDumper { 
    untracked int32 selfProcId = %(procId)i
  }

  sequence stuff = { patLayer0 & patLayer1 & heepSelection & simParticleCandidates }
  sequence ntpl = { muCandTK, seedMatchTKFS, seedMatchTKPR, bestMuons, ntupledump }

  # The filter that will determine the event output.
  path myfilter = { myHLTfilter & %(nocsa07muon)s %(weights2)s kine & TeVMuonReReco & stuff & ntpl }
}
'''

WemunuJets = [
    '/WmunuJets_Pt_80_120/CMSSW_1_6_7-CSA07-1192839569/RECO',
    '/WmunuJets_Pt_120_170/CMSSW_1_6_7-CSA07-1192839464/RECO',
    '/WmunuJets_Pt_170_230/CMSSW_1_6_7-CSA07-1195594220/RECO',
    '/WmunuJets_Pt_230_300/CMSSW_1_6_7-CSA07-1194806757/RECO',
    '/WmunuJets_Pt_300_380/CMSSW_1_6_7-CSA07-1195594273/RECO',
    '/WenuJets_Pt_80_120/CMSSW_1_6_7-CSA07-1192836181/RECO',
    '/WenuJets_Pt_120_170/CMSSW_1_6_7-CSA07-1194806653/RECO',
    '/WenuJets_Pt_170_230/CMSSW_1_6_7-CSA07-1192837699/RECO',
    '/WenuJets_Pt_230_300/CMSSW_1_6_7-CSA07-1192837907/RECO',
    '/WenuJets_Pt_300_380/CMSSW_1_6_7-CSA07-1192838062/RECO',
    '/WenuJets_Pt_380_470/CMSSW_1_6_7-CSA07-1192838218/RECO',
]

WemunuJetsExtra = [
    '/WenuJets_Pt_0_15/CMSSW_1_6_7-CSA07-1192837492/RECO',
    '/WenuJets_Pt_15_20/CMSSW_1_6_7-CSA07-1192837647/RECO',
    '/WenuJets_Pt_20_30/CMSSW_1_6_7-CSA07-1192837803/RECO',
    '/WenuJets_Pt_30_50/CMSSW_1_6_7-CSA07-1192838114/RECO',
    '/WenuJets_Pt_50_80/CMSSW_1_6_7-CSA07-1192838322/RECO',
    '/WenuJets_Pt_470_600/CMSSW_1_6_7-CSA07-1192838270/RECO',
    '/WenuJets_Pt_600_800/CMSSW_1_6_7-CSA07-1192838374/RECO',
    '/WenuJets_Pt_800_1000/CMSSW_1_6_7-CSA07-1192838425/RECO',
    '/WmunuJets_Pt_0_15/CMSSW_1_6_7-CSA07-1192839412/RECO',
    '/WmunuJets_Pt_15_20/CMSSW_1_6_7-CSA07-1192836446/RECO',
    '/WmunuJets_Pt_20_30/CMSSW_1_6_7-CSA07-1192836552/RECO',
    '/WmunuJets_Pt_30_50/CMSSW_1_6_7-CSA07-1195594167/RECO',
    '/WmunuJets_Pt_50_80/CMSSW_1_6_7-CSA07-1192839516/RECO',
    '/WmunuJets_Pt_380_470/CMSSW_1_6_7-CSA07-1195629943/RECO',
    '/WmunuJets_Pt_470_600/CMSSW_1_6_7-CSA07-1192836657/RECO',
    '/WmunuJets_Pt_600_800/CMSSW_1_6_7-CSA07-1193559937/RECO',
    '/WmunuJets_Pt_800_1000/CMSSW_1_6_7-CSA07-1203689988/RECO',
]

WemunuJetsExtraextra = [
    '/WenuJets_Pt_1000_1400/CMSSW_1_6_7-CSA07-1192837543/RECO',
    '/WenuJets_Pt_1400_1800/CMSSW_1_6_7-CSA07-1192837595/RECO',
    '/WenuJets_Pt_1800_2200/CMSSW_1_6_7-CSA07-1192837750/RECO',
    '/WenuJets_Pt_2200_2600/CMSSW_1_6_7-CSA07-1192837855/RECO',
    '/WenuJets_Pt_2600_3000/CMSSW_1_6_7-CSA07-1192837959/RECO',
    '/WenuJets_Pt_3000_3500/CMSSW_1_6_7-CSA07-1192838011/RECO',
    '/WenuJets_Pt_3500_-1/CMSSW_1_6_7-CSA07-1192838167/RECO',
    '/WmunuJets_Pt_1000_1400/CMSSW_1_6_7-CSA07-1193559885/RECO',
    '/WmunuJets_Pt_1400_1800/CMSSW_1_6_7-CSA07-1192836393/RECO',
    '/WmunuJets_Pt_1800_2200/CMSSW_1_6_7-CSA07-1192836498/RECO',
    '/WmunuJets_Pt_2200_2600/CMSSW_1_6_7-CSA07-1192836604/RECO',
    '/WmunuJets_Pt_2600_3000/CMSSW_1_6_7-CSA07-1192836554/RECO',
    '/WmunuJets_Pt_3000_3500/CMSSW_1_6_7-CSA07-1192836607/RECO',
    '/WmunuJets_Pt_3500_-1/CMSSW_1_6_7-CSA07-1194806808/RECO'
    ]

Soups_muonly = [
    ('/CSA07Muon/CMSSW_1_6_7-CSA07-Tier0-A1-Chowder/RECO', 'chowder_muonly'),
    ('/CSA07Muon/CMSSW_1_6_7-CSA07-Tier0-A1-Gumbo/RECO', 'gumbo_muonly'),
    ('/CSA07Muon/CMSSW_1_6_7-CSA07-Tier0-A2-Stew/RECO', 'stew_muonly')
]

datasets = [
    ('/WW_incl/CMSSW_1_6_7-CSA07-1196178448/RECO', 'ww', 71),
    ('/WZ_incl/CMSSW_1_6_7-CSA07-1195629996/RECO', 'wz', 72),
    ('/ZZ_incl/CMSSW_1_6_7-CSA07-1194964234/RECO', 'zz', 73),
    ('/DrellYan_ll_40/CMSSW_1_6_7-CSA07-1201885455/RECO', 'dy', 74),
    ('/tW_inclusive/CMSSW_1_6_7-CSA07-1195471738/RECO', 'tw', 76),
    ('/DYtautau_M10/CMSSW_1_6_7-CSA07-1193560401/RECO', 'dytt', 77),
    ]

datasets = [
#    ('/CSA07AllEvents/CMSSW_1_6_7-CSA07-Tier0-A3-Chowder/RECO', 'chowder', 0)
    ('/CSA07Muon/CMSSW_1_6_7-CSA07-Tier0-A2-Stew/RECO',     'mustew', 0),
#    ('/CSA07Electron/CMSSW_1_6_7-CSA07-Tier0-A2-Stew/RECO', 'elstew', 0),
#    ('/CSA07Electron/CMSSW_1_6_7-CSA07-Tier0-A1-Gumbo/RECO', 'elgumbo', 0),
#    ('/CSA07Muon/CMSSW_1_6_7-CSA07-Tier0-A1-Gumbo/RECO', 'mugumbo', 0)
    ]
    
crabFn = 'tmp_dilskim.cfg'
cmsswFn = 'tmp_dilskim_cmssw.cfg'

writeCfgsOnly = 'test' in sys.argv

for dataset in datasets:
    if type(dataset) == type(()):
        dataset, shortname, procId = dataset
    else:
        shortname = ''.join(dataset.split('/')[1].split('_')[:-1])
        raise 'set the pid for wemunujets!'
    print shortname

    # make and write out the tmp crab cfg
    if 'Chowder' in dataset:
        njobs = 2000
        whitelist = 'ufl'
#    elif 'Gumbo' in dataset:
#        njobs = 10
#        whitelist = 'unl'
#    elif 'Stew' in dataset:
#        njobs = 100
#        whitelist = 'unl'
    elif 'Gumbo' in dataset and 'CSA07Muon' in dataset:
        njobs = 50
        whitelist = 'unl'
    elif 'Gumbo' in dataset and 'CSA07Electron' in dataset:
        njobs = 250
        whitelist = 'ac.be'
    elif 'Stew' in dataset and 'CSA07Muon' in dataset:
        njobs = 250
        whitelist = 'unl'
    elif 'Stew' in dataset and 'CSA07Electron' in dataset:
        njobs = 20
        whitelist = 'in2p3.fr'
    elif 'Jets' in dataset:
        njobs = 15
        whitelist = 'wisc'
    elif 'WW' in dataset:
        njobs = 100
        whitelist = 'ultralight, purdue, polgrid'
    elif 'WZ' in dataset:
        njobs = 40
        whitelist = 'roma, ac.be'
    elif 'ZZ' in dataset:
        njobs = 20
        whitelist = 'roma, ac.be'
    elif 'tW' in dataset:
        njobs = 40
        whitelist = 'ultralight, purdue'
    elif 'tautau' in dataset:
        njobs = 50
        whitelist = 'wisc'
    elif 'DrellYan' in dataset:
        njobs = 500
        whitelist = 'ucsd, ifca'
    open(crabFn, 'wt').write(crabcfg % locals())

    # if we can expect the CSA07 process id to be there, do what we
    # need with it.
    has_CSA07_proc_id = False
    for x in ['CSA07Muon', 'CSA07Electron', 'CSA07AllEvents']:
        if x in dataset:
            has_CSA07_proc_id = True
    
    if has_CSA07_proc_id:
        weights = cmsswweights
        weights2 = 'csa07EventWeightProducer &'
    else:
        weights2 = ''

    # Filter Z+jets.
    if 'Chowder' in dataset:
        weights2 += ' csaids &'

    # if running on PDElectron, skip events that passed PDMuon so we
    # don't double count
    if 'CSA07Electron' in dataset:
        nocsa07muon = '!CSA07MuonFilter &'
    else:
        nocsa07muon = ''

    # make and write out the tmp cmssw cfg
    open(cmsswFn, 'wt').write(cmsswcfg % locals())

    if not writeCfgsOnly:
        # make the output dirs
        outdir = '/castor/cern.ch/user/t/tucker/crab/%s' % shortname
        os.system('rfmkdir -p %s' % outdir)
        os.system('nschmod 0777 %s' % outdir)

        # run crab
        os.system('crab -cfg %s -create -submit all' % crabFn)

        # mv the crab dir
        if not 'nomove' in sys.argv:
            os.system('mv crab_%s /scratchdisk1/tucker/crab' % shortname)

if not writeCfgsOnly:
    os.remove(crabFn)
    os.remove(cmsswFn)
    os.remove('glite.conf.CMS_CERN')
    os.remove('crab.history')
