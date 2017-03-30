#!/usr/bin/env python

miniAOD = True

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD
if miniAOD:
    from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
    electrons_miniAOD(process)
    from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
else:
    from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff import loose_cut, trigger_match, tight_cut, allDimuons
#from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff import loose_cut, trigger_match, tight_cut, allDimuons

#### if you run on data change HLT2 in
##Zprime2muAnalysis_cff
##(bits = cms.InputTag("TriggerResults","","HLT")) #### instead of HLT2 
#### METFilterMiniAOD_cfi.py
#### src = cms.InputTag("TriggerResults","","RECO"), #### instead of PAT


readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
#process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

process.source = cms.Source ("PoolSource",
                             fileNames =  cms.untracked.vstring(
# '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00312D7A-FEBD-E611-A713-002590DB923E.root',
# '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00355459-F4BD-E611-B5E8-D4AE526A11F3.root',
# '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00BCC036-F6BD-E611-92F5-0025905A6118.root',
# '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/028EE230-F7BD-E611-BCDB-0CC47AA9906E.root',
# '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/02C57965-0DBE-E611-91EF-70106F4A9254.root',


# '/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/06BD3929-EDC7-E611-8EC3-02163E019C96.root',
# '/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/88E1539B-0AC8-E611-B36F-02163E0140DC.root'

#                                                                 '/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/162AD1DB-1E98-E611-9893-008CFA56D58C.root',
                                                                #'/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/1A1F07FF-2698-E611-915C-0242AC130004.root'
#                                                                '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/18C80393-613A-E611-86DF-0090FAA573E0.root'
                                                                ),
                             secondaryFileNames = secFiles)



secFiles.extend( [
               ] )

process.maxEvents.input = -1
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v6'
#process.MessageLogger.cerr.FwkReport.reportEvery = 1 # default 1000

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('TauTauSelection',
                                       src = cms.InputTag('prunedGenParticles'),
                                       )

# Define the numerators and denominators, removing cuts from the
# allDimuons maker. "NoX" means remove cut X entirely (i.e. the
# loose_cut denominators), "TiX" means move cut X from the loose_cut
# to the tight_cut (meaning only one muon instead of two has to pass
# the cut).  "NoNo" means remove nothing (i.e. the numerator). This
# will break if loose_, tight_cut strings are changed upstream, so we
# try to check those with a simple string test below.

cuts = [
    ('Pt',      'pt > 53'),
    ('DB',      'abs(dB) < 0.2'),
    ('Iso',     'isolationR03.sumPt / innerTrack.pt < 0.10'),
    ('TkLayers','globalTrack.hitPattern.trackerLayersWithMeasurement > 5'),
    ('PxHits',  'globalTrack.hitPattern.numberOfValidPixelHits >= 1'),
    ('MuHits',  'globalTrack.hitPattern.numberOfValidMuonHits > 0'),
    ('MuMatch', ('numberOfMatchedStations > 1', 'isTrackerMuon')),
    #('MuMatch', ('( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))', 'isTrackerMuon')),
    ]

for name, cut in cuts:
    if type(cut) != tuple:
        cut = (cut,)
        
    lc = loose_cut
    for c in cut:
        if c not in lc:
            raise ValueError('cut "%s" not in cut string "%s"' % (c, lc))
        lc = lc.replace(' && ' + c, '') # Relies on none of the cuts above being first in the list.

    obj_no = allDimuons.clone(cut = 'daughter(0).pdgId() + daughter(1).pdgId() == 0', loose_cut = lc)
    #obj_no = allDimuons.clone(loose_cut = lc,tight_cut = tight_cut.replace(trigger_match, ''))#N-2
    setattr(process, 'allDimuonsNo' + name, obj_no)
    
    obj_ti = obj_no.clone(tight_cut = tight_cut + ' && ' + ' && '.join(cut))
    setattr(process, 'allDimuonsTi' + name, obj_ti)

process.allDimuonsNoNo      = allDimuons.clone(cut = 'daughter(0).pdgId() + daughter(1).pdgId() == 0')
#process.allDimuonsNoNo      = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))#N-2
process.allDimuonsNoTrgMtch = allDimuons.clone(cut = 'daughter(0).pdgId() + daughter(1).pdgId() == 0', tight_cut = tight_cut.replace(trigger_match, ''))

alldimus = [x for x in dir(process) if 'allDimuonsNo' in x or 'allDimuonsTi' in x]

# Sanity check that the replaces above did something.
for x in alldimus:
    if 'NoNo' in x:
        continue
    o = getattr(process, x)
    assert o.loose_cut.value() != loose_cut or o.tight_cut.value() != tight_cut

if miniAOD:
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DileptonPreselector_cfi')####?????
    process.leptons = process.leptonsMini.clone()
    process.p = cms.Path(process.egmGsfElectronIDSequence*process.dileptonPreseletor *  process.DYGenMassFilter * process.muonPhotonMatchMiniAOD * process.leptons * reduce(lambda x,y: x*y, [getattr(process, x) for x in alldimus]))
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
    for dataFilter in goodDataFiltersMiniAOD:
        #setattr(process,dataFilter
        process.p *= dataFilter
else:
    process.leptons = process.leptons.clone()
    process.p = cms.Path(process.goodDataFilter * process.muonPhotonMatch * process.leptons * reduce(lambda x,y: x*y, [getattr(process, x) for x in alldimus]))



# For all the allDimuons producers, make dimuons producers, and
# analyzers to make the histograms.
for alld in alldimus:
    dimu = process.dimuons.clone(src = alld)
    name = alld.replace('allD', 'd')
    setattr(process, name, dimu)
    hists = HistosFromPAT.clone(dilepton_src = name, leptonsFromDileptons = True)
    setattr(process, name.replace('dimuons', ''), hists)
    process.p *= dimu * hists

# Handle the cuts that have to be applied at the
# Zprime2muCompositeCandidatePicker level.
process.allDimuonsBASE = allDimuons.clone(cut = 'daughter(0).pdgId() + daughter(1).pdgId() == 0')
process.p *= process.allDimuonsBASE
process.dimuonsNoB2B     = process.dimuons.clone(src = 'allDimuonsBASE')
process.dimuonsNoVtxProb = process.dimuons.clone(src = 'allDimuonsBASE')
process.dimuonsNoDptPt   = process.dimuons.clone(src = 'allDimuonsBASE')
#process.dimuonsNoB2B     = process.dimuons.clone(src = 'allDimuonsN2')#N-2
#process.dimuonsNoVtxProb = process.dimuons.clone(src = 'allDimuonsN2')#N-2
#process.dimuonsNoDptPt   = process.dimuons.clone(src = 'allDimuonsN2')#N-2
delattr(process.dimuonsNoB2B,     'back_to_back_cos_angle_min')
delattr(process.dimuonsNoVtxProb, 'vertex_chi2_max')
delattr(process.dimuonsNoDptPt,   'dpt_over_pt_max')

process.p *= process.allDimuons
for dimu in ['dimuonsNoB2B', 'dimuonsNoVtxProb', 'dimuonsNoDptPt']:
    hists = HistosFromPAT.clone(dilepton_src = dimu, leptonsFromDileptons = True)
    setattr(process, dimu.replace('dimuons', ''), hists)
    process.p *= getattr(process, dimu) * hists

# Special case to remove |dB| and B2B cuts simultaneously, as they can
# be correlated (anti-cosmics).
process.allDimuonsNoCosm = process.allDimuons.clone(cut = 'daughter(0).pdgId() + daughter(1).pdgId() == 0', loose_cut = loose_cut.replace(' && abs(dB) < 0.2', ''))
#process.allDimuonsNoCosm = process.allDimuons.clone(loose_cut = loose_cut.replace(' && abs(dB) < 0.2', ''), tight_cut = tight_cut.replace(trigger_match, '')) #N-2
process.dimuonsNoCosm = process.dimuons.clone(src = 'allDimuonsNoCosm')
delattr(process.dimuonsNoCosm, 'back_to_back_cos_angle_min')
process.NoCosm = HistosFromPAT.clone(dilepton_src = 'dimuonsNoCosm', leptonsFromDileptons = True)
process.p *= process.allDimuonsNoCosm * process.dimuonsNoCosm * process.NoCosm

f = file('outfile', 'w')
f.write(process.dumpPython())
f.close()
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'ana_nminus1_%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nminus1effs_TauTau.py'
#config.JobType.priority = 1
#config.Data.inputDataset =  '%(ana_dataset)s' #for pattuples
config.Data.inputDataset =  '%(dataset)s' # for miniAOD
config.Data.inputDBS = 'global'
job_control
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/ferrico'
#config.Site.storageSite = 'T2_IT_Bari'
config.Site.storageSite = 'T2_IT_Bari'
'''

    just_testing = 'testing' in sys.argv
    if not 'no_data' in sys.argv:
        #running on miniaod we don't need of googlumis as it was
        #from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Run2016G_ll
        #Run2016G_ll.writeJSON('tmp.json')

        dataset_details = [
                           #                 ('SingleMuonRun2016F-ReReco-v1', '/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD'),
                           #                ('SingleMuonRun2016E-ReReco-v1', '/SingleMuon/Run2016E-23Sep2016-v1/MINIAOD'),
                           #               ('SingleMuonRun2016D-ReReco-v1', '/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD'),
                           #              ('SingleMuonRun2016C-ReReco-v1','/SingleMuon/Run2016C-23Sep2016-v1/MINIAOD')
                           #             ('SingleMuonRun2016B-ReReco-v3', '/SingleMuon/Run2016B-23Sep2016-v3/MINIAOD')
                           #            ('SingleMuonRun2016G-PromptReco-v1',  '/SingleMuon/Run2016G-PromptReco-v1/MINIAOD')
                           #              ('SingleMuonRun2016H-PromptReco-v3', '/SingleMuon/Run2016H-PromptReco-v3/MINIAOD'),
#                            ('SingleMuonRun2016H-PromptReco-v2', '/SingleMuon/Run2016H-PromptReco-v2/MINIAOD') 
                           #	    ('SingleMuonRun2016G-ReReco-v1', '/SingleMuon/Run2016G-23Sep2016-v1/MINIAOD')
            ]

#        for name, ana_dataset in dataset_details: # for pattuples
        for name, dataset in dataset_details:

            print name

            new_py = open('nminus1effs.py').read()
            new_py += "\nprocess.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v6'\n"
            open('nminus1effs_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()
            job_control = '''
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'tmp.json' #######
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-282037_13TeV_PromptReco_Collisions16_JSON_NoL1T_MuonPhys.txt'
'''
            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py ') #--dryrun

        if not just_testing:
            os.system('rm crabConfig.py nminus1effs_crab.py nminus1effs_crab.pyc tmp.json')

    if not 'no_mc' in sys.argv:
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 100000
''')

        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
        #samples =[DY120to200Powheg]
        samples =[ dyInclusive50 ]
        
        for sample in samples:
            #print sample.name
            open('crabConfig.py', 'wt').write(crab_cfg % sample)
            if not just_testing:
                os.system('crab submit -c crabConfig.py')
                #os.system('crab submit -c crabConfig.py --dryrun')
        if not just_testing:
            os.system('rm crabConfig.py')
