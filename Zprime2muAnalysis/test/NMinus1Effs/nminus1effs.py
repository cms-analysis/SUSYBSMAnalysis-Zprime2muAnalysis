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

#### if you run on data change HLT2 in
##Zprime2muAnalysis_cff
##(bits = cms.InputTag("TriggerResults","","HLT")) #### instead of HLT2 
#### METFilterMiniAOD_cfi.py
#### src = cms.InputTag("TriggerResults","","RECO"), #### instead of PAT


readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
#process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

process.source = cms.Source ("PoolSource",
                             fileNames =  cms.untracked.vstring("/store/data/Run2016G/SingleMuon/MINIAOD/PromptReco-v1/000/278/820/00000/0667AC34-2464-E611-84CE-02163E011979.root"),
                             secondaryFileNames = secFiles)



#readFiles.extend( ['store/data/Run2016G/SingleMuon/MINIAOD/PromptReco-v1/000/278/820/00000/0667AC34-2464-E611-84CE-02163E011979.root']);

#### m4500_6000
#readFiles.extend( ['/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/2E8F7144-A63A-E611-9DB8-001D0967DFF3.root',
#       '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/4A0CB50D-A63A-E611-A67A-0090FAA58D04.root',
#       '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/88049D17-A63A-E611-A789-0090FAA573F0.root' ] );

#### m120_200
#readFiles.extend( [
#        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/18C80393-613A-E611-86DF-0090FAA573E0.root',
#        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/AA745369-613A-E611-93A7-0025907B4EE6.root',
#        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/B2BFFF97-603A-E611-90A2-0CC47AB35D34.root' ] );

#### pattuples
#readFiles.extend( [
        #                   '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/18C80393-613A-E611-86DF-0090FAA573E0.root',
        #       '/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/4E45A151-D6FC-E511-9FE4-0090FAA58134.root',
        #'/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/18C80393-613A-E611-86DF-0090FAA573E0.root',
        #                   '/store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/datamc_dy200to400/160509_211417/0000/pat_1.root',
        
 #       ] );

#### m2300_3500
#readFiles.extend( [
#        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/8C0E76AA-D93C-E611-94E2-0CC47A4D76B6.root',
#        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/B449FDAB-D93C-E611-A15B-0025905B85FE.root',
#        '/store/mc/RunIISpring16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/F23DCBB0-D93C-E611-B5D6-0025905A606A.root' ] );


#### m50_120

secFiles.extend( [
               ] )

process.maxEvents.input = 100
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'
#process.MessageLogger.cerr.FwkReport.reportEvery = 1 # default 1000

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
    ]

for name, cut in cuts:
    if type(cut) != tuple:
        cut = (cut,)
        
    lc = loose_cut
    for c in cut:
        if c not in lc:
            raise ValueError('cut "%s" not in cut string "%s"' % (c, lc))
        lc = lc.replace(' && ' + c, '') # Relies on none of the cuts above being first in the list.

    obj_no = allDimuons.clone(loose_cut = lc)
    #obj_no = allDimuons.clone(loose_cut = lc,tight_cut = tight_cut.replace(trigger_match, ''))#N-2
    setattr(process, 'allDimuonsNo' + name, obj_no)
    
    obj_ti = obj_no.clone(tight_cut = tight_cut + ' && ' + ' && '.join(cut))
    setattr(process, 'allDimuonsTi' + name, obj_ti)

process.allDimuonsNoNo      = allDimuons.clone()
#process.allDimuonsNoNo      = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))#N-2
process.allDimuonsNoTrgMtch = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))

alldimus = [x for x in dir(process) if 'allDimuonsNo' in x or 'allDimuonsTi' in x]

# Sanity check that the replaces above did something.
for x in alldimus:
    if 'NoNo' in x:
        continue
    o = getattr(process, x)
    assert o.loose_cut.value() != loose_cut or o.tight_cut.value() != tight_cut

if miniAOD:
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DileptonPreselector_cfi')####?????
    process.leptons = process.leptons_mini.clone()
    process.p = cms.Path(process.egmGsfElectronIDSequence*process.dileptonPreseletor * process.muonPhotonMatchMiniAOD * process.leptons * reduce(lambda x,y: x*y, [getattr(process, x) for x in alldimus]))
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
#process.allDimuonsN2 = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))#N-2
#process.p *= process.allDimuonsN2#N-2
process.dimuonsNoB2B     = process.dimuons.clone()
process.dimuonsNoVtxProb = process.dimuons.clone()
process.dimuonsNoDptPt   = process.dimuons.clone()
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
process.allDimuonsNoCosm = process.allDimuons.clone(loose_cut = loose_cut.replace(' && abs(dB) < 0.2', ''))
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
config.JobType.psetName = 'nminus1effs.py'
#config.JobType.priority = 1
config.Data.inputDataset =  '%(ana_dataset)s'
config.Data.inputDBS = 'global'
job_control
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/federica'
#config.Site.storageSite = 'T2_IT_Bari'
config.Site.storageSite = 'T2_IT_Legnaro'
'''

    just_testing = 'testing' in sys.argv
    if not 'no_data' in sys.argv:
        #running on miniaod we don't need of googlumis as it was
        #from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Run2016G_ll
        #Run2016G_ll.writeJSON('tmp.json')

        dataset_details = [
           # ('SingleMuonRun2015B-Prompt_251162_251499',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251162_251499_20150713100409-3aa7688518cb1f1b044caf15b1a9ed05/USER'),
            ('SingleMuonRun2016G-PromptReco-v1',  '/SingleMuon/Run2016G-PromptReco-v1/MINIAOD')
            ]

        for name, ana_dataset in dataset_details:
            print name

            new_py = open('nminus1effs.py').read()
            new_py += "\nprocess.GlobalTag.globaltag = '80X_dataRun2_Prompt_v11'\n"
            open('nminus1effs_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()
            job_control = '''
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 10 ##100
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
config.Data.unitsPerJob  = 10000
''')

        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
        #samples =[DY120to200Powheg]
        samples =[DY50to120Powheg,DY120to200Powheg,DY200to400Powheg,DY400to800Powheg,DY800to1400Powheg,DY1400to2300Powheg,DY2300to3500Powheg,DY3500to4500Powheg,DY4500to6000Powheg,DY6000toInfPowheg]#,ttbar, wz, ww_incl, zz_incl, dy50to120]
        #samples =[dy50, dy120, dy200, dy400, dy800, dy1400, dy2300, dy3500, dy4500, dy6000, dy7500, dy8500, dy9500, zpsi5000, ttbar, inclmu15]
        for sample in samples:
            #print sample.name
            open('crabConfig.py', 'wt').write(crab_cfg % sample)
            if not just_testing:
                os.system('crab submit -c crabConfig.py')
                #os.system('crab submit -c crabConfig.py --dryrun')
        if not just_testing:
            os.system('rm crabConfig.py')
