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

from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2018_cff import loose_cut, trigger_match, tight_cut, allDimuons

#### METFilterMiniAOD_cfi.py
#### src = cms.InputTag("TriggerResults","","RECO"), #### instead of PAT

process.source = cms.Source ("PoolSource",
                             fileNames =  cms.untracked.vstring(
                                 '/store/data/Run2018C/SingleMuon/MINIAOD/17Sep2018-v1/120000/FAB77E55-E1DE-0D43-A907-BD709A4B2B1D.root'
                                 ),
                             )

process.maxEvents.input = -1
process.GlobalTag.globaltag ='dummy' # dummy value to be changed below
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
    ('MuHits', '( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) )'),
    ('MuMatch', ('( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))', 'isTrackerMuon')),
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
    process.load("SUSYBSMAnalysis.Zprime2muAnalysis.EventCounter_cfi")
    process.leptons = process.leptonsMini.clone()
    process.p = cms.Path(process.egmGsfElectronIDSequence*process.EventCounter * process.dileptonPreseletor * process.muonPhotonMatchMiniAOD * process.leptons * reduce(lambda x,y: x*y, [getattr(process, x) for x in alldimus]))
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
    for dataFilter in goodDataFiltersMiniAOD:
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

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = \
'''
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'ana_nminus1_%(name)s'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nminus1effs.py'
config.Data.inputDataset =  '%(ana_dataset)s' 
config.Data.inputDBS = 'global'
job_control
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_nminus1_%(name)s'
config.Data.outLFNDirBase = '/store/group/phys_exotica/dimuon/2018/nminus1effs/crab'
config.Site.storageSite = 'T2_CH_CERN'
'''

    just_testing = 'testing' in sys.argv
    if not 'no_data' in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import *
        dataset_details = [
            # Prompt Reco for ABC
            #('SingleMuonRun2018A-06June2018-v1', '/SingleMuon/Run2018A-06Jun2018-v1/MINIAOD'), 
            #('SingleMuonRun2018A-PromptReco-v3', '/SingleMuon/Run2018A-PromptReco-v3/MINIAOD'),
            #('SingleMuonRun2018B-PromptReco-v1', '/SingleMuon/Run2018B-PromptReco-v1/MINIAOD'),
            #('SingleMuonRun2018B-PromptReco-v2', '/SingleMuon/Run2018B-PromptReco-v2/MINIAOD'),
            #('SingleMuonRun2018C-PromptReco-v1', '/SingleMuon/Run2018C-PromptReco-v1/MINIAOD'),
            #('SingleMuonRun2018C-PromptReco-v2', '/SingleMuon/Run2018C-PromptReco-v2/MINIAOD'),
            #('SingleMuonRun2018C-PromptReco-v3', '/SingleMuon/Run2018C-PromptReco-v3/MINIAOD'),

            # PPD recommendation to use 17Sep2018 ReReco for ABC, Prompt for D
            ('SingleMuonRun2018A-17Sep2018-v2', '/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD'),
            ('SingleMuonRun2018B-17Sep2018-v1', '/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD'),
            ('SingleMuonRun2018C-17Sep2018-v1', '/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD'),
            ('SingleMuonRun2018D-PromptReco-v2', '/SingleMuon/Run2018D-PromptReco-v2/MINIAOD'),
        ]
        lumi_lists = [
			'Run2018MuonsOnly',
		]
        jobs = []
        for lumi_name in lumi_lists:
            ll = eval(lumi_name + '_ll') if lumi_name != 'NoLumiMask' else None
            for dd in dataset_details:
                jobs.append(dd + (lumi_name, ll))


        for dataset_name, ana_dataset, lumi_name, lumi_list in jobs:
            print lumi_name
            json_fn = 'tmp.json'
            lumi_list.writeJSON(json_fn)
            lumi_mask = json_fn

            name = '%s_%s' % (lumi_name, dataset_name)
            print name
            print lumi_mask

            new_py = open('nminus1effs.py').read()
            if '17Sept2018' in dataset_name:
                new_py += "\nprocess.GlobalTag.globaltag = '102X_dataRun2_Sep2018Rereco_v1'\n"
            else:
                new_py += "\nprocess.GlobalTag.globaltag = '102X_dataRun2_Prompt_v11'\n"
            open('nminus1effs_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()
            job_control = \
'''
config.Data.splitting = 'Automatic'
config.Data.lumiMask = '%(lumi_mask)s'
''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py') #--dryrun
                os.system('rm crabConfig.py nminus1effs_crab.py nminus1effs_crab.pyc tmp.json')

    if not 'no_mc' in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
        for sample in samples:
            print sample.name
            print sample.dataset
            name = sample.name
            ana_dataset = sample.dataset

            new_py = open('nminus1effs.py').read()
            new_py += "\nprocess.GlobalTag.globaltag = '102X_upgrade2018_realistic_v12'\n"  #### RunH
            open('nminus1effs_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()
            job_control = \
'''
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 100000
'''
            new_crab_cfg = new_crab_cfg.replace('job_control',job_control)
            open('crabConfig.py','wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py')
                os.system('rm crabConfig.py nminus1effs_crab.py nminus1effs_crab.pyc')
