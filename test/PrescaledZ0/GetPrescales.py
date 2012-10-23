import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.source.fileNames = ['/store/data/Run2012C/SingleMu/AOD/PromptReco-v2/000/202/272/94E42528-89F9-E111-BE44-BCAEC53296F8.root']
process.GlobalTag.globaltag = 'GR_P_V42_AN2::All'
process.maxEvents.input = 5000
process.options.wantSummary = True
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
process.CheckPrescale.dump_prescales = True

process.Mu17       = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu17_v3'))
process.Mu15eta2p1 = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu15_eta2p1_v3', 'HLT_Mu15_eta2p1_v4'))
process.Mu24eta2p1 = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu24_eta2p1_v3', 'HLT_Mu24_eta2p1_v4', 'HLT_Mu24_eta2p1_v5'))

process.MessageLogger.suppressWarning = cms.untracked.vstring('Mu17', 'Mu15eta2p1', 'Mu24eta2p1')

#process.p = cms.Path(process.Mu17 * process.Mu15eta2p1 * process.Mu24eta2p1)
process.p = cms.Path(process.Mu24eta2p1)

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = %(dataset)s
pset = GetPrescales.py
total_number_of_lumis = -1
number_of_jobs = 50
lumi_mask = /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-204601_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt

[USER]
ui_working_dir = crab/crab_getprescales_%(name)s
return_data = 1
'''

    just_testing = 'testing' in sys.argv

    dataset_details = [
        # ('SingleMu2012A_Prompt', '/SingleMu/Run2012A-PromptReco-v1/AOD'),
        # ('SingleMu2012B_Prompt', '/SingleMu/Run2012B-PromptReco-v1/AOD'),
        ('SingleMu2012C_24Aug2012', '/SingleMu/Run2012C-24Aug2012-v1/AOD'),
        ('SingleMu2012C_Prompt',    '/SingleMu/Run2012C-PromptReco-v2/AOD'),
        ('SingleMu2012D_Prompt',    '/SingleMu/Run2012D-PromptReco-v1/AOD'),
        ]

    for name, dataset in dataset_details:
        print name
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab -create -submit all')

    if not just_testing:
        os.system('rm -f crab.cfg tmp.json')
