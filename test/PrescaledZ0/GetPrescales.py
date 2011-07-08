import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.source.fileNames = ['/store/data/Run2011A/SingleMu/AOD/May10ReReco-v1/0000/E6E8249F-9E7D-E011-B106-002481A8A782.root']
process.GlobalTag.globaltag = 'GR_R_42_V13::All'
process.maxEvents.input = 5000
process.options.wantSummary = True
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

CheckPrescale = cms.EDAnalyzer('CheckPrescale',
                               hlt_process_name = cms.string('HLT'),
                               trigger_paths = cms.vstring(),
                               dump_prescales = cms.untracked.bool(True),
                               )

x = [
    (3,  (3,7)),
    (8,  (1,5)),
    (15, (2,6)),
    ]

objs = []
names = []
for pt, (v0, v1) in x:
    obj = CheckPrescale.clone(trigger_paths=cms.vstring(*['HLT_Mu%i_v%i' % (pt, ver) for ver in xrange(v0, v1+1)]))
    name = 'Mu%i' % pt
    names.append(name)
    setattr(process, name, obj)
    objs.append(obj)

process.MessageLogger.suppressWarning = cms.untracked.vstring(*names)

process.p = cms.Path(reduce(lambda x,y: x*y, objs))
    
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = glite

[CMSSW]
datasetpath = %(dataset)s
pset = GetPrescales.py
total_number_of_lumis = -1
number_of_jobs = 50
lumi_mask = tmp.json

[USER]
ui_working_dir = crab/crab_getprescales_%(name)s
return_data = 1

[GRID]
ce_white_list = T2_EE_Estonia,T2_CH_CSCS,T2_BR_UERJ,T2_BR_SPRACE
'''

    just_testing = 'testing' in sys.argv

    dataset_details = [
        ('SingleMu2011A_May10',  '/SingleMu/Run2011A-May10ReReco-v1/AOD'),
        ('SingleMu2011A_Prompt', '/SingleMu/Run2011A-PromptReco-v4/AOD'),
        ]

    json = ['"%i": [[1,26296]]' % r for r in xrange(160431, 167913+1)]
    open('tmp.json', 'wt').write('{' + ', '.join(json) + '}')

    for name, dataset in dataset_details:
        print name
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab -create -submit all')

    if not just_testing:
        os.system('rm crab.cfg tmp.json')
