#!/usr/bin/env python

import sys, os
from tuple_common import process, crab_cfg

process.source.fileNames = ['/store/data/Run2010A/Mu/RECO/Jul15thReReco-v1/0047/9C66C59B-9890-DF11-BA4A-003048678F02.root']
#process.source.fileNames = ['/store/data/Run2010A/Mu/RECO/v4/000/141/882/C453CB25-169B-DF11-B122-003048F1C58C.root']
process.maxEvents.input = 500
process.GlobalTag.globaltag = 'GR_R_37X_V6D::All'

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
removeMCUse(process)

if __name__ == '__main__' and 'submit' in sys.argv:
    scheduler = 'condor'
    job_control_ex = '''
lumi_mask = %(json)s
total_number_of_lumis = -1
lumis_per_job = 100
'''
    
    x = [
        ('jul15',          '/Mu/Run2010A-Jul15thReReco-v1/RECO', 'jsons/jul15.json',          'GR_R_37X_V6D'),
        ('prompt',         '/Mu/Run2010A-PromptReco-v4/RECO',    'jsons/prompt.json',         'GR10_P_V7'),
#        ('promptB',        '/Mu/Run2010B-PromptReco-v2/RECO',    'jsons/prompt.json',         'GR10_P_V10'),
#        ('prompt_dcsonly', '/Mu/Run2010B-PromptReco-v2/RECO',    'jsons/prompt_dcsonly.json', 'GR10_P_V10'),
        ]

    for name, dataset, json, tag in x:
        print name

        new_py = open('tuple_data.py').read()
        new_py += '\n\nprocess.GlobalTag.globaltag = "%s::All"\n' % tag

        pset = 'psets/tuple_data_crab_%(name)s.py' % locals()
        open(pset,'wt').write(new_py)

        job_control = job_control_ex % locals()
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        os.system('crab -create -submit all')

    os.system('rm crab.cfg')
