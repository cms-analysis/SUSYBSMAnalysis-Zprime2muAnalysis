#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName
from tuple_common import cms, process, crab_cfg

process.source.fileNames = ['/store/mc/Fall10/DYToMuMu_M-120_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/E8EE9B50-F1C8-DF11-B9B3-0018FE286F12.root']
process.maxEvents.input = 100
process.GlobalTag.globaltag = 'START38_V14::All'

if __name__ == '__main__' and 'submit' in sys.argv:
    job_control = '''
total_number_of_events = -1
events_per_job = 25000
'''

    just_testing = 'testing' in sys.argv
    from samples import samples
    for sample in samples:
        print sample.name

        new_py = open('tuple_mc.py').read()
        new_py += '\nswitchHLTProcessName(process, "%(hlt_process_name)s")\n' % sample.__dict__

        sample.pset = 'psets/tuple_mc_crab_%(name)s.py' % sample.__dict__
        open(sample.pset,'wt').write(new_py)

        sample.job_control = job_control % sample.__dict__
        open('crab.cfg', 'wt').write(crab_cfg % sample.__dict__)
        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm crab.cfg')
