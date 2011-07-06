#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName, AODOnly
from tuple_common import cms, process, crab_cfg

AODOnly(process)

process.source.fileNames = ['file:/uscms/home/tucker/scratch/FE434AF5-D188-E011-BDE0-0017A4771000.root']
process.maxEvents.input = 100
process.GlobalTag.globaltag = 'START42_V11::All'

#switchHLTProcessName(process, 'REDIGI311X')

if __name__ == '__main__' and 'submit' in sys.argv:
    job_control = '''
total_number_of_events = -1
events_per_job = 150000
'''

    just_testing = 'testing' in sys.argv
    create_only = 'create_only' in sys.argv
    from samples import samples
    for sample in samples:
        if 'Summer11' not in sample.dataset:
            continue
        print sample.name

        new_py = open('tuple_mc.py').read()
        new_py += '\nswitchHLTProcessName(process, "%(hlt_process_name)s")\n' % sample.__dict__

        sample.pset = 'psets/tuple_mc_crab_%(name)s.py' % sample.__dict__
        open(sample.pset,'wt').write(new_py)

        sample.job_control = job_control % sample.__dict__
        open('crab.cfg', 'wt').write(crab_cfg % sample.__dict__)
        if not just_testing:
            if create_only:
                os.system('crab -create')
            else:
                os.system('crab -create -submit all')
            os.system('rm crab.cfg')
