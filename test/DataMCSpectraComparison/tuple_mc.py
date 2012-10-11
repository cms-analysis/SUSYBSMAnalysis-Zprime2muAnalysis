#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName, AODOnly
from tuple_common import cms, process, crab_cfg

AODOnly(process)
#switchHLTProcessName(process, 'REDIGI311X')

process.source.fileNames = ['/store/mc/Summer12_DR53X/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/FE93A307-F2E2-E111-BC0D-E41F131815FC.root']
process.maxEvents.input = 1000
process.GlobalTag.globaltag = 'START53_V11::All'

if __name__ == '__main__' and hasattr(sys, 'argv') and 'submit' in sys.argv:
    job_control = '''
total_number_of_events = -1
events_per_job = 150000
'''

    just_testing = 'testing' in sys.argv
    create_only = 'create_only' in sys.argv

    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name

        new_py = open('tuple_mc.py').read()
        new_py += '\nswitchHLTProcessName(process, "%(hlt_process_name)s")\n' % sample.__dict__

        sample.pset = 'crab/psets/tuple_mc_crab_%(name)s.py' % sample.__dict__
        open(sample.pset,'wt').write(new_py)

        sample.job_control = job_control % sample.__dict__
        open('crab.cfg', 'wt').write(crab_cfg % sample.__dict__)
        if not just_testing:
            if create_only:
                os.system('crab -create')
            else:
                os.system('crab -create -submit all')
            os.system('rm crab.cfg')
