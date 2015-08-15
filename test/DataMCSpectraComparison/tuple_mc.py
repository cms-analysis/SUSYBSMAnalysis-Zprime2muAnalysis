#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName,AODOnly,removeMuonMCClassification,removeSimLeptons, pruneMCLeptons
from tuple_common import cms, process, crab_cfg

pruneMCLeptons(process, use_sim=True) # because of unscheduled I can't remove this for data.

AODOnly(process)# 
process.source.fileNames = ['/store/mc/RunIISpring15DR74/ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/FE091975-D534-E511-82AA-6C3BE5B5A4C8.root']

process.maxEvents.input = -1

process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

switchHLTProcessName(process, "HLT")


if __name__ == '__main__' and hasattr(sys, 'argv') and 'submit' in sys.argv:
    job_control = '''
total_number_of_events = -1
events_per_job = 10000
'''

    just_testing = 'testing' in sys.argv
    create_only = 'create_only' in sys.argv

    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name

        new_py = open('tuple_mc.py').read()
        new_py += '\nswitchHLTProcessName(process, "%(hlt_process_name)s")\n' % sample.__dict__

        sample.pset = 'tuple_mc_crab_%(name)s.py' % sample.__dict__
        open(sample.pset,'wt').write(new_py)

        #sample.job_control = job_control % sample.__dict__
        #print sample.__dict__

        sample.job = 'crab_%(name)s.py' % sample.__dict__
        open(sample.job, 'wt').write(crab_cfg % sample.__dict__)
        if not just_testing:
            if create_only:
                os.system('crab submit -c ' + sample.job)
            #else:
             #   os.system('crab -create -submit all')
           # os.system('rm crab.cfg')
