#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName
from tuple_common import cms, process, crab_cfg

process.source.fileNames = ['/store/mc/Spring10/TTbarJets-madgraph/GEN-SIM-RECO/START3X_V26_S09-v1/0016/6E7C4631-9D47-DF11-96CE-003048C69288.root']
#process.source.fileNames = ['/store/mc/Summer10/Zmumu_M20_CTEQ66-powheg/GEN-SIM-RECO/START36_V9_S09-v2/0032/FE2F8537-9A80-DF11-B80E-0026189438C9.root']
#process.maxEvents.input = 1000
process.GlobalTag.globaltag = 'START37_V6::All'

save_discriminatorSources = process.patJets.discriminatorSources[:]
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run36xOn35xInput
run36xOn35xInput(process, 'ak5GenJets')

if __name__ == '__main__' and 'submit' in sys.argv:
    job_control = '''
total_number_of_events = -1
events_per_job = 25000
'''

    from samples import samples
    for sample in samples:
        print sample.name

        new_py = open('tuple_mc.py').read()
        if not sample.is_35x:
            new_py += '\nprocess.patJets.discriminatorSources = save_discriminatorSources\n'
        new_py += '\nswitchHLTProcessName(process, "%(hlt_process_name)s")\n' % sample.__dict__

        sample.pset = 'psets/tuple_mc_crab_%(name)s.py' % sample.__dict__
        open(sample.pset,'wt').write(new_py)

        sample.job_control = job_control % sample.__dict__
        open('crab.cfg', 'wt').write(crab_cfg % sample.__dict__)
        os.system('crab -create -submit all')

    os.system('rm crab.cfg')
