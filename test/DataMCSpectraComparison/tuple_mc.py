#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName,AODOnly,removeMuonMCClassification,removeSimLeptons, pruneMCLeptons
from tuple_common import cms, process, crab_cfg

pruneMCLeptons(process, use_sim=True) # because of unscheduled I can't remove this for data.

AODOnly(process) #defined in PATTools and contains
#removeMuonMCClassification(process)#??? # throw the baby out with the
#removeSimLeptons(process)
#switchHLTProcessName(process, 'REDIGI311X')
#switchHLTProcessName(process, 'HLT') #default

process.source.fileNames = [
#
                            '/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/F6A03919-5232-E511-B0F4-002618943867.root'
]

process.maxEvents.input = -1

#process.GlobalTag.globaltag = 'MCRUN2_74_V9A'
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v10' #mc startup

if __name__ == '__main__' and hasattr(sys, 'argv') and 'submit' in sys.argv:
    job_control = '''
config.Data.splitting = 'EventAwareLumiBased'        
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 5000
'''

    just_testing = 'testing' in sys.argv
    #create_only = 'create_only' in sys.argv

    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name

        new_py = open('tuple_mc.py').read()
        new_py += '\nswitchHLTProcessName(process, "%(hlt_process_name)s")\n' % sample.__dict__

        sample.pset = 'crab/psets/tuple_mc_crab_%(name)s.py' % sample.__dict__
        open(sample.pset,'wt').write(new_py)

        sample.job_control = job_control % sample.__dict__
        #print sample.__dict__
        #sample.job = 'crab_%(name)s.py' % sample.__dict__
        open('crabConfig.py', 'wt').write(crab_cfg % sample.__dict__)
        if not just_testing:
            #if create_only:
                #os.system('crab submit -c ' + sample.job)
            #else:
            os.system('crab submit -c crabConfig.py')
            os.system('rm crabConfig.py')
