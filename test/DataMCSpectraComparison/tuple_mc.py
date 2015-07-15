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

process.source.fileNames = ['/store/relval/CMSSW_7_4_0/RelValZpMM_2250_13TeV_Tauola/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/AE8D58C2-14DB-E411-A038-002618943901.root',
 #                           '/store/relval/CMSSW_7_4_0/RelValZpMM_2250_13TeV_Tauola/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/B6A4F83A-14DB-E411-8A01-0025905B8596.root'
                            #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/06571208-C414-E511-B38A-003048FFCBA4.root',
]

process.maxEvents.input = 100

process.GlobalTag.globaltag = 'MCRUN2_74_V9A' # ok if you use Configuration.AlCa.GlobalTag
#process.GlobalTag.globaltag = 'MCRUN2_74_V7' # if you use condDBv2
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'
#process.GlobalTag.globaltag = 'auto:run2_mc'
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
