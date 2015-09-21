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
#                            '/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/025D8DD1-5532-E511-8B88-003048FFD730.root',
#                            '/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/5454A03F-5432-E511-B785-0025905B858A.root',
#                            '/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/5C0A93C3-5132-E511-AB6D-0025905A607E.root',
#                            '/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/BE773DF6-5432-E511-9C12-0025905B8562.root',
#                            '/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/DE2E4C8D-5A32-E511-B9FB-002618943922.root',
#                            '/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/E0ACB362-5232-E511-A9B5-00261894385A.root',
#                            '/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/F4423476-5A32-E511-9D17-0025905B8592.root',
                            
                            '/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/3C26ACDC-5032-E511-A912-00261894397A.root',
                            '/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/5CF3EBBB-4E32-E511-9E9C-002590596468.root',
                            '/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/9622B63D-5D32-E511-87DE-00261894389D.root',
                            '/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/BAB18E43-5D32-E511-A01C-0025905A6068.root',
                            '/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/DEF5E014-4F32-E511-9D63-0025905A48D8.root',
                            '/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/GEN-SIM-RECO/74X_mcRun2_startup_realistic50ns_v0_trackPog2015Jul24-v1/00000/F6A03919-5232-E511-B0F4-002618943867.root'
]

process.maxEvents.input = -1
#
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A' # ok if you use Configuration.AlCa.GlobalTag

#process.GlobalTag.globaltag = 'MCRUN2_74_V7' # if you use condDBv2
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'
#process.GlobalTag.globaltag = 'auto:run2_mc'
process.GlobalTag.globaltag = '74X_mcRun2_startup_realistic50ns_v0' #mc startup

if __name__ == '__main__' and hasattr(sys, 'argv') and 'submit' in sys.argv:
    job_control = '''
config.Data.splitting = 'EventAwareLumiBased'        
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 10000
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
