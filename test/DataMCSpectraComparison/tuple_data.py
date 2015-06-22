#!/usr/bin/env python

import sys, os, datetime, FWCore.ParameterSet.Config as cms
from tuple_common import process, crab_cfg

#process.source.fileNames = ['/store/data/Run2012A/SingleMu/AOD/13Jul2012-v1/00000/009C369E-85D0-E111-BD58-1CC1DE046FC0.root']
process.source.fileNames = [#'/store/data/Run2015A/SingleMuon/AOD/PromptReco-v1/000/248/036/00000/787C46EA-5714-E511-863E-02163E014609.root',
                            #'/store/data/Run2015A/EGamma/AOD/PromptReco-v1/000/246/958/00000/300FCE61-D00B-E511-9F20-02163E014686.root',
                            '/store/data/Run2015A/EGamma/RECO/PromptReco-v1/000/246/958/00000/0040F913-E80B-E511-9F86-02163E014220.root',
                            ]
#process.GlobalTag.globaltag = 'FT_53_V6C_AN4::All'
#process.GlobalTag.globaltag = 'GR_H_V33::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = 'GR_P_V56::All'
##process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:com10_2013', '')
process.maxEvents.input = 100

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
removeMCUse(process)

if __name__ == '__main__' and hasattr(sys, 'argv') and 'submit' in sys.argv:
    job_control_ex = '''
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = %(lumis_per_job)s
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange = '193093-193999' # '193093-194075'
config.Data.lumiMask = '%(lumi_mask)s'
'''

    lumis_per_job = 30
    lumi_mask = ''

    #create_only = 'create_only' in sys.argv
    just_testing = 'testing' in sys.argv
    scheduler = 'condor' if 'grid' not in sys.argv else 'glite'

    def submit(d):
        new_py = open('tuple_data.py').read()
        new_py += '\n\nprocess.GlobalTag.globaltag = "%(tag)s::All"\n' % d
        pset = 'crab/psets/tuple_data_crab_%(name)s.py' % d
        open(pset, 'wt').write(new_py)

        job_control = job_control_ex % d
        for k,v in locals().iteritems():
            d[k] = v
        open('crabConfig.py', 'wt').write(crab_cfg % d)
        if not just_testing:
            #if create_only:
                #os.system('crab -create')
            #else:
            os.system('crab submit -c crabConfig.py')
            os.system('rm -f crabConfig.py tmp.json')

    run_limits = []
    for x in sys.argv:
        try:
            run_limits.append(int(x))
        except ValueError:
            pass

    if run_limits:
        run1, run2 = run_limits
        if len(run_limits) != 2 or run1 > run2:
            raise RuntimeError('if any, must specify exactly two numeric arguments   min_run max_run  with max_run >= min_run')

        # Make up a fake lumi_mask that contains all lumis possible
        # for every run in the run range, since crab doesn't seem to
        # listen for a runselection parameter anymore.
        json = ['"%i": [[1,26296]]' % r for r in xrange(run_limits[0], run_limits[1] + 1)]
        open('tmp.json', 'wt').write('{' + ', '.join(json) + '}')
        lumi_mask = 'tmp.json'

        if run1 == 190782 and run2 == 190949:
            # Special settings for 6-Aug reprocessing of 5 runs
            dataset = '/SingleMu/Run2012A-recover-06Aug2012-v1/AOD'
            name    = 'SingleMuRun2012A-recover-06Aug2012'
            tag     = 'FT_53_V6C_AN4'
        elif run1 >= 190450 and run1 < 193752:
            dataset = '/SingleMu/Run2012A-13Jul2012-v1/AOD'
            name    = 'SingleMuRun2012A-13Jul2012'
            tag     = 'FT_53_V6C_AN4'
        elif run1 >= 193752 and run1 < 196532:
            dataset = '/SingleMu/Run2012B-13Jul2012-v1/AOD'
            name    = 'SingleMuRun2012B-13Jul2012'
            tag     = 'FT_53_V6C_AN4'
        elif run1 >= 197556 and run1 < 198914:
            dataset = '/SingleMu/Run2012C-24Aug2012-v1/AOD'
            name    = 'SingleMuRun2012C-24Aug2012'
            tag     = 'FT53_V10A_AN4'
        elif run1 == 201191 and run2 == 201191:
            # Special settings for Dec-11 reprocessing of 1 run
            # (Recovery of 134/pb for Golden JSON.)
            dataset = '/SingleMu/Run2012C-EcalRecover_11Dec2012-v1/AOD'
            name    = 'SingleMuRun2012C-EcalRecover_11Dec2012'
            tag     = 'GR_P_V42_AN2'            
        elif run1 >= 198934 and run1 < 203773:
            dataset = '/SingleMu/Run2012C-PromptReco-v2/AOD'
            name    = 'SingleMuRun2012C-Prompt'
            tag     = 'GR_P_V42_AN2'
        elif run1 == 206066 and run2 == 206066:
            dataset = '/SingleMu/Run2012D-PromptReco-v1/AOD'
            name    = 'SingleMuRun2012D-Prompt'
            tag     = 'GR_P_V42_AN2'
#        elif run1 >= 203773:
#            dataset = '/SingleMu/Run2012D-PromptReco-v1/AOD'
#            name    = 'SingleMuRun2012D-Prompt'
#            tag     = 'GR_P_V42_AN2'

        #### Run2015A ###
        if run1 >= 246908:
            dataset = '/EGamma/Run2015A-PromptReco-v1/RECO'
            name    = 'EGammaRun2015A-Prompt'
            tag     = 'GR_P_V56'
        else:
            raise ValueError("don't know how to do a run_limits production for run range [%i,%i]" % run_limits)

        name = '%s_%i_%i_%s' % (name, run_limits[0], run_limits[1], datetime.datetime.today().strftime('%Y%m%d%H%M%S'))
        print name, tag

        submit(locals())
    else:
        raise ValueError('must do a run-limits production until one dataset is closed')
        x = [
            ]
        for name, dataset, tag in x:
            submit(locals())

