#!/usr/bin/env python

import datetime, sys, os, FWCore.ParameterSet.Config as cms

process = cms.Process('MergePAT')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:pat.root'))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.out = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('merged.root'))
process.outp = cms.EndPath(process.out)

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = %(ana_dataset)s
dbs_url = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
pset = tuple_merge.py
get_edm_output = 1
number_of_jobs = 20
#total_number_of_lumis = -1
total_number_of_events = -1

[USER]
ui_working_dir = crab/crab_merge_%(name)s
copy_data = 1
storage_element = T3_US_FNALLPC
check_user_remote_dir = 0
publish_data = 1
publish_data_name = merge_%(uniq)s_%(name)s
dbs_url_for_publication = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
'''

    just_testing = 'testing' in sys.argv

    uniq = datetime.datetime.today().strftime('%Y%m%d%H%M%S')

    from samples import *
    to_merge = [wjets, inclmu15]

    for sample in to_merge:
        print sample.name
        sample.uniq = uniq
        open('crab.cfg', 'wt').write(crab_cfg % sample)
        if not just_testing:
            os.system('crab -create -submit all')

    if not just_testing:
        os.system('rm crab.cfg')
