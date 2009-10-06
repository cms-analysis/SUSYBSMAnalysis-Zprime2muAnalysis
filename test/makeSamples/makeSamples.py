#!/usr/bin/env python

# This script writes a CMSSW python config file to run
# GEN-SIM-DIGI-L1-DIGI2RAW-HLT and a crab.cfg for each of the samples
# in the list below, then submits the jobs. Actions are only taken if
# 'doit' is on the command line (useful to be able to inspect the
# generated scripts). It is currently CASTOR-oriented, i.e. you must
# have a directory /castor/cern.ch/user/u/username in which a
# directory for each sample will be made,
# e.g. crab/SAMPLE/GEN-...-HLT/. A python file with the appropriate
# PYTHIA parameters is looked for with name SAMPLE_cff.py in the
# fragments/ directory (these files are also compatible with running
# cmsDriver directly). The HLT.py, HLT8E29.py, and reco.py files are
# cleaned-up versions of what cmsDriver output at the time of their
# writing.

num_events = 20000
num_jobs = 25

samples = [
    'PYTHIA6_ZPSSMmumu_M1000_Mcut400_10TeV',
    'PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV',
    'PYTHIA6_ZPSSMmumu_M1300_Mcut600_10TeV',
    'PYTHIA6_DYmumu_Mcut200_10TeV',
    'PYTHIA6_DYmumu_Mcut500_10TeV',
    'PYTHIA6_Z0mumu_Mcut40_10TeV',
    'PYTHIA6_RS1Gmumu_M1250_Mcut600_c0m05_10TeV',
    ]

########################################################################

import os, sys

cmssw_version = os.environ['CMSSW_VERSION']
doit = 'doit' in sys.argv

# Make sure the fragments directory is python-importable.
os.system('touch fragments/__init__.py')

def cmd(cmd):
    print cmd
    if doit:
        os.system(cmd)
    else:
        print 'NOT running'
    
for sample in samples:
    print 'Sample:', sample
    
    # Prepare the output directory.
    output_dir = '/user/t/tucker/crab/%(cmssw_version)s/%(sample)s/GEN-SIM-DIGI-L1-DIGI2RAW-HLT/' % locals()
    output_dir_abs = '/castor/cern.ch' + output_dir
    cmd('nsmkdir -p %(output_dir_abs)s ; nschmod 0777 %(output_dir_abs)s' % locals())

    # Prepare the CMSSW cfg.
    cmssw_cfg = open('fragments/HLT.py').read()
    cmssw_cfg += '''
from fragments.%(sample)s_cff import generator
process.generator = generator
'''
    open('HLT.py', 'wt').write(cmssw_cfg % locals())
    
    # Prepare the CRAB cfg.
    crab_cfg = open('fragments/crab.cfg').read()
    open('crab.cfg', 'wt').write(crab_cfg % locals())

    # Create and submit the jobs.
    cmd('crab -create -submit all')

    # Clean up.
    cmd('rm crab.cfg ; rm HLT.py')
