#!/usr/bin/env python

import os

samples = [
    'ZPSSMmumu_M1000_Mcut400',
    'ZPSSMmumu_M1500_Mcut600',
    'ZPSSMmumu_M2000_Mcut1000',
    'ZPSSMmumu_M3000_Mcut1500',
    'DYmumu_M200_filter',
    'DYmumu_M500_filter',
    'DYmumu_M1000_filter',
    'DYmumu_M2000_filter'
    ]

#samples = samples[0:1]

throughGEN = False
throughHLT = True
HLT2RECO = False

events = 1
tag = 'IDEAL_V9::All'

for sample in samples:
    if 'DY' in sample:
        extra = ':ProductionFilterSequence'
    else:
        extra = ''

    if throughGEN:
        cmd = 'cmsDriver.py Configuration/GenProduction/python/PYTHIA6_%(sample)s_10TeV_cff.py -s GEN%(extra)s --eventcontent RAWSIM --datatier GEN --conditions FrontierConditions_GlobalTag,%(tag)s -n %(events)i --no_exec'
    elif throughHLT:
        cmd = 'cmsDriver.py Configuration/GenProduction/python/PYTHIA6_%(sample)s_10TeV_cff.py -s GEN%(extra)s,SIM,DIGI,L1,DIGI2RAW,HLT --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions FrontierConditions_GlobalTag,%(tag)s -n %(events)i --no_exec'
    elif HLT2RECO:
        cmd = 'cmsDriver.py Configuration/GenProduction/python/PYTHIA6_%(sample)s_10TeV_cff.py -s RAW2DIGI,RECO --filein file:hlt/PYTHIA6_%(sample)s_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.root --eventcontent RECO --datatier RECO --conditions FrontierConditions_GlobalTag,%(tag)s -n %(events)i --no_exec'
    os.system(cmd % locals())
