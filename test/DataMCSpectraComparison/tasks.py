#!/usr/bin/env python

import sys, os
from samples import samples

just_testing = 'testing' in sys.argv
if just_testing:
    sys.argv.remove('testing')
    
try:
    cmd, extra = sys.argv[1].lower(), sys.argv[2:]
except IndexError:
    print 'usage: utils.py command [extra]'
    sys.exit(1)

def do(cmd):
    print cmd
    if not just_testing:
        os.system(cmd)

if cmd == 'hadd':
    if extra:
        extra = '_' + extra[0]
    for sample in samples:
        name = sample.name
        do('hadd ana_datamc_%(name)s.root crab/crab_ana%(extra)s_datamc_%(name)s/res/zp2mu_histos*root' % locals())

