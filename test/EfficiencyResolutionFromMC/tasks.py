#!/usr/bin/env python

import sys, os, glob
import das_client
import json
from itertools import combinations
from FWCore.PythonUtilities.LumiList import LumiList
from SUSYBSMAnalysis.Zprime2muAnalysis.hadd import hadd
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn

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
    ret = []
    if not just_testing:
        cmds = cmd.split('\n') if '\n' in cmd else [cmd]
        for cmd in cmds:
            if cmd != '' and not cmd.startswith('#'):
                ret.append(os.system(cmd))
    if len(ret) == 1:
        ret = ret[0]
    return ret

if cmd == 'setdirs':
    crab_dirs_location = extra[0]
    do('mkdir -p ' + os.path.join(crab_dirs_location, 'psets'))
    do('ln -s %s crab' % crab_dirs_location)
    do('ln -s . crab/crab') # this is so crab -publish won't complain about being unable to find the pset if you launch it from crab/
    do('mkdir crab/publish_logs')

elif cmd == 'maketagdirs':
    extra = extra[0]
    do('rm data mc plots')
    for which in ['data', 'mc', 'plots']:
        #d = '~/nobackup/zp2mu_ana_effres_%s/%s' % (which,extra)
        d = './zp2mu_ana_effres_%s/%s' % (which,extra)
        do('mkdir -p %s' % d)
        do('ln -s %s %s' % (d, which))

elif cmd == 'checkevents':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('grep TrigReport crab/crab_effres_%s/res/*stdout | grep \' p$\' | sed -e "s/ +/ /g" | awk \'{ s += $4; t += $5; u += $6; } END { print "summary: total: ", s, "passed: ", t, "failed: ", u }\'' % sample.name)

elif cmd == 'checkstatus':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('crab status -d crab/crab_ana_effres_%(name)s*' %sample)
        
elif cmd == 'report':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('crab report -d crab/crab_ana_effres_%(name)s* ' % sample)
        
elif cmd == 'resubmit':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('crab resubmit -d crab/crab_ana_effres_%(name)s* ' % sample)

elif cmd == 'kill':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('crab kill -d crab/crab_ana_effres_%(name)s* ' % sample)

elif cmd == 'getoutput':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('crab getoutput -d crab/crab_ana_effres_%(name)s* --checksum=no ' % sample)

elif cmd == 'gathermc':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    extra = '_' + extra[0] if extra else ''
    for sample in samples:
        name = sample.name
        pattern = 'crab/crab_ana%(extra)s_effres_%(name)s*/results/zp2mu_histos*root' % locals()
        fn = 'ana_effres_%(name)s.root' % locals()
        n = len(glob.glob(pattern))
        if n == 0:
            big_warn('no files matching %s' % pattern)
        else:
            files = glob.glob('crab/crab_ana%(extra)s_effres_%(name)s*/results/zp2mu_histos*root' % locals())
            hadd('mc/ana_effres_%s.root' % name, files)


