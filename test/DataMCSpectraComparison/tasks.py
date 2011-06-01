#!/usr/bin/env python

import sys, os, glob
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
        cmds = cmd.split('\n') if '\n' in cmd else [cmd]
        for c in cmds:
            if cmd != '' and not cmd.startswith('#'):
                os.system(c)

if cmd == 'setdirs':
    do('''ln -s /uscms/home/tucker/nobackup/crab_dirs crab
mkdir -p psets
ln -s `pwd`/psets crab/psets''')

elif cmd == 'checkevents':
    for sample in samples:
        print sample.name
        do('grep TrigReport crab/crab_datamc_%s/res/*stdout | grep \' p$\' | sed -e "s/ +/ /g" | awk \'{ s += $4; t += $5; u += $6; } END { print "summary: total: ", s, "passed: ", t, "failed: ", u }\'' % sample.name)

elif cmd == 'publishmc':
    for sample in samples:
        do('crab -c crab/crab_datamc_${x} -publish >&! crab/publish_logs/publish.${x} &')

elif cmd == 'hadd':
    extra = '_' + extra[0] if extra else ''
    for sample in samples:
        name = sample.name
        pattern = 'crab/crab_ana%(extra)s_datamc_%(name)s/res/zp2mu_histos*root' % locals()
        fn = 'ana_datamc_%(name)s.root' % locals()
        n = len(glob.glob(pattern))
        if n == 0:
            raise RuntimeError('no files matching %s' % pattern)
        elif n == 1:
            do('cp %s %s' % (pattern, fn))
        else:
            do('hadd ana_datamc_%(name)s.root crab/crab_ana%(extra)s_datamc_%(name)s/res/zp2mu_histos*root' % locals())

elif cmd == 'gatherhistos':
    extra = extra[0] if extra else 'renameme'
    which = 'Run2011APlusDCSOnlyMuonsOnly'
    which = 'Run2011A'
    dirs = glob.glob('crab/crab_ana_datamc_%s_SingleMu2011A_*' % which)
    files_glob = ' '.join([os.path.join(x, 'res/*.root') for x in dirs])
    do('''
mkdir -p ana_datamc_%(extra)s
ln -s /uscms_data/d2/tucker/zp2mu_ana_datamc_mc/V00-10-12 ana_datamc_%(extra)s/mc
hadd ana_datamc_%(extra)s/ana_datamc_data.root %(files_glob)s
''' % locals())
    for dir in dirs:
        do('crab -c %s -report' % dir)
    jsons = [os.path.join(dir, 'res/lumiSummary.json') for dir in dirs]
    for i,json1 in enumerate(jsons):
        for json2 in jsons[i+1:]:
            # probably a better way to do this using LumiList but oh well
            print 'checking overlap between', json1, json2
            if os.popen('compareJSON.py --and %s %s' % (json1, json2)).read() != '{}\n':
                raise RuntimeError('overlap between %s and %s lumisections' % (json1, json2))
    jsons = ' '.join(jsons)
    do('''
mergeJSON.py %(jsons)s --output ana_datamc_%(extra)s/ana_datamc_data.forlumi.json
lumiCalc.py -i ana_datamc_%(extra)s/ana_datamc_data.forlumi.json overview > ana_datamc_%(extra)s/ana_datamc_data.lumi
''' % locals())

elif cmd == 'mclinks':
    extra = extra[0] if extra else 'renameme'
    for x in ['muonsonly', 'allgood']:
        print('ln -s ~/nobackup/ana_datamc_mc/%s %s' % (extra, x))
