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
    for x in ['MuonsOnly', 'AllGood']:
        xl = x.lower()
        do('''
mkdir -p ana_datamc_%(extra)s/%(xl)s
hadd ana_datamc_%(extra)s/%(xl)s/ana_datamc_data.root crab/crab_ana_datamc_Run2010?%(x)s/res/*root
crab -c crab/crab_ana_datamc_Run2010A%(x)s -report
crab -c crab/crab_ana_datamc_Run2010B%(x)s -report
compareJSON.py --and crab/crab_ana_datamc_Run2010?%(x)s/res/lumiSummary.json
mergeJSON.py crab/crab_ana_datamc_Run2010?%(x)s/res/lumiSummary.json --output=ana_datamc_%(extra)s/%(xl)s/ana_datamc_data.forlumi.json
# next two if comparing to an existing copy with little differences
diff -s ana_datamc_%(extra)s/%(xl)s/ana_datamc_data.forlumi.json ana_datamc_nov4/%(xl)s/
cp ana_datamc_nov4/%(xl)s/ana_datamc_data.lumi ana_datamc_%(extra)s/%(xl)s/
# else this
#lumiCalc.py -c frontier://LumiProd/CMS_LUMI_PROD -i ana_datamc_%(extra)s/%(xl)s/ana_datamc_data.forlumi.json overview > ana_datamc_%(extra)s/%(xl)s/ana_datamc_data.lumi
''' % locals())



'''
setenv XXXDIR2 ~/nobackup/ana_datamc_mc/renameme
mkdir -p $XXXDIR2
# copy the mc
foreach y (muonsonly allgood)
  cd ana_datamc_renameme/${y}
  pwd
  foreach x (${XXXDIR2}/*.root)
    ln -s $x
  end
  cd ../..
end


Workflow to run over data without using crab:

setenv XXXDIR withNewStuff
setenv XXXTODO data # or nov4

cmsRun -j ana_datamc_data.fjr.xml histos.py $XXXTODO >&! /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out
gzip /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out
mv /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out.gz .
fjr2json.py --output=ana_datamc_data.forlumi.json ana_datamc_data.fjr.xml
lumiCalc.py -i ana_datamc_data.forlumi.json overview > ana_datamc_data.lumi
mkdir ana_datamc_${XXXDIR}/muonsonly
mv ana_datamc_data.* ana_datamc_${XXXDIR}/muonsonly/

cmsRun -j ana_datamc_data.fjr.xml histos.py $XXXTODO all_good >&! /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out
gzip /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out
mv /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out.gz .
fjr2json.py --output=ana_datamc_data.forlumi.json ana_datamc_data.fjr.xml
lumiCalc.py -i ana_datamc_data.forlumi.json overview > ana_datamc_data.lumi
mkdir ana_datamc_${XXXDIR}/allgood
mv ana_datamc_data.* ana_datamc_${XXXDIR}/allgood/

cd ana_datamc_${XXXDIR}
foreach x (`ls -1 --color=no mc/*.root`)
  cd muonsonly
  ln -sf ../${x}
  cd ../allgood
  ln -sf ../${x}
  cd ..
end
cd ..
'''
