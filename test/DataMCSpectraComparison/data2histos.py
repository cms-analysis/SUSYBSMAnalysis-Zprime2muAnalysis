import os, sys, glob

just_testing = 'testing' in sys.argv

def do(cmd):
    print cmd
    if not just_testing:
        os.system(cmd)

#in_fns = glob.glob('crab/crab_datamc_promptB*/res/merged.root')
in_fns = ['crab/crab_datamc_promptB_all/res/merged.root']
out_fn_base = 'ana_datamc/ana_datamc_data_promptB_%s.rout'
open(out_fn_base.replace('_%s.rout', '.whichcrabs'), 'wt').write(' '.join(in_fns))

#for x in ['allgood', 'zp2mu']:
for x in ['allgood']:
    lumis_fn = '/afs/cern.ch/user/t/tucker/runreg/output/full_%s.cmssw' % x
    json_fn = lumis_fn.replace('.cmssw', '.json')
    out_fn = out_fn_base % x
    fjr_fn = out_fn.replace('.rout', '.fjr.xml')
    for_lumi_fn = out_fn.replace('.rout', '.forlumi.json')
    do('cp %s %s' % (lumis_fn, out_fn.replace('.rout', '.cmssw')))
    do('cp %s %s' % (json_fn,  out_fn.replace('.rout', '.runreg.json')))
    do('cmsRun -j %s histos.py data %s %s %s' % (fjr_fn, ' '.join(in_fns), lumis_fn, out_fn))
    do('fjr2json.py --output=%s %s ' % (for_lumi_fn, fjr_fn))
    do('lumiCalc.py -i %s overview > %s' % (for_lumi_fn, out_fn.replace('.rout', '.lumi')))
